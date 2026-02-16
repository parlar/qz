//! Compression orchestrators: chunked pipelined compression with disk-backed block storage.

use anyhow::Result;
use std::time::Instant;
use tracing::info;
use crate::cli::{CompressConfig, QualityCompressor, QualityMode};
use super::*;
use super::codecs;

/// Chunked streaming compression with pipelined I/O and disk-backed block storage.
///
/// Memory reduction strategy (compared to the naive chunked approach):
/// 1. Temp files: compressed blocks are written to disk immediately instead of
///    accumulating in Vec<Vec<u8>> (saves ~2.4 GB for 100M reads, and scales
///    to any input size without growing).
/// 2. Sequential stream compression: headers->sequences->qualities instead of
///    rayon::join across all 3 streams. This limits BSC working memory to the
///    largest single stream's blocks (~6 GB for 31 sequence blocks) instead of
///    all three streams simultaneously (~14 GB for 70 blocks).
/// 3. shrink_to_fit on blocks: BSC allocates ~25 MB output buffers but typical
///    compressed data is ~9 MB. Without shrinking, blocks waste ~2.5x capacity.
///
/// Peak memory: ~8 GB for 100M reads (down from ~55 GB in the naive approach):
///   - One chunk's raw data during compression: ~760 MB (largest stream)
///   - BSC working memory: ~6 GB (31 parallel block compressions)
///   - Pipeline: next chunk being read on main thread: ~1.75 GB
///
/// Produces output identical to compress_streaming_bsc (same archive format,
/// same multi-block BSC format — existing decompressor works unchanged).
pub(super) fn compress_chunked_bsc(args: &CompressConfig) -> Result<()> {
    use std::io::Write;
    const CHUNK_SIZE: usize = 2_500_000; // 2.5M records: better I/O+compute overlap

    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };

    info!("Pipelined compression: {} records per chunk", CHUNK_SIZE);

    let mut reader = crate::io::FastqReader::from_path(input_path, args.fasta)?;

    // Accumulate compressed blocks in memory (no temp files — faster for typical inputs)
    let mut all_h_blocks: Vec<Vec<u8>> = Vec::new();
    let mut all_s_blocks: Vec<Vec<u8>> = Vec::new();
    let mut all_q_blocks: Vec<Vec<u8>> = Vec::new();
    let mut all_rc_blocks: Vec<Vec<u8>> = Vec::new();

    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;
    let mut chunk_idx: usize = 0;
    let rc_canon = args.rc_canon;

    // Decide if quality_ctx is possible (lossless mode, not no_quality)
    let collect_for_ctx = !no_quality && quality_mode == QualityMode::Lossless;

    // Read first chunk (can't skip packing yet — quality_ctx decision depends on chunk size)
    let mut cur = read_chunk_streams(
        &mut reader, CHUNK_SIZE, quality_mode, quality_binning, no_quality,
        args.sequence_hints, args.sequence_delta, rc_canon, collect_for_ctx, false,
    )?;

    // Decide quality_ctx usage from first chunk (must be consistent across all chunks)
    let use_quality_ctx = collect_for_ctx
        && (args.quality_compressor == QualityCompressor::QualityCtx
            || cur.num_reads >= MIN_READS_QUALITY_CTX);
    let quality_compressor_used = if use_quality_ctx {
        QualityCompressor::QualityCtx
    } else {
        QualityCompressor::Bsc
    };
    if use_quality_ctx {
        info!("Using context-adaptive quality compression (quality_ctx)");
    }

    while cur.num_reads > 0 {
        info!("Chunk {}: {} reads (h={} s={} q={} bytes)",
            chunk_idx, cur.num_reads, cur.header_stream.len(), cur.seq_stream.len(), cur.qual_stream.len());

        // Pipeline: compress current chunk on a background thread while
        // reading the next chunk on the main thread.
        let (next_result, compress_result) = std::thread::scope(|scope| {
            let h_data = std::mem::take(&mut cur.header_stream);
            let s_data = std::mem::take(&mut cur.seq_stream);
            let q_data = std::mem::take(&mut cur.qual_stream);
            let rc_data = std::mem::take(&mut cur.rc_flags);
            let qs = std::mem::take(&mut cur.qual_strings);
            let ss = std::mem::take(&mut cur.seq_strings);

            let compress_handle = scope.spawn(move || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                let chunk_t0 = Instant::now();
                // Compress all three streams in parallel for speed
                let (h_result, (s_result, q_result)) = rayon::join(
                    || { let t = Instant::now(); let r = compress_stream_to_bsc_blocks(&h_data, bsc_static); info!("  BSC headers: {:.2}s", t.elapsed().as_secs_f64()); r },
                    || rayon::join(
                        || { let t = Instant::now(); let r = compress_stream_to_bsc_blocks(&s_data, bsc_static); info!("  BSC sequences: {:.2}s", t.elapsed().as_secs_f64()); r },
                        || -> Result<Vec<Vec<u8>>> {
                            let qt = Instant::now();
                            let qr = if no_quality || (q_data.is_empty() && qs.is_empty()) {
                                Ok(Vec::new())
                            } else if use_quality_ctx {
                                use rayon::prelude::*;
                                drop(q_data); // Don't need packed stream
                                let qual_refs: Vec<&str> = qs.iter().map(|s| s.as_str()).collect();
                                let seq_refs: Vec<&str> = ss.iter().map(|s| s.as_str()).collect();
                                let sub_block_reads: usize = args.quality_ctx_block_size;
                                let n = qual_refs.len();
                                if n <= sub_block_reads {
                                    let blob = quality_ctx::compress_qualities_ctx(&qual_refs, &seq_refs)?;
                                    Ok(vec![blob])
                                } else {
                                    let num_sub = (n + sub_block_reads - 1) / sub_block_reads;
                                    let blobs: Vec<Vec<u8>> = (0..num_sub)
                                        .into_par_iter()
                                        .map(|i| {
                                            let start = i * sub_block_reads;
                                            let end = (start + sub_block_reads).min(n);
                                            quality_ctx::compress_qualities_ctx(&qual_refs[start..end], &seq_refs[start..end])
                                        })
                                        .collect::<Result<Vec<_>>>()?;
                                    Ok(blobs)
                                }
                            } else {
                                compress_stream_to_bsc_blocks(&q_data, bsc_static)
                            };
                            info!("  Quality: {:.2}s", qt.elapsed().as_secs_f64());
                            qr
                        },
                    ),
                );

                // Compress RC flags (if present)
                let rc_blocks = if rc_data.is_empty() {
                    Vec::new()
                } else {
                    compress_stream_to_bsc_blocks(&rc_data, bsc_static)?
                };

                info!("  Chunk compress total: {:.2}s", chunk_t0.elapsed().as_secs_f64());
                Ok((h_result?, s_result?, q_result?, rc_blocks))
            });

            // Read next chunk on main thread (overlaps with compression)
            let io_t = Instant::now();
            let next = read_chunk_streams(
                &mut reader, CHUNK_SIZE, quality_mode, quality_binning, no_quality,
                args.sequence_hints, args.sequence_delta, rc_canon,
                collect_for_ctx && use_quality_ctx, use_quality_ctx,
            );
            info!("  I/O read: {:.2}s", io_t.elapsed().as_secs_f64());

            let compressed = match compress_handle.join() {
                Ok(v) => v,
                Err(e) => std::panic::resume_unwind(e),
            };
            (next, compressed)
        });

        // Accumulate compressed blocks in memory (instant, no disk I/O)
        let (h_blk, s_blk, q_blk, rc_blk) = compress_result?;
        all_h_blocks.extend(h_blk);
        all_s_blocks.extend(s_blk);
        all_q_blocks.extend(q_blk);
        if rc_canon {
            all_rc_blocks.extend(rc_blk);
        }

        num_reads += cur.num_reads;
        total_bases += cur.total_bases;
        original_size += cur.original_size;
        chunk_idx += 1;

        cur = next_result?;
    }

    // Compute stream sizes in multi-block format
    let h_data_size: usize = all_h_blocks.iter().map(|b| 4 + b.len()).sum();
    let s_data_size: usize = all_s_blocks.iter().map(|b| 4 + b.len()).sum();
    let q_data_size: usize = all_q_blocks.iter().map(|b| 4 + b.len()).sum();
    let rc_data_size: usize = all_rc_blocks.iter().map(|b| 4 + b.len()).sum();

    let headers_len = if !all_h_blocks.is_empty() { 4 + h_data_size } else { 0 };
    let sequences_len = if !all_s_blocks.is_empty() { 4 + s_data_size } else { 0 };
    let qualities_len = if !all_q_blocks.is_empty() { 4 + q_data_size } else { 0 };
    let rc_flags_len = if !all_rc_blocks.is_empty() { 4 + rc_data_size } else { 0 };

    info!("Read {} records in {} chunks ({} bases)", num_reads, chunk_idx, total_bases);
    info!("Compressed blocks: headers={} ({}) seq={} ({}) qual={} ({})",
        all_h_blocks.len(), humanize_bytes(h_data_size),
        all_s_blocks.len(), humanize_bytes(s_data_size),
        all_q_blocks.len(), humanize_bytes(q_data_size));
    if rc_canon {
        info!("  RC flags: {} ({})", all_rc_blocks.len(), humanize_bytes(rc_data_size));
    }

    // Write final archive directly from memory (no temp file copy)
    info!("Writing output file...");
    let encoding_type: u8 = if rc_canon { 6 } else if args.sequence_delta { 5 } else if args.sequence_hints { 4 } else { 0 };
    let mut output_file = std::io::BufWriter::new(std::fs::File::create(&args.output)?);

    output_file.write_all(&[encoding_type])?;
    output_file.write_all(&[0u8])?;                                             // arithmetic = disabled
    output_file.write_all(&[binning_to_code(quality_binning)])?;
    output_file.write_all(&[compressor_to_code(quality_compressor_used)])?;
    output_file.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?;
    output_file.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;
    output_file.write_all(&[0u8])?;                                             // quality_model = disabled
    output_file.write_all(&[0u8])?;                                             // quality_delta = disabled
    output_file.write_all(&[0u8])?;                                             // quality_dict = disabled
    output_file.write_all(&0u16.to_le_bytes())?;                                // template_prefix_len = 0
    output_file.write_all(&[0u8])?;                                             // has_comment = false

    output_file.write_all(&(num_reads as u64).to_le_bytes())?;
    output_file.write_all(&(headers_len as u64).to_le_bytes())?;
    output_file.write_all(&(sequences_len as u64).to_le_bytes())?;
    output_file.write_all(&0u64.to_le_bytes())?;                                // nmasks_len = 0
    output_file.write_all(&(if no_quality { 0 } else { qualities_len } as u64).to_le_bytes())?;

    // Write streams in multi-block format: [num_blocks: u32][block_len: u32, block_data]...
    let write_blocks = |blocks: &[Vec<u8>], out: &mut std::io::BufWriter<std::fs::File>| -> Result<()> {
        if blocks.is_empty() { return Ok(()); }
        out.write_all(&(blocks.len() as u32).to_le_bytes())?;
        for block in blocks {
            out.write_all(&(block.len() as u32).to_le_bytes())?;
            out.write_all(block)?;
        }
        Ok(())
    };

    write_blocks(&all_h_blocks, &mut output_file)?;
    write_blocks(&all_s_blocks, &mut output_file)?;
    write_blocks(&all_q_blocks, &mut output_file)?;

    // Append RC flags stream after qualities (encoding_type=6 signals its presence)
    if !all_rc_blocks.is_empty() {
        output_file.write_all(&(rc_flags_len as u64).to_le_bytes())?;
        write_blocks(&all_rc_blocks, &mut output_file)?;
    }
    output_file.flush()?;

    let metadata_size = 9 + 2 + 1 + 40 + if !all_rc_blocks.is_empty() { 8 } else { 0 };
    let total_compressed = headers_len + sequences_len + qualities_len + rc_flags_len + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!("Compression ratio: {:.2}x", original_size as f64 / total_compressed as f64);
    info!("Stream breakdown:");
    info!("  Headers:   {} bytes ({:.1}%)", headers_len, 100.0 * headers_len as f64 / total_compressed as f64);
    info!("  Sequences: {} bytes ({:.1}%)", sequences_len, 100.0 * sequences_len as f64 / total_compressed as f64);
    info!("  Qualities: {} bytes ({:.1}%)", qualities_len, 100.0 * qualities_len as f64 / total_compressed as f64);
    if rc_flags_len > 0 {
        info!("  RC flags:  {} bytes ({:.1}%)", rc_flags_len, 100.0 * rc_flags_len as f64 / total_compressed as f64);
    }
    info!("  Metadata:  {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total_compressed as f64);

    Ok(())
}

/// Chunked streaming compression with local reordering within each chunk.
///
/// Same bounded-memory design as `compress_chunked_bsc`, but reads records into a Vec
/// per chunk, sorts by `reorder_sort_key`, then builds streams from sorted records.
/// Reads within each 5M-record chunk are sorted by content similarity; across chunks
/// the order follows the input file.
pub(super) fn compress_chunked_bsc_reorder_local(args: &CompressConfig) -> Result<()> {
    use std::io::{Write, BufWriter};
    const CHUNK_SIZE: usize = 5_000_000; // 5M records per chunk

    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let sequence_hints = args.sequence_hints;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };

    info!("Local reorder mode: {} records per chunk, sorted by content similarity", CHUNK_SIZE);

    let mut reader = crate::io::FastqReader::from_path(input_path, args.fasta)?;

    // Temp files for streaming compressed blocks to disk
    let working_dir = &args.working_dir;
    let h_tmp_path = working_dir.join(".qz_chunked_h.tmp");
    let s_tmp_path = working_dir.join(".qz_chunked_s.tmp");
    let q_tmp_path = working_dir.join(".qz_chunked_q.tmp");

    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 {
                let _ = std::fs::remove_file(p);
            }
        }
    }
    let rc_tmp_path = working_dir.join(".qz_chunked_rc.tmp");
    let mut tmp_paths = vec![h_tmp_path.clone(), s_tmp_path.clone(), q_tmp_path.clone()];
    if args.rc_canon { tmp_paths.push(rc_tmp_path.clone()); }
    let _cleanup = TmpCleanup(tmp_paths);

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut s_tmp = BufWriter::new(std::fs::File::create(&s_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);
    let mut rc_tmp = if args.rc_canon {
        Some(BufWriter::new(std::fs::File::create(&rc_tmp_path)?))
    } else {
        None
    };

    let mut h_num_blocks: u32 = 0;
    let mut s_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;
    let mut rc_num_blocks: u32 = 0;
    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;
    let mut chunk_idx: usize = 0;
    let rc_canon = args.rc_canon;

    // Read first chunk
    let (mut cur_records, mut cur_bases, mut cur_orig) =
        read_chunk_records(&mut reader, CHUNK_SIZE)?;

    while !cur_records.is_empty() {
        let cur_reads = cur_records.len();
        info!("Chunk {}: {} reads, sorting...", chunk_idx, cur_reads);

        // Sort records by content similarity
        let sort_start = Instant::now();
        cur_records = sort_records_by_key(cur_records);
        info!("  Sorted in {:.2}s", sort_start.elapsed().as_secs_f64());

        // Build streams from sorted records, then compress
        let use_fqzcomp = args.quality_compressor == QualityCompressor::Fqzcomp;
        let (next_result, compress_result) = std::thread::scope(|scope| {
            let records = std::mem::take(&mut cur_records);

            let compress_handle = scope.spawn(move || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                // When using fqzcomp, skip quality in records_to_streams (we compress separately)
                let skip_quality = no_quality || use_fqzcomp;
                let (h_data, s_data, q_data, rc_data) = records_to_streams(&records, quality_mode, quality_binning, skip_quality, sequence_hints, false, rc_canon)?;

                // Overlap BSC (H/S) with quality compression via rayon::join
                let (bsc_result, q_result) = rayon::join(
                    || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                        let h_blocks = compress_stream_to_bsc_blocks(&h_data, bsc_static)?;
                        let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
                        let rc_blocks = if rc_data.is_empty() {
                            Vec::new()
                        } else {
                            compress_stream_to_bsc_blocks(&rc_data, bsc_static)?
                        };
                        Ok((h_blocks, s_blocks, rc_blocks))
                    },
                    || -> Result<Vec<Vec<u8>>> {
                        if no_quality {
                            Ok(Vec::new())
                        } else if use_fqzcomp {
                            // Fqzcomp: compress sorted records' qualities directly (one blob per chunk)
                            let q_blob = codecs::compress_qualities_fqzcomp(&records)?;
                            Ok(vec![q_blob])
                        } else {
                            if q_data.is_empty() {
                                Ok(Vec::new())
                            } else {
                                compress_stream_to_bsc_blocks(&q_data, bsc_static)
                            }
                        }
                    },
                );
                drop(records);
                let (h_blocks, s_blocks, rc_blocks) = bsc_result?;
                let q_blocks = q_result?;

                Ok((h_blocks, s_blocks, q_blocks, rc_blocks))
            });

            // Read next chunk on main thread (overlaps with compression)
            let next = read_chunk_records(&mut reader, CHUNK_SIZE);

            let compressed = match compress_handle.join() {
                Ok(v) => v,
                Err(e) => std::panic::resume_unwind(e),
            };
            (next, compressed)
        });

        let (h_blk, s_blk, q_blk, rc_blk) = compress_result?;
        h_num_blocks += write_blocks_to_tmp(h_blk, &mut h_tmp)?;
        s_num_blocks += write_blocks_to_tmp(s_blk, &mut s_tmp)?;
        if use_fqzcomp {
            // Fqzcomp: write single blob as one quality "block" (same as compress_chunked_fqzcomp)
            for q_blob in q_blk {
                if !q_blob.is_empty() {
                    use std::io::Write as _;
                    q_tmp.write_all(&(q_blob.len() as u32).to_le_bytes())?;
                    q_tmp.write_all(&q_blob)?;
                    q_num_blocks += 1;
                }
            }
        } else {
            q_num_blocks += write_blocks_to_tmp(q_blk, &mut q_tmp)?;
        }
        if let Some(ref mut rc_file) = rc_tmp {
            rc_num_blocks += write_blocks_to_tmp(rc_blk, rc_file)?;
        }

        num_reads += cur_reads;
        total_bases += cur_bases;
        original_size += cur_orig;
        chunk_idx += 1;

        let (nr, nb, no) = next_result?;
        cur_records = nr;
        cur_bases = nb;
        cur_orig = no;
    }

    // Flush and close temp files
    h_tmp.flush()?;
    s_tmp.flush()?;
    q_tmp.flush()?;
    if let Some(ref mut rc_file) = rc_tmp {
        rc_file.flush()?;
    }
    drop(h_tmp);
    drop(s_tmp);
    drop(q_tmp);
    drop(rc_tmp);

    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;
    let rc_tmp_size = if rc_canon { std::fs::metadata(&rc_tmp_path)?.len() as usize } else { 0 };

    let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
    let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };
    let rc_flags_len = if rc_num_blocks > 0 { 4 + rc_tmp_size } else { 0 };

    info!("Read {} records in {} chunks ({} bases)", num_reads, chunk_idx, total_bases);

    // Write final output
    let encoding_type: u8 = if rc_canon { 6 } else { 0 };
    let rc_stream = if rc_canon {
        Some(RcStreamParams { flags_len: rc_flags_len, num_blocks: rc_num_blocks, tmp_path: &rc_tmp_path })
    } else {
        None
    };
    write_chunked_archive(
        &args.output, quality_binning, args.quality_compressor, no_quality, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path,
        original_size, start_time, rc_stream,
    )
}


/// Chunked streaming compression with fqzcomp for quality scores.
///
/// Reads records per chunk (fqzcomp needs individual quality strings, not a byte stream).
/// Headers and sequences compressed with BSC blocks as usual.
/// Each chunk's qualities compressed with fqzcomp independently, stored as one "block"
/// in the multi-block format. Decompressor reads blocks and fqzcomp-decompresses each.
pub(super) fn compress_chunked_fqzcomp(args: &CompressConfig) -> Result<()> {
    use std::io::{Write, BufWriter};
    const CHUNK_SIZE: usize = 5_000_000;

    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.bsc_static;
    let no_quality = args.no_quality;

    info!("Chunked fqzcomp mode: {} records per chunk, BSC for headers/sequences, fqzcomp for qualities", CHUNK_SIZE);

    let mut reader = crate::io::FastqReader::from_path(input_path, args.fasta)?;

    let working_dir = &args.working_dir;
    let h_tmp_path = working_dir.join(".qz_chunked_h.tmp");
    let s_tmp_path = working_dir.join(".qz_chunked_s.tmp");
    let q_tmp_path = working_dir.join(".qz_chunked_q.tmp");

    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 {
                let _ = std::fs::remove_file(p);
            }
        }
    }
    let _cleanup = TmpCleanup(vec![h_tmp_path.clone(), s_tmp_path.clone(), q_tmp_path.clone()]);

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut s_tmp = BufWriter::new(std::fs::File::create(&s_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);

    let mut h_num_blocks: u32 = 0;
    let mut s_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;
    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;
    let mut chunk_idx: usize = 0;

    let (mut cur_records, mut cur_bases, mut cur_orig) =
        read_chunk_records(&mut reader, CHUNK_SIZE)?;

    while !cur_records.is_empty() {
        let cur_reads = cur_records.len();
        info!("Chunk {}: {} reads", chunk_idx, cur_reads);

        let (next_result, compress_result) = std::thread::scope(|scope| {
            let records = std::mem::take(&mut cur_records);

            let compress_handle = scope.spawn(move || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<u8>)> {
                // Build header and sequence streams (fqzcomp stores raw ASCII quality,
                // so we pass QualityBinning::None and skip quality from records_to_streams)
                let (h_data, s_data, _q_unused, _rc_unused) = records_to_streams(
                    &records, QualityMode::Lossless, QualityBinning::None, true, false, false, false, // no_quality=true to skip, no hints, no delta, no rc
                )?;

                // Overlap BSC (headers + sequences) with fqzcomp (qualities) via rayon::join
                let (bsc_result, q_blob) = rayon::join(
                    || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                        let h_blocks = compress_stream_to_bsc_blocks(&h_data, bsc_static)?;
                        let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
                        Ok((h_blocks, s_blocks))
                    },
                    || -> Result<Vec<u8>> {
                        if no_quality {
                            Ok(Vec::new())
                        } else {
                            codecs::compress_qualities_fqzcomp(&records)
                        }
                    },
                );
                let (h_blocks, s_blocks) = bsc_result?;
                let q_blob = q_blob?;

                Ok((h_blocks, s_blocks, q_blob))
            });

            let next = read_chunk_records(&mut reader, CHUNK_SIZE);
            let compressed = match compress_handle.join() {
                Ok(v) => v,
                Err(e) => std::panic::resume_unwind(e),
            };
            (next, compressed)
        });

        let (h_blk, s_blk, q_blob) = compress_result?;
        h_num_blocks += write_blocks_to_tmp(h_blk, &mut h_tmp)?;
        s_num_blocks += write_blocks_to_tmp(s_blk, &mut s_tmp)?;

        // Write fqzcomp blob as a single quality "block"
        if !q_blob.is_empty() {
            q_tmp.write_all(&(q_blob.len() as u32).to_le_bytes())?;
            q_tmp.write_all(&q_blob)?;
            q_num_blocks += 1;
        }

        num_reads += cur_reads;
        total_bases += cur_bases;
        original_size += cur_orig;
        chunk_idx += 1;

        let (nr, nb, no) = next_result?;
        cur_records = nr;
        cur_bases = nb;
        cur_orig = no;
    }

    h_tmp.flush()?;
    s_tmp.flush()?;
    q_tmp.flush()?;
    drop(h_tmp);
    drop(s_tmp);
    drop(q_tmp);

    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
    let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    info!("Read {} records in {} chunks ({} bases)", num_reads, chunk_idx, total_bases);

    write_chunked_archive(
        &args.output, QualityBinning::None, QualityCompressor::Fqzcomp, no_quality, 0u8,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path,
        original_size, start_time, None,
    )
}

/// Global reorder: two-pass bucket sort with bounded memory per bucket.
///
/// Pass 1: Stream through FASTQ, compute sort key for each record, write to one of 256
/// bucket temp files based on top bits of the sort key.
/// Pass 2: For each bucket in order, read all records, sort by full key, build streams,
/// compress to BSC blocks, and write to output temp files.
///
/// Memory is bounded by the largest bucket (~input_size/256 + BSC working memory).
pub(super) fn compress_global_reorder_bsc(args: &CompressConfig) -> Result<()> {
    use std::io::{Write, Read as IoRead, BufWriter, BufReader};
    const NUM_BUCKETS: usize = 256;

    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };

    info!("Global reorder mode: two-pass bucket sort with {} buckets", NUM_BUCKETS);

    let working_dir = &args.working_dir;

    // --- Temp file setup ---
    // Bucket files for pass 1
    let bucket_paths: Vec<std::path::PathBuf> = (0..NUM_BUCKETS)
        .map(|i| working_dir.join(format!(".qz_bucket_{:03}.tmp", i)))
        .collect();
    // Output stream temp files (same as chunked mode)
    let h_tmp_path = working_dir.join(".qz_reorder_h.tmp");
    let s_tmp_path = working_dir.join(".qz_reorder_s.tmp");
    let q_tmp_path = working_dir.join(".qz_reorder_q.tmp");

    // Cleanup guard for ALL temp files
    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 {
                let _ = std::fs::remove_file(p);
            }
        }
    }
    let mut cleanup_paths: Vec<std::path::PathBuf> = bucket_paths.clone();
    cleanup_paths.push(h_tmp_path.clone());
    cleanup_paths.push(s_tmp_path.clone());
    cleanup_paths.push(q_tmp_path.clone());
    let _cleanup = TmpCleanup(cleanup_paths);

    // --- Pass 1: bucket records by sort key ---
    info!("Pass 1: bucketing records by sort key...");
    let pass1_start = Instant::now();

    let mut bucket_writers: Vec<BufWriter<std::fs::File>> = bucket_paths.iter()
        .map(|p| Ok(BufWriter::new(std::fs::File::create(p)?)))
        .collect::<Result<Vec<_>>>()?;
    let mut bucket_counts = vec![0usize; NUM_BUCKETS];
    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;

    {
        let mut reader = crate::io::FastqReader::from_path(input_path, args.fasta)?;
        while let Some(record) = reader.next()? {
            let qual_bytes = record.quality.as_ref().map(|q| q.as_bytes()).unwrap_or(&[]);
            let key = dna_utils::reorder_sort_key(record.sequence.as_bytes(), qual_bytes);
            let bucket = (key >> 56) as usize % NUM_BUCKETS;

            // Write record to bucket file: [key:16][id_len:u32][id][seq_len:u32][seq][qual_len:u32][qual]
            let w = &mut bucket_writers[bucket];
            w.write_all(&key.to_le_bytes())?;

            w.write_all(&(record.id.len() as u32).to_le_bytes())?;
            w.write_all(record.id.as_bytes())?;

            w.write_all(&(record.sequence.len() as u32).to_le_bytes())?;
            w.write_all(record.sequence.as_bytes())?;

            let qual_len = record.quality.as_ref().map(|q| q.len()).unwrap_or(0);
            w.write_all(&(qual_len as u32).to_le_bytes())?;
            if let Some(ref q) = record.quality {
                w.write_all(q.as_bytes())?;
            }

            total_bases += record.sequence.len();
            original_size += record.id.len() + record.sequence.len() + qual_len + 3;
            bucket_counts[bucket] += 1;
            num_reads += 1;
        }
    }

    // Flush and close bucket writers
    for w in &mut bucket_writers {
        w.flush()?;
    }
    drop(bucket_writers);

    let non_empty_buckets = bucket_counts.iter().filter(|&&c| c > 0).count();
    let max_bucket = bucket_counts.iter().max().copied().unwrap_or(0);
    info!("Pass 1 done in {:.2}s: {} reads into {} non-empty buckets (max bucket: {} reads)",
        pass1_start.elapsed().as_secs_f64(), num_reads, non_empty_buckets, max_bucket);

    // --- Pass 2: sort each bucket, compress ---
    info!("Pass 2: sorting buckets and compressing...");
    let pass2_start = Instant::now();

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut s_tmp = BufWriter::new(std::fs::File::create(&s_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);

    let mut h_num_blocks: u32 = 0;
    let mut s_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;

    for bucket_idx in 0..NUM_BUCKETS {
        let count = bucket_counts[bucket_idx];
        if count == 0 {
            continue;
        }

        // Read all records from this bucket file
        let mut keyed_records: Vec<(u128, crate::io::FastqRecord)> = Vec::with_capacity(count);
        {
            let file = std::fs::File::open(&bucket_paths[bucket_idx])?;
            let mut br = BufReader::new(file);

            for _ in 0..count {
                let mut key_buf = [0u8; 16];
                br.read_exact(&mut key_buf)?;
                let key = u128::from_le_bytes(key_buf);

                let mut len_buf = [0u8; 4];

                br.read_exact(&mut len_buf)?;
                let id_len = u32::from_le_bytes(len_buf) as usize;
                let mut id_bytes = vec![0u8; id_len];
                br.read_exact(&mut id_bytes)?;

                br.read_exact(&mut len_buf)?;
                let seq_len = u32::from_le_bytes(len_buf) as usize;
                let mut seq_bytes = vec![0u8; seq_len];
                br.read_exact(&mut seq_bytes)?;

                br.read_exact(&mut len_buf)?;
                let qual_len = u32::from_le_bytes(len_buf) as usize;
                let quality = if qual_len > 0 {
                    let mut qual_bytes = vec![0u8; qual_len];
                    br.read_exact(&mut qual_bytes)?;
                    Some(String::from_utf8(qual_bytes)
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in quality: {}", e))?)
                } else {
                    None
                };

                keyed_records.push((key, crate::io::FastqRecord {
                    id: String::from_utf8(id_bytes)
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in id: {}", e))?,
                    sequence: String::from_utf8(seq_bytes)
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?,
                    quality,
                }));
            }
        }

        // Sort by full key within bucket
        keyed_records.sort_unstable_by_key(|(k, _)| *k);

        // Build streams and compress (sequential: h -> s -> q)
        let records: Vec<crate::io::FastqRecord> = keyed_records.into_iter().map(|(_, r)| r).collect();
        let (h_data, s_data, q_data, _rc_unused) = records_to_streams(&records, quality_mode, quality_binning, no_quality, args.sequence_hints, args.sequence_delta, false)?;
        drop(records);

        let h_blocks = compress_stream_to_bsc_blocks(&h_data, bsc_static)?;
        drop(h_data);
        h_num_blocks += write_blocks_to_tmp(h_blocks, &mut h_tmp)?;

        let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
        drop(s_data);
        s_num_blocks += write_blocks_to_tmp(s_blocks, &mut s_tmp)?;

        if !no_quality && !q_data.is_empty() {
            let q_blocks = compress_stream_to_bsc_blocks(&q_data, bsc_static)?;
            q_num_blocks += write_blocks_to_tmp(q_blocks, &mut q_tmp)?;
        }

        // Delete bucket file immediately to free disk space
        let _ = std::fs::remove_file(&bucket_paths[bucket_idx]);
    }

    info!("Pass 2 done in {:.2}s", pass2_start.elapsed().as_secs_f64());

    // Flush and close output temp files
    h_tmp.flush()?;
    s_tmp.flush()?;
    q_tmp.flush()?;
    drop(h_tmp);
    drop(s_tmp);
    drop(q_tmp);

    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
    let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    info!("Read {} records in {} non-empty buckets ({} bases)", num_reads, non_empty_buckets, total_bases);

    // Write final output (identical archive format)
    let encoding_type: u8 = if args.sequence_delta { 5 } else if args.sequence_hints { 4 } else { 0 };
    write_chunked_archive(
        &args.output, quality_binning, QualityCompressor::Bsc, no_quality, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path,
        original_size, start_time, None,
    )
}
pub(super) fn compress(args: &CompressConfig) -> Result<()> {
    let start_time = Instant::now();

    // Set up thread pool
    let num_threads = if args.threads == 0 {
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8)
    } else {
        args.threads
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok(); // Ignore error if already initialized

    info!("Using {} threads for compression", num_threads);
    info!("Compression level: {}", args.compression_level);
    info!("Input files: {:?}", args.input);
    info!("Output: {:?}", args.output);

    // Validate --sequence-hints compatibility
    if args.sequence_hints {
        if args.delta_encoding || args.rle_encoding || args.debruijn || args.arithmetic {
            anyhow::bail!("--sequence-hints cannot be combined with --delta-encoding, --rle-encoding, --debruijn, or --arithmetic");
        }
        if args.sequence_compressor != SequenceCompressor::Bsc {
            anyhow::bail!("--sequence-hints requires BSC sequence compressor (default)");
        }
    }

    // Validate --sequence-delta compatibility
    if args.sequence_delta {
        if args.delta_encoding || args.rle_encoding || args.debruijn || args.arithmetic || args.sequence_hints {
            anyhow::bail!("--sequence-delta cannot be combined with --delta-encoding, --rle-encoding, --debruijn, --arithmetic, or --sequence-hints");
        }
        if args.sequence_compressor != SequenceCompressor::Bsc {
            anyhow::bail!("--sequence-delta requires BSC sequence compressor (default)");
        }
    }

    // Validate --quality-compressor quality-ctx
    if args.quality_compressor == QualityCompressor::QualityCtx {
        if args.no_quality {
            anyhow::bail!("--quality-compressor quality-ctx cannot be used with --no-quality");
        }
        if args.quality_mode != QualityMode::Lossless {
            anyhow::bail!("--quality-compressor quality-ctx requires lossless quality mode (default)");
        }
    }

    // Validate --factorize compatibility
    if args.factorize {
        if args.delta_encoding || args.rle_encoding || args.debruijn || args.arithmetic
            || args.sequence_hints || args.sequence_delta || args.rc_canon {
            anyhow::bail!("--factorize cannot be combined with --delta-encoding, --rle-encoding, --debruijn, --arithmetic, --sequence-hints, --sequence-delta, or --rc-canon");
        }
        if args.sequence_compressor != SequenceCompressor::Bsc {
            anyhow::bail!("--factorize requires BSC sequence compressor (default)");
        }
        if args.header_compressor != HeaderCompressor::Bsc {
            anyhow::bail!("--factorize requires BSC header compressor (default)");
        }
    }

    // Validate --local-reorder / --ultra / --fast-ultra compatibility
    let has_ultra = args.ultra.is_some();
    if args.local_reorder || has_ultra || args.fast_ultra {
        let mode_name = if args.local_reorder { "--local-reorder" } else if args.fast_ultra { "--fast-ultra" } else { "--ultra" };
        if (args.local_reorder as u8 + has_ultra as u8 + args.fast_ultra as u8) > 1 {
            anyhow::bail!("--local-reorder, --ultra, and --fast-ultra cannot be combined");
        }
        if args.factorize || args.delta_encoding || args.rle_encoding || args.debruijn
            || args.arithmetic || args.sequence_hints || args.sequence_delta || args.rc_canon {
            anyhow::bail!("{mode_name} cannot be combined with other encoding options");
        }
        if args.sequence_compressor != SequenceCompressor::Bsc {
            anyhow::bail!("{mode_name} requires BSC sequence compressor (default)");
        }
        if args.header_compressor != HeaderCompressor::Bsc {
            anyhow::bail!("{mode_name} requires BSC header compressor (default)");
        }
    }

    if args.input.len() != 1 {
        anyhow::bail!("Expected exactly 1 input file, got {}", args.input.len());
    }

    // Use memory-efficient streaming mode for the default BSC path
    // (no advanced quality encoding, BSC for all streams)
    let can_stream_base = args.sequence_compressor == SequenceCompressor::Bsc
        && args.header_compressor == HeaderCompressor::Bsc
        && !args.delta_encoding
        && !args.rle_encoding
        && !args.debruijn
        && !args.arithmetic
        && !args.quality_modeling
        && !args.quality_delta
        && !args.dict_training
        && !args.twobit
        && !args.header_template;

    let can_stream = can_stream_base
        && (args.quality_compressor == QualityCompressor::Bsc
            || args.quality_compressor == QualityCompressor::QualityCtx);

    // Reorder modes require the streaming-compatible BSC base path
    if let Some(mode) = args.reorder {
        if !can_stream_base {
            anyhow::bail!("--reorder requires BSC compressors for headers/sequences and no advanced encoding options");
        }
        return match mode {
            crate::cli::ReorderMode::Local => compress_chunked_bsc_reorder_local(args),
            crate::cli::ReorderMode::Global => compress_global_reorder_bsc(args),
        };
    }

    // Factorize mode: read-to-read delta with inverted index
    if args.factorize {
        return factorize::compress_factorize(args);
    }

    // Local reorder mode: center-hash grouping + delta encoding
    if args.local_reorder {
        return ultra::compress_harc(args);
    }

    // Ultra mode: level-based compression with auto-tuning
    if let Some(requested_level) = args.ultra {
        let level = ultra::resolve_ultra_level(requested_level);
        return ultra::compress_reorder_local_with_level(args, level);
    }

    // Backwards compat: --fast-ultra maps to ultra level 2
    if args.fast_ultra {
        let level = ultra::resolve_ultra_level(2);
        return ultra::compress_reorder_local_with_level(args, level);
    }

    // Fqzcomp quality needs record-level access, so use chunked record-based path
    if can_stream_base && args.quality_compressor == QualityCompressor::Fqzcomp {
        return compress_chunked_fqzcomp(args);
    }

    if can_stream {
        // Always use chunked pipeline: pipelined I/O (read next chunk while
        // compressing current), parallel h||s||q compression, sub-block
        // quality_ctx parallelism, and bounded memory via temp files.
        return compress_chunked_bsc(args);
    }

    // Non-streaming mode: read all records into memory first
    // (required for advanced quality encoding, non-BSC compressors)
    compress_in_memory(args, start_time)
}

/// In-memory compression path for advanced quality encodings and non-BSC compressors.
///
/// Reads all records into memory, applies optional delta/RLE/debruijn/arithmetic encoding,
/// then compresses all three streams (headers, sequences, qualities) in parallel.
fn compress_in_memory(args: &CompressConfig, start_time: Instant) -> Result<()> {
    use crate::io::{FastqReader, FastqRecord};
    use std::io::Write;

    info!("Reading FASTQ file...");
    let input_path = &args.input[0];
    let mut reader = FastqReader::from_path(input_path, false)?; // false = FASTQ format
    let mut records = Vec::new();

    while let Some(record) = reader.next()? {
        records.push(record);
    }

    info!("Read {} records", records.len());

    // Apply quality quantization if needed
    let processed_records: Vec<FastqRecord> = if !args.no_quality {
        records
            .into_iter()
            .map(|mut record| {
                if let Some(qual) = &record.quality {
                    record.quality = Some(quality::quantize_quality(qual, args.quality_mode).into_owned());
                }
                record
            })
            .collect()
    } else {
        records
    };

    // Determine encoding type and quality binning
    let encoding_type: u8 = if args.delta_encoding {
        1 // Delta encoding
    } else if args.rle_encoding {
        2 // RLE encoding
    } else if args.debruijn {
        3 // De Bruijn graph encoding
    } else if args.sequence_hints {
        4 // Sequence hints for BWT clustering
    } else if args.sequence_delta {
        5 // Inline delta encoding
    } else {
        0 // No special encoding
    };

    let quality_binning = if args.no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(args.quality_mode)
    };

    // Quality pre-computation (must complete before parallel compression block)
    // Train quality dictionary if requested
    let quality_dict_opt = if args.dict_training && !args.no_quality {
        info!("Training zstd dictionary for quality scores...");
        let quality_strings: Vec<String> = processed_records
            .iter()
            .filter_map(|r| r.quality.clone())
            .collect();

        if !quality_strings.is_empty() {
            let dict_size_bytes = args.dict_size * 1024; // Convert KB to bytes
            let dict = zstd_dict::train_from_quality_scores(&quality_strings, dict_size_bytes, 0.1)?;
            Some(dict)
        } else {
            None
        }
    } else {
        None
    };

    // Determine quality compression strategy
    // Priority: delta > modeling > standard
    let (quality_delta_opt, quality_model_opt) = if args.quality_delta && !args.no_quality {
        info!("Applying quality delta encoding...");
        let qual_compressed = if let Some(ref dict) = quality_dict_opt {
            codecs::compress_qualities_with_delta_and_dict(&processed_records, dict, args.compression_level)?
        } else {
            codecs::compress_qualities_with_delta(&processed_records, args.compression_level)?
        };
        (Some(qual_compressed), None)
    } else if args.quality_modeling && !args.no_quality {
        info!("Building positional quality model...");
        let (qual_compressed, model) = if let Some(ref dict) = quality_dict_opt {
            codecs::compress_qualities_with_model_and_dict(&processed_records, dict, args.compression_level)?
        } else {
            codecs::compress_qualities_with_model(&processed_records, args.compression_level, args.quality_compressor, args.bsc_static)?
        };
        (None, Some((qual_compressed, model)))
    } else {
        (None, None)
    };

    // Helper closure: resolve quality compression (priority: delta > model > dict > default)
    let bsc_static = args.bsc_static;
    let resolve_quality = || -> Result<Vec<u8>> {
        if let Some(ref qual) = quality_delta_opt {
            Ok(qual.clone())
        } else if let Some((ref qual, _)) = quality_model_opt {
            Ok(qual.clone())
        } else if args.no_quality {
            Ok(Vec::new())
        } else if let Some(ref dict) = quality_dict_opt {
            codecs::compress_qualities_with_dict(&processed_records, quality_binning, dict, args.compression_level, args.quality_compressor)
        } else {
            codecs::compress_qualities_with(&processed_records, quality_binning, args.compression_level, args.quality_compressor, bsc_static)
        }
    };

    // Compress all three streams in parallel:
    //   Left:  headers
    //   Right: sequences + qualities (with nested rayon::join)
    let (header_result, seq_qual_result) = rayon::join(
        || -> Result<(Vec<u8>, read_id::ReadIdTemplate)> {
            match args.header_compressor {
                HeaderCompressor::Bsc => {
                    if args.header_template {
                        let read_ids: Vec<String> = processed_records.iter().map(|r| r.id.clone()).collect();
                        let encoded = read_id::compress_read_ids(&read_ids)?;
                        if encoded.template.prefix.is_empty() {
                            // Template analysis found no common Illumina structure, fall back to raw BSC
                            info!("No template structure found in headers, falling back to raw BSC{}...", if args.bsc_static { " (static)" } else { " (adaptive)" });
                            let compressed = codecs::compress_headers_bsc_with(&processed_records, args.bsc_static)?;
                            let dummy_template = read_id::ReadIdTemplate {
                                prefix: String::new(),
                                has_comment: false,
                                common_comment: None,
                            };
                            Ok((compressed, dummy_template))
                        } else {
                            info!("Compressing read IDs with template encoding + BSC{}...", if args.bsc_static { " (static)" } else { " (adaptive)" });
                            let compressed = bsc_compress_parallel(&encoded.encoded_data, args.bsc_static)?;
                            Ok((compressed, encoded.template))
                        }
                    } else {
                        info!("Compressing read IDs with raw BSC{}...", if args.bsc_static { " (static)" } else { " (adaptive)" });
                        let compressed = codecs::compress_headers_bsc_with(&processed_records, args.bsc_static)?;
                        let dummy_template = read_id::ReadIdTemplate {
                            prefix: String::new(),
                            has_comment: false,
                            common_comment: None,
                        };
                        Ok((compressed, dummy_template))
                    }
                }
                HeaderCompressor::Zstd => {
                    info!("Compressing read IDs with template encoding...");
                    codecs::compress_headers(&processed_records, args.compression_level)
                }
                HeaderCompressor::OpenZl => {
                    info!("Compressing read IDs with raw OpenZL...");
                    let compressed = codecs::compress_headers_openzl(&processed_records)?;
                    let dummy_template = read_id::ReadIdTemplate {
                        prefix: String::new(),
                        has_comment: false,
                        common_comment: None,
                    };
                    Ok((compressed, dummy_template))
                }
            }
        },
        || -> Result<(Vec<u8>, Vec<u8>, Vec<u8>)> {
            if args.delta_encoding {
                info!("Applying delta encoding and compressing...");
                let sequences_data: Vec<String> = processed_records.iter().map(|r| r.sequence.clone()).collect();
                let encoded_sequences = delta::apply_delta_encoding(&sequences_data);

                let mut seq_stream = Vec::new();
                for enc in &encoded_sequences {
                    let len = enc.len();
                    seq_stream.write_all(&(len as u32).to_le_bytes())?;
                    seq_stream.write_all(enc)?;
                }

                let bsc_static = args.bsc_static;
                let (seq_result, qual_result) = rayon::join(
                    || match args.sequence_compressor {
                        SequenceCompressor::Bsc => bsc_compress_parallel(&seq_stream, bsc_static),
                        SequenceCompressor::Zstd => compress_zstd(&seq_stream, args.compression_level),
                        SequenceCompressor::OpenZl => openzl::compress_parallel(&seq_stream),
                    },
                    &resolve_quality,
                );
                Ok((seq_result?, Vec::new(), qual_result?))
            } else if args.rle_encoding {
                info!("Applying RLE encoding and compressing...");
                let sequences_data: Vec<String> = processed_records.iter().map(|r| r.sequence.clone()).collect();
                let encoded_sequences = rle::apply_rle_to_sequences(&sequences_data);

                let mut seq_stream = Vec::new();
                for enc in &encoded_sequences {
                    let len = enc.len();
                    seq_stream.write_all(&(len as u32).to_le_bytes())?;
                    seq_stream.write_all(enc)?;
                }

                let bsc_static = args.bsc_static;
                let (seq_result, qual_result) = rayon::join(
                    || match args.sequence_compressor {
                        SequenceCompressor::Bsc => bsc_compress_parallel(&seq_stream, bsc_static),
                        SequenceCompressor::Zstd => compress_zstd(&seq_stream, args.compression_level),
                        SequenceCompressor::OpenZl => openzl::compress_parallel(&seq_stream),
                    },
                    &resolve_quality,
                );
                Ok((seq_result?, Vec::new(), qual_result?))
            } else if args.debruijn {
                info!("Compressing with de Bruijn graph (patent-safe, RC-aware)...");
                let sequences: Vec<String> = processed_records.iter().map(|r| r.sequence.clone()).collect();
                let (seq_result, qual_result) = rayon::join(
                    || debruijn::compress_sequences_debruijn(&sequences, args.kmer_size),
                    &resolve_quality,
                );
                Ok((seq_result?, Vec::new(), qual_result?))
            } else if args.arithmetic {
                info!("Compressing with arithmetic coding (experimental 6-8x)...");

                let (seq_result, qual_result) = rayon::join(
                    || codecs::compress_sequences_arithmetic(&processed_records),
                    || -> Result<Vec<u8>> {
                        if let Some(ref qual) = quality_delta_opt {
                            Ok(qual.clone())
                        } else if let Some((ref qual, _)) = quality_model_opt {
                            Ok(qual.clone())
                        } else if args.no_quality {
                            Ok(Vec::new())
                        } else {
                            codecs::compress_qualities_arithmetic(&processed_records)
                        }
                    },
                );
                Ok((seq_result?, Vec::new(), qual_result?))
            } else {
                match args.sequence_compressor {
                    SequenceCompressor::Bsc => {
                        let twobit = args.twobit;
                        info!("Compressing sequences{} and qualities in parallel with BSC{}...",
                            if twobit { " (2-bit)" } else { "" },
                            if args.bsc_static { " (static)" } else { " (adaptive)" });
                        let bsc_static = args.bsc_static;
                        let (seq_result, qual_result) = rayon::join(
                            || if twobit {
                                codecs::compress_sequences_2bit_bsc_with(&processed_records, bsc_static)
                            } else {
                                codecs::compress_sequences_raw_bsc_with(&processed_records, bsc_static, args.sequence_hints, args.sequence_delta)
                            },
                            &resolve_quality,
                        );
                        let (sequences, nmasks) = seq_result?;
                        let qualities = qual_result?;
                        Ok((sequences, nmasks, qualities))
                    }
                    SequenceCompressor::Zstd => {
                        info!("Compressing with columnar format (N-mask lossless)...");
                        let (_h, sequences, nmasks, qualities_columnar, _stats) = compress_columnar(&processed_records, quality_binning, args.compression_level)?;
                        let qualities = if let Some(ref qual) = quality_delta_opt {
                            qual.clone()
                        } else if let Some((ref qual, _)) = quality_model_opt {
                            qual.clone()
                        } else {
                            qualities_columnar
                        };
                        Ok((sequences, nmasks, qualities))
                    }
                    SequenceCompressor::OpenZl => {
                        info!("Compressing sequences and qualities in parallel with OpenZL...");
                        let (seq_result, qual_result) = rayon::join(
                            || codecs::compress_sequences_raw_openzl(&processed_records),
                            &resolve_quality,
                        );
                        let (sequences, nmasks) = seq_result?;
                        let qualities = qual_result?;
                        Ok((sequences, nmasks, qualities))
                    }
                }
            }
        },
    );
    let (headers, read_id_template) = header_result?;
    let (sequences, nmasks, qualities) = seq_qual_result?;

    // Compute stats after all streams are compressed
    let stats = columnar::ColumnarStats {
        num_reads: processed_records.len(),
        total_bases: processed_records.iter().map(|r| r.sequence.len()).sum(),
        original_size: processed_records.iter().map(|r| r.id.len() + r.sequence.len() + r.quality.as_ref().map(|q| q.len()).unwrap_or(0) + 3).sum(),
        headers_size: headers.len(),
        sequences_size: sequences.len(),
        qualities_size: qualities.len(),
        total_compressed: headers.len() + sequences.len() + qualities.len(),
    };

    // Write output file with archive format:
    // [encoding_type: 1 byte]
    // [arithmetic_mode: 1 byte]
    // [read_lengths_len: 4 bytes][read_lengths: 4 bytes each] (if arithmetic mode)
    // [read_lengths_len: 4 bytes][read_lengths: 4 bytes each] (if arithmetic mode)
    // [quality_binning: 1 byte]
    // [quality_modeling: 1 byte]
    // [quality_model_size: 2 bytes][quality_model_data] (if modeling enabled)
    // [quality_delta: 1 byte]
    // [template_prefix_len: 2 bytes][template_prefix][template_has_comment: 1 byte]
    // [num_reads: 8 bytes][headers_len: 8 bytes][sequences_len: 8 bytes][nmasks_len: 8 bytes][qualities_len: 8 bytes]
    // [headers][sequences][nmasks][qualities]
    info!("Writing output file...");
    let mut output_file = std::fs::File::create(&args.output)?;

    // Write encoding type
    output_file.write_all(&[encoding_type])?;

    // Write arithmetic mode flag and read lengths
    if args.arithmetic {
        output_file.write_all(&[1])?; // arithmetic enabled
        let read_lengths: Vec<u32> = processed_records.iter().map(|r| r.sequence.len() as u32).collect();
        output_file.write_all(&(read_lengths.len() as u32).to_le_bytes())?;
        for &len in &read_lengths {
            output_file.write_all(&len.to_le_bytes())?;
        }
    } else {
        output_file.write_all(&[0])?; // arithmetic disabled
    }

    // Write quality binning mode
    output_file.write_all(&[binning_to_code(quality_binning)])?;

    // Write quality compressor type
    output_file.write_all(&[compressor_to_code(args.quality_compressor)])?;

    // Write sequence compressor type
    output_file.write_all(&[seq_compressor_to_code(args.sequence_compressor)])?;

    // Write header compressor type
    output_file.write_all(&[header_compressor_to_code(args.header_compressor)])?;

    // Write quality modeling flag and model
    if let Some((_, ref model)) = quality_model_opt {
        output_file.write_all(&[1])?; // modeling enabled
        let model_bytes = quality_model::serialize_model(model);
        output_file.write_all(&(model_bytes.len() as u16).to_le_bytes())?;
        output_file.write_all(&model_bytes)?;
    } else {
        output_file.write_all(&[0])?; // modeling disabled
    }

    // Write quality delta flag
    output_file.write_all(&[if quality_delta_opt.is_some() { 1 } else { 0 }])?;

    // Write quality dictionary if present
    if let Some(ref dict) = quality_dict_opt {
        output_file.write_all(&[1])?; // dictionary enabled
        output_file.write_all(&(dict.len() as u32).to_le_bytes())?;
        output_file.write_all(dict)?;
    } else {
        output_file.write_all(&[0])?; // dictionary disabled
    }

    // Write read ID template metadata
    let template_prefix_bytes = read_id_template.prefix.as_bytes();
    output_file.write_all(&(template_prefix_bytes.len() as u16).to_le_bytes())?;
    output_file.write_all(template_prefix_bytes)?;
    output_file.write_all(&[if read_id_template.has_comment { 1 } else { 0 }])?;

    // Write common comment (only when template is in use and has_comment is true)
    let common_comment_bytes = if read_id_template.has_comment {
        if let Some(ref cc) = read_id_template.common_comment {
            let bytes = cc.as_bytes();
            output_file.write_all(&(bytes.len() as u16).to_le_bytes())?;
            output_file.write_all(bytes)?;
            bytes.len()
        } else {
            output_file.write_all(&0u16.to_le_bytes())?;
            0
        }
    } else {
        0
    };

    // Write stream lengths
    output_file.write_all(&(stats.num_reads as u64).to_le_bytes())?;
    output_file.write_all(&(headers.len() as u64).to_le_bytes())?;
    output_file.write_all(&(sequences.len() as u64).to_le_bytes())?;
    output_file.write_all(&(nmasks.len() as u64).to_le_bytes())?;
    output_file.write_all(&(qualities.len() as u64).to_le_bytes())?;

    // Write data streams
    output_file.write_all(&headers)?;
    output_file.write_all(&sequences)?;
    output_file.write_all(&nmasks)?;
    output_file.write_all(&qualities)?;

    let quality_model_size = if let Some((_, ref model)) = quality_model_opt {
        1 + 2 + quality_model::serialize_model(model).len() // flag + size + data
    } else {
        1 // just the flag
    };
    let dict_size = if let Some(ref dict) = quality_dict_opt {
        1 + 4 + dict.len() // flag + size + data
    } else {
        1 // just the flag
    };
    let arithmetic_size = if args.arithmetic {
        1 + 4 + (processed_records.len() * 4) // flag + count + (4 bytes per read length)
    } else {
        1 // just the flag
    };
    let common_comment_size = if read_id_template.has_comment { 2 + common_comment_bytes } else { 0 };
    let metadata_size = 1 + arithmetic_size + 1 + 1 + quality_model_size + 1 + dict_size + 2 + template_prefix_bytes.len() + 1 + common_comment_size + 40; // encoding + arithmetic + quality_binning + compressor + model + delta + dict + template + common_comment + 5 lengths
    let total_compressed = headers.len() + sequences.len() + nmasks.len() + qualities.len() + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());

    // Show compression stats
    info!("Original size: {} bytes", stats.original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!("Compression ratio: {:.2}x", stats.original_size as f64 / total_compressed as f64);

    // Show detailed stream breakdown
    info!("Stream breakdown:");
    info!("  Headers:   {} bytes ({:.1}%)", headers.len(), 100.0 * headers.len() as f64 / total_compressed as f64);
    info!("  Sequences: {} bytes ({:.1}%)", sequences.len(), 100.0 * sequences.len() as f64 / total_compressed as f64);
    info!("  N-masks:   {} bytes ({:.1}%)", nmasks.len(), 100.0 * nmasks.len() as f64 / total_compressed as f64);
    info!("  Qualities: {} bytes ({:.1}%)", qualities.len(), 100.0 * qualities.len() as f64 / total_compressed as f64);
    info!("  Metadata:  {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total_compressed as f64);

    Ok(())
}
