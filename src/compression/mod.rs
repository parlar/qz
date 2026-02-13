pub mod benchmark;
mod columnar;
mod delta;
mod n_mask;
pub mod paired_end;
mod quality;
mod quality_delta;
mod quality_model;
mod read_id;
mod reorder;
pub mod reorder_safe;
mod rle;
mod zstd_dict;
pub mod arithmetic;
pub mod arithmetic_quality;
pub mod arithmetic_sequence;
pub mod bsc;
pub mod debruijn;
pub mod dna_utils;

use crate::cli::{CompressArgs, DecompressArgs, QualityMode, QualityCompressor, SequenceCompressor, HeaderCompressor};
use anyhow::Result;
use std::time::Instant;
use tracing::info;

/// Convert QualityMode to QualityBinning
fn quality_mode_to_binning(mode: QualityMode) -> QualityBinning {
    match mode {
        QualityMode::Lossless => QualityBinning::None, // None = lossless (keep original)
        QualityMode::IlluminaBin => QualityBinning::Illumina8,
        QualityMode::Illumina4 => QualityBinning::FourLevel,
        QualityMode::Binary => QualityBinning::Binary { threshold: 20 },
        QualityMode::Qvz => QualityBinning::None, // TODO: Implement QVZ
        QualityMode::Discard => QualityBinning::Binary { threshold: 255 }, // All qualities -> same value
    }
}

/// Convert QualityBinning to a byte code for storage
fn binning_to_code(binning: QualityBinning) -> u8 {
    match binning {
        QualityBinning::None => 0, // Lossless
        QualityBinning::Illumina8 => 1,
        QualityBinning::Binary { .. } => 2,
        QualityBinning::FourLevel => 3,
    }
}

/// Convert byte code back to QualityBinning
fn code_to_binning(code: u8) -> Result<QualityBinning> {
    match code {
        0 => Ok(QualityBinning::None), // Lossless
        1 => Ok(QualityBinning::Illumina8),
        2 => Ok(QualityBinning::Binary { threshold: 20 }),
        3 => Ok(QualityBinning::FourLevel),
        _ => anyhow::bail!("Invalid quality binning code: {}", code),
    }
}

/// Convert QualityCompressor to a byte code for storage
fn compressor_to_code(compressor: QualityCompressor) -> u8 {
    match compressor {
        QualityCompressor::Zstd => 0,
        QualityCompressor::Bsc => 1,
    }
}

/// Convert byte code back to QualityCompressor
fn code_to_compressor(code: u8) -> Result<QualityCompressor> {
    match code {
        0 => Ok(QualityCompressor::Zstd),
        1 => Ok(QualityCompressor::Bsc),
        _ => anyhow::bail!("Invalid quality compressor code: {}", code),
    }
}

/// Convert SequenceCompressor to a byte code for storage
fn seq_compressor_to_code(compressor: SequenceCompressor) -> u8 {
    match compressor {
        SequenceCompressor::Zstd => 0,
        SequenceCompressor::Bsc => 1,
    }
}

/// Convert byte code back to SequenceCompressor
fn code_to_seq_compressor(code: u8) -> Result<SequenceCompressor> {
    match code {
        0 => Ok(SequenceCompressor::Zstd),
        1 => Ok(SequenceCompressor::Bsc),
        _ => anyhow::bail!("Invalid sequence compressor code: {}", code),
    }
}

/// Convert HeaderCompressor to a byte code for storage
fn header_compressor_to_code(compressor: HeaderCompressor) -> u8 {
    match compressor {
        HeaderCompressor::Zstd => 0,
        HeaderCompressor::Bsc => 1,
    }
}

/// Convert byte code back to HeaderCompressor
fn code_to_header_compressor(code: u8) -> Result<HeaderCompressor> {
    match code {
        0 => Ok(HeaderCompressor::Zstd),
        1 => Ok(HeaderCompressor::Bsc),
        _ => anyhow::bail!("Invalid header compressor code: {}", code),
    }
}

pub use columnar::{
    compress_columnar,
    QualityBinning,
};

// Advanced and tier2 features not currently exposed in public API

/// Dispatch BSC parallel compression: adaptive (default) or static (--bsc-static).
fn bsc_compress_parallel(data: &[u8], use_static: bool) -> Result<Vec<u8>> {
    if use_static {
        bsc::compress_parallel(data)
    } else {
        bsc::compress_parallel_adaptive(data)
    }
}

/// Memory-efficient streaming compression for the default BSC path.
/// Builds compression streams during FASTQ reading instead of storing all records.
/// Reduces peak memory from ~2x input size (records + streams) to ~1x (streams only).
fn compress_streaming_bsc(args: &CompressArgs) -> Result<()> {
    use crate::io::FastqReader;
    use std::io::Write;

    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.bsc_static;
    let quality_binning = if args.no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(args.quality_mode)
    };

    info!("Streaming mode: building compression streams during read");
    info!("Reading FASTQ and building streams...");

    let mut reader = FastqReader::from_path(input_path, args.fasta)?;
    let mut header_stream: Vec<u8> = Vec::new();
    let mut seq_stream: Vec<u8> = Vec::new();
    let mut qual_stream: Vec<u8> = Vec::new();
    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;

    while let Some(record) = reader.next()? {
        // Header stream: varint(len) + raw bytes
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.extend_from_slice(record.id.as_bytes());

        // Sequence stream: varint(len) + raw bytes
        write_varint(&mut seq_stream, record.sequence.len())?;
        seq_stream.extend_from_slice(record.sequence.as_bytes());

        // Quality stream: varint(len) + packed bytes
        if !args.no_quality {
            if let Some(ref qual) = record.quality {
                let quantized = quality::quantize_quality(qual, args.quality_mode);
                write_varint(&mut qual_stream, quantized.len())?;
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                qual_stream.extend_from_slice(&packed);
            }
        }

        total_bases += record.sequence.len();
        original_size += record.id.len() + record.sequence.len()
            + record.quality.as_ref().map(|q| q.len()).unwrap_or(0) + 3;
        num_reads += 1;

        if num_reads % 10_000_000 == 0 {
            info!("  {} million records read...", num_reads / 1_000_000);
        }
    }

    info!("Read {} records ({} bases)", num_reads, total_bases);
    info!("Stream sizes: headers={} seq={} qual={}",
        header_stream.len(), seq_stream.len(), qual_stream.len());

    // Compress all three streams in parallel with BSC
    info!("Compressing streams with BSC{}...",
        if bsc_static { " (static)" } else { " (adaptive)" });

    let (header_result, (seq_result, qual_result)) = rayon::join(
        || bsc_compress_parallel(&header_stream, bsc_static),
        || rayon::join(
            || bsc_compress_parallel(&seq_stream, bsc_static),
            || if args.no_quality {
                Ok(Vec::new())
            } else {
                bsc_compress_parallel(&qual_stream, bsc_static)
            },
        ),
    );

    drop(header_stream);
    drop(seq_stream);
    drop(qual_stream);

    let headers = header_result?;
    let sequences = seq_result?;
    let qualities = qual_result?;

    // Write output file (same archive format as non-streaming path)
    info!("Writing output file...");
    let mut output_file = std::fs::File::create(&args.output)?;

    output_file.write_all(&[0u8])?;                                             // encoding_type = 0
    output_file.write_all(&[0u8])?;                                             // arithmetic = disabled
    output_file.write_all(&[binning_to_code(quality_binning)])?;                // quality_binning
    output_file.write_all(&[compressor_to_code(QualityCompressor::Bsc)])?;      // quality_compressor
    output_file.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?; // seq_compressor
    output_file.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;// header_compressor
    output_file.write_all(&[0u8])?;                                             // quality_model = disabled
    output_file.write_all(&[0u8])?;                                             // quality_delta = disabled
    output_file.write_all(&[0u8])?;                                             // quality_dict = disabled
    output_file.write_all(&0u16.to_le_bytes())?;                                // template_prefix_len = 0
    output_file.write_all(&[0u8])?;                                             // has_comment = false

    // Stream lengths
    output_file.write_all(&(num_reads as u64).to_le_bytes())?;
    output_file.write_all(&(headers.len() as u64).to_le_bytes())?;
    output_file.write_all(&(sequences.len() as u64).to_le_bytes())?;
    output_file.write_all(&0u64.to_le_bytes())?;                                // nmasks_len = 0
    output_file.write_all(&(qualities.len() as u64).to_le_bytes())?;

    // Data streams
    output_file.write_all(&headers)?;
    output_file.write_all(&sequences)?;
    output_file.write_all(&qualities)?;

    let metadata_size = 9 + 2 + 1 + 40; // 9 flag bytes + template header (2+1) + 5 u64 lengths
    let total_compressed = headers.len() + sequences.len() + qualities.len() + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!("Compression ratio: {:.2}x", original_size as f64 / total_compressed as f64);
    info!("Stream breakdown:");
    info!("  Headers:   {} bytes ({:.1}%)", headers.len(), 100.0 * headers.len() as f64 / total_compressed as f64);
    info!("  Sequences: {} bytes ({:.1}%)", sequences.len(), 100.0 * sequences.len() as f64 / total_compressed as f64);
    info!("  Qualities: {} bytes ({:.1}%)", qualities.len(), 100.0 * qualities.len() as f64 / total_compressed as f64);
    info!("  Metadata:  {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total_compressed as f64);

    Ok(())
}

/// Read a chunk of FASTQ records into raw compression streams.
/// Returns (header_stream, seq_stream, qual_stream, num_reads, total_bases, original_size).
fn read_chunk_streams<R: std::io::BufRead>(
    reader: &mut crate::io::FastqReader<R>,
    chunk_size: usize,
    quality_mode: QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>, usize, usize, usize)> {
    let mut header_stream = Vec::new();
    let mut seq_stream = Vec::new();
    let mut qual_stream = Vec::new();
    let mut chunk_reads = 0usize;
    let mut chunk_bases = 0usize;
    let mut chunk_orig = 0usize;

    for _ in 0..chunk_size {
        match reader.next()? {
            Some(record) => {
                write_varint(&mut header_stream, record.id.len())?;
                header_stream.extend_from_slice(record.id.as_bytes());

                write_varint(&mut seq_stream, record.sequence.len())?;
                seq_stream.extend_from_slice(record.sequence.as_bytes());

                if !no_quality {
                    if let Some(ref qual) = record.quality {
                        let quantized = quality::quantize_quality(qual, quality_mode);
                        write_varint(&mut qual_stream, quantized.len())?;
                        let packed = columnar::pack_qualities(&quantized, quality_binning);
                        qual_stream.extend_from_slice(&packed);
                    }
                }

                chunk_bases += record.sequence.len();
                chunk_orig += record.id.len() + record.sequence.len()
                    + record.quality.as_ref().map(|q| q.len()).unwrap_or(0) + 3;
                chunk_reads += 1;
            }
            None => break,
        }
    }

    Ok((header_stream, seq_stream, qual_stream, chunk_reads, chunk_bases, chunk_orig))
}

/// Format bytes as human-readable string (e.g., "1.23 GiB").
fn humanize_bytes(bytes: usize) -> String {
    const KIB: f64 = 1024.0;
    const MIB: f64 = KIB * 1024.0;
    const GIB: f64 = MIB * 1024.0;
    let b = bytes as f64;
    if b >= GIB { format!("{:.2} GiB", b / GIB) }
    else if b >= MIB { format!("{:.1} MiB", b / MIB) }
    else if b >= KIB { format!("{:.0} KiB", b / KIB) }
    else { format!("{} B", bytes) }
}

/// Compress a single stream into 25 MB BSC blocks using rayon par_iter.
/// Blocks are shrunk to actual compressed size to avoid excess capacity
/// (BSC allocates ~25 MB output buffers but compressed data is typically ~9 MB).
fn compress_stream_to_bsc_blocks(data: &[u8], bsc_static: bool) -> Result<Vec<Vec<u8>>> {
    use rayon::prelude::*;
    const BSC_BLOCK_SIZE: usize = 25 * 1024 * 1024;

    if data.is_empty() {
        return Ok(Vec::new());
    }

    let compress_fn: fn(&[u8]) -> Result<Vec<u8>> = if bsc_static {
        bsc::compress
    } else {
        bsc::compress_adaptive
    };

    let chunks: Vec<&[u8]> = data.chunks(BSC_BLOCK_SIZE).collect();
    let mut blocks: Vec<Vec<u8>> = chunks.par_iter()
        .map(|chunk| compress_fn(chunk))
        .collect::<Result<Vec<_>>>()?;

    // Shrink each block to actual compressed size. Without this, 1400
    // accumulated blocks x 25 MB excess capacity = ~35 GB wasted.
    for b in &mut blocks {
        b.shrink_to_fit();
    }

    Ok(blocks)
}

/// Write compressed blocks to a temp file in BSC multi-block body format:
/// [block_len: u32, block_data]... (the num_blocks header is written separately).
fn write_blocks_to_tmp(blocks: Vec<Vec<u8>>, tmp: &mut std::io::BufWriter<std::fs::File>) -> Result<u32> {
    use std::io::Write;
    let mut count = 0u32;
    for block in blocks {
        tmp.write_all(&(block.len() as u32).to_le_bytes())?;
        tmp.write_all(&block)?;
        count += 1;
    }
    Ok(count)
}

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
fn compress_chunked_bsc(args: &CompressArgs) -> Result<()> {
    use std::io::{Write, BufWriter, BufReader};
    const CHUNK_SIZE: usize = 5_000_000; // 5M records per chunk

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

    info!("Chunked streaming mode: {} records per chunk, pipelined I/O", CHUNK_SIZE);

    let mut reader = crate::io::FastqReader::from_path(input_path, args.fasta)?;

    // Temp files for streaming compressed blocks to disk (constant memory)
    let working_dir = &args.working_dir;
    let h_tmp_path = working_dir.join(".fqz_chunked_h.tmp");
    let s_tmp_path = working_dir.join(".fqz_chunked_s.tmp");
    let q_tmp_path = working_dir.join(".fqz_chunked_q.tmp");

    // Cleanup guard: remove temp files on drop (handles success, error, and panic)
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

    // Read first chunk
    let (mut cur_h, mut cur_s, mut cur_q, mut cur_reads, mut cur_bases, mut cur_orig) =
        read_chunk_streams(&mut reader, CHUNK_SIZE, quality_mode, quality_binning, no_quality)?;

    while cur_reads > 0 {
        info!("Chunk {}: {} reads (h={} s={} q={} bytes)",
            chunk_idx, cur_reads, cur_h.len(), cur_s.len(), cur_q.len());

        // Pipeline: compress current chunk on a background thread while
        // reading the next chunk on the main thread.
        //
        // The compression thread processes streams SEQUENTIALLY (headers ->
        // sequences -> qualities) to limit BSC working memory. Raw data for
        // each stream is dropped after compression before starting the next.
        let (next_result, compress_result) = std::thread::scope(|scope| {
            let h_data = std::mem::take(&mut cur_h);
            let s_data = std::mem::take(&mut cur_s);
            let q_data = std::mem::take(&mut cur_q);

            let compress_handle = scope.spawn(move || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                // Compress headers, then free raw header data before starting sequences
                let h_blocks = compress_stream_to_bsc_blocks(&h_data, bsc_static)?;
                drop(h_data);

                // Compress sequences, then free raw sequence data before starting qualities
                let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
                drop(s_data);

                // Compress qualities (if present)
                let q_blocks = if no_quality || q_data.is_empty() {
                    Vec::new()
                } else {
                    compress_stream_to_bsc_blocks(&q_data, bsc_static)?
                };

                Ok((h_blocks, s_blocks, q_blocks))
            });

            // Read next chunk on main thread (overlaps with compression)
            let next = read_chunk_streams(
                &mut reader, CHUNK_SIZE, quality_mode, quality_binning, no_quality,
            );

            let compressed = compress_handle.join().unwrap();
            (next, compressed)
        });

        // Write compressed blocks to temp files immediately, freeing block memory
        let (h_blk, s_blk, q_blk) = compress_result?;
        h_num_blocks += write_blocks_to_tmp(h_blk, &mut h_tmp)?;
        s_num_blocks += write_blocks_to_tmp(s_blk, &mut s_tmp)?;
        q_num_blocks += write_blocks_to_tmp(q_blk, &mut q_tmp)?;

        num_reads += cur_reads;
        total_bases += cur_bases;
        original_size += cur_orig;
        chunk_idx += 1;

        let (nh, ns, nq, nr, nb, no) = next_result?;
        cur_h = nh;
        cur_s = ns;
        cur_q = nq;
        cur_reads = nr;
        cur_bases = nb;
        cur_orig = no;
    }

    // Flush and close temp file writers
    h_tmp.flush()?;
    s_tmp.flush()?;
    q_tmp.flush()?;
    drop(h_tmp);
    drop(s_tmp);
    drop(q_tmp);

    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    // Stream sizes in multi-block format: [num_blocks: u32] + block data from temp file
    let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
    let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    info!("Read {} records in {} chunks ({} bases)", num_reads, chunk_idx, total_bases);
    info!("Compressed blocks: headers={} ({}) seq={} ({}) qual={} ({})",
        h_num_blocks, humanize_bytes(h_tmp_size),
        s_num_blocks, humanize_bytes(s_tmp_size),
        q_num_blocks, humanize_bytes(q_tmp_size));

    // Write final output file (identical archive format to compress_streaming_bsc)
    info!("Writing output file...");
    let mut output_file = std::fs::File::create(&args.output)?;

    output_file.write_all(&[0u8])?;                                             // encoding_type = 0
    output_file.write_all(&[0u8])?;                                             // arithmetic = disabled
    output_file.write_all(&[binning_to_code(quality_binning)])?;                // quality_binning
    output_file.write_all(&[compressor_to_code(QualityCompressor::Bsc)])?;      // quality_compressor
    output_file.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?; // seq_compressor
    output_file.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;// header_compressor
    output_file.write_all(&[0u8])?;                                             // quality_model = disabled
    output_file.write_all(&[0u8])?;                                             // quality_delta = disabled
    output_file.write_all(&[0u8])?;                                             // quality_dict = disabled
    output_file.write_all(&0u16.to_le_bytes())?;                                // template_prefix_len = 0
    output_file.write_all(&[0u8])?;                                             // has_comment = false

    // Stream lengths
    output_file.write_all(&(num_reads as u64).to_le_bytes())?;
    output_file.write_all(&(headers_len as u64).to_le_bytes())?;
    output_file.write_all(&(sequences_len as u64).to_le_bytes())?;
    output_file.write_all(&0u64.to_le_bytes())?;                                // nmasks_len = 0
    output_file.write_all(&(qualities_len as u64).to_le_bytes())?;

    // Copy streams from temp files into output (streaming copy, ~zero memory)
    let copy_stream = |num_blocks: u32, tmp_path: &std::path::Path, out: &mut std::fs::File| -> Result<()> {
        if num_blocks == 0 { return Ok(()); }
        out.write_all(&num_blocks.to_le_bytes())?;
        let mut tmp_reader = BufReader::new(std::fs::File::open(tmp_path)?);
        std::io::copy(&mut tmp_reader, out)?;
        Ok(())
    };

    copy_stream(h_num_blocks, &h_tmp_path, &mut output_file)?;
    copy_stream(s_num_blocks, &s_tmp_path, &mut output_file)?;
    copy_stream(q_num_blocks, &q_tmp_path, &mut output_file)?;

    let metadata_size = 9 + 2 + 1 + 40;
    let total_compressed = headers_len + sequences_len + qualities_len + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!("Compression ratio: {:.2}x", original_size as f64 / total_compressed as f64);
    info!("Stream breakdown:");
    info!("  Headers:   {} bytes ({:.1}%)", headers_len, 100.0 * headers_len as f64 / total_compressed as f64);
    info!("  Sequences: {} bytes ({:.1}%)", sequences_len, 100.0 * sequences_len as f64 / total_compressed as f64);
    info!("  Qualities: {} bytes ({:.1}%)", qualities_len, 100.0 * qualities_len as f64 / total_compressed as f64);
    info!("  Metadata:  {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total_compressed as f64);

    // Temp files cleaned up by TmpCleanup drop guard
    Ok(())
}

pub fn compress(args: &CompressArgs) -> Result<()> {
    use crate::io::{FastqReader, FastqRecord};
    use std::io::Write;

    let start_time = Instant::now();

    // Set up thread pool
    let num_threads = if args.threads == 0 {
        num_cpus::get()
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

    // Detect paired-end mode
    if args.input.len() == 2 {
        info!("Paired-end mode detected - using correlation compression");
        return compress_paired_end(&args.input[0], &args.input[1], args);
    } else if args.input.len() > 2 {
        anyhow::bail!("More than 2 input files not supported");
    }

    // Use memory-efficient streaming mode for the default BSC path
    // (no reordering, no advanced quality encoding, BSC for all streams)
    let can_stream = args.sequence_compressor == SequenceCompressor::Bsc
        && args.header_compressor == HeaderCompressor::Bsc
        && args.quality_compressor == QualityCompressor::Bsc
        && !args.delta_encoding
        && !args.rle_encoding
        && !args.debruijn
        && !args.arithmetic
        && !args.quality_modeling
        && !args.quality_delta
        && !args.dict_training
        && args.reorder_by == crate::cli::ReorderMode::None
        && !args.allow_reordering;

    if can_stream {
        if args.chunked {
            return compress_chunked_bsc(args);
        }
        return compress_streaming_bsc(args);
    }

    // Non-streaming mode: read all records into memory first
    // (required for reordering, advanced quality encoding, non-BSC compressors)
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
                    record.quality = Some(quality::quantize_quality(qual, args.quality_mode));
                }
                record
            })
            .collect()
    } else {
        records
    };

    // Apply patent-safe reordering if requested
    let mut processed_records = processed_records;
    if args.reorder_by != crate::cli::ReorderMode::None {
        info!("Reordering reads using {:?} strategy (patent-safe)...", args.reorder_by);
        let reorder_mode = match args.reorder_by {
            crate::cli::ReorderMode::None => reorder_safe::ReorderMode::None,
            crate::cli::ReorderMode::Flowcell => reorder_safe::ReorderMode::Flowcell,
            crate::cli::ReorderMode::Gc => reorder_safe::ReorderMode::GcContent,
            crate::cli::ReorderMode::Length => reorder_safe::ReorderMode::Length,
            crate::cli::ReorderMode::Lexicographic => reorder_safe::ReorderMode::Lexicographic,
            crate::cli::ReorderMode::Smart => reorder_safe::ReorderMode::Smart,
        };
        reorder_safe::reorder_reads(&mut processed_records, reorder_mode)?;
    }

    // Reorder if enabled (homology-based - may have patent implications)
    let processed_records = if args.allow_reordering {
        info!("Reordering reads for better compression (homology-based)...");
        reorder::reorder_auto(processed_records)?
    } else {
        processed_records
    };

    // Determine encoding type and quality binning
    let encoding_type: u8 = if args.delta_encoding {
        1 // Delta encoding
    } else if args.rle_encoding {
        2 // RLE encoding
    } else if args.debruijn {
        3 // De Bruijn graph encoding
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
        if !args.allow_reordering {
            info!("WARNING: Quality delta encoding works best with --allow-reordering");
        }
        info!("Applying quality delta encoding...");
        let qual_compressed = if let Some(ref dict) = quality_dict_opt {
            compress_qualities_with_delta_and_dict(&processed_records, dict, args.compression_level)?
        } else {
            compress_qualities_with_delta(&processed_records, args.compression_level)?
        };
        (Some(qual_compressed), None)
    } else if args.quality_modeling && !args.no_quality {
        info!("Building positional quality model...");
        let (qual_compressed, model) = if let Some(ref dict) = quality_dict_opt {
            compress_qualities_with_model_and_dict(&processed_records, dict, args.compression_level)?
        } else {
            compress_qualities_with_model(&processed_records, args.compression_level)?
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
            compress_qualities_with_dict(&processed_records, quality_binning, dict, args.compression_level, args.quality_compressor)
        } else {
            compress_qualities_with(&processed_records, quality_binning, args.compression_level, args.quality_compressor, bsc_static)
        }
    };

    // Compress all three streams in parallel:
    //   Left:  headers
    //   Right: sequences + qualities (with nested rayon::join)
    let (header_result, seq_qual_result) = rayon::join(
        || -> Result<(Vec<u8>, read_id::ReadIdTemplate)> {
            match args.header_compressor {
                HeaderCompressor::Bsc => {
                    info!("Compressing read IDs with raw BSC{}...", if args.bsc_static { " (static)" } else { " (adaptive)" });
                    let compressed = compress_headers_bsc_with(&processed_records, args.bsc_static)?;
                    let dummy_template = read_id::ReadIdTemplate {
                        prefix: String::new(),
                        has_comment: false,
                        common_comment: None,
                    };
                    Ok((compressed, dummy_template))
                }
                HeaderCompressor::Zstd => {
                    info!("Compressing read IDs with template encoding...");
                    compress_headers(&processed_records, args.compression_level)
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
                    || compress_sequences_arithmetic(&processed_records),
                    || -> Result<Vec<u8>> {
                        if let Some(ref qual) = quality_delta_opt {
                            Ok(qual.clone())
                        } else if let Some((ref qual, _)) = quality_model_opt {
                            Ok(qual.clone())
                        } else if args.no_quality {
                            Ok(Vec::new())
                        } else {
                            compress_qualities_arithmetic(&processed_records)
                        }
                    },
                );
                Ok((seq_result?, Vec::new(), qual_result?))
            } else {
                match args.sequence_compressor {
                    SequenceCompressor::Bsc => {
                        info!("Compressing sequences and qualities in parallel with BSC{}...", if args.bsc_static { " (static)" } else { " (adaptive)" });
                        let bsc_static = args.bsc_static;
                        let (seq_result, qual_result) = rayon::join(
                            || compress_sequences_raw_bsc_with(&processed_records, bsc_static),
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
    let metadata_size = 1 + arithmetic_size + 1 + 1 + quality_model_size + 1 + dict_size + 2 + template_prefix_bytes.len() + 1 + 40; // encoding + arithmetic + quality_binning + compressor + model + delta + dict + template + 5 lengths
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

/// Block-by-block BSC stream decoder for memory-efficient decompression.
///
/// Reads one 25 MB BSC block at a time from the archive file, handling
/// records that span block boundaries via an overlap buffer.
/// Peak memory per decoder: ~25 MB (one decompressed block + small overlap).
/// Reads compressed BSC blocks from a stream region of the archive file,
/// decompresses them in parallel batches (using rayon), and sends
/// decompressed blocks through a bounded channel to the consumer thread.
fn stream_decompressor(
    path: &std::path::Path,
    offset: u64,
    stream_len: usize,
    tx: std::sync::mpsc::SyncSender<std::result::Result<Vec<u8>, String>>,
) {
    if let Err(e) = stream_decompressor_inner(path, offset, stream_len, &tx) {
        let _ = tx.send(Err(e.to_string()));
    }
}

fn stream_decompressor_inner(
    path: &std::path::Path,
    offset: u64,
    stream_len: usize,
    tx: &std::sync::mpsc::SyncSender<std::result::Result<Vec<u8>, String>>,
) -> Result<()> {
    use std::io::{Read, Seek, SeekFrom};
    use rayon::prelude::*;

    if stream_len == 0 {
        return Ok(());
    }

    let mut file = std::io::BufReader::with_capacity(
        4 * 1024 * 1024,
        std::fs::File::open(path)?,
    );
    file.seek(SeekFrom::Start(offset))?;

    let mut buf4 = [0u8; 4];
    file.read_exact(&mut buf4)?;
    let num_blocks = u32::from_le_bytes(buf4) as usize;

    const BATCH: usize = 8;
    let mut blocks_done = 0;

    while blocks_done < num_blocks {
        let batch_size = BATCH.min(num_blocks - blocks_done);

        // Read compressed blocks from file (sequential I/O)
        let mut compressed_blocks = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            file.read_exact(&mut buf4)?;
            let block_len = u32::from_le_bytes(buf4) as usize;
            let mut data = vec![0u8; block_len];
            file.read_exact(&mut data)?;
            compressed_blocks.push(data);
        }

        // Decompress batch in parallel using rayon
        let decompressed: Result<Vec<Vec<u8>>> = compressed_blocks
            .into_par_iter()
            .map(|block| bsc::decompress(&block))
            .collect();
        let decompressed = decompressed?;

        // Send decompressed blocks through channel (may block if full)
        for block in decompressed {
            if tx.send(Ok(block)).is_err() {
                return Ok(()); // Receiver dropped
            }
        }

        blocks_done += batch_size;
    }

    Ok(())
}

/// Buffer that receives decompressed blocks from a background decompressor
/// thread via a bounded channel. Provides varint/bytes reading interface
/// for streaming record reconstruction.
struct ChannelStreamBuffer {
    rx: std::sync::mpsc::Receiver<std::result::Result<Vec<u8>, String>>,
    buf: Vec<u8>,
    pos: usize,
}

impl ChannelStreamBuffer {
    fn new(rx: std::sync::mpsc::Receiver<std::result::Result<Vec<u8>, String>>) -> Self {
        Self {
            rx,
            buf: Vec::new(),
            pos: 0,
        }
    }

    /// Ensure at least `min` unread bytes are available.
    /// Receives decompressed blocks from the channel as needed.
    fn fill(&mut self, min: usize) -> Result<bool> {
        while self.buf.len() - self.pos < min {
            // Compact consumed data
            if self.pos > 0 {
                self.buf.drain(..self.pos);
                self.pos = 0;
            }
            match self.rx.recv() {
                Ok(Ok(block)) => self.buf.extend_from_slice(&block),
                Ok(Err(e)) => anyhow::bail!("Stream decompressor error: {}", e),
                Err(_) => return Ok(self.buf.len() - self.pos >= min),
            }
        }
        Ok(true)
    }

    fn read_varint(&mut self) -> Result<usize> {
        let mut value = 0usize;
        let mut shift = 0;
        loop {
            if !self.fill(1)? {
                anyhow::bail!("Unexpected end of stream reading varint");
            }
            let byte = self.buf[self.pos];
            self.pos += 1;
            value |= ((byte & 0x7F) as usize) << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
            shift += 7;
        }
    }

    fn read_bytes(&mut self, n: usize) -> Result<&[u8]> {
        if !self.fill(n)? {
            anyhow::bail!(
                "Unexpected end of stream: needed {} bytes, have {}",
                n,
                self.buf.len() - self.pos
            );
        }
        let start = self.pos;
        self.pos += n;
        Ok(&self.buf[start..start + n])
    }
}

/// Check if archive can use the streaming BSC decompression path.
/// Reads only the first 12 bytes of the header (no bulk load).
fn can_stream_decompress(args: &DecompressArgs) -> Result<bool> {
    use std::io::Read;

    if args.range.is_some() {
        return Ok(false); // Range queries need random access
    }

    let mut file = std::fs::File::open(&args.input)?;
    let mut header = [0u8; 12];
    if file.read_exact(&mut header).is_err() {
        return Ok(false);
    }

    let encoding_type = header[0];
    let arithmetic_enabled = header[1] != 0;
    let quality_compressor = header[3];
    let sequence_compressor = header[4];
    let header_compressor = header[5];
    let quality_model_enabled = header[6] != 0;
    let quality_delta_enabled = header[7] != 0;
    let quality_dict_present = header[8] != 0;
    let template_prefix_len = u16::from_le_bytes([header[9], header[10]]) as usize;

    Ok(encoding_type == 0
        && !arithmetic_enabled
        && !quality_model_enabled
        && !quality_delta_enabled
        && !quality_dict_present
        && template_prefix_len == 0
        && quality_compressor == 1  // BSC
        && sequence_compressor == 1 // BSC
        && header_compressor == 1)  // BSC
}

/// Memory-efficient parallel streaming decompression for BSC archives.
///
/// Spawns 3 decompressor threads (one per stream: headers, sequences, qualities)
/// that read and decompress BSC blocks in parallel batches using rayon.
/// Decompressed blocks are fed through bounded channels to the main thread,
/// which reconstructs FASTQ records and writes them to the output file.
///
/// Peak memory: ~300 MB (3 streams × ~4 decompressed blocks × 25 MB + output buffer)
/// regardless of input size. Uses all available CPU cores via rayon.
fn decompress_streaming_bsc(args: &DecompressArgs) -> Result<()> {
    use std::io::{Read, Write};

    let start_time = Instant::now();

    info!("Input file: {:?}", args.input);
    info!("Output files: {:?}", args.output);
    info!("Streaming decompression mode (parallel BSC)");

    // Read archive header (52 bytes for BSC path with no template/model/dict)
    let mut file = std::fs::File::open(&args.input)?;
    let mut header = [0u8; 52];
    file.read_exact(&mut header)?;
    drop(file);

    let quality_binning = code_to_binning(header[2])?;
    let _has_comment = header[11] != 0;

    let num_reads = u64::from_le_bytes(header[12..20].try_into().unwrap()) as usize;
    let headers_len = u64::from_le_bytes(header[20..28].try_into().unwrap()) as usize;
    let sequences_len = u64::from_le_bytes(header[28..36].try_into().unwrap()) as usize;
    let nmasks_len = u64::from_le_bytes(header[36..44].try_into().unwrap()) as usize;
    let qualities_len = u64::from_le_bytes(header[44..52].try_into().unwrap()) as usize;

    info!("Archive: {} reads, headers={} seq={} qual={}",
        num_reads, humanize_bytes(headers_len), humanize_bytes(sequences_len), humanize_bytes(qualities_len));

    // Compute stream positions in file
    let data_start = 52u64;
    let h_offset = data_start;
    let s_offset = h_offset + headers_len as u64;
    let q_offset = s_offset + sequences_len as u64 + nmasks_len as u64;
    let has_quality = qualities_len > 0;
    let bits_per_qual = quality_binning.bits_per_quality();

    // Open output file
    if args.output.is_empty() {
        anyhow::bail!("No output file specified");
    }
    let output_path = &args.output[0];
    let mut output: Box<dyn Write> = if args.gzipped {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        Box::new(std::io::BufWriter::with_capacity(
            8 * 1024 * 1024,
            GzEncoder::new(
                std::fs::File::create(output_path)?,
                Compression::new(args.gzip_level),
            ),
        ))
    } else {
        Box::new(std::io::BufWriter::with_capacity(
            8 * 1024 * 1024,
            std::fs::File::create(output_path)?,
        ))
    };

    // Spawn parallel decompressor threads and stream records to output
    info!("Decompressing {} records with parallel BSC...", num_reads);
    let archive_path = &args.input;

    std::thread::scope(|scope| -> Result<()> {
        // Bounded channels: capacity 2 = max 2 pre-decompressed blocks buffered per stream
        let (h_tx, h_rx) = std::sync::mpsc::sync_channel(2);
        let (s_tx, s_rx) = std::sync::mpsc::sync_channel(2);
        let h_path = archive_path.clone();
        let s_path = archive_path.clone();

        scope.spawn(move || stream_decompressor(&h_path, h_offset, headers_len, h_tx));
        scope.spawn(move || stream_decompressor(&s_path, s_offset, sequences_len, s_tx));

        let q_rx = if has_quality {
            let (q_tx, q_rx) = std::sync::mpsc::sync_channel(2);
            let q_path = archive_path.clone();
            scope.spawn(move || stream_decompressor(&q_path, q_offset, qualities_len, q_tx));
            Some(q_rx)
        } else {
            None
        };

        let mut h_buf = ChannelStreamBuffer::new(h_rx);
        let mut s_buf = ChannelStreamBuffer::new(s_rx);
        let mut q_buf = q_rx.map(ChannelStreamBuffer::new);

        for i in 0..num_reads {
            // Header: varint(len) + raw bytes
            let h_len = h_buf.read_varint()?;
            output.write_all(h_buf.read_bytes(h_len)?)?;
            output.write_all(b"\n")?;

            // Sequence: varint(len) + raw bytes
            let s_len = s_buf.read_varint()?;
            output.write_all(s_buf.read_bytes(s_len)?)?;
            output.write_all(b"\n+\n")?;

            // Quality: varint(orig_len) + packed bytes
            if let Some(ref mut q) = q_buf {
                let q_len = q.read_varint()?;
                let packed_len = (q_len * bits_per_qual + 7) / 8;
                let packed = q.read_bytes(packed_len)?;
                let quality_str = columnar::unpack_qualities(packed, q_len, quality_binning);
                output.write_all(quality_str.as_bytes())?;
            }
            output.write_all(b"\n")?;

            if (i + 1) % 10_000_000 == 0 {
                info!("  {} million records written...", (i + 1) / 1_000_000);
            }
        }

        output.flush()?;
        Ok(())
    })?;

    info!("Decompressed {} records", num_reads);
    let elapsed = start_time.elapsed();
    info!("Decompression completed in {:.2}s", elapsed.as_secs_f64());

    Ok(())
}

pub fn decompress(args: &DecompressArgs) -> Result<()> {
    use std::io::{Read, Write};

    let start_time = Instant::now();

    // Try streaming decompression for standard BSC archives (>100x less memory)
    if can_stream_decompress(args)? {
        return decompress_streaming_bsc(args);
    }

    info!("Input file: {:?}", args.input);
    info!("Output files: {:?}", args.output);

    // Read compressed archive
    info!("Reading compressed file...");
    let mut input_file = std::fs::File::open(&args.input)?;
    let mut archive_data = Vec::new();
    input_file.read_to_end(&mut archive_data)?;

    // Parse archive format
    if archive_data.len() < 50 {
        anyhow::bail!("Invalid archive: too small");
    }

    let mut offset = 0;

    // Read encoding type
    let encoding_type = archive_data[offset];
    offset += 1;

    // Read arithmetic mode flag and read lengths
    let arithmetic_enabled = archive_data[offset] != 0;
    offset += 1;
    let read_lengths_opt = if arithmetic_enabled {
        let num_lengths = u32::from_le_bytes(archive_data[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;
        let mut lengths = Vec::with_capacity(num_lengths);
        for _ in 0..num_lengths {
            let len = u32::from_le_bytes(archive_data[offset..offset + 4].try_into().unwrap()) as usize;
            offset += 4;
            lengths.push(len);
        }
        Some(lengths)
    } else {
        None
    };

    // Read quality binning mode
    let quality_binning = code_to_binning(archive_data[offset])?;
    offset += 1;

    // Read quality compressor type (default to zstd for backward compatibility)
    let quality_compressor = if offset < archive_data.len() && archive_data[offset] <= 1 {
        let compressor = code_to_compressor(archive_data[offset])?;
        offset += 1;
        compressor
    } else {
        // Old format without compressor field - assume zstd
        QualityCompressor::Zstd
    };

    // Read sequence compressor type
    let sequence_compressor = code_to_seq_compressor(archive_data[offset])?;
    offset += 1;

    // Read header compressor type
    let header_compressor = code_to_header_compressor(archive_data[offset])?;
    offset += 1;

    // Read quality modeling flag and model
    let quality_model_enabled = archive_data[offset] != 0;
    offset += 1;
    let quality_model_opt = if quality_model_enabled {
        let model_size = u16::from_le_bytes(archive_data[offset..offset + 2].try_into().unwrap()) as usize;
        offset += 2;
        let model_bytes = &archive_data[offset..offset + model_size];
        offset += model_size;
        Some(quality_model::deserialize_model(model_bytes)?)
    } else {
        None
    };

    // Read quality delta flag
    let quality_delta_enabled = archive_data[offset] != 0;
    offset += 1;

    // Read quality dictionary if present
    let quality_dict_opt = if archive_data[offset] != 0 {
        offset += 1;
        let dict_size = u32::from_le_bytes(archive_data[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;
        let dict = archive_data[offset..offset + dict_size].to_vec();
        offset += dict_size;
        Some(dict)
    } else {
        offset += 1;
        None
    };

    // Read read ID template metadata
    let template_prefix_len = u16::from_le_bytes(archive_data[offset..offset + 2].try_into().unwrap()) as usize;
    offset += 2;
    let template_prefix = String::from_utf8_lossy(&archive_data[offset..offset + template_prefix_len]).to_string();
    offset += template_prefix_len;
    let template_has_comment = archive_data[offset] != 0;
    offset += 1;

    let read_id_template = read_id::ReadIdTemplate {
        prefix: template_prefix,
        has_comment: template_has_comment,
        common_comment: None, // Will be read from archive if present
    };

    let read_u64 = |data: &[u8], off: &mut usize| -> u64 {
        let val = u64::from_le_bytes(data[*off..*off + 8].try_into().unwrap());
        *off += 8;
        val
    };

    let num_reads = read_u64(&archive_data, &mut offset) as usize;
    let headers_len = read_u64(&archive_data, &mut offset) as usize;
    let sequences_len = read_u64(&archive_data, &mut offset) as usize;
    let nmasks_len = read_u64(&archive_data, &mut offset) as usize;
    let qualities_len = read_u64(&archive_data, &mut offset) as usize;

    if archive_data.len() < offset + headers_len + sequences_len + nmasks_len + qualities_len {
        anyhow::bail!("Invalid archive: data truncated");
    }

    let headers = &archive_data[offset..offset + headers_len];
    offset += headers_len;
    let sequences = &archive_data[offset..offset + sequences_len];
    offset += sequences_len;
    let nmasks = &archive_data[offset..offset + nmasks_len];
    offset += nmasks_len;
    let qualities = &archive_data[offset..offset + qualities_len];

    // Decompress based on encoding mode
    info!("Decompressing...");
    let records = if arithmetic_enabled {
        info!("Decompressing arithmetic-coded data...");

        // Decompress headers
        let read_ids = match header_compressor {
            HeaderCompressor::Bsc => decompress_headers_bsc(headers, num_reads)?,
            HeaderCompressor::Zstd => {
                let header_data = decompress_zstd(headers)?;
                read_id::decode_read_ids(&header_data, &read_id_template, num_reads)?
            }
        };

        // Decompress sequences with arithmetic coding
        let read_lengths = read_lengths_opt.ok_or_else(|| anyhow::anyhow!("Arithmetic mode requires read lengths"))?;
        let decoded_sequences = arithmetic_sequence::decode_sequences_arithmetic(sequences, &read_lengths, num_reads)?;

        // Decompress qualities with arithmetic coding
        let decoded_qualities = if qualities_len > 0 {
            let read_length = if !read_lengths.is_empty() { read_lengths[0] } else { 0 };
            arithmetic_quality::decode_qualities_arithmetic(qualities, &decoded_sequences, read_length, num_reads)?
        } else {
            vec![String::new(); num_reads]
        };

        // Reconstruct records
        let mut records = Vec::new();
        for i in 0..num_reads {
            let id = read_ids.get(i).cloned().unwrap_or_else(|| format!("@UNKNOWN_{}", i));
            let seq = decoded_sequences.get(i).cloned().unwrap_or_default();
            let qual = if qualities_len > 0 {
                Some(decoded_qualities.get(i).cloned().unwrap_or_default())
            } else {
                None
            };
            records.push(crate::io::FastqRecord::new(id, seq, qual));
        }

        records
    } else if encoding_type == 3 {
        // De Bruijn graph mode
        info!("Decompressing de Bruijn graph-coded sequences...");

        // Decompress headers
        let read_ids = match header_compressor {
            HeaderCompressor::Bsc => decompress_headers_bsc(headers, num_reads)?,
            HeaderCompressor::Zstd => {
                let header_data = decompress_zstd(headers)?;
                read_id::decode_read_ids(&header_data, &read_id_template, num_reads)?
            }
        };

        // Decompress sequences with de Bruijn graph
        let decoded_sequences = debruijn::decompress_sequences_debruijn(sequences, num_reads)?;

        // Decompress qualities
        let qualities_data = if qualities_len == 0 {
            Vec::new()
        } else if let Some(ref dict) = quality_dict_opt {
            zstd_dict::decompress_with_dict(qualities, dict)?
        } else {
            decompress_qualities_data(qualities, quality_compressor)?
        };

        // Reconstruct records
        let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
            crate::io::FastqRecord::new(id, seq, None)
        }).collect();

        // Decode quality strings
        if !qualities_data.is_empty() {
            if quality_delta_enabled {
                let mut encoded_deltas = Vec::new();
                let mut qual_offset = 0;
                while qual_offset < qualities_data.len() {
                    let q_len = read_varint(&qualities_data, &mut qual_offset).ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                    if qual_offset + q_len > qualities_data.len() {
                        break;
                    }
                    let deltas = quality_delta::unpack_deltas(&qualities_data[qual_offset..qual_offset + q_len]);
                    qual_offset += q_len;
                    encoded_deltas.push(deltas);
                }
                let decoded_qualities = quality_delta::decode_quality_deltas(&encoded_deltas)?;
                for (i, quality) in decoded_qualities.into_iter().enumerate() {
                    if i < records.len() {
                        records[i].quality = Some(quality);
                    }
                }
            } else {
                let mut qual_offset = 0;
                for record in &mut records {
                    if qual_offset < qualities_data.len() {
                        let q_len = read_varint(&qualities_data, &mut qual_offset)
                            .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                        let quality_str = if let Some(ref model) = quality_model_opt {
                            if qual_offset + q_len <= qualities_data.len() {
                                let deltas = quality_model::unpack_deltas(&qualities_data[qual_offset..qual_offset + q_len]);
                                qual_offset += q_len;
                                quality_model::decode_with_model(&deltas, model)
                            } else {
                                break;
                            }
                        } else {
                            let bits_per_qual = quality_binning.bits_per_quality();
                            let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                            if qual_offset + q_encoded_len <= qualities_data.len() {
                                let quality_str = columnar::unpack_qualities(&qualities_data[qual_offset..qual_offset + q_encoded_len], q_len, quality_binning);
                                qual_offset += q_encoded_len;
                                quality_str
                            } else {
                                break;
                            }
                        };
                        record.quality = Some(quality_str);
                    }
                }
            }
        }

        records
    } else if encoding_type == 1 || encoding_type == 2 {
        // Delta or RLE mode: decompress all streams in parallel
        let ((header_result, seq_result), qual_result) = rayon::join(
            || rayon::join(
                || -> Result<Vec<String>> {
                    match header_compressor {
                        HeaderCompressor::Bsc => decompress_headers_bsc(headers, num_reads),
                        HeaderCompressor::Zstd => {
                            let header_data = decompress_zstd(headers)?;
                            read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                        }
                    }
                },
                || match sequence_compressor {
                    SequenceCompressor::Bsc => bsc::decompress_parallel(sequences),
                    SequenceCompressor::Zstd => decompress_zstd(sequences),
                },
            ),
            || if let Some(ref dict) = quality_dict_opt {
                zstd_dict::decompress_with_dict(qualities, dict)
            } else {
                decompress_qualities_data(qualities, quality_compressor)
            },
        );
        let read_ids = header_result?;
        let seq_stream = seq_result?;
        let quality_stream = qual_result?;

        // Parse encoded sequences
        let mut encoded_sequences = Vec::new();
        let mut s_offset = 0;
        while s_offset + 4 <= seq_stream.len() {
            let len = u32::from_le_bytes(seq_stream[s_offset..s_offset + 4].try_into().unwrap()) as usize;
            s_offset += 4;
            if s_offset + len > seq_stream.len() {
                break;
            }
            encoded_sequences.push(seq_stream[s_offset..s_offset + len].to_vec());
            s_offset += len;
        }

        // Decode sequences
        let decoded_sequences = if encoding_type == 1 {
            info!("Decoding delta-encoded sequences...");
            delta::decode_delta_encoding(&encoded_sequences)?
        } else {
            info!("Decoding RLE-encoded sequences...");
            rle::decode_rle_from_sequences(&encoded_sequences)?
        };

        // Reconstruct records
        let mut records = Vec::new();
        let mut _q_offset = 0;

        for (i, seq) in decoded_sequences.iter().enumerate() {
            let id = read_ids.get(i).cloned().unwrap_or_else(|| format!("@UNKNOWN_{}", i));

            records.push(crate::io::FastqRecord::new(id, seq.clone(), None));
        }

        // Decode qualities based on encoding mode
        if quality_delta_enabled {
            // Quality delta encoding: decode all at once
            let mut encoded_deltas = Vec::new();
            let mut q_offset = 0;
            while q_offset < quality_stream.len() {
                let q_len = read_varint(&quality_stream, &mut q_offset).ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                if q_offset + q_len > quality_stream.len() {
                    break;
                }
                let deltas = quality_delta::unpack_deltas(&quality_stream[q_offset..q_offset + q_len]);
                q_offset += q_len;
                encoded_deltas.push(deltas);
            }
            let decoded_qualities = quality_delta::decode_quality_deltas(&encoded_deltas)?;
            for (i, quality) in decoded_qualities.into_iter().enumerate() {
                if i < records.len() {
                    records[i].quality = Some(quality);
                }
            }
        } else {
            // Standard or model-based quality decoding (per-read)
            let mut q_offset = 0;
            for record in &mut records {
                if q_offset < quality_stream.len() {
                    let q_len = read_varint(&quality_stream, &mut q_offset).ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;

                    let quality_str = if let Some(ref model) = quality_model_opt {
                        // Quality modeling mode: decode delta values
                        let deltas = quality_model::unpack_deltas(&quality_stream[q_offset..q_offset + q_len]);
                        q_offset += q_len;
                        quality_model::decode_with_model(&deltas, model)
                    } else {
                        // Standard binning mode
                        let bits_per_qual = quality_binning.bits_per_quality();
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        let quality_str = columnar::unpack_qualities(&quality_stream[q_offset..q_offset + q_encoded_len], q_len, quality_binning);
                        q_offset += q_encoded_len;
                        quality_str
                    };
                    record.quality = Some(quality_str);
                }
            }
        }

        records
    } else {
        // Normal mode: decompress all streams in parallel
        info!("Decompressing streams in parallel...");

        // Step 1: Parallel decompression of headers, sequences, qualities
        let (header_result, (seq_result, qual_result)) = rayon::join(
            || -> Result<Vec<String>> {
                match header_compressor {
                    HeaderCompressor::Bsc => decompress_headers_bsc(headers, num_reads),
                    HeaderCompressor::Zstd => {
                        let header_data = decompress_zstd(headers)?;
                        read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                    }
                }
            },
            || rayon::join(
                || -> Result<Vec<String>> {
                    match sequence_compressor {
                        SequenceCompressor::Bsc => {
                            decompress_sequences_raw_bsc(sequences, num_reads)
                        }
                        SequenceCompressor::Zstd => {
                            let sequences_data = decompress_zstd(sequences)?;
                            let nmasks_data = decompress_zstd(nmasks)?;
                            let mut decoded = Vec::with_capacity(num_reads);
                            let mut seq_offset = 0;
                            let mut nmask_offset = 0;
                            for _ in 0..num_reads {
                                let seq_len = read_varint(&sequences_data, &mut seq_offset)
                                    .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length"))?;
                                let seq_2bit_len = (seq_len + 3) / 4;
                                if seq_offset + seq_2bit_len > sequences_data.len() {
                                    anyhow::bail!("Truncated sequence data");
                                }
                                let sequence_2bit = &sequences_data[seq_offset..seq_offset + seq_2bit_len];
                                seq_offset += seq_2bit_len;
                                let nmask_len = (seq_len + 7) / 8;
                                let n_mask = if nmask_offset + nmask_len <= nmasks_data.len() {
                                    let mask = &nmasks_data[nmask_offset..nmask_offset + nmask_len];
                                    nmask_offset += nmask_len;
                                    mask
                                } else {
                                    &[]
                                };
                                let encoding = n_mask::NMaskEncoding {
                                    sequence_2bit: sequence_2bit.to_vec(),
                                    n_mask: n_mask.to_vec(),
                                    length: seq_len,
                                };
                                decoded.push(n_mask::decode_with_n_mask(&encoding));
                            }
                            Ok(decoded)
                        }
                    }
                },
                || -> Result<Vec<u8>> {
                    if qualities_len == 0 {
                        Ok(Vec::new())
                    } else if let Some(ref dict) = quality_dict_opt {
                        zstd_dict::decompress_with_dict(qualities, dict)
                    } else {
                        decompress_qualities_data(qualities, quality_compressor)
                    }
                },
            ),
        );

        let read_ids = header_result?;
        let decoded_sequences = seq_result?;
        let qualities_data = qual_result?;

        // Step 2: Construct records from decoded headers and sequences
        let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
            crate::io::FastqRecord::new(id, seq, None)
        }).collect();

        // Step 3: Decode quality strings and attach to records
        if !qualities_data.is_empty() {
            if quality_delta_enabled {
                let mut encoded_deltas = Vec::new();
                let mut qual_offset = 0;
                while qual_offset < qualities_data.len() {
                    let q_len = read_varint(&qualities_data, &mut qual_offset).ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                    if qual_offset + q_len > qualities_data.len() {
                        break;
                    }
                    let deltas = quality_delta::unpack_deltas(&qualities_data[qual_offset..qual_offset + q_len]);
                    qual_offset += q_len;
                    encoded_deltas.push(deltas);
                }
                let decoded_qualities = quality_delta::decode_quality_deltas(&encoded_deltas)?;
                for (i, quality) in decoded_qualities.into_iter().enumerate() {
                    if i < records.len() {
                        records[i].quality = Some(quality);
                    }
                }
            } else {
                let mut qual_offset = 0;
                for record in &mut records {
                    if qual_offset < qualities_data.len() {
                        let q_len = read_varint(&qualities_data, &mut qual_offset)
                            .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;

                        let quality_str = if let Some(ref model) = quality_model_opt {
                            if qual_offset + q_len <= qualities_data.len() {
                                let deltas = quality_model::unpack_deltas(&qualities_data[qual_offset..qual_offset + q_len]);
                                qual_offset += q_len;
                                quality_model::decode_with_model(&deltas, model)
                            } else {
                                break;
                            }
                        } else {
                            let bits_per_qual = quality_binning.bits_per_quality();
                            let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                            if qual_offset + q_encoded_len <= qualities_data.len() {
                                let quality_str = columnar::unpack_qualities(&qualities_data[qual_offset..qual_offset + q_encoded_len], q_len, quality_binning);
                                qual_offset += q_encoded_len;
                                quality_str
                            } else {
                                break;
                            }
                        };
                        record.quality = Some(quality_str);
                    }
                }
            }
        }

        records
    };

    info!("Decompressed {} records", records.len());

    // Write output
    info!("Writing output file...");
    if args.output.is_empty() {
        anyhow::bail!("No output file specified");
    }

    let output_path = &args.output[0];
    let mut output = if args.gzipped {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        Box::new(std::io::BufWriter::with_capacity(
            8 * 1024 * 1024,
            GzEncoder::new(
                std::fs::File::create(output_path)?,
                Compression::new(args.gzip_level),
            ),
        )) as Box<dyn Write>
    } else {
        Box::new(std::io::BufWriter::with_capacity(
            8 * 1024 * 1024,
            std::fs::File::create(output_path)?,
        )) as Box<dyn Write>
    };

    for record in &records {
        output.write_all(record.id.as_bytes())?;
        output.write_all(b"\n")?;
        output.write_all(record.sequence.as_bytes())?;
        output.write_all(b"\n+\n")?;
        if let Some(qual) = &record.quality {
            output.write_all(qual.as_bytes())?;
        }
        output.write_all(b"\n")?;
    }

    let elapsed = start_time.elapsed();
    info!("Decompression completed in {:.2}s", elapsed.as_secs_f64());

    Ok(())
}

/// Helper: Compress headers stream with template encoding
pub(crate) fn compress_headers(records: &[crate::io::FastqRecord], level: i32) -> Result<(Vec<u8>, read_id::ReadIdTemplate)> {
    // Extract read IDs
    let read_ids: Vec<String> = records.iter().map(|r| r.id.clone()).collect();

    // Apply read ID template compression
    let encoded = read_id::compress_read_ids(&read_ids)?;

    // Compress the encoded data with zstd
    let compressed = compress_zstd(&encoded.encoded_data, level)?;

    Ok((compressed, encoded.template))
}

/// Helper: Compress headers as raw ASCII + BSC (no template encoding)
/// Benchmark showed Raw + BSC (6.90x) dramatically beats Template + Zstd (3.64x)
pub(crate) fn compress_headers_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool) -> Result<Vec<u8>> {
    let mut header_stream = Vec::new();
    for record in records {
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.extend_from_slice(record.id.as_bytes());
    }
    bsc_compress_parallel(&header_stream, bsc_static)
}

pub(crate) fn compress_headers_bsc(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    compress_headers_bsc_with(records, false)
}

/// Helper: Decompress raw BSC-compressed headers
fn decompress_headers_bsc(compressed: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let decompressed = bsc::decompress_parallel(compressed)?;
    let mut headers = Vec::with_capacity(num_reads);
    let mut offset = 0;

    for _ in 0..num_reads {
        let hdr_len = read_varint(&decompressed, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Failed to read header length varint"))?;
        if offset + hdr_len > decompressed.len() {
            anyhow::bail!("Truncated header data at offset {}", offset);
        }
        let header = String::from_utf8(decompressed[offset..offset + hdr_len].to_vec())
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in header: {}", e))?;
        offset += hdr_len;
        headers.push(header);
    }

    Ok(headers)
}

/// Helper: Compress sequences as raw ASCII + BSC (no 2-bit encoding, no N-mask)
pub(crate) fn compress_sequences_raw_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut seq_stream = Vec::new();
    for record in records {
        write_varint(&mut seq_stream, record.sequence.len())?;
        seq_stream.extend_from_slice(record.sequence.as_bytes());
    }

    let compressed = bsc_compress_parallel(&seq_stream, bsc_static)?;
    Ok((compressed, Vec::new())) // empty nmasks
}

pub(crate) fn compress_sequences_raw_bsc(records: &[crate::io::FastqRecord]) -> Result<(Vec<u8>, Vec<u8>)> {
    compress_sequences_raw_bsc_with(records, false)
}

/// Helper: Decompress raw ASCII BSC-compressed sequences
fn decompress_sequences_raw_bsc(compressed: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let decompressed = bsc::decompress_parallel(compressed)?;
    let mut sequences = Vec::with_capacity(num_reads);
    let mut offset = 0;

    for _ in 0..num_reads {
        let seq_len = read_varint(&decompressed, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length varint"))?;
        if offset + seq_len > decompressed.len() {
            anyhow::bail!("Truncated sequence data at offset {}", offset);
        }
        let sequence = String::from_utf8(decompressed[offset..offset + seq_len].to_vec())
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?;
        offset += seq_len;
        sequences.push(sequence);
    }

    Ok(sequences)
}

/// Helper: Compress qualities stream (standard mode)
pub(crate) fn compress_qualities_with(records: &[crate::io::FastqRecord], binning: columnar::QualityBinning, level: i32, compressor: QualityCompressor, bsc_static: bool) -> Result<Vec<u8>> {
    use std::io::Write;

    let mut quality_stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut quality_stream, qual.len())?;
            let packed = columnar::pack_qualities(qual, binning);
            quality_stream.write_all(&packed)?;
        }
    }

    match compressor {
        QualityCompressor::Zstd => compress_zstd(&quality_stream, level),
        QualityCompressor::Bsc => bsc_compress_parallel(&quality_stream, bsc_static),
    }
}

pub(crate) fn compress_qualities(records: &[crate::io::FastqRecord], binning: columnar::QualityBinning, level: i32, compressor: QualityCompressor) -> Result<Vec<u8>> {
    compress_qualities_with(records, binning, level, compressor, false)
}

/// Helper: Compress qualities stream with dictionary
fn compress_qualities_with_dict(
    records: &[crate::io::FastqRecord],
    binning: columnar::QualityBinning,
    dictionary: &[u8],
    level: i32,
    compressor: QualityCompressor,
) -> Result<Vec<u8>> {
    use std::io::Write;

    let mut quality_stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut quality_stream, qual.len())?;
            let packed = columnar::pack_qualities(qual, binning);
            quality_stream.write_all(&packed)?;
        }
    }

    match compressor {
        QualityCompressor::Zstd => zstd_dict::compress_with_dict(&quality_stream, dictionary, level),
        QualityCompressor::Bsc => {
            // BSC doesn't use dictionary - just use regular BSC compression
            info!("Note: BSC doesn't support dictionary mode, using standard BSC");
            bsc::compress_parallel_adaptive(&quality_stream) // dict path always adaptive
        }
    }
}

/// Helper: Compress qualities stream with positional modeling
fn compress_qualities_with_model(
    records: &[crate::io::FastqRecord],
    level: i32,
) -> Result<(Vec<u8>, quality_model::QualityModel)> {
    use std::io::Write;

    // Extract quality strings
    let quality_strings: Vec<String> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores to model");
    }

    // Build quality model
    let model = quality_model::build_quality_model(&quality_strings)?;
    info!("Built quality model: {} positions, median quality {:.1}",
        model.positional_medians.len(),
        model.positional_medians.iter().map(|&q| q as f64).sum::<f64>() / model.positional_medians.len() as f64
    );

    // Encode qualities using model
    let mut quality_stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut quality_stream, qual.len())?;
            let deltas = quality_model::encode_with_model(qual, &model);
            let packed = quality_model::pack_deltas(&deltas);
            quality_stream.write_all(&packed)?;
        }
    }

    // Compress delta stream with zstd
    let compressed = compress_zstd(&quality_stream, level)?;

    Ok((compressed, model))
}

/// Helper: Compress qualities stream with delta encoding between adjacent reads
fn compress_qualities_with_delta(records: &[crate::io::FastqRecord], level: i32) -> Result<Vec<u8>> {
    use std::io::Write;

    // Extract quality strings
    let quality_strings: Vec<String> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores for delta encoding");
    }

    // Encode qualities as deltas from previous read
    let encoded_deltas = quality_delta::encode_quality_deltas(&quality_strings);

    // Analyze delta statistics for logging
    let stats = quality_delta::analyze_deltas(&encoded_deltas);
    info!("Quality delta encoding: {:.1}% zeros, {:.1}% small (|Δ|≤3), max delta: {}",
        stats.zero_percent, stats.small_percent, stats.max_delta);

    // Pack and compress delta stream with zstd
    let mut quality_stream = Vec::new();
    for deltas in &encoded_deltas {
        write_varint(&mut quality_stream, deltas.len())?;
        let packed = quality_delta::pack_deltas(deltas);
        quality_stream.write_all(&packed)?;
    }

    compress_zstd(&quality_stream, level)
}

/// Write variable-length integer
pub(crate) fn write_varint<W: std::io::Write>(writer: &mut W, mut value: usize) -> std::io::Result<()> {
    while value >= 0x80 {
        writer.write_all(&[((value & 0x7F) | 0x80) as u8])?;
        value >>= 7;
    }
    writer.write_all(&[value as u8])
}

/// Read variable-length integer
fn read_varint(data: &[u8], offset: &mut usize) -> Option<usize> {
    let mut value = 0usize;
    let mut shift = 0;

    loop {
        if *offset >= data.len() {
            return None;
        }

        let byte = data[*offset];
        *offset += 1;

        value |= ((byte & 0x7F) as usize) << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
    }
}

/// Compress data with Zstd (alternative to gzip)
fn compress_zstd(data: &[u8], level: i32) -> Result<Vec<u8>> {
    zstd::bulk::compress(data, level)
        .map_err(|e| anyhow::anyhow!("Zstd compression failed: {}", e))
}

/// Decompress data with Zstd
fn decompress_zstd(compressed: &[u8]) -> Result<Vec<u8>> {
    zstd::bulk::decompress(compressed, 100_000_000) // 100MB limit
        .map_err(|e| anyhow::anyhow!("Zstd decompression failed: {}", e))
}

fn decompress_qualities_data(compressed: &[u8], compressor: QualityCompressor) -> Result<Vec<u8>> {
    match compressor {
        QualityCompressor::Zstd => decompress_zstd(compressed),
        QualityCompressor::Bsc => bsc::decompress_parallel(compressed),
    }
}

/// Compress sequences with arithmetic coding
fn compress_sequences_arithmetic(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    let sequences: Vec<String> = records
        .iter()
        .map(|r| r.sequence.clone())
        .collect();

    arithmetic_sequence::encode_sequences_arithmetic(&sequences)
}

/// Compress qualities with arithmetic coding
fn compress_qualities_arithmetic(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    let sequences: Vec<String> = records
        .iter()
        .map(|r| r.sequence.clone())
        .collect();

    let qualities: Vec<String> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if qualities.is_empty() {
        anyhow::bail!("No quality scores for arithmetic coding");
    }

    // Get read length (assume all reads same length for now)
    let read_length = if !records.is_empty() {
        records[0].sequence.len()
    } else {
        anyhow::bail!("No records to compress");
    };

    arithmetic_quality::encode_qualities_arithmetic(&qualities, &sequences, read_length)
}

/// Helper: Compress qualities with model and dictionary
fn compress_qualities_with_model_and_dict(
    records: &[crate::io::FastqRecord],
    dictionary: &[u8],
    level: i32,
) -> Result<(Vec<u8>, quality_model::QualityModel)> {
    use std::io::Write;

    // Extract quality strings
    let quality_strings: Vec<String> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores to model");
    }

    // Build quality model
    let model = quality_model::build_quality_model(&quality_strings)?;
    info!("Built quality model: {} positions, median quality {:.1}",
        model.positional_medians.len(),
        model.positional_medians.iter().map(|&q| q as f64).sum::<f64>() / model.positional_medians.len() as f64
    );

    // Encode qualities using model
    let mut quality_stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut quality_stream, qual.len())?;
            let deltas = quality_model::encode_with_model(qual, &model);
            let packed = quality_model::pack_deltas(&deltas);
            quality_stream.write_all(&packed)?;
        }
    }

    // Compress delta stream with dictionary
    let compressed = zstd_dict::compress_with_dict(&quality_stream, dictionary, level)?;

    Ok((compressed, model))
}

/// Helper: Compress qualities with delta and dictionary
fn compress_qualities_with_delta_and_dict(
    records: &[crate::io::FastqRecord],
    dictionary: &[u8],
    level: i32,
) -> Result<Vec<u8>> {
    use std::io::Write;

    // Extract quality strings
    let quality_strings: Vec<String> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores for delta encoding");
    }

    // Encode qualities as deltas from previous read
    let encoded_deltas = quality_delta::encode_quality_deltas(&quality_strings);

    // Analyze delta statistics for logging
    let stats = quality_delta::analyze_deltas(&encoded_deltas);
    info!("Quality delta encoding: {:.1}% zeros, {:.1}% small (|Δ|≤3), max delta: {}",
        stats.zero_percent, stats.small_percent, stats.max_delta);

    // Pack and compress delta stream with dictionary
    let mut quality_stream = Vec::new();
    for deltas in &encoded_deltas {
        write_varint(&mut quality_stream, deltas.len())?;
        let packed = quality_delta::pack_deltas(deltas);
        quality_stream.write_all(&packed)?;
    }

    zstd_dict::compress_with_dict(&quality_stream, dictionary, level)
}


/// Compress paired-end FASTQ files with R2 correlation compression
pub fn compress_paired_end(
    r1_path: &std::path::Path,
    r2_path: &std::path::Path,
    args: &CompressArgs,
) -> Result<()> {
    use crate::io::FastqReader;
    use std::io::Write;

    let start_time = Instant::now();

    // Set up thread pool
    let num_threads = if args.threads == 0 {
        num_cpus::get()
    } else {
        args.threads
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok();

    info!("Using {} threads for compression", num_threads);
    info!("Compression level: {}", args.compression_level);
    info!("Paired-end mode: R1 = {:?}, R2 = {:?}", r1_path, r2_path);
    info!("Output: {:?}", args.output);

    // Read R1 records
    info!("Reading R1 file...");
    let mut r1_reader = FastqReader::from_path(r1_path, false)?;
    let mut r1_records = Vec::new();
    while let Some(record) = r1_reader.next()? {
        r1_records.push(record);
    }
    info!("Read {} R1 records", r1_records.len());

    // Read R2 records
    info!("Reading R2 file...");
    let mut r2_reader = FastqReader::from_path(r2_path, false)?;
    let mut r2_records = Vec::new();
    while let Some(record) = r2_reader.next()? {
        r2_records.push(record);
    }
    info!("Read {} R2 records", r2_records.len());

    if r1_records.len() != r2_records.len() {
        anyhow::bail!("R1 and R2 have different number of reads: {} vs {}", 
            r1_records.len(), r2_records.len());
    }

    // Apply quality quantization to both R1 and R2
    let quality_binning = if args.no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(args.quality_mode)
    };

    if !args.no_quality {
        for record in &mut r1_records {
            if let Some(qual) = &record.quality {
                record.quality = Some(quality::quantize_quality(qual, args.quality_mode));
            }
        }
        for record in &mut r2_records {
            if let Some(qual) = &record.quality {
                record.quality = Some(quality::quantize_quality(qual, args.quality_mode));
            }
        }
    }

    // Calculate paired-end correlation statistics
    let r1_sequences: Vec<String> = r1_records.iter().map(|r| r.sequence.clone()).collect();
    let r2_sequences: Vec<String> = r2_records.iter().map(|r| r.sequence.clone()).collect();
    
    let diff_rate = paired_end::calculate_difference_rate(&r1_sequences, &r2_sequences);
    info!("Paired-end correlation: {:.2}% difference rate", diff_rate);

    // Encode R2 as differences from R1
    info!("Encoding R2 with paired-end correlation...");
    let (r2_diff_counts, r2_diff_data) = paired_end::encode_paired_differences(&r1_sequences, &r2_sequences)?;

    // Compress R1 normally (sequences + qualities + headers)
    info!("Compressing R1 data...");
    let (r1_headers, r1_template) = compress_headers(&r1_records, args.compression_level)?;
    let (r1_sequences_enc, r1_nmasks_enc) = compress_seqs_nmask(&r1_records, args.compression_level)?;
    let r1_qualities = if !args.no_quality {
        Some(compress_qualities(&r1_records, columnar::QualityBinning::None, args.compression_level, args.quality_compressor)?)
    } else {
        None
    };

    // Compress R2 data (differences + qualities + headers)
    info!("Compressing R2 data...");
    let (r2_headers, r2_template) = compress_headers(&r2_records, args.compression_level)?;
    let r2_diff_counts_compressed = compress_zstd(&bincode::serialize(&r2_diff_counts)?, args.compression_level)?;
    let r2_diff_data_compressed = compress_zstd(&r2_diff_data, args.compression_level)?;
    let r2_qualities = if !args.no_quality {
        Some(compress_qualities(&r2_records, columnar::QualityBinning::None, args.compression_level, args.quality_compressor)?)
    } else {
        None
    };

    // Write paired-end archive format
    info!("Writing paired-end archive...");
    let mut output_file = std::fs::File::create(&args.output)?;

    // Magic + version + flags
    output_file.write_all(b"FQZ\0")?;  // Magic
    output_file.write_all(&[2u8])?;   // Version 2 = paired-end
    output_file.write_all(&[0u8])?;   // Encoding type
    output_file.write_all(&[binning_to_code(quality_binning)])?;  // Quality binning
    output_file.write_all(&[compressor_to_code(args.quality_compressor)])?;  // Quality compressor

    // R1 template
    write_paired_end_template(&mut output_file, &r1_template)?;

    // R2 template
    write_paired_end_template(&mut output_file, &r2_template)?;

    // Number of reads
    output_file.write_all(&(r1_records.len() as u64).to_le_bytes())?;

    // R1 stream lengths
    output_file.write_all(&(r1_headers.len() as u64).to_le_bytes())?;
    output_file.write_all(&(r1_sequences_enc.len() as u64).to_le_bytes())?;
    output_file.write_all(&(r1_nmasks_enc.len() as u64).to_le_bytes())?;
    output_file.write_all(&(r1_qualities.as_ref().map(|q| q.len()).unwrap_or(0) as u64).to_le_bytes())?;

    // R2 stream lengths
    output_file.write_all(&(r2_headers.len() as u64).to_le_bytes())?;
    output_file.write_all(&(r2_diff_counts_compressed.len() as u64).to_le_bytes())?;
    output_file.write_all(&(r2_diff_data_compressed.len() as u64).to_le_bytes())?;
    output_file.write_all(&(r2_qualities.as_ref().map(|q| q.len()).unwrap_or(0) as u64).to_le_bytes())?;

    // Write data streams
    output_file.write_all(&r1_headers)?;
    output_file.write_all(&r1_sequences_enc)?;
    output_file.write_all(&r1_nmasks_enc)?;
    if let Some(ref qual) = r1_qualities {
        output_file.write_all(qual)?;
    }

    output_file.write_all(&r2_headers)?;
    output_file.write_all(&r2_diff_counts_compressed)?;
    output_file.write_all(&r2_diff_data_compressed)?;
    if let Some(ref qual) = r2_qualities {
        output_file.write_all(qual)?;
    }

    let elapsed = start_time.elapsed();
    let output_size = std::fs::metadata(&args.output)?.len();
    let input_size = std::fs::metadata(r1_path)?.len() + std::fs::metadata(r2_path)?.len();
    let ratio = input_size as f64 / output_size as f64;

    info!("Paired-end compression complete!");
    info!("Input size: {:.2} MB (R1 + R2)", input_size as f64 / (1024.0 * 1024.0));
    info!("Output size: {:.2} MB", output_size as f64 / (1024.0 * 1024.0));
    info!("Compression ratio: {:.2}x", ratio);
    info!("Time: {:.2}s", elapsed.as_secs_f64());

    Ok(())
}

// Helper function to write read ID template for paired-end archives
fn write_paired_end_template<W: std::io::Write>(writer: &mut W, template: &read_id::ReadIdTemplate) -> Result<()> {
    let prefix_bytes = template.prefix.as_bytes();
    writer.write_all(&(prefix_bytes.len() as u16).to_le_bytes())?;
    writer.write_all(prefix_bytes)?;
    writer.write_all(&[template.has_comment as u8])?;
    
    if let Some(ref comment) = template.common_comment {
        let comment_bytes = comment.as_bytes();
        writer.write_all(&(comment_bytes.len() as u16).to_le_bytes())?;
        writer.write_all(comment_bytes)?;
    } else {
        writer.write_all(&[0u8, 0u8])?;  // No common comment
    }
    
    Ok(())
}

// Helper function for simple sequence compression with n-mask (for paired-end)
fn compress_seqs_nmask(records: &[crate::io::FastqRecord], level: i32) -> Result<(Vec<u8>, Vec<u8>)> {
    let sequences: Vec<String> = records.iter().map(|r| r.sequence.clone()).collect();
    let encodings: Vec<n_mask::NMaskEncoding> = sequences
        .iter()
        .map(|seq| n_mask::encode_with_n_mask(seq))
        .collect();
    
    // Concatenate all 2-bit sequences
    let mut all_sequences = Vec::new();
    for enc in &encodings {
        all_sequences.extend_from_slice(&enc.sequence_2bit);
    }
    
    // Concatenate all N-masks
    let mut all_nmasks = Vec::new();
    for enc in &encodings {
        all_nmasks.extend_from_slice(&enc.n_mask);
    }
    
    let sequences_compressed = compress_zstd(&all_sequences, level)?;
    let nmasks_compressed = compress_zstd(&all_nmasks, level)?;
    
    Ok((sequences_compressed, nmasks_compressed))
}
