/// Two-pass pattern routing for sequence compression.
///
/// Pass 1: Scan first N reads to learn globally frequent syncmer hash patterns.
/// Pass 2: Route each read's sequence to main or pattern stream based on hash match.
///         The pattern stream gets all reads matching ANY learned pattern — BWT
///         compresses one larger block more efficiently than 100 tiny blocks.
///
/// Archive format (encoding_type=7, version=2):
///   sequences region = [version:1B][num_patterns:1B]
///     [routing_len:8B][routing_data]          // BSC-compressed, 1 byte per read (0=main, 1=pattern)
///     [main_seq_len:8B][main_seq_blocks]      // multi-block BSC
///     [pattern_seq_len:8B][pattern_seq_blocks] // single combined stream, multi-block BSC

use anyhow::{Context, Result};
use rustc_hash::FxHashMap;
use tracing::info;
use std::time::Instant;

use super::*;

const SCAN_LIMIT: usize = 5_000_000;
const MAX_PATTERNS: usize = 100;
const MIN_PATTERN_COUNT: usize = 50;
const CHUNK_SIZE: usize = 5_000_000;
/// Version byte for the new pattern routing format (distinguishes from old factorize v1).
const FORMAT_VERSION: u8 = 2;

// ── Pass 1: Learn patterns ──────────────────────────────────────────────

/// Scan first `scan_limit` reads and return the top frequent syncmer hashes.
fn learn_patterns<R: std::io::BufRead>(
    reader: &mut crate::io::FastqReader<R>,
    scan_limit: usize,
    min_count: usize,
    max_patterns: usize,
) -> Result<Vec<u64>> {
    let mut hash_counts: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(2_000_000, Default::default());
    let mut scanned = 0usize;

    for _ in 0..scan_limit {
        match reader.next()? {
            Some(record) => {
                let hash = dna_utils::compute_min_syncmer_hash(record.sequence.as_bytes());
                if hash != u64::MAX {
                    *hash_counts.entry(hash).or_insert(0) += 1;
                }
                scanned += 1;
            }
            None => break,
        }
    }

    // Sort by frequency descending, take top patterns above threshold
    let mut counts: Vec<(u64, u32)> = hash_counts.into_iter()
        .filter(|&(_, c)| c >= min_count as u32)
        .collect();
    counts.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    counts.truncate(max_patterns);

    let patterns: Vec<u64> = counts.iter().map(|&(h, _)| h).collect();

    info!("Pass 1: scanned {} reads, found {} patterns (min count ≥ {})",
        scanned, patterns.len(), min_count);
    for (i, &(hash, count)) in counts.iter().take(10).enumerate() {
        info!("  pattern {}: hash={:#018x} count={}", i + 1, hash, count);
    }

    Ok(patterns)
}

// ── Pass 2: Route and compress ──────────────────────────────────────────

pub(super) fn compress_factorize(args: &crate::cli::CompressConfig) -> Result<()> {
    use std::io::{Write, BufWriter};

    let start_time = Instant::now();
    let bsc_static = args.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };

    info!("Factorized compression v2: two-pass pattern routing");

    // ── Pass 1: Learn patterns ──────────────────────────────────────────
    let mut reader = crate::io::FastqReader::from_path(&args.input[0], args.fasta)?;
    let patterns = learn_patterns(&mut reader, SCAN_LIMIT, MIN_PATTERN_COUNT, MAX_PATTERNS)?;
    let num_patterns = patterns.len();
    drop(reader);

    // Build routing set: any hash that matched a pattern → route to pattern stream
    let routing_set: rustc_hash::FxHashSet<u64> = patterns.into_iter().collect();

    // ── Pass 2: Route reads to streams ──────────────────────────────────
    info!("Pass 2: routing reads to 1 pattern stream + main stream ({} learned patterns)", num_patterns);

    let working_dir = &args.working_dir;

    // Temp files: headers, main sequences, combined pattern sequences, qualities
    let h_tmp_path = working_dir.join(".qz_fact2_h.tmp");
    let ms_tmp_path = working_dir.join(".qz_fact2_ms.tmp");
    let ps_tmp_path = working_dir.join(".qz_fact2_ps.tmp");
    let q_tmp_path = working_dir.join(".qz_fact2_q.tmp");

    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 { let _ = std::fs::remove_file(p); }
        }
    }
    let _cleanup = TmpCleanup(vec![
        h_tmp_path.clone(), ms_tmp_path.clone(), ps_tmp_path.clone(), q_tmp_path.clone(),
    ]);

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut ms_tmp = BufWriter::new(std::fs::File::create(&ms_tmp_path)?);
    let mut ps_tmp = BufWriter::new(std::fs::File::create(&ps_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);

    let mut h_num_blocks: u32 = 0;
    let mut ms_num_blocks: u32 = 0;
    let mut ps_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;
    let mut num_reads: usize = 0;
    let mut original_size: usize = 0;
    let mut routing: Vec<u8> = Vec::new();
    let mut pattern_count: usize = 0;
    let mut main_count: usize = 0;

    let mut reader = crate::io::FastqReader::from_path(&args.input[0], args.fasta)?;
    let (mut cur_records, _, mut cur_orig) = read_chunk_records(&mut reader, CHUNK_SIZE)?;
    let mut chunk_idx = 0usize;

    while !cur_records.is_empty() {
        let cur_reads = cur_records.len();
        info!("Chunk {}: {} reads", chunk_idx, cur_reads);

        let (next_result, compress_result) = std::thread::scope(|scope| {
            let records = std::mem::take(&mut cur_records);
            let routing_set_ref = &routing_set;

            let compress_handle = scope.spawn(move || -> Result<ChunkResult> {
                compress_chunk_v2(
                    records, quality_mode, quality_binning, no_quality, bsc_static,
                    routing_set_ref,
                )
            });

            let next = read_chunk_records(&mut reader, CHUNK_SIZE);
            let compressed = compress_handle.join().unwrap();
            (next, compressed)
        });

        let chunk = compress_result?;

        h_num_blocks += write_blocks_to_tmp(chunk.h_blocks, &mut h_tmp)?;
        ms_num_blocks += write_blocks_to_tmp(chunk.ms_blocks, &mut ms_tmp)?;
        if !chunk.ps_blocks.is_empty() {
            ps_num_blocks += write_blocks_to_tmp(chunk.ps_blocks, &mut ps_tmp)?;
        }
        if !chunk.q_blocks.is_empty() {
            q_num_blocks += write_blocks_to_tmp(chunk.q_blocks, &mut q_tmp)?;
        }

        routing.extend_from_slice(&chunk.routing);
        main_count += chunk.main_count;
        pattern_count += chunk.pattern_count;

        num_reads += cur_reads;
        original_size += cur_orig;
        chunk_idx += 1;

        let (nr, _, no) = next_result?;
        cur_records = nr;
        cur_orig = no;
    }

    h_tmp.flush()?; ms_tmp.flush()?; ps_tmp.flush()?; q_tmp.flush()?;
    drop(h_tmp); drop(ms_tmp); drop(ps_tmp); drop(q_tmp);

    info!("Read {} records in {} chunks", num_reads, chunk_idx);
    info!("  Main stream: {} reads ({:.1}%)", main_count, 100.0 * main_count as f64 / num_reads.max(1) as f64);
    info!("  Pattern stream: {} reads ({:.1}%)", pattern_count, 100.0 * pattern_count as f64 / num_reads.max(1) as f64);

    // Compress routing metadata
    info!("Compressing routing metadata ({} bytes)...", routing.len());
    let routing_compressed = bsc::compress_parallel(&routing)?;
    drop(routing);
    info!("  Routing compressed: {} bytes", routing_compressed.len());

    // Get temp file sizes
    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let ms_tmp_size = std::fs::metadata(&ms_tmp_path)?.len() as usize;
    let ps_tmp_size = std::fs::metadata(&ps_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    let h_size = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let ms_size = if ms_num_blocks > 0 { 4 + ms_tmp_size } else { 0 };
    let ps_size = if ps_num_blocks > 0 { 4 + ps_tmp_size } else { 0 };
    let q_size = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    // Sequences region: version(1) + num_patterns(1)
    //   + routing_len(8) + routing_data
    //   + main_seq_len(8) + main_seq_blocks
    //   + pattern_seq_len(8) + pattern_seq_blocks
    let routing_compressed_len = routing_compressed.len();
    let seq_region_size = 1 + 1
        + 8 + routing_compressed_len
        + 8 + ms_size
        + 8 + ps_size;

    // ── Write archive ───────────────────────────────────────────────────
    info!("Writing output file...");
    let mut out = BufWriter::new(std::fs::File::create(&args.output)?);

    // Archive header (encoding_type=7)
    out.write_all(&[7u8])?;                                           // encoding_type
    out.write_all(&[0u8])?;                                           // arithmetic_mode
    out.write_all(&[binning_to_code(quality_binning)])?;              // quality_binning
    out.write_all(&[compressor_to_code(QualityCompressor::Bsc)])?;    // quality_compressor
    out.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?; // seq_compressor
    out.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?; // header_compressor
    out.write_all(&[0u8; 3])?;                                       // padding
    out.write_all(&0u16.to_le_bytes())?;                              // template_prefix
    out.write_all(&[0u8])?;                                           // has_comment

    out.write_all(&(num_reads as u64).to_le_bytes())?;
    out.write_all(&(h_size as u64).to_le_bytes())?;
    out.write_all(&(seq_region_size as u64).to_le_bytes())?;
    out.write_all(&0u64.to_le_bytes())?;                              // nmasks_len
    out.write_all(&(q_size as u64).to_le_bytes())?;

    let copy_stream = |num_blocks: u32, tmp_path: &std::path::Path, out: &mut BufWriter<std::fs::File>| -> Result<()> {
        if num_blocks == 0 { return Ok(()); }
        out.write_all(&num_blocks.to_le_bytes())?;
        let mut f = std::fs::File::open(tmp_path)?;
        std::io::copy(&mut f, out)?;
        Ok(())
    };

    // Headers
    copy_stream(h_num_blocks, &h_tmp_path, &mut out)?;

    // Sequences region
    out.write_all(&[FORMAT_VERSION])?;
    out.write_all(&[num_patterns as u8])?;

    out.write_all(&(routing_compressed_len as u64).to_le_bytes())?;
    out.write_all(&routing_compressed)?;
    drop(routing_compressed);

    out.write_all(&(ms_size as u64).to_le_bytes())?;
    copy_stream(ms_num_blocks, &ms_tmp_path, &mut out)?;

    out.write_all(&(ps_size as u64).to_le_bytes())?;
    copy_stream(ps_num_blocks, &ps_tmp_path, &mut out)?;

    // Qualities
    copy_stream(q_num_blocks, &q_tmp_path, &mut out)?;

    out.flush()?;

    let metadata_size = 9 + 2 + 1 + 40;
    let total = h_size + seq_region_size + q_size + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total);
    info!("Compression ratio: {:.2}x", original_size as f64 / total as f64);
    info!("Stream breakdown:");
    info!("  Headers:    {} bytes ({:.1}%)", h_size, 100.0 * h_size as f64 / total as f64);
    info!("  Routing:    {} bytes ({:.1}%)", routing_compressed_len, 100.0 * routing_compressed_len as f64 / total as f64);
    info!("  Main seqs:  {} bytes ({:.1}%)", ms_size, 100.0 * ms_size as f64 / total as f64);
    info!("  Pattern seqs: {} bytes ({:.1}%)", ps_size, 100.0 * ps_size as f64 / total as f64);
    info!("  Qualities:  {} bytes ({:.1}%)", q_size, 100.0 * q_size as f64 / total as f64);

    Ok(())
}

// ── Per-chunk compression ───────────────────────────────────────────────

struct ChunkResult {
    h_blocks: Vec<Vec<u8>>,
    ms_blocks: Vec<Vec<u8>>,
    ps_blocks: Vec<Vec<u8>>,  // combined pattern stream blocks
    q_blocks: Vec<Vec<u8>>,
    routing: Vec<u8>,
    main_count: usize,
    pattern_count: usize,
}

fn compress_chunk_v2(
    records: Vec<crate::io::FastqRecord>,
    quality_mode: crate::cli::QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    bsc_static: bool,
    routing_set: &rustc_hash::FxHashSet<u64>,
) -> Result<ChunkResult> {
    let n = records.len();

    let mut header_stream = Vec::with_capacity(n * 64);
    let mut main_seq_stream = Vec::new();
    let mut pattern_seq_stream = Vec::new();
    let mut qual_stream = Vec::with_capacity(if no_quality { 0 } else { n * 136 });
    let mut routing = Vec::with_capacity(n);
    let mut main_count = 0usize;
    let mut pattern_count = 0usize;

    for record in &records {
        // Header
        dna_utils::write_varint(&mut header_stream, record.id.len() as u64);
        header_stream.extend_from_slice(record.id.as_bytes());

        // Quality
        if !no_quality {
            if let Some(ref qual) = record.quality {
                let quantized = quality::quantize_quality(qual, quality_mode);
                dna_utils::write_varint(&mut qual_stream, quantized.len() as u64);
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                qual_stream.extend_from_slice(&packed);
            }
        }

        // Route sequence: 0=main, 1=pattern
        let seq = record.sequence.as_bytes();
        let hash = dna_utils::compute_min_syncmer_hash(seq);
        let is_pattern = routing_set.contains(&hash);
        routing.push(if is_pattern { 1 } else { 0 });

        let target = if is_pattern {
            pattern_count += 1;
            &mut pattern_seq_stream
        } else {
            main_count += 1;
            &mut main_seq_stream
        };
        dna_utils::write_varint(target, seq.len() as u64);
        target.extend_from_slice(seq);
    }

    // BSC-compress streams sequentially to limit memory
    let h_blocks = compress_stream_to_bsc_blocks(&header_stream, bsc_static)?;
    drop(header_stream);

    let ms_blocks = compress_stream_to_bsc_blocks(&main_seq_stream, bsc_static)?;
    drop(main_seq_stream);

    let ps_blocks = if pattern_seq_stream.is_empty() {
        Vec::new()
    } else {
        compress_stream_to_bsc_blocks(&pattern_seq_stream, bsc_static)?
    };
    drop(pattern_seq_stream);

    let q_blocks = if no_quality || qual_stream.is_empty() {
        Vec::new()
    } else {
        compress_stream_to_bsc_blocks(&qual_stream, bsc_static)?
    };
    drop(qual_stream);

    Ok(ChunkResult {
        h_blocks, ms_blocks, ps_blocks, q_blocks,
        routing, main_count, pattern_count,
    })
}

// ── Decoder ─────────────────────────────────────────────────────────────

/// Decode pattern-routed sequences from the sequences region (version=2 format).
///
/// Format: [version:1B][num_patterns:1B][routing_len:8B][routing_data]
///         [main_len:8B][main_blocks][pattern_len:8B][pattern_blocks]
pub(super) fn decode_factorized_sequences_v2(
    seq_region: &[u8],
    num_reads: usize,
) -> Result<Vec<String>> {
    if seq_region.len() < 2 {
        anyhow::bail!("Factorize v2: sequences region too small");
    }

    let version = seq_region[0];
    if version != FORMAT_VERSION {
        anyhow::bail!("Factorize v2: unexpected version byte {}", version);
    }
    let num_patterns = seq_region[1] as usize;
    let mut off = 2usize;

    info!("Factorize v2: {} learned patterns, 2 streams (main + pattern)", num_patterns);

    // Read routing
    if off + 8 > seq_region.len() {
        anyhow::bail!("Factorize v2: truncated routing length");
    }
    let routing_len = u64::from_le_bytes(seq_region[off..off + 8].try_into().unwrap()) as usize;
    off += 8;
    if off + routing_len > seq_region.len() {
        anyhow::bail!("Factorize v2: truncated routing data");
    }
    let routing_compressed = &seq_region[off..off + routing_len];
    off += routing_len;

    // Read main sequence stream
    if off + 8 > seq_region.len() {
        anyhow::bail!("Factorize v2: truncated main seq length");
    }
    let main_len = u64::from_le_bytes(seq_region[off..off + 8].try_into().unwrap()) as usize;
    off += 8;
    if off + main_len > seq_region.len() {
        anyhow::bail!("Factorize v2: truncated main seq data");
    }
    let main_compressed = &seq_region[off..off + main_len];
    off += main_len;

    // Read combined pattern stream
    if off + 8 > seq_region.len() {
        anyhow::bail!("Factorize v2: truncated pattern seq length");
    }
    let pattern_len = u64::from_le_bytes(seq_region[off..off + 8].try_into().unwrap()) as usize;
    off += 8;
    if off + pattern_len > seq_region.len() {
        anyhow::bail!("Factorize v2: truncated pattern seq data");
    }
    let pattern_compressed = &seq_region[off..off + pattern_len];

    // Decompress routing
    let routing_data = bsc::decompress_parallel(routing_compressed)
        .context("Failed to decompress routing")?;

    if routing_data.len() != num_reads {
        anyhow::bail!("Factorize v2: routing length {} != num_reads {}",
            routing_data.len(), num_reads);
    }

    // Decompress main + pattern streams in parallel
    let (main_result, pattern_result) = rayon::join(
        || if main_compressed.is_empty() { Ok(Vec::new()) }
           else { bsc::decompress_parallel(main_compressed)
               .context("Failed to decompress main sequence stream") },
        || if pattern_compressed.is_empty() { Ok(Vec::new()) }
           else { bsc::decompress_parallel(pattern_compressed)
               .context("Failed to decompress pattern sequence stream") },
    );
    let main_data = main_result?;
    let pattern_data = pattern_result?;

    // Reconstruct sequences: walk routing, pull from main or pattern cursor
    let pattern_reads: usize = routing_data.iter().filter(|&&r| r != 0).count();
    info!("Reconstructing {} sequences ({} main + {} pattern)...",
        num_reads, num_reads - pattern_reads, pattern_reads);

    let mut main_off = 0usize;
    let mut pattern_off = 0usize;
    let mut sequences = Vec::with_capacity(num_reads);

    for i in 0..num_reads {
        let (data, cursor) = if routing_data[i] == 0 {
            (&main_data, &mut main_off)
        } else {
            (&pattern_data, &mut pattern_off)
        };

        let seq_len = dna_utils::read_varint(data, cursor)
            .ok_or_else(|| anyhow::anyhow!("Factorize v2: truncated seq_len for read {}", i))?
            as usize;

        if *cursor + seq_len > data.len() {
            anyhow::bail!("Factorize v2: truncated sequence data for read {} (need {}, have {})",
                i, seq_len, data.len() - *cursor);
        }

        let seq = String::from_utf8(data[*cursor..*cursor + seq_len].to_vec())
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence {}: {}", i, e))?;
        *cursor += seq_len;
        sequences.push(seq);
    }

    Ok(sequences)
}

/// Legacy decoder for old 3-sub-stream format (factorize v1).
pub(super) fn decode_factorized_sequences(
    meta_data: &[u8],
    raw_seq_data: &[u8],
    delta_seq_data: &[u8],
    num_reads: usize,
) -> Result<Vec<String>> {
    let mut sequences: Vec<Vec<u8>> = Vec::with_capacity(num_reads);
    let mut meta_off = 0usize;
    let mut raw_off = 0usize;
    let mut delta_off = 0usize;
    let mut delta_count = 0usize;

    for i in 0..num_reads {
        if meta_off >= meta_data.len() {
            anyhow::bail!("Truncated meta stream at read {}", i);
        }
        let class = meta_data[meta_off];
        meta_off += 1;

        if class == 0 {
            let seq_len = dna_utils::read_varint(raw_seq_data, &mut raw_off)
                .ok_or_else(|| anyhow::anyhow!("Failed to read raw seq_len for read {}", i))?
                as usize;
            if raw_off + seq_len > raw_seq_data.len() {
                anyhow::bail!("Truncated raw seq data for read {}", i);
            }
            let seq = raw_seq_data[raw_off..raw_off + seq_len].to_vec();
            raw_off += seq_len;
            sequences.push(seq);
        } else {
            let ref_idx = dna_utils::read_varint(meta_data, &mut meta_off)
                .ok_or_else(|| anyhow::anyhow!("Failed to read ref_index for read {}", i))?
                as usize;
            if ref_idx >= sequences.len() {
                anyhow::bail!("Invalid ref_index {} for read {} (only {} decoded)",
                    ref_idx, i, sequences.len());
            }

            let seq_len = dna_utils::read_varint(delta_seq_data, &mut delta_off)
                .ok_or_else(|| anyhow::anyhow!("Failed to read delta seq_len for read {}", i))?
                as usize;
            if delta_off + seq_len > delta_seq_data.len() {
                anyhow::bail!("Truncated delta seq data for read {}", i);
            }
            let delta_bytes = &delta_seq_data[delta_off..delta_off + seq_len];
            delta_off += seq_len;

            let ref_seq = &sequences[ref_idx];
            let mut seq = Vec::with_capacity(seq_len);

            if class == 1 {
                for (j, &db) in delta_bytes.iter().enumerate() {
                    if db == 0x00 {
                        seq.push(if j < ref_seq.len() { ref_seq[j] } else { b'N' });
                    } else {
                        seq.push(db);
                    }
                }
            } else if class == 2 {
                let rn = ref_seq.len();
                for (j, &db) in delta_bytes.iter().enumerate() {
                    if db == 0x00 {
                        seq.push(if j < rn { complement(ref_seq[rn - 1 - j]) } else { b'N' });
                    } else {
                        seq.push(db);
                    }
                }
            } else {
                anyhow::bail!("Unknown class {} for read {}", class, i);
            }

            delta_count += 1;
            sequences.push(seq);
        }
    }

    info!("  Decoded {} delta-encoded reads out of {}", delta_count, num_reads);

    sequences.into_iter().enumerate().map(|(i, bytes)| {
        String::from_utf8(bytes)
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence {}: {}", i, e))
    }).collect()
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        x => x,
    }
}
