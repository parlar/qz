mod columnar;
mod delta;
mod n_mask;
pub mod paired_end;
mod quality;
mod quality_delta;
mod quality_model;
mod read_id;
mod rle;
mod zstd_dict;
pub mod arithmetic_quality;
pub mod arithmetic_sequence;
pub mod bsc;
pub mod debruijn;
pub mod dna_utils;
pub mod fqzcomp;
pub mod greedy_contig;
pub mod quality_context;
pub mod openzl;
pub mod template;
mod factorize;
mod harc;
mod quality_ctx;

use crate::cli::{CompressArgs, DecompressArgs, QualityMode, QualityCompressor, SequenceCompressor, HeaderCompressor};
use anyhow::{Context, Result};
use std::time::Instant;
use tracing::info;

/// Number of BSC blocks to decompress in each parallel batch
const DECOMPRESS_BATCH_SIZE: usize = 8;

/// Minimum reads for quality_ctx to outperform BSC (below this, models don't converge)
const MIN_READS_QUALITY_CTX: usize = 100_000;
/// I/O buffer size for reading archive files during decompression
const IO_BUFFER_SIZE: usize = 8 * 1024 * 1024; // 8 MB

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
        QualityCompressor::OpenZl => 2,
        QualityCompressor::Fqzcomp => 3,
        QualityCompressor::QualityCtx => 4,
    }
}

/// Convert byte code back to QualityCompressor
fn code_to_compressor(code: u8) -> Result<QualityCompressor> {
    match code {
        0 => Ok(QualityCompressor::Zstd),
        1 => Ok(QualityCompressor::Bsc),
        2 => Ok(QualityCompressor::OpenZl),
        3 => Ok(QualityCompressor::Fqzcomp),
        4 => Ok(QualityCompressor::QualityCtx),
        _ => anyhow::bail!("Invalid quality compressor code: {}", code),
    }
}

/// Convert SequenceCompressor to a byte code for storage
fn seq_compressor_to_code(compressor: SequenceCompressor) -> u8 {
    match compressor {
        SequenceCompressor::Zstd => 0,
        SequenceCompressor::Bsc => 1,
        SequenceCompressor::OpenZl => 2,
    }
}

/// Convert byte code back to SequenceCompressor
fn code_to_seq_compressor(code: u8) -> Result<SequenceCompressor> {
    match code {
        0 => Ok(SequenceCompressor::Zstd),
        1 => Ok(SequenceCompressor::Bsc),
        2 => Ok(SequenceCompressor::OpenZl),
        _ => anyhow::bail!("Invalid sequence compressor code: {}", code),
    }
}

/// Convert HeaderCompressor to a byte code for storage
fn header_compressor_to_code(compressor: HeaderCompressor) -> u8 {
    match compressor {
        HeaderCompressor::Zstd => 0,
        HeaderCompressor::Bsc => 1,
        HeaderCompressor::OpenZl => 2,
    }
}

/// Convert byte code back to HeaderCompressor
fn code_to_header_compressor(code: u8) -> Result<HeaderCompressor> {
    match code {
        0 => Ok(HeaderCompressor::Zstd),
        1 => Ok(HeaderCompressor::Bsc),
        2 => Ok(HeaderCompressor::OpenZl),
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

/// Encode a sequence into seq_stream using inline delta against a cache of previously seen reads.
/// If a similar read exists in the cache (same min syncmer hash, >60% positional identity,
/// same length), encodes as delta (0x00=match, actual base=mismatch). Otherwise emits raw.
/// Stream format: [varint(len)] [flag: 0x00=raw, 0x01=delta] [ref_idx varint if delta] [bytes]
fn encode_sequence_delta(
    seq: &[u8],
    seq_stream: &mut Vec<u8>,
    cache: &mut std::collections::HashMap<u64, usize>,
    reads: &mut Vec<Vec<u8>>,
) -> Result<bool> {
    let hash = dna_utils::compute_min_syncmer_hash(seq);
    let current_idx = reads.len();

    let mut used_delta = false;

    if hash != u64::MAX {
        if let Some(&ref_idx) = cache.get(&hash) {
            let ref_seq = &reads[ref_idx];
            if ref_seq.len() == seq.len() {
                // Count matching bases
                let matches = seq.iter().zip(ref_seq.iter()).filter(|(a, b)| a == b).count();
                let identity = matches as f64 / seq.len() as f64;

                if identity > 0.6 {
                    // Delta encode: [varint(len)] [0x01] [ref_idx varint] [delta bytes]
                    write_varint(seq_stream, seq.len())?;
                    seq_stream.push(0x01);
                    write_varint(seq_stream, ref_idx)?;
                    for (i, &base) in seq.iter().enumerate() {
                        if base == ref_seq[i] {
                            seq_stream.push(0x00); // match
                        } else {
                            seq_stream.push(base); // mismatch: actual base
                        }
                    }
                    used_delta = true;
                }
            }
        }
    }

    if !used_delta {
        // Raw encode: [varint(len)] [0x00] [raw bases]
        write_varint(seq_stream, seq.len())?;
        seq_stream.push(0x00);
        seq_stream.extend_from_slice(seq);
    }

    // Always update cache and store the read
    if hash != u64::MAX {
        cache.insert(hash, current_idx);
    }
    reads.push(seq.to_vec());

    Ok(used_delta)
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
    if args.sequence_delta {
        info!("Inline delta encoding enabled (syncmer-based cache)");
    }
    info!("Reading FASTQ and building streams...");

    let mut reader = FastqReader::from_path(input_path, args.fasta)?;
    let est_reads = 1_000_000;
    let mut header_stream: Vec<u8> = Vec::with_capacity(est_reads * 64);
    let mut seq_stream: Vec<u8> = Vec::with_capacity(est_reads * 152);
    let mut qual_stream: Vec<u8> = Vec::with_capacity(est_reads * 136);
    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;

    // Collect raw strings for quality_ctx (if applicable)
    let collect_for_ctx = !args.no_quality && args.quality_mode == QualityMode::Lossless;
    let mut qual_strings: Vec<String> = if collect_for_ctx { Vec::with_capacity(est_reads) } else { Vec::new() };
    let mut seq_strings: Vec<String> = if collect_for_ctx { Vec::with_capacity(est_reads) } else { Vec::new() };
    let mut fixed_read_len: Option<usize> = None;
    let mut variable_lengths = false;

    // Delta encoding state (only used when --sequence-delta)
    let mut delta_cache: std::collections::HashMap<u64, usize> = if args.sequence_delta {
        std::collections::HashMap::with_capacity(est_reads)
    } else {
        std::collections::HashMap::new()
    };
    let mut delta_reads: Vec<Vec<u8>> = if args.sequence_delta {
        Vec::with_capacity(est_reads)
    } else {
        Vec::new()
    };
    let mut delta_count: usize = 0;

    while let Some(record) = reader.next()? {
        // Header stream: varint(len) + raw bytes
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.extend_from_slice(record.id.as_bytes());

        // Sequence stream
        if args.sequence_delta {
            if encode_sequence_delta(record.sequence.as_bytes(), &mut seq_stream, &mut delta_cache, &mut delta_reads)? {
                delta_count += 1;
            }
        } else {
            write_varint(&mut seq_stream, record.sequence.len())?;
            if args.sequence_hints {
                seq_stream.push(dna_utils::compute_sequence_hint(record.sequence.as_bytes()));
            }
            seq_stream.extend_from_slice(record.sequence.as_bytes());
        }

        // Quality stream: varint(len) + packed bytes
        if !args.no_quality {
            if let Some(ref qual) = record.quality {
                let quantized = quality::quantize_quality(qual, args.quality_mode);
                write_varint(&mut qual_stream, quantized.len())?;
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                qual_stream.extend_from_slice(&packed);
            }
        }

        // Track sizes before potential move
        let seq_len = record.sequence.len();
        let qual_len = record.quality.as_ref().map(|q| q.len()).unwrap_or(0);
        total_bases += seq_len;
        original_size += record.id.len() + seq_len + qual_len + 3;

        // Track read length consistency and collect strings for quality_ctx
        if collect_for_ctx {
            match fixed_read_len {
                None => fixed_read_len = Some(seq_len),
                Some(prev) if prev != seq_len => variable_lengths = true,
                _ => {}
            }
            seq_strings.push(record.sequence);
            qual_strings.push(record.quality.unwrap_or_default());
        }

        num_reads += 1;

        if num_reads % 10_000_000 == 0 {
            info!("  {} million records read...", num_reads / 1_000_000);
        }
    }

    info!("Read {} records ({} bases)", num_reads, total_bases);
    if args.sequence_delta {
        info!("Delta stats: {}/{} reads delta-encoded ({:.1}%)",
            delta_count, num_reads, 100.0 * delta_count as f64 / num_reads as f64);
        drop(delta_cache);
        drop(delta_reads);
    }
    info!("Stream sizes: headers={} seq={} qual={}",
        header_stream.len(), seq_stream.len(), qual_stream.len());

    // Decide quality compression method
    // Auto-select quality_ctx for large fixed-length lossless datasets,
    // or respect an explicit QualityCompressor::QualityCtx setting
    let use_quality_ctx = collect_for_ctx
        && !variable_lengths
        && !args.no_quality
        && (args.quality_compressor == QualityCompressor::QualityCtx
            || num_reads >= MIN_READS_QUALITY_CTX);

    if args.quality_compressor == QualityCompressor::QualityCtx && !use_quality_ctx && variable_lengths {
        info!("Warning: quality-ctx requires fixed-length reads, falling back to BSC");
    }

    let quality_compressor_used = if use_quality_ctx {
        QualityCompressor::QualityCtx
    } else {
        QualityCompressor::Bsc
    };

    // Compress streams
    let (headers, sequences, qualities) = if use_quality_ctx {
        info!("Compressing qualities with context-adaptive range coder ({} reads)...", num_reads);
        drop(qual_stream); // Don't need packed stream for quality_ctx

        let qual_refs: Vec<&str> = qual_strings.iter().map(|s| s.as_str()).collect();
        let seq_refs: Vec<&str> = seq_strings.iter().map(|s| s.as_str()).collect();

        // Compress headers + sequences in parallel, then quality_ctx (needs sequence context)
        let (header_result, seq_result) = rayon::join(
            || bsc_compress_parallel(&header_stream, bsc_static),
            || bsc_compress_parallel(&seq_stream, bsc_static),
        );
        drop(header_stream);
        drop(seq_stream);

        let ctx_blob = quality_ctx::compress_qualities_ctx(&qual_refs, &seq_refs)?;
        drop(qual_strings);
        drop(seq_strings);
        let qualities = quality_ctx::wrap_as_multiblock(ctx_blob);

        info!("Quality ctx: {} bytes", qualities.len());
        (header_result?, seq_result?, qualities)
    } else {
        drop(qual_strings);
        drop(seq_strings);

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

        (header_result?, seq_result?, qual_result?)
    };

    // Write output file (same archive format as non-streaming path)
    info!("Writing output file...");
    let mut output_file = std::fs::File::create(&args.output)?;

    let encoding_type: u8 = if args.sequence_delta { 5 } else if args.sequence_hints { 4 } else { 0 };
    output_file.write_all(&[encoding_type])?;                                   // encoding_type
    output_file.write_all(&[0u8])?;                                             // arithmetic = disabled
    output_file.write_all(&[binning_to_code(quality_binning)])?;                // quality_binning
    output_file.write_all(&[compressor_to_code(quality_compressor_used)])?;     // quality_compressor
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

/// Result of reading a chunk of FASTQ records into compression streams.
struct ChunkStreams {
    header_stream: Vec<u8>,
    seq_stream: Vec<u8>,
    qual_stream: Vec<u8>,
    rc_flags: Vec<u8>,
    num_reads: usize,
    total_bases: usize,
    original_size: usize,
    /// Raw quality strings (only populated when collect_strings=true)
    qual_strings: Vec<String>,
    /// Raw sequence strings (only populated when collect_strings=true)
    seq_strings: Vec<String>,
    /// Whether reads have variable lengths (quality_ctx requires fixed length)
    variable_lengths: bool,
}

/// Read a chunk of FASTQ records into raw compression streams.
/// rc_flags is empty unless rc_canon=true, in which case it has 1 byte per read (0=fwd, 1=revcomp).
/// When collect_strings=true, also keeps raw quality and sequence strings for quality_ctx.
fn read_chunk_streams<R: std::io::BufRead>(
    reader: &mut crate::io::FastqReader<R>,
    chunk_size: usize,
    quality_mode: QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    sequence_hints: bool,
    sequence_delta: bool,
    rc_canon: bool,
    collect_strings: bool,
) -> Result<ChunkStreams> {
    let mut header_stream = Vec::with_capacity(chunk_size * 64);
    let mut seq_stream = Vec::with_capacity(chunk_size * 152);
    let mut qual_stream = Vec::with_capacity(chunk_size * 136);
    let mut rc_flags: Vec<u8> = if rc_canon { Vec::with_capacity(chunk_size) } else { Vec::new() };
    let mut qual_strings: Vec<String> = if collect_strings { Vec::with_capacity(chunk_size) } else { Vec::new() };
    let mut seq_strings: Vec<String> = if collect_strings { Vec::with_capacity(chunk_size) } else { Vec::new() };
    let mut chunk_reads = 0usize;
    let mut chunk_bases = 0usize;
    let mut chunk_orig = 0usize;
    let mut fixed_read_len: Option<usize> = None;
    let mut variable_lengths = false;

    let mut delta_cache = std::collections::HashMap::new();
    let mut delta_reads = Vec::new();

    for _ in 0..chunk_size {
        match reader.next()? {
            Some(record) => {
                write_varint(&mut header_stream, record.id.len())?;
                header_stream.extend_from_slice(record.id.as_bytes());

                if sequence_delta {
                    encode_sequence_delta(record.sequence.as_bytes(), &mut seq_stream, &mut delta_cache, &mut delta_reads)?;
                } else {
                    write_varint(&mut seq_stream, record.sequence.len())?;
                    if sequence_hints {
                        seq_stream.push(dna_utils::compute_sequence_hint(record.sequence.as_bytes()));
                    }
                    if rc_canon {
                        let (canon_seq, was_reversed) = dna_utils::canonicalize_sequence(record.sequence.as_bytes());
                        seq_stream.extend_from_slice(&canon_seq);
                        rc_flags.push(was_reversed as u8);
                    } else {
                        seq_stream.extend_from_slice(record.sequence.as_bytes());
                    }
                }

                if !no_quality {
                    if let Some(ref qual) = record.quality {
                        let quantized = quality::quantize_quality(qual, quality_mode);
                        write_varint(&mut qual_stream, quantized.len())?;
                        let packed = columnar::pack_qualities(&quantized, quality_binning);
                        qual_stream.extend_from_slice(&packed);
                    }
                }

                let seq_len = record.sequence.len();
                let qual_len = record.quality.as_ref().map(|q| q.len()).unwrap_or(0);
                chunk_bases += seq_len;
                chunk_orig += record.id.len() + seq_len + qual_len + 3;

                if collect_strings {
                    match fixed_read_len {
                        None => fixed_read_len = Some(seq_len),
                        Some(prev) if prev != seq_len => variable_lengths = true,
                        _ => {}
                    }
                    seq_strings.push(record.sequence);
                    qual_strings.push(record.quality.unwrap_or_default());
                }

                chunk_reads += 1;
            }
            None => break,
        }
    }

    Ok(ChunkStreams {
        header_stream,
        seq_stream,
        qual_stream,
        rc_flags,
        num_reads: chunk_reads,
        total_bases: chunk_bases,
        original_size: chunk_orig,
        qual_strings,
        seq_strings,
        variable_lengths,
    })
}

/// Read a chunk of FASTQ records into a Vec (for reorder modes that need to sort before streaming).
/// Returns (records, total_bases, original_size).
fn read_chunk_records<R: std::io::BufRead>(
    reader: &mut crate::io::FastqReader<R>,
    chunk_size: usize,
) -> Result<(Vec<crate::io::FastqRecord>, usize, usize)> {
    let mut records = Vec::with_capacity(chunk_size.min(1_000_000));
    let mut chunk_bases = 0usize;
    let mut chunk_orig = 0usize;

    for _ in 0..chunk_size {
        match reader.next()? {
            Some(record) => {
                chunk_bases += record.sequence.len();
                chunk_orig += record.id.len() + record.sequence.len()
                    + record.quality.as_ref().map(|q| q.len()).unwrap_or(0) + 3;
                records.push(record);
            }
            None => break,
        }
    }

    Ok((records, chunk_bases, chunk_orig))
}

/// Convert a slice of FastqRecords into raw compression streams.
/// Returns (header_stream, seq_stream, qual_stream).
/// Convert records to compression streams.
/// Returns (header_stream, seq_stream, qual_stream, rc_flags).
/// rc_flags is empty unless rc_canon=true.
fn records_to_streams(
    records: &[crate::io::FastqRecord],
    quality_mode: QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    sequence_hints: bool,
    sequence_delta: bool,
    rc_canon: bool,
) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)> {
    let n = records.len();
    let mut header_stream = Vec::with_capacity(n * 64);
    let mut seq_stream = Vec::with_capacity(n * 152);
    let mut qual_stream = Vec::with_capacity(n * 136);
    let mut rc_flags: Vec<u8> = if rc_canon { Vec::with_capacity(n) } else { Vec::new() };

    let mut delta_cache = std::collections::HashMap::new();
    let mut delta_reads = Vec::new();

    for record in records {
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.extend_from_slice(record.id.as_bytes());

        if sequence_delta {
            encode_sequence_delta(record.sequence.as_bytes(), &mut seq_stream, &mut delta_cache, &mut delta_reads)?;
        } else {
            write_varint(&mut seq_stream, record.sequence.len())?;
            if sequence_hints {
                seq_stream.push(dna_utils::compute_sequence_hint(record.sequence.as_bytes()));
            }
            if rc_canon {
                let (canon_seq, was_reversed) = dna_utils::canonicalize_sequence(record.sequence.as_bytes());
                seq_stream.extend_from_slice(&canon_seq);
                rc_flags.push(was_reversed as u8);
            } else {
                seq_stream.extend_from_slice(record.sequence.as_bytes());
            }
        }

        if !no_quality {
            if let Some(ref qual) = record.quality {
                let quantized = quality::quantize_quality(qual, quality_mode);
                write_varint(&mut qual_stream, quantized.len())?;
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                qual_stream.extend_from_slice(&packed);
            }
        }
    }

    Ok((header_stream, seq_stream, qual_stream, rc_flags))
}

/// Sort records by reorder sort key, returning a new sorted Vec.
fn sort_records_by_key(records: Vec<crate::io::FastqRecord>) -> Vec<crate::io::FastqRecord> {
    use rayon::prelude::*;
    let keys: Vec<u128> = records
        .par_iter()
        .map(|r| {
            let qual_bytes = r.quality.as_ref().map(|q| q.as_bytes()).unwrap_or(&[]);
            dna_utils::reorder_sort_key(r.sequence.as_bytes(), qual_bytes)
        })
        .collect();
    let mut keyed: Vec<(u128, crate::io::FastqRecord)> = keys.into_iter()
        .zip(records.into_iter())
        .collect();
    keyed.sort_unstable_by_key(|(k, _)| *k);
    keyed.into_iter().map(|(_, r)| r).collect()
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
    use std::io::{Write, BufWriter};
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
    let h_tmp_path = working_dir.join(".qz_chunked_h.tmp");
    let s_tmp_path = working_dir.join(".qz_chunked_s.tmp");
    let q_tmp_path = working_dir.join(".qz_chunked_q.tmp");

    // Cleanup guard: remove temp files on drop (handles success, error, and panic)
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

    // Decide if quality_ctx is possible (lossless mode, not no_quality)
    let collect_for_ctx = !no_quality && quality_mode == QualityMode::Lossless;

    // Read first chunk
    let mut cur = read_chunk_streams(
        &mut reader, CHUNK_SIZE, quality_mode, quality_binning, no_quality,
        args.sequence_hints, args.sequence_delta, rc_canon, collect_for_ctx,
    )?;

    // Decide quality_ctx usage from first chunk (must be consistent across all chunks)
    let use_quality_ctx = collect_for_ctx
        && !cur.variable_lengths
        && (args.quality_compressor == QualityCompressor::QualityCtx
            || cur.num_reads >= MIN_READS_QUALITY_CTX);
    let quality_compressor_used = if use_quality_ctx {
        QualityCompressor::QualityCtx
    } else {
        QualityCompressor::Bsc
    };
    if use_quality_ctx {
        info!("Using context-adaptive quality compression (quality_ctx)");
    } else if args.quality_compressor == QualityCompressor::QualityCtx && cur.variable_lengths {
        info!("Warning: quality-ctx requires fixed-length reads, falling back to BSC");
    }

    while cur.num_reads > 0 {
        info!("Chunk {}: {} reads (h={} s={} q={} bytes)",
            chunk_idx, cur.num_reads, cur.header_stream.len(), cur.seq_stream.len(), cur.qual_stream.len());

        // Pipeline: compress current chunk on a background thread while
        // reading the next chunk on the main thread.
        //
        // The compression thread processes streams SEQUENTIALLY (headers ->
        // sequences -> qualities) to limit BSC working memory. Raw data for
        // each stream is dropped after compression before starting the next.
        let (next_result, compress_result) = std::thread::scope(|scope| {
            let h_data = std::mem::take(&mut cur.header_stream);
            let s_data = std::mem::take(&mut cur.seq_stream);
            let q_data = std::mem::take(&mut cur.qual_stream);
            let rc_data = std::mem::take(&mut cur.rc_flags);
            let qs = std::mem::take(&mut cur.qual_strings);
            let ss = std::mem::take(&mut cur.seq_strings);

            let compress_handle = scope.spawn(move || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                // Compress headers, then free raw header data before starting sequences
                let h_blocks = compress_stream_to_bsc_blocks(&h_data, bsc_static)?;
                drop(h_data);

                // Compress sequences, then free raw sequence data before starting qualities
                let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
                drop(s_data);

                // Compress qualities
                let q_blocks = if no_quality || (q_data.is_empty() && qs.is_empty()) {
                    Vec::new()
                } else if use_quality_ctx {
                    drop(q_data); // Don't need packed stream
                    let qual_refs: Vec<&str> = qs.iter().map(|s| s.as_str()).collect();
                    let seq_refs: Vec<&str> = ss.iter().map(|s| s.as_str()).collect();
                    let blob = quality_ctx::compress_qualities_ctx(&qual_refs, &seq_refs)?;
                    drop(qs);
                    drop(ss);
                    vec![blob] // Single blob as one "block"
                } else {
                    drop(qs);
                    drop(ss);
                    compress_stream_to_bsc_blocks(&q_data, bsc_static)?
                };

                // Compress RC flags (if present)
                let rc_blocks = if rc_data.is_empty() {
                    Vec::new()
                } else {
                    compress_stream_to_bsc_blocks(&rc_data, bsc_static)?
                };

                Ok((h_blocks, s_blocks, q_blocks, rc_blocks))
            });

            // Read next chunk on main thread (overlaps with compression)
            let next = read_chunk_streams(
                &mut reader, CHUNK_SIZE, quality_mode, quality_binning, no_quality,
                args.sequence_hints, args.sequence_delta, rc_canon,
                collect_for_ctx && use_quality_ctx,
            );

            let compressed = compress_handle.join().unwrap();
            (next, compressed)
        });

        // Write compressed blocks to temp files immediately, freeing block memory
        let (h_blk, s_blk, q_blk, rc_blk) = compress_result?;
        h_num_blocks += write_blocks_to_tmp(h_blk, &mut h_tmp)?;
        s_num_blocks += write_blocks_to_tmp(s_blk, &mut s_tmp)?;
        q_num_blocks += write_blocks_to_tmp(q_blk, &mut q_tmp)?;
        if let Some(ref mut rc_file) = rc_tmp {
            rc_num_blocks += write_blocks_to_tmp(rc_blk, rc_file)?;
        }

        num_reads += cur.num_reads;
        total_bases += cur.total_bases;
        original_size += cur.original_size;
        chunk_idx += 1;

        cur = next_result?;
    }

    // Flush and close temp file writers
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

    // Stream sizes in multi-block format: [num_blocks: u32] + block data from temp file
    let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
    let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };
    let rc_flags_len = if rc_num_blocks > 0 { 4 + rc_tmp_size } else { 0 };

    info!("Read {} records in {} chunks ({} bases)", num_reads, chunk_idx, total_bases);
    info!("Compressed blocks: headers={} ({}) seq={} ({}) qual={} ({})",
        h_num_blocks, humanize_bytes(h_tmp_size),
        s_num_blocks, humanize_bytes(s_tmp_size),
        q_num_blocks, humanize_bytes(q_tmp_size));
    if rc_canon {
        info!("  RC flags: {} ({})", rc_num_blocks, humanize_bytes(rc_tmp_size));
    }

    // Write final output
    let encoding_type: u8 = if rc_canon { 6 } else if args.sequence_delta { 5 } else if args.sequence_hints { 4 } else { 0 };
    write_chunked_archive_rc(
        &args.output, quality_binning, quality_compressor_used, no_quality, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len, rc_flags_len,
        h_num_blocks, s_num_blocks, q_num_blocks, rc_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path, &rc_tmp_path,
        original_size, start_time,
    )
    // Temp files cleaned up by TmpCleanup drop guard
}

/// Chunked streaming compression with local reordering within each chunk.
///
/// Same bounded-memory design as `compress_chunked_bsc`, but reads records into a Vec
/// per chunk, sorts by `reorder_sort_key`, then builds streams from sorted records.
/// Reads within each 5M-record chunk are sorted by content similarity; across chunks
/// the order follows the input file.
fn compress_chunked_bsc_reorder_local(args: &CompressArgs) -> Result<()> {
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
                            let q_blob = compress_qualities_fqzcomp(&records)?;
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

            let compressed = compress_handle.join().unwrap();
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
    write_chunked_archive_rc(
        &args.output, quality_binning, args.quality_compressor, no_quality, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len, rc_flags_len,
        h_num_blocks, s_num_blocks, q_num_blocks, rc_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path, &rc_tmp_path,
        original_size, start_time,
    )
}


/// Chunked streaming compression with fqzcomp for quality scores.
///
/// Reads records per chunk (fqzcomp needs individual quality strings, not a byte stream).
/// Headers and sequences compressed with BSC blocks as usual.
/// Each chunk's qualities compressed with fqzcomp independently, stored as one "block"
/// in the multi-block format. Decompressor reads blocks and fqzcomp-decompresses each.
fn compress_chunked_fqzcomp(args: &CompressArgs) -> Result<()> {
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
                            compress_qualities_fqzcomp(&records)
                        }
                    },
                );
                let (h_blocks, s_blocks) = bsc_result?;
                let q_blob = q_blob?;

                Ok((h_blocks, s_blocks, q_blob))
            });

            let next = read_chunk_records(&mut reader, CHUNK_SIZE);
            let compressed = compress_handle.join().unwrap();
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

    write_chunked_archive_with_compressor(
        &args.output, QualityBinning::None, no_quality,
        QualityCompressor::Fqzcomp, 0u8,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path,
        original_size, start_time,
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
fn compress_global_reorder_bsc(args: &CompressArgs) -> Result<()> {
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
        &args.output, quality_binning, no_quality, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path,
        original_size, start_time,
    )
}

/// Write the final archive from temp files (shared by chunked, local reorder, and global reorder).
fn write_chunked_archive(
    output_path: &std::path::Path,
    quality_binning: QualityBinning,
    no_quality: bool,
    encoding_type: u8,
    num_reads: usize,
    headers_len: usize,
    sequences_len: usize,
    qualities_len: usize,
    h_num_blocks: u32,
    s_num_blocks: u32,
    q_num_blocks: u32,
    h_tmp_path: &std::path::Path,
    s_tmp_path: &std::path::Path,
    q_tmp_path: &std::path::Path,
    original_size: usize,
    start_time: Instant,
) -> Result<()> {
    write_chunked_archive_with_compressor(
        output_path, quality_binning, no_quality,
        QualityCompressor::Bsc, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        h_tmp_path, s_tmp_path, q_tmp_path,
        original_size, start_time,
    )
}

fn write_chunked_archive_with_compressor(
    output_path: &std::path::Path,
    quality_binning: QualityBinning,
    no_quality: bool,
    quality_compressor: QualityCompressor,
    encoding_type: u8,
    num_reads: usize,
    headers_len: usize,
    sequences_len: usize,
    qualities_len: usize,
    h_num_blocks: u32,
    s_num_blocks: u32,
    q_num_blocks: u32,
    h_tmp_path: &std::path::Path,
    s_tmp_path: &std::path::Path,
    q_tmp_path: &std::path::Path,
    original_size: usize,
    start_time: Instant,
) -> Result<()> {
    use std::io::{Write, BufReader};

    info!("Writing output file...");
    let mut output_file = std::fs::File::create(output_path)?;

    output_file.write_all(&[encoding_type])?;                                    // encoding_type
    output_file.write_all(&[0u8])?;                                             // arithmetic = disabled
    output_file.write_all(&[binning_to_code(quality_binning)])?;                // quality_binning
    output_file.write_all(&[compressor_to_code(quality_compressor)])?;           // quality_compressor
    output_file.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?; // seq_compressor
    output_file.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;// header_compressor
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

    let copy_stream = |num_blocks: u32, tmp_path: &std::path::Path, out: &mut std::fs::File| -> Result<()> {
        if num_blocks == 0 { return Ok(()); }
        out.write_all(&num_blocks.to_le_bytes())?;
        let mut tmp_reader = BufReader::new(std::fs::File::open(tmp_path)?);
        std::io::copy(&mut tmp_reader, out)?;
        Ok(())
    };

    copy_stream(h_num_blocks, h_tmp_path, &mut output_file)?;
    copy_stream(s_num_blocks, s_tmp_path, &mut output_file)?;
    copy_stream(q_num_blocks, q_tmp_path, &mut output_file)?;

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

    Ok(())
}

/// Write chunked archive with optional RC flags stream appended after qualities.
/// RC flags stream layout in archive: [rc_flags_len: u64 in header] ... [num_blocks: u32][block_data]
fn write_chunked_archive_rc(
    output_path: &std::path::Path,
    quality_binning: QualityBinning,
    quality_compressor: QualityCompressor,
    no_quality: bool,
    encoding_type: u8,
    num_reads: usize,
    headers_len: usize,
    sequences_len: usize,
    qualities_len: usize,
    rc_flags_len: usize,
    h_num_blocks: u32,
    s_num_blocks: u32,
    q_num_blocks: u32,
    rc_num_blocks: u32,
    h_tmp_path: &std::path::Path,
    s_tmp_path: &std::path::Path,
    q_tmp_path: &std::path::Path,
    rc_tmp_path: &std::path::Path,
    original_size: usize,
    start_time: Instant,
) -> Result<()> {
    use std::io::{Write, BufReader};

    info!("Writing output file...");
    let mut output_file = std::fs::File::create(output_path)?;

    output_file.write_all(&[encoding_type])?;                                    // encoding_type (6=rc_canon)
    output_file.write_all(&[0u8])?;                                             // arithmetic = disabled
    output_file.write_all(&[binning_to_code(quality_binning)])?;                // quality_binning
    output_file.write_all(&[compressor_to_code(quality_compressor)])?;           // quality_compressor
    output_file.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?;  // seq_compressor
    output_file.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?; // header_compressor
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

    let copy_stream = |num_blocks: u32, tmp_path: &std::path::Path, out: &mut std::fs::File| -> Result<()> {
        if num_blocks == 0 { return Ok(()); }
        out.write_all(&num_blocks.to_le_bytes())?;
        let mut tmp_reader = BufReader::new(std::fs::File::open(tmp_path)?);
        std::io::copy(&mut tmp_reader, out)?;
        Ok(())
    };

    copy_stream(h_num_blocks, h_tmp_path, &mut output_file)?;
    copy_stream(s_num_blocks, s_tmp_path, &mut output_file)?;
    copy_stream(q_num_blocks, q_tmp_path, &mut output_file)?;

    // Append RC flags stream after qualities (encoding_type=6 signals its presence)
    if rc_num_blocks > 0 {
        output_file.write_all(&(rc_flags_len as u64).to_le_bytes())?;
        copy_stream(rc_num_blocks, rc_tmp_path, &mut output_file)?;
    }

    let metadata_size = 9 + 2 + 1 + 40 + if rc_num_blocks > 0 { 8 } else { 0 };
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

pub fn compress(args: &CompressArgs) -> Result<()> {
    use crate::io::{FastqReader, FastqRecord};
    use std::io::Write;

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

    // Validate --local-reorder / --ultra compatibility
    if args.local_reorder || args.ultra {
        let mode_name = if args.local_reorder { "--local-reorder" } else { "--ultra" };
        if args.local_reorder && args.ultra {
            anyhow::bail!("--local-reorder and --ultra cannot be combined");
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

    // Detect paired-end mode
    if args.input.len() == 2 {
        info!("Paired-end mode detected - using correlation compression");
        return compress_paired_end(&args.input[0], &args.input[1], args);
    } else if args.input.len() > 2 {
        anyhow::bail!("More than 2 input files not supported");
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
        return harc::compress_harc(args);
    }

    // Ultra mode: single large BSC block with parallel BWT (best ratio)
    if args.ultra {
        return harc::compress_reorder_local(args);
    }

    // Fqzcomp quality needs record-level access, so use chunked record-based path
    if can_stream_base && args.quality_compressor == QualityCompressor::Fqzcomp {
        return compress_chunked_fqzcomp(args);
    }

    if can_stream {
        // RC canon always uses chunked path (needs temp file for flags stream)
        if args.chunked || args.rc_canon {
            return compress_chunked_bsc(args);
        }
        return compress_streaming_bsc(args);
    }

    // Non-streaming mode: read all records into memory first
    // (required for advanced quality encoding, non-BSC compressors)
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
            compress_qualities_with_model(&processed_records, args.compression_level, args.quality_compressor, args.bsc_static)?
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
                    if args.header_template {
                        let read_ids: Vec<String> = processed_records.iter().map(|r| r.id.clone()).collect();
                        let encoded = read_id::compress_read_ids(&read_ids)?;
                        if encoded.template.prefix.is_empty() {
                            // Template analysis found no common Illumina structure, fall back to raw BSC
                            info!("No template structure found in headers, falling back to raw BSC{}...", if args.bsc_static { " (static)" } else { " (adaptive)" });
                            let compressed = compress_headers_bsc_with(&processed_records, args.bsc_static)?;
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
                        let compressed = compress_headers_bsc_with(&processed_records, args.bsc_static)?;
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
                    compress_headers(&processed_records, args.compression_level)
                }
                HeaderCompressor::OpenZl => {
                    info!("Compressing read IDs with raw OpenZL...");
                    let compressed = compress_headers_openzl(&processed_records)?;
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
                        let twobit = args.twobit;
                        info!("Compressing sequences{} and qualities in parallel with BSC{}...",
                            if twobit { " (2-bit)" } else { "" },
                            if args.bsc_static { " (static)" } else { " (adaptive)" });
                        let bsc_static = args.bsc_static;
                        let (seq_result, qual_result) = rayon::join(
                            || if twobit {
                                compress_sequences_2bit_bsc_with(&processed_records, bsc_static)
                            } else {
                                compress_sequences_raw_bsc_with(&processed_records, bsc_static, args.sequence_hints, args.sequence_delta)
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
                            || compress_sequences_raw_openzl(&processed_records),
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

    let mut blocks_done = 0;

    while blocks_done < num_blocks {
        let batch_size = DECOMPRESS_BATCH_SIZE.min(num_blocks - blocks_done);

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
                anyhow::bail!("Decompressor channel closed: receiver dropped before all blocks were consumed");
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
            // Compact consumed data: copy_within + truncate avoids drain's iterator overhead
            if self.pos > 0 {
                let remaining = self.buf.len() - self.pos;
                if remaining > 0 {
                    self.buf.copy_within(self.pos.., 0);
                }
                self.buf.truncate(remaining);
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
        let mut bytes_read = 0;
        loop {
            if !self.fill(1)? {
                anyhow::bail!("Unexpected end of stream reading varint");
            }
            let byte = self.buf[self.pos];
            self.pos += 1;
            bytes_read += 1;
            value |= ((byte & 0x7F) as usize) << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
            shift += 7;
            if bytes_read >= MAX_VARINT_BYTES {
                anyhow::bail!("Malformed varint: too many continuation bytes");
            }
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
/// Reads the first 52 bytes of the header to check all flags.
fn can_stream_decompress(args: &DecompressArgs) -> Result<bool> {
    use std::io::Read;

    let mut file = std::fs::File::open(&args.input)?;
    let mut header = [0u8; 52];
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
    let nmasks_len = u64::from_le_bytes(header[36..44].try_into().unwrap());

    Ok((encoding_type == 0 || encoding_type == 4 || encoding_type == 6) // 0=raw, 4=raw+hints, 6=rc_canon
        && !arithmetic_enabled
        && !quality_model_enabled
        && !quality_delta_enabled
        && !quality_dict_present
        && template_prefix_len == 0
        && nmasks_len == 0
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
    let mut file = std::fs::File::open(&args.input)
        .with_context(|| format!("Failed to open archive: {:?}", args.input))?;
    let mut header = [0u8; 52];
    file.read_exact(&mut header)
        .context("Failed to read archive header")?;
    drop(file);

    let encoding_type = header[0];
    let has_sequence_hints = encoding_type == 4;
    let has_rc_canon = encoding_type == 6;
    let quality_binning = code_to_binning(header[2])?;
    let _has_comment = header[11] != 0;

    let num_reads = u64::from_le_bytes(header[12..20].try_into().unwrap()) as usize;
    let headers_len = u64::from_le_bytes(header[20..28].try_into().unwrap()) as usize;
    let sequences_len = u64::from_le_bytes(header[28..36].try_into().unwrap()) as usize;
    let nmasks_len = u64::from_le_bytes(header[36..44].try_into().unwrap()) as usize;
    let qualities_len = u64::from_le_bytes(header[44..52].try_into().unwrap()) as usize;

    // For RC canon archives, read the rc_flags_len from after the standard streams
    let (rc_flags_len, rc_flags_offset) = if has_rc_canon {
        let q_end = 52u64 + headers_len as u64 + sequences_len as u64 + nmasks_len as u64 + qualities_len as u64;
        let mut file = std::fs::File::open(&args.input)?;
        use std::io::Seek;
        file.seek(std::io::SeekFrom::Start(q_end))?;
        let mut len_buf = [0u8; 8];
        file.read_exact(&mut len_buf)?;
        let rc_len = u64::from_le_bytes(len_buf) as usize;
        let rc_off = q_end + 8;
        (rc_len, rc_off)
    } else {
        (0, 0u64)
    };

    info!("Archive: {} reads, headers={} seq={} qual={}{}",
        num_reads, humanize_bytes(headers_len), humanize_bytes(sequences_len), humanize_bytes(qualities_len),
        if has_rc_canon { format!(" rc_flags={}", humanize_bytes(rc_flags_len)) } else { String::new() });

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
            IO_BUFFER_SIZE,
            GzEncoder::new(
                std::fs::File::create(output_path)
                    .with_context(|| format!("Failed to create output file: {:?}", output_path))?,
                Compression::new(args.gzip_level),
            ),
        ))
    } else {
        Box::new(std::io::BufWriter::with_capacity(
            IO_BUFFER_SIZE,
            std::fs::File::create(output_path)
                .with_context(|| format!("Failed to create output file: {:?}", output_path))?,
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

        // RC flags stream (encoding_type=6)
        let rc_rx = if has_rc_canon && rc_flags_len > 0 {
            let (rc_tx, rc_rx) = std::sync::mpsc::sync_channel(2);
            let rc_path = archive_path.clone();
            scope.spawn(move || stream_decompressor(&rc_path, rc_flags_offset, rc_flags_len, rc_tx));
            Some(rc_rx)
        } else {
            None
        };

        let mut h_buf = ChannelStreamBuffer::new(h_rx);
        let mut s_buf = ChannelStreamBuffer::new(s_rx);
        let mut q_buf = q_rx.map(ChannelStreamBuffer::new);
        let mut rc_buf = rc_rx.map(ChannelStreamBuffer::new);

        for i in 0..num_reads {
            // Header: varint(len) + raw bytes
            let h_len = h_buf.read_varint()?;
            output.write_all(h_buf.read_bytes(h_len)?)?;
            output.write_all(b"\n")?;

            // Sequence: varint(len) [+ hint byte] + raw bytes
            let s_len = s_buf.read_varint()?;
            if has_sequence_hints {
                s_buf.read_bytes(1)?; // skip hint byte
            }
            let seq_bytes = s_buf.read_bytes(s_len)?;

            // If RC-canonicalized, check flag and reverse-complement if needed
            if let Some(ref mut rc) = rc_buf {
                let flag = rc.read_bytes(1)?[0];
                if flag != 0 {
                    let revcomp = dna_utils::reverse_complement(seq_bytes);
                    output.write_all(&revcomp)?;
                } else {
                    output.write_all(seq_bytes)?;
                }
            } else {
                output.write_all(seq_bytes)?;
            }
            output.write_all(b"\n+\n")?;

            // Quality: varint(orig_len) + packed bytes
            if let Some(ref mut q) = q_buf {
                let q_len = q.read_varint()?;
                let packed_len = (q_len * bits_per_qual + 7) / 8;
                let packed = q.read_bytes(packed_len)?;
                columnar::unpack_qualities_to_writer(packed, q_len, quality_binning, &mut output)?;
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
    let mut input_file = std::fs::File::open(&args.input)
        .with_context(|| format!("Failed to open archive: {:?}", args.input))?;
    let mut archive_data = Vec::new();
    input_file.read_to_end(&mut archive_data)
        .context("Failed to read archive data")?;

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
    let quality_compressor = if offset < archive_data.len() && archive_data[offset] <= 4 {
        let compressor = code_to_compressor(archive_data[offset])?;
        offset += 1;
        compressor
    } else {
        // Old format without compressor field - assume zstd
        QualityCompressor::Zstd
    };

    // Fqzcomp stores raw ASCII quality bytes (not bit-packed), so override binning
    let quality_binning = if quality_compressor == QualityCompressor::Fqzcomp {
        QualityBinning::None
    } else {
        quality_binning
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

    // Read common comment (only present when template is in use and has_comment is true)
    let common_comment = if template_prefix_len > 0 && template_has_comment {
        let cc_len = u16::from_le_bytes(archive_data[offset..offset + 2].try_into().unwrap()) as usize;
        offset += 2;
        if cc_len > 0 {
            let cc = String::from_utf8_lossy(&archive_data[offset..offset + cc_len]).to_string();
            offset += cc_len;
            Some(cc)
        } else {
            None
        }
    } else {
        None
    };

    let read_id_template = read_id::ReadIdTemplate {
        prefix: template_prefix,
        has_comment: template_has_comment,
        common_comment,
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
            HeaderCompressor::Bsc => {
                if template_prefix_len > 0 {
                    let header_data = bsc::decompress_parallel(headers)
                        .context("Failed to decompress headers (BSC)")?;
                    read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                        .context("Failed to decode template-encoded read IDs")?
                } else {
                    decompress_headers_bsc(headers, num_reads)
                        .context("Failed to decompress headers (BSC)")?
                }
            }
            HeaderCompressor::Zstd => {
                let header_data = decompress_zstd(headers)
                    .context("Failed to decompress headers (zstd)")?;
                read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                    .context("Failed to decode read IDs")?
            }
            HeaderCompressor::OpenZl => {
                if template_prefix_len > 0 {
                    let header_data = openzl::decompress_parallel(headers)
                        .context("Failed to decompress headers (OpenZL)")?;
                    read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                        .context("Failed to decode template-encoded read IDs")?
                } else {
                    decompress_headers_openzl(headers, num_reads)
                        .context("Failed to decompress headers (OpenZL)")?
                }
            }
        };

        // Decompress sequences with arithmetic coding
        let read_lengths = read_lengths_opt.ok_or_else(|| anyhow::anyhow!("Arithmetic mode requires read lengths"))?;
        let decoded_sequences = arithmetic_sequence::decode_sequences_arithmetic(sequences, &read_lengths, num_reads)
            .context("Failed to decode arithmetic-coded sequences")?;

        // Decompress qualities with arithmetic coding
        let decoded_qualities = if qualities_len > 0 {
            let read_length = if !read_lengths.is_empty() { read_lengths[0] } else { 0 };
            arithmetic_quality::decode_qualities_arithmetic(qualities, &decoded_sequences, read_length, num_reads)
                .context("Failed to decode arithmetic-coded quality scores")?
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
            HeaderCompressor::Bsc => {
                if template_prefix_len > 0 {
                    let header_data = bsc::decompress_parallel(headers)
                        .context("Failed to decompress headers (BSC)")?;
                    read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                        .context("Failed to decode template-encoded read IDs")?
                } else {
                    decompress_headers_bsc(headers, num_reads)
                        .context("Failed to decompress headers (BSC)")?
                }
            }
            HeaderCompressor::Zstd => {
                let header_data = decompress_zstd(headers)
                    .context("Failed to decompress headers (zstd)")?;
                read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                    .context("Failed to decode read IDs")?
            }
            HeaderCompressor::OpenZl => {
                if template_prefix_len > 0 {
                    let header_data = openzl::decompress_parallel(headers)
                        .context("Failed to decompress headers (OpenZL)")?;
                    read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                        .context("Failed to decode template-encoded read IDs")?
                } else {
                    decompress_headers_openzl(headers, num_reads)
                        .context("Failed to decompress headers (OpenZL)")?
                }
            }
        };

        // Decompress sequences with de Bruijn graph
        let decoded_sequences = debruijn::decompress_sequences_debruijn(sequences, num_reads)
            .context("Failed to decompress de Bruijn-coded sequences")?;

        // Decompress qualities
        let qualities_data = if qualities_len == 0 {
            Vec::new()
        } else if let Some(ref dict) = quality_dict_opt {
            zstd_dict::decompress_with_dict(qualities, dict)
                .context("Failed to decompress quality scores (dictionary mode)")?
        } else {
            decompress_qualities_data(qualities, quality_compressor)
                .context("Failed to decompress quality scores")?
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
                        let quality_str = if quality_compressor == QualityCompressor::Fqzcomp {
                            if qual_offset + q_len <= qualities_data.len() {
                                let raw = &qualities_data[qual_offset..qual_offset + q_len];
                                qual_offset += q_len;
                                unsafe { String::from_utf8_unchecked(raw.to_vec()) }
                            } else {
                                break;
                            }
                        } else if let Some(ref model) = quality_model_opt {
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
    } else if encoding_type == 7 {
        info!("Decompressing factorized sequences (encoding_type=7)...");

        // Detect format version: v2 starts with version byte 0x02,
        // v1 starts with u64 meta_size whose LSB is always ≥ 32.
        let is_v2 = !sequences.is_empty() && sequences[0] == 2;

        if is_v2 {
            // ── Factorize v2: pattern-routed streams ──
            // Decompress headers + qualities in parallel with sequence decoding
            let (header_result, qual_result) = rayon::join(
                || decompress_headers_bsc(headers, num_reads)
                    .context("Failed to decompress headers"),
                || -> Result<Vec<u8>> {
                    if qualities_len == 0 {
                        Ok(Vec::new())
                    } else {
                        decompress_qualities_data(qualities, quality_compressor)
                            .context("Failed to decompress qualities")
                    }
                },
            );

            let read_ids = header_result?;
            let qualities_data = qual_result?;

            let decoded_sequences = factorize::decode_factorized_sequences_v2(
                sequences, num_reads)?;

            // Build records
            let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
                crate::io::FastqRecord::new(id, seq, None)
            }).collect();

            // Decode qualities
            if !qualities_data.is_empty() {
                let mut qual_offset = 0;
                for record in &mut records {
                    if qual_offset < qualities_data.len() {
                        let q_len = read_varint(&qualities_data, &mut qual_offset)
                            .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                        let bits_per_qual = quality_binning.bits_per_quality();
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        if qual_offset + q_encoded_len <= qualities_data.len() {
                            let quality_str = columnar::unpack_qualities(&qualities_data[qual_offset..qual_offset + q_encoded_len], q_len, quality_binning);
                            qual_offset += q_encoded_len;
                            record.quality = Some(quality_str);
                        } else {
                            break;
                        }
                    }
                }
            }

            records
        } else {
            // ── Factorize v1: 3 sub-streams (meta, raw_seq, delta_seq) ──
            let mut soff = 0usize;

            let mut read_sub = |label: &str| -> Result<&[u8]> {
                if soff + 8 > sequences.len() {
                    anyhow::bail!("Truncated factorized {} length", label);
                }
                let len = u64::from_le_bytes(sequences[soff..soff + 8].try_into().unwrap()) as usize;
                soff += 8;
                if soff + len > sequences.len() {
                    anyhow::bail!("Truncated factorized {} data (need {}, have {})", label, len, sequences.len() - soff);
                }
                let slice = &sequences[soff..soff + len];
                soff += len;
                Ok(slice)
            };

            let meta_compressed = read_sub("meta")?;
            let raw_seq_compressed = read_sub("raw_seq")?;
            let delta_seq_compressed = read_sub("delta_seq")?;

            let ((header_result, qual_result), (meta_result, (raw_seq_result, delta_seq_result))) = rayon::join(
                || rayon::join(
                    || decompress_headers_bsc(headers, num_reads)
                        .context("Failed to decompress headers"),
                    || -> Result<Vec<u8>> {
                        if qualities_len == 0 {
                            Ok(Vec::new())
                        } else {
                            decompress_qualities_data(qualities, quality_compressor)
                                .context("Failed to decompress qualities")
                        }
                    },
                ),
                || rayon::join(
                    || bsc::decompress_parallel(meta_compressed)
                        .context("Failed to decompress meta stream"),
                    || rayon::join(
                        || bsc::decompress_parallel(raw_seq_compressed)
                            .context("Failed to decompress raw_seq stream"),
                        || if delta_seq_compressed.is_empty() { Ok(Vec::new()) }
                           else { bsc::decompress_parallel(delta_seq_compressed)
                               .context("Failed to decompress delta_seq stream") },
                    ),
                ),
            );

            let read_ids = header_result?;
            let meta_data = meta_result?;
            let raw_seq_data = raw_seq_result?;
            let delta_seq_data = delta_seq_result?;
            let qualities_data = qual_result?;

            info!("Reconstructing factorized sequences...");
            let decoded_sequences = factorize::decode_factorized_sequences(
                &meta_data, &raw_seq_data, &delta_seq_data, num_reads)?;

            let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
                crate::io::FastqRecord::new(id, seq, None)
            }).collect();

            if !qualities_data.is_empty() {
                let mut qual_offset = 0;
                for record in &mut records {
                    if qual_offset < qualities_data.len() {
                        let q_len = read_varint(&qualities_data, &mut qual_offset)
                            .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                        let bits_per_qual = quality_binning.bits_per_quality();
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        if qual_offset + q_encoded_len <= qualities_data.len() {
                            let quality_str = columnar::unpack_qualities(&qualities_data[qual_offset..qual_offset + q_encoded_len], q_len, quality_binning);
                            qual_offset += q_encoded_len;
                            record.quality = Some(quality_str);
                        } else {
                            break;
                        }
                    }
                }
            }

            records
        }
    } else if encoding_type == 8 || encoding_type == 9 {
        // Local reorder + delta (8) or ultra/big-block (9): both store permutation for original order restore
        let mode_name = if encoding_type == 8 { "local-reorder" } else { "ultra" };
        let use_quality_ctx = quality_compressor == QualityCompressor::QualityCtx;
        info!("Decompressing {} sequences (encoding_type={}, quality={})...",
            mode_name, encoding_type, if use_quality_ctx { "quality_ctx" } else { "BSC" });

        // Decode headers and sequences in parallel
        // (quality_ctx needs sequences as context, so can't decompress qualities in parallel)
        let (header_result, seq_result) = rayon::join(
            || decompress_headers_bsc(headers, num_reads)
                .context("Failed to decompress headers"),
            || -> Result<harc::HarcDecoded> {
                if encoding_type == 8 {
                    harc::decode_harc_sequences(sequences, num_reads)
                } else {
                    harc::decode_reorder_local(sequences, num_reads)
                }
            },
        );

        let read_ids = header_result?;
        let decoded = seq_result?;

        // Decode qualities (may need sequences for quality_ctx)
        let mut qual_strings: Vec<Option<String>> = vec![None; num_reads];
        if qualities_len > 0 {
            if use_quality_ctx {
                // Quality_ctx: decompress using sequences as context (in reordered order)
                let ctx_qualities = quality_ctx::decompress_quality_ctx_multiblock(
                    qualities, &decoded.sequences,
                ).context("Failed to decompress qualities (quality_ctx)")?;
                for (i, q) in ctx_qualities.into_iter().enumerate() {
                    qual_strings[i] = Some(q);
                }
            } else {
                // BSC: decompress and parse packed quality bytes
                let qualities_data = decompress_qualities_data(qualities, quality_compressor)
                    .context("Failed to decompress qualities")?;
                let mut qual_offset = 0;
                for i in 0..num_reads {
                    if qual_offset < qualities_data.len() {
                        let q_len = read_varint(&qualities_data, &mut qual_offset)
                            .ok_or_else(|| anyhow::anyhow!("Failed to read quality length at read {}", i))?;
                        let bits_per_qual = quality_binning.bits_per_quality();
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        if qual_offset + q_encoded_len <= qualities_data.len() {
                            let quality_str = columnar::unpack_qualities(
                                &qualities_data[qual_offset..qual_offset + q_encoded_len],
                                q_len, quality_binning,
                            );
                            qual_offset += q_encoded_len;
                            qual_strings[i] = Some(quality_str);
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        // Build records in reordered order, then un-permute to original order
        let mut records: Vec<Option<crate::io::FastqRecord>> = (0..num_reads).map(|_| None).collect();
        for (reorder_pos, (id, seq)) in read_ids.into_iter()
            .zip(decoded.sequences.into_iter())
            .enumerate()
        {
            let orig_idx = decoded.order[reorder_pos] as usize;
            let qual = std::mem::take(&mut qual_strings[reorder_pos]);
            records[orig_idx] = Some(crate::io::FastqRecord { id, sequence: seq, quality: qual });
        }

        records.into_iter().map(|r| r.unwrap()).collect()
    } else if encoding_type == 1 || encoding_type == 2 {
        // Delta or RLE mode: decompress all streams in parallel
        let ((header_result, seq_result), qual_result) = rayon::join(
            || rayon::join(
                || -> Result<Vec<String>> {
                    match header_compressor {
                        HeaderCompressor::Bsc => {
                            if template_prefix_len > 0 {
                                let header_data = bsc::decompress_parallel(headers)
                                    .context("Failed to decompress headers (BSC)")?;
                                read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                                    .context("Failed to decode template-encoded read IDs")
                            } else {
                                decompress_headers_bsc(headers, num_reads)
                                    .context("Failed to decompress headers (BSC)")
                            }
                        }
                        HeaderCompressor::Zstd => {
                            let header_data = decompress_zstd(headers)
                                .context("Failed to decompress headers (zstd)")?;
                            read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                                .context("Failed to decode read IDs")
                        }
                        HeaderCompressor::OpenZl => {
                            if template_prefix_len > 0 {
                                let header_data = openzl::decompress_parallel(headers)
                                    .context("Failed to decompress headers (OpenZL)")?;
                                read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                                    .context("Failed to decode template-encoded read IDs")
                            } else {
                                decompress_headers_openzl(headers, num_reads)
                                    .context("Failed to decompress headers (OpenZL)")
                            }
                        }
                    }
                },
                || match sequence_compressor {
                    SequenceCompressor::Bsc => bsc::decompress_parallel(sequences)
                        .context("Failed to decompress sequences (BSC)"),
                    SequenceCompressor::Zstd => decompress_zstd(sequences)
                        .context("Failed to decompress sequences (zstd)"),
                    SequenceCompressor::OpenZl => openzl::decompress_parallel(sequences)
                        .context("Failed to decompress sequences (OpenZL)"),
                },
            ),
            || if let Some(ref dict) = quality_dict_opt {
                zstd_dict::decompress_with_dict(qualities, dict)
                    .context("Failed to decompress quality scores (dictionary mode)")
            } else {
                decompress_qualities_data(qualities, quality_compressor)
                    .context("Failed to decompress quality scores")
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

                    let quality_str = if quality_compressor == QualityCompressor::Fqzcomp {
                        if q_offset + q_len <= quality_stream.len() {
                            let raw = &quality_stream[q_offset..q_offset + q_len];
                            q_offset += q_len;
                            unsafe { String::from_utf8_unchecked(raw.to_vec()) }
                        } else {
                            break;
                        }
                    } else if let Some(ref model) = quality_model_opt {
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
        let use_quality_ctx = quality_compressor == QualityCompressor::QualityCtx;

        // Step 1: Parallel decompression of headers, sequences, qualities
        // For quality_ctx: skip quality decompression here (needs sequences first)
        let (header_result, (seq_result, qual_result)) = rayon::join(
            || -> Result<Vec<String>> {
                match header_compressor {
                    HeaderCompressor::Bsc => {
                        if template_prefix_len > 0 {
                            let header_data = bsc::decompress_parallel(headers)
                                .context("Failed to decompress headers (BSC)")?;
                            read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                                .context("Failed to decode template-encoded read IDs")
                        } else {
                            decompress_headers_bsc(headers, num_reads)
                                .context("Failed to decompress headers (BSC)")
                        }
                    }
                    HeaderCompressor::Zstd => {
                        let header_data = decompress_zstd(headers)
                            .context("Failed to decompress headers (zstd)")?;
                        read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                            .context("Failed to decode read IDs")
                    }
                    HeaderCompressor::OpenZl => {
                        if template_prefix_len > 0 {
                            let header_data = openzl::decompress_parallel(headers)
                                .context("Failed to decompress headers (OpenZL)")?;
                            read_id::decode_read_ids(&header_data, &read_id_template, num_reads)
                                .context("Failed to decode template-encoded read IDs")
                        } else {
                            decompress_headers_openzl(headers, num_reads)
                                .context("Failed to decompress headers (OpenZL)")
                        }
                    }
                }
            },
            || rayon::join(
                || -> Result<Vec<String>> {
                    match sequence_compressor {
                        SequenceCompressor::Bsc => {
                            if nmasks_len > 0 {
                                decompress_sequences_2bit_bsc(sequences, nmasks, num_reads)
                                    .context("Failed to decompress 2-bit sequences (BSC)")
                            } else {
                                decompress_sequences_raw_bsc(sequences, num_reads, encoding_type)
                                    .context("Failed to decompress sequences (BSC)")
                            }
                        }
                        SequenceCompressor::OpenZl => {
                            decompress_sequences_raw_openzl(sequences, num_reads)
                                .context("Failed to decompress sequences (OpenZL)")
                        }
                        SequenceCompressor::Zstd => {
                            let sequences_data = decompress_zstd(sequences)
                                .context("Failed to decompress sequences (zstd)")?;
                            let nmasks_data = decompress_zstd(nmasks)
                                .context("Failed to decompress N-masks (zstd)")?;
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
                    if use_quality_ctx || qualities_len == 0 {
                        Ok(Vec::new())
                    } else if let Some(ref dict) = quality_dict_opt {
                        zstd_dict::decompress_with_dict(qualities, dict)
                            .context("Failed to decompress quality scores (dictionary mode)")
                    } else {
                        decompress_qualities_data(qualities, quality_compressor)
                            .context("Failed to decompress quality scores")
                    }
                },
            ),
        );

        let read_ids = header_result?;
        let decoded_sequences = seq_result?;
        let qualities_data = qual_result?;

        // For quality_ctx: decompress qualities using sequences as context
        // Must happen after sequences are available (breaks seq||qual parallelism)
        let quality_ctx_strings = if use_quality_ctx {
            info!("Decompressing quality scores (quality_ctx)...");
            Some(quality_ctx::decompress_quality_ctx_multiblock(qualities, &decoded_sequences)
                .context("Failed to decompress quality scores (quality_ctx)")?)
        } else {
            None
        };

        // Step 2: Construct records from decoded headers and sequences
        let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
            crate::io::FastqRecord::new(id, seq, None)
        }).collect();

        // Step 3: Decode quality strings and attach to records
        if let Some(quality_strings) = quality_ctx_strings {
            for (i, qual) in quality_strings.into_iter().enumerate() {
                if i < records.len() {
                    records[i].quality = Some(qual);
                }
            }
        } else if !qualities_data.is_empty() {
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

                        let quality_str = if quality_compressor == QualityCompressor::Fqzcomp {
                            if qual_offset + q_len <= qualities_data.len() {
                                let raw = &qualities_data[qual_offset..qual_offset + q_len];
                                qual_offset += q_len;
                                unsafe { String::from_utf8_unchecked(raw.to_vec()) }
                            } else {
                                break;
                            }
                        } else if let Some(ref model) = quality_model_opt {
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
            IO_BUFFER_SIZE,
            GzEncoder::new(
                std::fs::File::create(output_path)
                    .with_context(|| format!("Failed to create output file: {:?}", output_path))?,
                Compression::new(args.gzip_level),
            ),
        )) as Box<dyn Write>
    } else {
        Box::new(std::io::BufWriter::with_capacity(
            IO_BUFFER_SIZE,
            std::fs::File::create(output_path)
                .with_context(|| format!("Failed to create output file: {:?}", output_path))?,
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
pub(crate) fn compress_sequences_raw_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool, sequence_hints: bool, sequence_delta: bool) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut seq_stream = Vec::new();

    if sequence_delta {
        let mut delta_cache = std::collections::HashMap::new();
        let mut delta_reads = Vec::new();
        let mut delta_count = 0usize;
        for record in records {
            if encode_sequence_delta(record.sequence.as_bytes(), &mut seq_stream, &mut delta_cache, &mut delta_reads)? {
                delta_count += 1;
            }
        }
        info!("Inline delta: {}/{} reads delta-encoded ({:.1}%)",
            delta_count, records.len(), delta_count as f64 / records.len().max(1) as f64 * 100.0);
    } else {
        for record in records {
            write_varint(&mut seq_stream, record.sequence.len())?;
            if sequence_hints {
                seq_stream.push(dna_utils::compute_sequence_hint(record.sequence.as_bytes()));
            }
            seq_stream.extend_from_slice(record.sequence.as_bytes());
        }
    }

    let compressed = bsc_compress_parallel(&seq_stream, bsc_static)?;
    Ok((compressed, Vec::new())) // empty nmasks
}


/// Helper: Decompress raw ASCII BSC-compressed sequences
fn decompress_sequences_raw_bsc(compressed: &[u8], num_reads: usize, encoding_type: u8) -> Result<Vec<String>> {
    let decompressed = bsc::decompress_parallel(compressed)?;
    let mut sequences = Vec::with_capacity(num_reads);
    let mut offset = 0;
    let skip_hints = encoding_type == 4;
    let decode_delta = encoding_type == 5;

    // Delta decoding needs a cache of all previously decoded sequences
    let mut delta_reads: Vec<Vec<u8>> = if decode_delta { Vec::with_capacity(num_reads) } else { Vec::new() };

    for _ in 0..num_reads {
        let seq_len = read_varint(&decompressed, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length varint"))?;

        if decode_delta {
            // Read flag byte
            if offset >= decompressed.len() {
                anyhow::bail!("Truncated delta flag at offset {}", offset);
            }
            let flag = decompressed[offset];
            offset += 1;

            if flag == 0x01 {
                // Delta-encoded: read ref_idx varint, then delta bytes
                let ref_idx = read_varint(&decompressed, &mut offset)
                    .ok_or_else(|| anyhow::anyhow!("Failed to read delta ref_idx varint"))?;
                if ref_idx >= delta_reads.len() {
                    anyhow::bail!("Invalid delta ref_idx {} (only {} reads decoded so far)", ref_idx, delta_reads.len());
                }
                if offset + seq_len > decompressed.len() {
                    anyhow::bail!("Truncated delta sequence data at offset {}", offset);
                }
                let ref_seq = &delta_reads[ref_idx];
                let mut seq = Vec::with_capacity(seq_len);
                for i in 0..seq_len {
                    let byte = decompressed[offset + i];
                    if byte == 0x00 {
                        // Match: copy from reference
                        if i < ref_seq.len() {
                            seq.push(ref_seq[i]);
                        } else {
                            seq.push(b'N');
                        }
                    } else {
                        seq.push(byte); // Mismatch: actual base
                    }
                }
                offset += seq_len;
                delta_reads.push(seq.clone());
                let sequence = String::from_utf8(seq)
                    .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in delta sequence: {}", e))?;
                sequences.push(sequence);
            } else {
                // Raw: read seq_len bytes directly
                if offset + seq_len > decompressed.len() {
                    anyhow::bail!("Truncated sequence data at offset {}", offset);
                }
                let seq_bytes = decompressed[offset..offset + seq_len].to_vec();
                offset += seq_len;
                delta_reads.push(seq_bytes.clone());
                let sequence = String::from_utf8(seq_bytes)
                    .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?;
                sequences.push(sequence);
            }
        } else {
            if skip_hints {
                offset += 1; // skip hint byte
            }
            if offset + seq_len > decompressed.len() {
                anyhow::bail!("Truncated sequence data at offset {}", offset);
            }
            let sequence = String::from_utf8(decompressed[offset..offset + seq_len].to_vec())
                .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?;
            offset += seq_len;
            sequences.push(sequence);
        }
    }

    Ok(sequences)
}

/// Helper: Compress headers with OpenZL
fn compress_headers_openzl(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    let mut header_stream = Vec::new();
    for record in records {
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.extend_from_slice(record.id.as_bytes());
    }
    openzl::compress_parallel(&header_stream)
}

/// Helper: Decompress raw OpenZL-compressed headers
fn decompress_headers_openzl(compressed: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let decompressed = openzl::decompress_parallel(compressed)?;
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

/// Helper: Compress sequences as raw ASCII + OpenZL
fn compress_sequences_raw_openzl(records: &[crate::io::FastqRecord]) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut seq_stream = Vec::new();
    for record in records {
        write_varint(&mut seq_stream, record.sequence.len())?;
        seq_stream.extend_from_slice(record.sequence.as_bytes());
    }

    let compressed = openzl::compress_parallel(&seq_stream)?;
    Ok((compressed, Vec::new())) // empty nmasks
}

/// Helper: Decompress raw ASCII OpenZL-compressed sequences
fn decompress_sequences_raw_openzl(compressed: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let decompressed = openzl::decompress_parallel(compressed)?;
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

/// Helper: Compress sequences as 2-bit encoded + N-mask + BSC
fn compress_sequences_2bit_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool) -> Result<(Vec<u8>, Vec<u8>)> {
    use std::io::Write;

    let mut seq_stream = Vec::new();
    let mut nmask_stream = Vec::new();
    for record in records {
        let encoding = n_mask::encode_with_n_mask(&record.sequence);
        write_varint(&mut seq_stream, encoding.length)?;
        seq_stream.write_all(&encoding.sequence_2bit)?;
        nmask_stream.write_all(&encoding.n_mask)?;
    }

    let (seq_result, nmask_result) = rayon::join(
        || bsc_compress_parallel(&seq_stream, bsc_static),
        || bsc_compress_parallel(&nmask_stream, bsc_static),
    );
    Ok((seq_result?, nmask_result?))
}

/// Helper: Decompress 2-bit BSC-compressed sequences with N-mask
fn decompress_sequences_2bit_bsc(sequences: &[u8], nmasks: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let (seq_data, nmask_data) = rayon::join(
        || bsc::decompress_parallel(sequences),
        || bsc::decompress_parallel(nmasks),
    );
    let seq_data = seq_data?;
    let nmask_data = nmask_data?;

    let mut decoded = Vec::with_capacity(num_reads);
    let mut seq_offset = 0;
    let mut nmask_offset = 0;

    for _ in 0..num_reads {
        let seq_len = read_varint(&seq_data, &mut seq_offset)
            .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length varint"))?;
        let seq_2bit_len = (seq_len + 3) / 4;
        if seq_offset + seq_2bit_len > seq_data.len() {
            anyhow::bail!("Truncated 2-bit sequence data at offset {}", seq_offset);
        }
        let sequence_2bit = &seq_data[seq_offset..seq_offset + seq_2bit_len];
        seq_offset += seq_2bit_len;

        let nmask_len = (seq_len + 7) / 8;
        let n_mask = if nmask_offset + nmask_len <= nmask_data.len() {
            let mask = &nmask_data[nmask_offset..nmask_offset + nmask_len];
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

/// Helper: Compress qualities stream (standard mode)
pub(crate) fn compress_qualities_with(records: &[crate::io::FastqRecord], binning: columnar::QualityBinning, level: i32, compressor: QualityCompressor, bsc_static: bool) -> Result<Vec<u8>> {
    use std::io::Write;

    if compressor == QualityCompressor::Fqzcomp {
        // Wrap in multi-block format (1 block) for consistency with chunked path
        let blob = compress_qualities_fqzcomp(records)?;
        let mut out = Vec::with_capacity(8 + blob.len());
        out.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
        out.extend_from_slice(&(blob.len() as u32).to_le_bytes()); // block_len
        out.extend_from_slice(&blob);
        return Ok(out);
    }

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
        QualityCompressor::OpenZl => openzl::compress_parallel(&quality_stream),
        QualityCompressor::Fqzcomp => unreachable!(),
        QualityCompressor::QualityCtx => unreachable!("quality_ctx handled separately"),
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
        QualityCompressor::OpenZl => {
            info!("Note: OpenZL doesn't support dictionary mode, using standard OpenZL");
            openzl::compress_parallel(&quality_stream)
        }
        QualityCompressor::Fqzcomp => {
            info!("Note: fqzcomp doesn't support dictionary mode, using standard fqzcomp");
            compress_qualities_fqzcomp(records)
        }
        QualityCompressor::QualityCtx => unreachable!("quality_ctx handled separately"),
    }
}

/// Helper: Compress qualities stream with positional modeling
fn compress_qualities_with_model(
    records: &[crate::io::FastqRecord],
    level: i32,
    quality_compressor: QualityCompressor,
    bsc_static: bool,
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

    // Compress delta stream
    let compressed = match quality_compressor {
        QualityCompressor::Bsc => bsc_compress_parallel(&quality_stream, bsc_static)?,
        QualityCompressor::Zstd => compress_zstd(&quality_stream, level)?,
        QualityCompressor::OpenZl => openzl::compress_parallel(&quality_stream)?,
        QualityCompressor::Fqzcomp => {
            info!("Note: fqzcomp doesn't support quality modeling, using BSC for model deltas");
            bsc_compress_parallel(&quality_stream, bsc_static)?
        }
        QualityCompressor::QualityCtx => unreachable!("quality_ctx handled separately"),
    };

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

/// Maximum varint size in bytes (10 bytes = up to 70 bits, enough for usize on 64-bit)
const MAX_VARINT_BYTES: usize = 10;

/// Read variable-length integer
fn read_varint(data: &[u8], offset: &mut usize) -> Option<usize> {
    let mut value = 0usize;
    let mut shift = 0;
    let mut bytes_read = 0;

    loop {
        if *offset >= data.len() {
            return None;
        }

        let byte = data[*offset];
        *offset += 1;
        bytes_read += 1;

        value |= ((byte & 0x7F) as usize) << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
        if bytes_read >= MAX_VARINT_BYTES {
            return None; // Malformed varint: too many continuation bytes
        }
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
        QualityCompressor::OpenZl => openzl::decompress_parallel(compressed),
        QualityCompressor::Fqzcomp => decompress_qualities_fqzcomp_multiblock(compressed),
        QualityCompressor::QualityCtx => unreachable!("quality_ctx decompression handled separately"),
    }
}

/// Decompress multi-block fqzcomp quality data.
///
/// Format: [num_blocks: u32][block_len: u32, block_data]...
/// Each block is an independent fqzcomp blob (from compress_qualities_fqzcomp).
/// Format: [num_blocks: u32][block_len: u32, block_data]...
fn decompress_qualities_fqzcomp_multiblock(compressed: &[u8]) -> Result<Vec<u8>> {
    if compressed.len() < 4 {
        anyhow::bail!("fqzcomp multi-block: data too short");
    }

    let num_blocks = u32::from_le_bytes(compressed[0..4].try_into().unwrap()) as usize;

    if num_blocks == 0 {
        return Ok(Vec::new());
    }

    decompress_fqzcomp_blocks(compressed, num_blocks)
}

fn decompress_fqzcomp_blocks(compressed: &[u8], num_blocks: usize) -> Result<Vec<u8>> {
    let mut offset = 4; // skip num_blocks
    let mut output = Vec::new();

    for i in 0..num_blocks {
        if offset + 4 > compressed.len() {
            anyhow::bail!("fqzcomp multi-block: truncated block {} header", i);
        }
        let block_len = u32::from_le_bytes(compressed[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;
        if offset + block_len > compressed.len() {
            anyhow::bail!("fqzcomp multi-block: truncated block {} data ({} bytes, {} available)", i, block_len, compressed.len() - offset);
        }
        let block_data = &compressed[offset..offset + block_len];
        offset += block_len;

        let decompressed = decompress_qualities_fqzcomp(block_data)?;
        output.extend_from_slice(&decompressed);
    }

    Ok(output)
}

/// Compress quality scores using mean-quality reordering + fqzcomp.
///
/// Strategy: sort reads by mean quality score (1 byte key per read), then
/// compress reordered qualities with fqzcomp strat=0. Grouping reads by
/// similar quality profiles improves fqzcomp's context model adaptation.
///
/// Archive format: [num_reads: u32][sort_keys_bsc_len: u32][sort_keys_bsc]
///                 [num_sub_chunks: u32][sub_chunk_len: u32, sub_chunk_data]...
///
/// Qualities are sorted by mean quality, then split into sub-chunks that are
/// compressed with fqzcomp in parallel via rayon.
fn compress_qualities_fqzcomp(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    const FQZCOMP_SUB_CHUNK: usize = 500_000;

    let qual_strs: Vec<&str> = records.iter()
        .filter_map(|r| r.quality.as_deref())
        .collect();

    if qual_strs.is_empty() {
        return Ok(Vec::new());
    }

    // Compute mean quality per read (1 byte each)
    let sort_keys: Vec<u8> = qual_strs.iter()
        .map(|q| {
            let sum: u64 = q.bytes().map(|b| b as u64).sum();
            (sum / q.len().max(1) as u64) as u8
        })
        .collect();

    // Sort by mean quality (stable sort)
    let mut indices: Vec<usize> = (0..qual_strs.len()).collect();
    indices.sort_by_key(|&i| sort_keys[i]);

    // Reorder quality strings
    let sorted_quals: Vec<&str> = indices.iter().map(|&i| qual_strs[i]).collect();

    let num_sub_chunks = (sorted_quals.len() + FQZCOMP_SUB_CHUNK - 1) / FQZCOMP_SUB_CHUNK;
    info!("fqzcomp: {} reads, sorting by mean quality, {} sub-chunks, compressing with strat=0",
        qual_strs.len(), num_sub_chunks);

    // Compress sub-chunks in parallel with rayon, and BSC-compress sort keys concurrently
    let (sub_chunk_blobs, sort_keys_bsc) = rayon::join(
        || -> Result<Vec<Vec<u8>>> {
            sorted_quals.par_chunks(FQZCOMP_SUB_CHUNK)
                .map(|chunk| fqzcomp::compress(chunk, 0))
                .collect::<Result<Vec<_>>>()
        },
        || bsc::compress_parallel_adaptive(&sort_keys),
    );
    let sub_chunk_blobs = sub_chunk_blobs?;
    let sort_keys_bsc = sort_keys_bsc?;

    let fqzcomp_total: usize = sub_chunk_blobs.iter().map(|b| b.len()).sum();
    info!("fqzcomp: sort keys {} B (BSC'd), fqzcomp {} B in {} sub-chunks",
        sort_keys_bsc.len(), fqzcomp_total, sub_chunk_blobs.len());

    // Pack: [num_reads][sort_keys_bsc_len][sort_keys_bsc][num_sub_chunks][sub_chunk_len, data]...
    let mut output = Vec::with_capacity(12 + sort_keys_bsc.len() + fqzcomp_total + sub_chunk_blobs.len() * 4);
    output.extend_from_slice(&(qual_strs.len() as u32).to_le_bytes());
    output.extend_from_slice(&(sort_keys_bsc.len() as u32).to_le_bytes());
    output.extend_from_slice(&sort_keys_bsc);
    output.extend_from_slice(&(sub_chunk_blobs.len() as u32).to_le_bytes());
    for blob in &sub_chunk_blobs {
        output.extend_from_slice(&(blob.len() as u32).to_le_bytes());
        output.extend_from_slice(blob);
    }

    Ok(output)
}

/// Decompress fqzcomp quality data with mean-quality unsort.
///
/// Format: [num_reads: u32][sort_keys_bsc_len: u32][sort_keys_bsc]
///         [num_sub_chunks: u32][sub_chunk_len: u32, sub_chunk_data]...
///
/// Returns varint-prefixed raw ASCII quality stream (same format as other compressors
/// with QualityBinning::None) so the existing quality parsing code works unchanged.
fn decompress_qualities_fqzcomp(compressed: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    const FQZCOMP_SUB_CHUNK: usize = 500_000;

    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    if compressed.len() < 12 {
        anyhow::bail!("fqzcomp quality data too short");
    }

    let mut offset = 0;

    // Parse header
    let num_reads = u32::from_le_bytes(compressed[offset..offset + 4].try_into().unwrap()) as usize;
    offset += 4;
    let sort_keys_len = u32::from_le_bytes(compressed[offset..offset + 4].try_into().unwrap()) as usize;
    offset += 4;

    info!("fqzcomp decompress: {} reads, sort_keys_bsc {} bytes, total blob {} bytes",
        num_reads, sort_keys_len, compressed.len());

    if offset + sort_keys_len > compressed.len() {
        anyhow::bail!("fqzcomp: truncated sort keys");
    }

    // BSC-decompress sort keys
    let sort_keys_data = &compressed[offset..offset + sort_keys_len];
    offset += sort_keys_len;

    // Parse sub-chunk table
    if offset + 4 > compressed.len() {
        anyhow::bail!("fqzcomp: truncated sub-chunk header");
    }
    let num_sub_chunks = u32::from_le_bytes(compressed[offset..offset + 4].try_into().unwrap()) as usize;
    offset += 4;

    // Collect sub-chunk slices
    let mut sub_chunk_slices: Vec<&[u8]> = Vec::with_capacity(num_sub_chunks);
    for i in 0..num_sub_chunks {
        if offset + 4 > compressed.len() {
            anyhow::bail!("fqzcomp: truncated sub-chunk {} header", i);
        }
        let chunk_len = u32::from_le_bytes(compressed[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;
        if offset + chunk_len > compressed.len() {
            anyhow::bail!("fqzcomp: truncated sub-chunk {} data", i);
        }
        sub_chunk_slices.push(&compressed[offset..offset + chunk_len]);
        offset += chunk_len;
    }

    // Compute reads per sub-chunk (must match compression side)
    let reads_per_sub: Vec<usize> = (0..num_sub_chunks).map(|i| {
        let start = i * FQZCOMP_SUB_CHUNK;
        FQZCOMP_SUB_CHUNK.min(num_reads - start)
    }).collect();

    info!("fqzcomp decompress: {} sub-chunks, reads per sub: {:?}", num_sub_chunks, &reads_per_sub);

    // Decompress sort keys and fqzcomp sub-chunks in parallel
    let (sort_keys_result, sub_results): (Result<Vec<u8>>, Vec<Result<(Vec<u8>, Vec<i32>)>>) = rayon::join(
        || bsc::decompress_parallel(sort_keys_data),
        || sub_chunk_slices.par_iter().zip(reads_per_sub.par_iter())
            .map(|(data, &n_reads)| fqzcomp::decompress(data, n_reads))
            .collect(),
    );
    let sort_keys = sort_keys_result?;

    // Concatenate decoded qualities and lengths from all sub-chunks
    let mut decoded_quals = Vec::new();
    let mut lengths: Vec<i32> = Vec::with_capacity(num_reads);
    for res in sub_results {
        let (quals, lens) = res?;
        decoded_quals.extend_from_slice(&quals);
        lengths.extend_from_slice(&lens);
    }

    // Reconstruct permutation (same stable sort as compression)
    let mut indices: Vec<usize> = (0..num_reads).collect();
    indices.sort_by_key(|&i| sort_keys[i]);

    // Split decoded qualities into individual reads
    let mut sorted_quals: Vec<&[u8]> = Vec::with_capacity(num_reads);
    let mut q_offset = 0;
    for i in 0..num_reads {
        let len = lengths[i] as usize;
        if q_offset + len > decoded_quals.len() {
            anyhow::bail!("fqzcomp: truncated quality data at read {}", i);
        }
        sorted_quals.push(&decoded_quals[q_offset..q_offset + len]);
        q_offset += len;
    }

    // Inverse permute: unsorted[indices[i]] = sorted[i]
    let mut unsorted_quals: Vec<&[u8]> = vec![&[]; num_reads];
    for (sorted_idx, &original_idx) in indices.iter().enumerate() {
        unsorted_quals[original_idx] = sorted_quals[sorted_idx];
    }

    // Build varint-prefixed stream (raw ASCII, compatible with QualityBinning::None)
    let total_est: usize = unsorted_quals.iter().map(|q| q.len() + 2).sum();
    let mut output = Vec::with_capacity(total_est);
    for qual in &unsorted_quals {
        write_varint(&mut output, qual.len())?;
        output.extend_from_slice(qual);
    }

    Ok(output)
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
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8)
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
                record.quality = Some(quality::quantize_quality(qual, args.quality_mode).into_owned());
            }
        }
        for record in &mut r2_records {
            if let Some(qual) = &record.quality {
                record.quality = Some(quality::quantize_quality(qual, args.quality_mode).into_owned());
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
    let r2_diff_counts_bytes: Vec<u8> = r2_diff_counts.iter()
        .flat_map(|&c| c.to_le_bytes())
        .collect();
    let r2_diff_counts_compressed = compress_zstd(&r2_diff_counts_bytes, args.compression_level)?;
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
    output_file.write_all(b"QZ\0\0")?;  // Magic
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
