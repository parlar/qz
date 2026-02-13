//! Stream-level compression strategy benchmarking
//!
//! Systematically tests every compression strategy and variant on each
//! FASTQ stream type (headers, sequences, qualities), producing ranked
//! results for comparison.
//!
//! Usage via CLI: `fqz benchmark -i input.fastq`

use crate::cli::QualityCompressor;
use crate::io::FastqRecord;
use anyhow::Result;
use std::io::Write;
use std::time::Instant;
use super::arithmetic_quality;
use super::arithmetic_sequence;
use super::bsc;
use super::columnar::{self, QualityBinning};
use super::delta;
use super::n_mask;
use super::quality_delta;
use super::quality_model;
use super::read_id;
use super::rle;
use super::debruijn;
use super::zstd_dict;

/// Result of benchmarking a single strategy on a single stream.
pub struct StrategyResult {
    #[allow(dead_code)]
    pub stream: &'static str,
    pub strategy: String,
    pub original_bytes: usize,
    pub compressed_bytes: usize,
    pub ratio: f64,
    pub bits_per_element: f64,
    pub compress_time_ms: f64,
    pub lossy: bool,
}

/// Result of compressing the full FASTQ file with a specific tool.
pub struct WholeFileResult {
    pub tool: String,
    pub original_bytes: usize,
    pub compressed_bytes: usize,
    pub ratio: f64,
    pub compress_time_ms: f64,
    pub notes: String,
}

/// Whole-file benchmark comparing FQZ against other compressors.
pub struct WholeFileBenchmark {
    pub num_reads: usize,
    pub original_fastq_bytes: usize,
    pub results: Vec<WholeFileResult>,
}

/// Aggregated benchmark report for all three streams.
pub struct BenchmarkReport {
    pub num_reads: usize,
    pub total_bases: usize,
    pub total_quality_scores: usize,
    pub header_results: Vec<StrategyResult>,
    pub sequence_results: Vec<StrategyResult>,
    pub quality_results: Vec<StrategyResult>,
    pub whole_file: Option<WholeFileBenchmark>,
}

/// Run all benchmarks on the given FASTQ records.
/// If `stream_filter` is Some, only benchmark that stream (e.g. "sequences").
pub fn run_benchmark(records: &[FastqRecord], stream_filter: Option<&str>) -> Result<BenchmarkReport> {
    let sequences: Vec<String> = records.iter().map(|r| r.sequence.clone()).collect();
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();

    let run_headers = stream_filter.is_none_or(|s| s == "headers");
    let run_sequences = stream_filter.is_none_or(|s| s == "sequences");
    let run_qualities = stream_filter.is_none_or(|s| s == "qualities");
    let run_whole = stream_filter.is_none();

    let headers: Vec<String> = if run_headers {
        records.iter().map(|r| r.id.clone()).collect()
    } else {
        Vec::new()
    };
    let qualities: Vec<String> = if run_qualities {
        records.iter().filter_map(|r| r.quality.clone()).collect()
    } else {
        Vec::new()
    };
    let total_quality_scores: usize = qualities.iter().map(|q| q.len()).sum();

    // Run requested stream benchmarks in parallel
    let ((header_results, sequence_results), (quality_results, whole_file)) = rayon::join(
        || rayon::join(
            || if run_headers { benchmark_headers(&headers) } else { Vec::new() },
            || if run_sequences { benchmark_sequences(&sequences) } else { Vec::new() },
        ),
        || rayon::join(
            || if run_qualities { benchmark_qualities(&sequences, &qualities) } else { Vec::new() },
            || if run_whole {
                match run_whole_file_benchmark(records) {
                    Ok(wf) => Some(wf),
                    Err(e) => {
                        eprintln!("Warning: whole-file benchmark failed: {}", e);
                        None
                    }
                }
            } else { None },
        ),
    );

    Ok(BenchmarkReport {
        num_reads: records.len(),
        total_bases,
        total_quality_scores,
        header_results,
        sequence_results,
        quality_results,
        whole_file,
    })
}

// ---------------------------------------------------------------------------
// Header benchmarks
// ---------------------------------------------------------------------------

fn benchmark_headers(headers: &[String]) -> Vec<StrategyResult> {
    let raw_data: Vec<u8> = headers.join("\n").into_bytes();
    let original_bytes = raw_data.len();
    let num_elements = headers.len();
    let mut results = Vec::new();

    // Raw + zstd at various levels
    for level in [1, 3, 9, 19] {
        if let Some(r) = bench_zstd("headers", &format!("Raw + zstd {}", level), &raw_data, original_bytes, num_elements, level, false) {
            results.push(r);
        }
    }

    // Raw + gzip
    if let Some(r) = bench_gzip("headers", "Raw + gzip", &raw_data, original_bytes, num_elements, false) {
        results.push(r);
    }

    // Raw + bsc
    if let Some(r) = bench_bsc("headers", "Raw + bsc", &raw_data, original_bytes, num_elements, false) {
        results.push(r);
    }

    // Raw + bsc adaptive
    {
        let start = Instant::now();
        if let Ok(compressed) = bsc::compress_adaptive(&raw_data) {
            let elapsed = start.elapsed().as_secs_f64() * 1000.0;
            results.push(StrategyResult {
                stream: "headers",
                strategy: "Raw + bsc adaptive".to_string(),
                original_bytes,
                compressed_bytes: compressed.len(),
                ratio: original_bytes as f64 / compressed.len() as f64,
                bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                compress_time_ms: elapsed,
                lossy: false,
            });
        }
    }

    // Template encoding + zstd at various levels
    if let Ok(encoded) = read_id::compress_read_ids(headers) {
        let template_data = &encoded.encoded_data;
        // Template metadata overhead: prefix bytes + flags
        let template_overhead = encoded.template.prefix.len() + 2;

        for level in [1, 3, 9, 19] {
            let start = Instant::now();
            if let Ok(compressed) = zstd::bulk::compress(template_data, level) {
                let total = compressed.len() + template_overhead;
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.push(StrategyResult {
                    stream: "headers",
                    strategy: format!("Template + zstd {}", level),
                    original_bytes,
                    compressed_bytes: total,
                    ratio: original_bytes as f64 / total as f64,
                    bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        }

        // Template + zstd with dictionary
        let header_bytes: Vec<&[u8]> = headers.iter().map(|h| h.as_bytes()).collect();
        if let Ok(dict) = zstd_dict::train_dictionary(&header_bytes, 8192) {
            let start = Instant::now();
            if let Ok(compressed) = zstd_dict::compress_with_dict(template_data, &dict, 19) {
                let total = compressed.len() + template_overhead + dict.len();
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.push(StrategyResult {
                    stream: "headers",
                    strategy: "Template + zstd 19 + dict".to_string(),
                    original_bytes,
                    compressed_bytes: total,
                    ratio: original_bytes as f64 / total as f64,
                    bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        }
    }

    results.sort_by(|a, b| b.ratio.partial_cmp(&a.ratio).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// ---------------------------------------------------------------------------
// Sequence benchmarks
// ---------------------------------------------------------------------------

fn benchmark_sequences(sequences: &[String]) -> Vec<StrategyResult> {
    let raw_data: Vec<u8> = sequences.concat().into_bytes();
    let original_bytes = raw_data.len();
    let num_elements = original_bytes; // bits per base
    let results = std::sync::Mutex::new(Vec::new());
    let res = &results;

    // Pre-compute shared encodings outside the parallel scope
    let two_bit = encode_2bit(sequences);
    let three_bit = encode_3bit(sequences);

    rayon::scope(|s| {
        // Raw + zstd 3
        s.spawn(|_| {
            if let Some(r) = bench_zstd("sequences", "Raw + zstd 3", &raw_data, original_bytes, num_elements, 3, false) {
                res.lock().unwrap().push(r);
            }
        });
        // Raw + zstd 19
        s.spawn(|_| {
            if let Some(r) = bench_zstd("sequences", "Raw + zstd 19", &raw_data, original_bytes, num_elements, 19, false) {
                res.lock().unwrap().push(r);
            }
        });
        // Raw + gzip
        s.spawn(|_| {
            if let Some(r) = bench_gzip("sequences", "Raw + gzip", &raw_data, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });
        // Raw + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("sequences", "Raw + bsc", &raw_data, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });
        // 2-bit + zstd levels
        for level in [3, 9, 19] {
            let two_bit = &two_bit;
            s.spawn(move |_| {
                let start = Instant::now();
                if let Ok(compressed) = zstd::bulk::compress(two_bit, level) {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: format!("2-bit + zstd {}", level),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            });
        }
        // 2-bit + N-mask + zstd
        for level in [3, 19] {
            s.spawn(move |_| {
                let start = Instant::now();
                let mut seq_stream = Vec::new();
                let mut nmask_stream = Vec::new();
                for seq in sequences {
                    let enc = n_mask::encode_with_n_mask(seq);
                    seq_stream.extend_from_slice(&enc.sequence_2bit);
                    nmask_stream.extend_from_slice(&enc.n_mask);
                }
                let encoding_time = start.elapsed().as_secs_f64() * 1000.0;
                let start2 = Instant::now();
                let seq_compressed = zstd::bulk::compress(&seq_stream, level);
                let nmask_compressed = zstd::bulk::compress(&nmask_stream, 19);
                if let (Ok(sc), Ok(nc)) = (seq_compressed, nmask_compressed) {
                    let total = sc.len() + nc.len();
                    let elapsed = encoding_time + start2.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: format!("2-bit + N-mask + zstd {}", level),
                        original_bytes,
                        compressed_bytes: total,
                        ratio: original_bytes as f64 / total as f64,
                        bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            });
        }
        // Delta encoding + zstd
        s.spawn(|_| {
            let start = Instant::now();
            let encoded = delta::apply_delta_encoding(sequences);
            let mut stream = Vec::new();
            for enc in &encoded {
                let _ = stream.write_all(&(enc.len() as u32).to_le_bytes());
                let _ = stream.write_all(enc);
            }
            if let Ok(compressed) = zstd::bulk::compress(&stream, 3) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Delta + zstd 3".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // RLE encoding + zstd
        s.spawn(|_| {
            let start = Instant::now();
            let encoded = rle::apply_rle_to_sequences(sequences);
            let mut stream = Vec::new();
            for enc in &encoded {
                let _ = stream.write_all(&(enc.len() as u32).to_le_bytes());
                let _ = stream.write_all(enc);
            }
            if let Ok(compressed) = zstd::bulk::compress(&stream, 3) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "RLE + zstd 3".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Arithmetic (ANS) coding
        s.spawn(|_| {
            let start = Instant::now();
            match arithmetic_sequence::encode_sequences_arithmetic(sequences) {
                Ok(compressed) => {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: "Arithmetic (ANS)".to_string(),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
                Err(e) => {
                    eprintln!("  [skip] Arithmetic sequence coding failed: {}", e);
                }
            }
        });
        // Neural context mixing (GeCo3-style) — skipped, too slow for large datasets
        // To re-enable: super::neural_sequence::compress_sequences(sequences)
        // 3-bit bitpacked + zstd levels
        for level in [3, 19] {
            let three_bit = &three_bit;
            s.spawn(move |_| {
                let start = Instant::now();
                if let Ok(compressed) = zstd::bulk::compress(three_bit, level) {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: format!("3-bit bitpack + zstd {}", level),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            });
        }
        // 3-bit bitpacked + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("sequences", "3-bit bitpack + bsc", &three_bit, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });
        // 2-bit + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("sequences", "2-bit + bsc", &two_bit, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });
        // 2-bit + N-mask + bsc
        s.spawn(|_| {
            let start = Instant::now();
            let mut seq_stream = Vec::new();
            let mut nmask_stream = Vec::new();
            for seq in sequences {
                let enc = n_mask::encode_with_n_mask(seq);
                seq_stream.extend_from_slice(&enc.sequence_2bit);
                nmask_stream.extend_from_slice(&enc.n_mask);
            }
            let encoding_time = start.elapsed().as_secs_f64() * 1000.0;
            let start2 = Instant::now();
            let seq_compressed = bsc::compress(&seq_stream);
            let nmask_compressed = bsc::compress(&nmask_stream);
            if let (Ok(sc), Ok(nc)) = (seq_compressed, nmask_compressed) {
                let total = sc.len() + nc.len();
                let elapsed = encoding_time + start2.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "2-bit + N-mask + bsc".to_string(),
                    original_bytes,
                    compressed_bytes: total,
                    ratio: original_bytes as f64 / total as f64,
                    bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Context MTF (order 8) → 2-bit bitpack → BSC
        s.spawn(|_| {
            let start = Instant::now();
            let mtf = sequence_context_mtf(sequences, 8);
            let bitpacked = bitpack_values(&mtf, 2);
            if let Ok(compressed) = bsc::compress(&bitpacked) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Context MTF-8 + 2-bit → BSC".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Context MTF (order 6) → 2-bit bitpack → BSC
        s.spawn(|_| {
            let start = Instant::now();
            let mtf = sequence_context_mtf(sequences, 6);
            let bitpacked = bitpack_values(&mtf, 2);
            if let Ok(compressed) = bsc::compress(&bitpacked) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Context MTF-6 + 2-bit → BSC".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Context MTF (order 8) → 2-bit bitpack → zstd 19
        s.spawn(|_| {
            let start = Instant::now();
            let mtf = sequence_context_mtf(sequences, 8);
            let bitpacked = bitpack_values(&mtf, 2);
            if let Ok(compressed) = zstd::bulk::compress(&bitpacked, 19) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Context MTF-8 + 2-bit + zstd 19".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Context MTF (order 8) raw → BSC
        s.spawn(|_| {
            let start = Instant::now();
            let mtf = sequence_context_mtf(sequences, 8);
            if let Ok(compressed) = bsc::compress(&mtf) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                results.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Context MTF-8 raw → BSC".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Raw + BSC adaptive (LZP + adaptive QLFC coder)
        s.spawn(|_| {
            let start = Instant::now();
            if let Ok(compressed) = bsc::compress_adaptive(&raw_data) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Raw + bsc adaptive".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Raw + xz (LZMA)
        s.spawn(|_| {
            use xz2::write::XzEncoder;
            let start = Instant::now();
            let mut encoder = XzEncoder::new(Vec::new(), 9);
            if encoder.write_all(&raw_data).is_ok() {
                if let Ok(compressed) = encoder.finish() {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: "Raw + xz 9 (LZMA)".to_string(),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            }
        });
        // 2-bit + xz (LZMA)
        s.spawn(|_| {
            use xz2::write::XzEncoder;
            let start = Instant::now();
            let mut encoder = XzEncoder::new(Vec::new(), 9);
            if encoder.write_all(&two_bit).is_ok() {
                if let Ok(compressed) = encoder.finish() {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: "2-bit + xz 9 (LZMA)".to_string(),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            }
        });
        // Raw + zstd 19 with trained dictionary
        s.spawn(|_| {
            // Train dictionary on individual sequences
            let seq_bytes: Vec<&[u8]> = sequences.iter().map(|s| s.as_bytes()).collect();
            if let Ok(dict) = zstd_dict::train_dictionary(&seq_bytes, 16384) {
                let start = Instant::now();
                if let Ok(compressed) = zstd_dict::compress_with_dict(&raw_data, &dict, 19) {
                    let total = compressed.len() + dict.len();
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: "Raw + zstd 19 + dict".to_string(),
                        original_bytes,
                        compressed_bytes: total,
                        ratio: original_bytes as f64 / total as f64,
                        bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            }
        });
        // 2-bit + BSC adaptive
        s.spawn(|_| {
            let start = Instant::now();
            if let Ok(compressed) = bsc::compress_adaptive(&two_bit) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "2-bit + bsc adaptive".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Column-major transpose + BSC (exploits positional base bias)
        s.spawn(|_| {
            let start = Instant::now();
            let transposed = transpose_sequences_column_major(sequences);
            if let Ok(compressed) = bsc::compress(&transposed) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Column-major + BSC".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Column-major transpose + zstd 19
        s.spawn(|_| {
            let start = Instant::now();
            let transposed = transpose_sequences_column_major(sequences);
            if let Ok(compressed) = zstd::bulk::compress(&transposed, 19) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Column-major + zstd 19".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // Column-major 2-bit transpose + BSC
        s.spawn(|_| {
            let start = Instant::now();
            let transposed = transpose_sequences_column_major(sequences);
            // 2-bit encode the transposed data
            let two_bit_transposed = encode_2bit_raw(&transposed);
            if let Ok(compressed) = bsc::compress(&two_bit_transposed) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "sequences",
                    strategy: "Column-major 2-bit + BSC".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });
        // De Bruijn graph — auto k (adaptive)
        s.spawn(|_| {
            let start = Instant::now();
            match debruijn::compress_sequences_debruijn(sequences, 0) {
                Ok(compressed) => {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "sequences",
                        strategy: "De Bruijn (auto-k) + BSC".to_string(),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
                Err(e) => {
                    eprintln!("  [skip] De Bruijn (auto-k) failed: {}", e);
                }
            }
        });
        // De Bruijn graph — sweep k values for comparison
        for k in [15, 19, 23, 31] {
            s.spawn(move |_| {
                let start = Instant::now();
                match debruijn::compress_sequences_debruijn(sequences, k) {
                    Ok(compressed) => {
                        let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                        res.lock().unwrap().push(StrategyResult {
                            stream: "sequences",
                            strategy: format!("De Bruijn (k={}) + BSC", k),
                            original_bytes,
                            compressed_bytes: compressed.len(),
                            ratio: original_bytes as f64 / compressed.len() as f64,
                            bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                            compress_time_ms: elapsed,
                            lossy: false,
                        });
                    }
                    Err(e) => {
                        eprintln!("  [skip] De Bruijn (k={}) failed: {}", k, e);
                    }
                }
            });
        }
    });

    let mut results = results.into_inner().unwrap();
    results.sort_by(|a, b| b.ratio.partial_cmp(&a.ratio).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// ---------------------------------------------------------------------------
// Quality benchmarks
// ---------------------------------------------------------------------------

fn benchmark_qualities(sequences: &[String], qualities: &[String]) -> Vec<StrategyResult> {
    let raw_data: Vec<u8> = qualities.concat().into_bytes();
    let original_bytes = raw_data.len();
    let num_elements = original_bytes; // bits per quality score
    let results = std::sync::Mutex::new(Vec::new());
    let res = &results;

    // Pre-compute shared data outside the parallel scope
    let packed_none = pack_all_qualities(qualities, QualityBinning::None);
    let (bitpacked_auto, auto_bits) = bitpack_qualities_auto(qualities);
    let packed_illumina8 = pack_all_qualities(qualities, QualityBinning::Illumina8);
    let (bitpacked_illumina8, illumina8_bits) = bitpack_binned_qualities(qualities, QualityBinning::Illumina8);
    let packed_fourlevel = pack_all_qualities(qualities, QualityBinning::FourLevel);
    let (bitpacked_fourlevel, fourlevel_bits) = bitpack_binned_qualities(qualities, QualityBinning::FourLevel);
    let packed_binary = pack_all_qualities(qualities, QualityBinning::Binary { threshold: 20 });
    let (bitpacked_binary, binary_bits) = bitpack_binned_qualities(qualities, QualityBinning::Binary { threshold: 20 });

    // Create references for use in move closures (prevents ownership transfer)
    let raw_ref = &raw_data;
    let packed_none_ref = &packed_none;
    let bitpacked_auto_ref = &bitpacked_auto;
    let packed_illumina8_ref = &packed_illumina8;
    let bitpacked_illumina8_ref = &bitpacked_illumina8;
    let packed_fourlevel_ref = &packed_fourlevel;
    let bitpacked_fourlevel_ref = &bitpacked_fourlevel;
    let packed_binary_ref = &packed_binary;
    let bitpacked_binary_ref = &bitpacked_binary;

    rayon::scope(|s| {
        // --- Lossless strategies ---

        // Raw + zstd
        for level in [3, 9, 19] {
            s.spawn(move |_| {
                if let Some(r) = bench_zstd("qualities", &format!("Raw + zstd {}", level), raw_ref, original_bytes, num_elements, level, false) {
                    res.lock().unwrap().push(r);
                }
            });
        }

        // Raw + gzip
        s.spawn(|_| {
            if let Some(r) = bench_gzip("qualities", "Raw + gzip", raw_ref, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });

        // Raw + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("qualities", "Raw + bsc", raw_ref, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });

        // Raw + bsc adaptive
        s.spawn(|_| {
            let start = Instant::now();
            if let Ok(compressed) = bsc::compress_adaptive(raw_ref) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "qualities",
                    strategy: "Raw + bsc adaptive".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });

        // Packed 7-bit (lossless) + zstd
        for level in [3, 19] {
            s.spawn(move |_| {
                if let Some(r) = bench_zstd("qualities", &format!("Packed 7-bit + zstd {}", level), packed_none_ref, original_bytes, num_elements, level, false) {
                    res.lock().unwrap().push(r);
                }
            });
        }

        // Positional model + zstd
        s.spawn(|_| {
            if let Ok(model) = quality_model::build_quality_model(qualities) {
                let start = Instant::now();
                let mut model_stream = Vec::new();
                for qual in qualities {
                    let deltas = quality_model::encode_with_model(qual, &model);
                    let packed = quality_model::pack_deltas(&deltas);
                    model_stream.extend_from_slice(&packed);
                }
                let encoding_time = start.elapsed().as_secs_f64() * 1000.0;
                let model_overhead = quality_model::serialize_model(&model).len();

                for level in [3, 19] {
                    let start = Instant::now();
                    if let Ok(compressed) = zstd::bulk::compress(&model_stream, level) {
                        let total = compressed.len() + model_overhead;
                        let elapsed = encoding_time + start.elapsed().as_secs_f64() * 1000.0;
                        res.lock().unwrap().push(StrategyResult {
                            stream: "qualities",
                            strategy: format!("Positional model + zstd {}", level),
                            original_bytes,
                            compressed_bytes: total,
                            ratio: original_bytes as f64 / total as f64,
                            bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                            compress_time_ms: elapsed,
                            lossy: false,
                        });
                    }
                }
            }
        });

        // Inter-read delta + zstd
        s.spawn(|_| {
            let start = Instant::now();
            let encoded_deltas = quality_delta::encode_quality_deltas(qualities);
            let mut delta_stream = Vec::new();
            for deltas in &encoded_deltas {
                let packed = quality_delta::pack_deltas(deltas);
                delta_stream.extend_from_slice(&packed);
            }
            if let Ok(compressed) = zstd::bulk::compress(&delta_stream, 3) {
                let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "qualities",
                    strategy: "Inter-read delta + zstd 3".to_string(),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });

        // Arithmetic (ANS) coding — skipped, too slow for large datasets (per-symbol model creation)

        // Context-adaptive rANS — skipped, too slow for large datasets (per-symbol model creation)

        // rANS → BSC — skipped, too slow for large datasets


        // --- Bitpacked lossless strategies ---

        // Auto-width bitpacked + zstd
        for level in [3, 19] {
            s.spawn(move |_| {
                let start = Instant::now();
                if let Ok(compressed) = zstd::bulk::compress(bitpacked_auto_ref, level) {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "qualities",
                        strategy: format!("Bitpack {}b + zstd {}", auto_bits, level),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: false,
                    });
                }
            });
        }

        // Auto-width bitpacked + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("qualities", &format!("Bitpack {}b + bsc", auto_bits), bitpacked_auto_ref, original_bytes, num_elements, false) {
                res.lock().unwrap().push(r);
            }
        });

        // Bitpacked positional model deltas + zstd
        s.spawn(|_| {
            if let Ok(model) = quality_model::build_quality_model(qualities) {
                let start = Instant::now();
                let mut all_deltas: Vec<u8> = Vec::new();
                for qual in qualities {
                    let deltas = quality_model::encode_with_model(qual, &model);
                    for d in &deltas {
                        all_deltas.push((*d as i16 + 128) as u8);
                    }
                }
                let max_delta = all_deltas.iter().copied().max().unwrap_or(0);
                let bits = min_bits_for_range(max_delta);
                let bitpacked = bitpack_values(&all_deltas, bits);
                let model_overhead = quality_model::serialize_model(&model).len();
                let encoding_time = start.elapsed().as_secs_f64() * 1000.0;

                for level in [3, 19] {
                    let start = Instant::now();
                    if let Ok(compressed) = zstd::bulk::compress(&bitpacked, level) {
                        let total = compressed.len() + model_overhead;
                        let elapsed = encoding_time + start.elapsed().as_secs_f64() * 1000.0;
                        res.lock().unwrap().push(StrategyResult {
                            stream: "qualities",
                            strategy: format!("Model + bitpack {}b + zstd {}", bits, level),
                            original_bytes,
                            compressed_bytes: total,
                            ratio: original_bytes as f64 / total as f64,
                            bits_per_element: (total as f64 * 8.0) / num_elements as f64,
                            compress_time_ms: elapsed,
                            lossy: false,
                        });
                    }
                }
            }
        });

        // Bitpacked inter-read deltas + zstd
        s.spawn(|_| {
            let start = Instant::now();
            let encoded_deltas = quality_delta::encode_quality_deltas(qualities);
            let mut all_deltas: Vec<u8> = Vec::new();
            for deltas in &encoded_deltas {
                for d in deltas {
                    all_deltas.push((*d as i16 + 128) as u8);
                }
            }
            let max_delta = all_deltas.iter().copied().max().unwrap_or(0);
            let bits = min_bits_for_range(max_delta);
            let bitpacked = bitpack_values(&all_deltas, bits);
            let encoding_time = start.elapsed().as_secs_f64() * 1000.0;

            let start = Instant::now();
            if let Ok(compressed) = zstd::bulk::compress(&bitpacked, 3) {
                let elapsed = encoding_time + start.elapsed().as_secs_f64() * 1000.0;
                res.lock().unwrap().push(StrategyResult {
                    stream: "qualities",
                    strategy: format!("Inter-read delta + bitpack {}b + zstd 3", bits),
                    original_bytes,
                    compressed_bytes: compressed.len(),
                    ratio: original_bytes as f64 / compressed.len() as f64,
                    bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                    compress_time_ms: elapsed,
                    lossy: false,
                });
            }
        });

        // --- Lossy strategies ---

        // Illumina8 binned + zstd
        s.spawn(|_| {
            if let Some(r) = bench_zstd("qualities", "Illumina8 + zstd 3", packed_illumina8_ref, original_bytes, num_elements, 3, true) {
                res.lock().unwrap().push(r);
            }
        });

        // Illumina8 binned + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("qualities", "Illumina8 + bsc", packed_illumina8_ref, original_bytes, num_elements, true) {
                res.lock().unwrap().push(r);
            }
        });

        // Illumina8 bitpacked + zstd
        for level in [3, 19] {
            s.spawn(move |_| {
                let start = Instant::now();
                if let Ok(compressed) = zstd::bulk::compress(bitpacked_illumina8_ref, level) {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "qualities",
                        strategy: format!("Illumina8 bitpack {}b + zstd {}", illumina8_bits, level),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: true,
                    });
                }
            });
        }

        // Illumina8 bitpacked + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("qualities", &format!("Illumina8 bitpack {}b + bsc", illumina8_bits), bitpacked_illumina8_ref, original_bytes, num_elements, true) {
                res.lock().unwrap().push(r);
            }
        });

        // FourLevel binned + zstd
        s.spawn(|_| {
            if let Some(r) = bench_zstd("qualities", "FourLevel + zstd 3", packed_fourlevel_ref, original_bytes, num_elements, 3, true) {
                res.lock().unwrap().push(r);
            }
        });

        // FourLevel bitpacked + zstd
        for level in [3, 19] {
            s.spawn(move |_| {
                let start = Instant::now();
                if let Ok(compressed) = zstd::bulk::compress(bitpacked_fourlevel_ref, level) {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "qualities",
                        strategy: format!("FourLevel bitpack {}b + zstd {}", fourlevel_bits, level),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: true,
                    });
                }
            });
        }

        // FourLevel bitpacked + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("qualities", &format!("FourLevel bitpack {}b + bsc", fourlevel_bits), bitpacked_fourlevel_ref, original_bytes, num_elements, true) {
                res.lock().unwrap().push(r);
            }
        });

        // Binary binned + zstd
        s.spawn(|_| {
            if let Some(r) = bench_zstd("qualities", "Binary + zstd 3", packed_binary_ref, original_bytes, num_elements, 3, true) {
                res.lock().unwrap().push(r);
            }
        });

        // Binary bitpacked + zstd
        for level in [3, 19] {
            s.spawn(move |_| {
                let start = Instant::now();
                if let Ok(compressed) = zstd::bulk::compress(bitpacked_binary_ref, level) {
                    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
                    res.lock().unwrap().push(StrategyResult {
                        stream: "qualities",
                        strategy: format!("Binary bitpack {}b + zstd {}", binary_bits, level),
                        original_bytes,
                        compressed_bytes: compressed.len(),
                        ratio: original_bytes as f64 / compressed.len() as f64,
                        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
                        compress_time_ms: elapsed,
                        lossy: true,
                    });
                }
            });
        }

        // Binary bitpacked + bsc
        s.spawn(|_| {
            if let Some(r) = bench_bsc("qualities", &format!("Binary bitpack {}b + bsc", binary_bits), bitpacked_binary_ref, original_bytes, num_elements, true) {
                res.lock().unwrap().push(r);
            }
        });
    }); // end rayon::scope

    let mut results = results.into_inner().unwrap();
    results.sort_by(|a, b| b.ratio.partial_cmp(&a.ratio).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// ---------------------------------------------------------------------------
// Shared compression helpers
// ---------------------------------------------------------------------------

fn bench_zstd(
    stream: &'static str,
    name: &str,
    data: &[u8],
    original_bytes: usize,
    num_elements: usize,
    level: i32,
    lossy: bool,
) -> Option<StrategyResult> {
    let start = Instant::now();
    let compressed = zstd::bulk::compress(data, level).ok()?;
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    Some(StrategyResult {
        stream,
        strategy: name.to_string(),
        original_bytes,
        compressed_bytes: compressed.len(),
        ratio: original_bytes as f64 / compressed.len() as f64,
        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
        compress_time_ms: elapsed,
        lossy,
    })
}

fn bench_gzip(
    stream: &'static str,
    name: &str,
    data: &[u8],
    original_bytes: usize,
    num_elements: usize,
    lossy: bool,
) -> Option<StrategyResult> {
    let start = Instant::now();
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(data).ok()?;
    let compressed = encoder.finish().ok()?;
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    Some(StrategyResult {
        stream,
        strategy: name.to_string(),
        original_bytes,
        compressed_bytes: compressed.len(),
        ratio: original_bytes as f64 / compressed.len() as f64,
        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
        compress_time_ms: elapsed,
        lossy,
    })
}

fn bench_bsc(
    stream: &'static str,
    name: &str,
    data: &[u8],
    original_bytes: usize,
    num_elements: usize,
    lossy: bool,
) -> Option<StrategyResult> {
    let start = Instant::now();
    let compressed = bsc::compress(data).ok()?;
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    Some(StrategyResult {
        stream,
        strategy: name.to_string(),
        original_bytes,
        compressed_bytes: compressed.len(),
        ratio: original_bytes as f64 / compressed.len() as f64,
        bits_per_element: (compressed.len() as f64 * 8.0) / num_elements as f64,
        compress_time_ms: elapsed,
        lossy,
    })
}

// ---------------------------------------------------------------------------
// Bitpacking helpers
// ---------------------------------------------------------------------------

/// Pack an array of values using a fixed number of bits per value.
///
/// Values are packed MSB-first into a byte stream. The caller must ensure
/// that every value fits within `bits_per_value` bits (i.e., value < 2^bits_per_value).
fn bitpack_values(values: &[u8], bits_per_value: u8) -> Vec<u8> {
    if bits_per_value == 0 || values.is_empty() {
        return Vec::new();
    }
    if bits_per_value >= 8 {
        return values.to_vec();
    }

    let total_bits = values.len() * bits_per_value as usize;
    let mut output = Vec::with_capacity((total_bits + 7) / 8);
    let mut current_byte: u8 = 0;
    let mut bits_used: u8 = 0;

    for &val in values {
        // Insert bits_per_value bits of val into the output stream
        let mut remaining_bits = bits_per_value;
        let masked_val = val & ((1u8 << bits_per_value) - 1);

        while remaining_bits > 0 {
            let space = 8 - bits_used;
            if remaining_bits <= space {
                current_byte |= masked_val.checked_shl((space - remaining_bits) as u32).unwrap_or(0);
                bits_used += remaining_bits;
                remaining_bits = 0;
            } else {
                current_byte |= masked_val.checked_shr((remaining_bits - space) as u32).unwrap_or(0);
                remaining_bits -= space;
                bits_used = 8;
            }

            if bits_used == 8 {
                output.push(current_byte);
                current_byte = 0;
                bits_used = 0;
            }
        }
    }
    if bits_used > 0 {
        output.push(current_byte);
    }
    output
}

/// Determine the minimum number of bits needed to represent the max value in a dataset.
fn min_bits_for_range(max_value: u8) -> u8 {
    if max_value == 0 {
        return 1;
    }
    8 - max_value.leading_zeros() as u8
}

/// Encode sequences to 3-bit representation (A=0, C=1, G=2, T=3, N=4).
/// Unlike 2-bit, this preserves N bases without a separate mask.
fn encode_3bit(sequences: &[String]) -> Vec<u8> {
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let mut values = Vec::with_capacity(total_bases);

    for seq in sequences {
        for base in seq.bytes() {
            let val = match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4, // N and other ambiguous bases
            };
            values.push(val);
        }
    }
    bitpack_values(&values, 3)
}

/// Convert quality strings to raw Phred scores (subtract 33 from ASCII).
fn qualities_to_phred(qualities: &[String]) -> Vec<u8> {
    let total: usize = qualities.iter().map(|q| q.len()).sum();
    let mut phred = Vec::with_capacity(total);
    for qual in qualities {
        for &byte in qual.as_bytes() {
            phred.push(byte.saturating_sub(33));
        }
    }
    phred
}

/// Bitpack quality scores at their natural minimum bit width.
fn bitpack_qualities_auto(qualities: &[String]) -> (Vec<u8>, u8) {
    let phred = qualities_to_phred(qualities);
    let max_val = phred.iter().copied().max().unwrap_or(0);
    let bits = min_bits_for_range(max_val);
    (bitpack_values(&phred, bits), bits)
}

/// Bitpack quality scores after applying a binning scheme.
/// Returns (packed_bytes, bits_per_value).
fn bitpack_binned_qualities(qualities: &[String], binning: QualityBinning) -> (Vec<u8>, u8) {
    let mut values = Vec::new();
    for qual in qualities {
        for &byte in qual.as_bytes() {
            let phred = byte.saturating_sub(33);
            values.push(binning.encode(phred));
        }
    }
    let max_val = values.iter().copied().max().unwrap_or(0);
    let bits = min_bits_for_range(max_val);
    (bitpack_values(&values, bits), bits)
}

/// Encode sequences to 2-bit representation (no N-mask; N mapped to A).
fn encode_2bit(sequences: &[String]) -> Vec<u8> {
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let mut output = Vec::with_capacity((total_bases + 3) / 4);
    let mut current_byte = 0u8;
    let mut bits_used = 0;

    for seq in sequences {
        for base in seq.bytes() {
            let bits = match base {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _ => 0b00,
            };
            current_byte = (current_byte << 2) | bits;
            bits_used += 2;
            if bits_used == 8 {
                output.push(current_byte);
                current_byte = 0;
                bits_used = 0;
            }
        }
    }
    if bits_used > 0 {
        current_byte <<= 8 - bits_used;
        output.push(current_byte);
    }
    output
}

/// Apply move-to-front encoding using an adaptive order-N context model for DNA.
///
/// For each base, the context model predicts a probability ranking of {A,C,G,T}.
/// The output is the rank of the actual base (0=most probable, 3=least probable).
/// This concentrates mass on symbol 0, making the output highly compressible.
fn sequence_context_mtf(sequences: &[String], order: usize) -> Vec<u8> {
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let mut output = Vec::with_capacity(total_bases);

    let table_size: usize = 1 << (order * 2).min(20); // 4^order, capped at 1M
    let mask = table_size - 1;
    // Flat array: [context][4 bases] counts
    let mut counts = vec![[1u32; 4]; table_size];

    for seq in sequences {
        let mut ctx_hash: usize = 0;
        let mut bases_in_ctx: usize = 0;

        for &byte in seq.as_bytes() {
            let base = match byte {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0, // N → A for context model
            };

            // Get current context's counts and rank bases
            let ctx = ctx_hash & mask;
            let ctx_counts = &counts[ctx];
            let mut ranked: [(u32, u8); 4] = [
                (ctx_counts[0], 0),
                (ctx_counts[1], 1),
                (ctx_counts[2], 2),
                (ctx_counts[3], 3),
            ];
            ranked.sort_by(|a, b| b.0.cmp(&a.0)); // Most probable first

            // Find rank of actual base
            let rank = ranked.iter().position(|&(_, b)| b == base).unwrap_or(0) as u8;
            output.push(rank);

            // Update counts
            counts[ctx][base as usize] += 1;
            // Rescale to prevent overflow
            let total: u32 = counts[ctx].iter().sum();
            if total >= 1 << 15 {
                for c in counts[ctx].iter_mut() {
                    *c = (*c + 1) / 2;
                }
            }

            // Update context hash (shift in 2 bits for this base)
            if bases_in_ctx < order {
                ctx_hash = (ctx_hash << 2) | base as usize;
                bases_in_ctx += 1;
            } else {
                ctx_hash = ((ctx_hash << 2) | base as usize) & mask;
            }
        }
    }

    output
}

/// Transpose sequences from row-major (read-by-read) to column-major (position-by-position).
/// For each position p (0..max_len), outputs all bases at position p across all reads.
/// Reads shorter than max_len contribute no byte at that position (variable-length aware).
fn transpose_sequences_column_major(sequences: &[String]) -> Vec<u8> {
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    let num_reads = sequences.len();

    // Pre-convert to bytes for fast indexing
    let seq_bytes: Vec<&[u8]> = sequences.iter().map(|s| s.as_bytes()).collect();

    let mut output = Vec::with_capacity(num_reads * max_len);
    for pos in 0..max_len {
        for seq in &seq_bytes {
            if pos < seq.len() {
                output.push(seq[pos]);
            }
            // Skip reads shorter than this position (no padding needed)
        }
    }
    output
}

/// Encode raw DNA bytes (ASCII A/C/G/T/N) to 2-bit (4 bases per byte, N→A).
fn encode_2bit_raw(data: &[u8]) -> Vec<u8> {
    let mut output = Vec::with_capacity((data.len() + 3) / 4);
    let mut current_byte = 0u8;
    let mut bits_used = 0;

    for &base in data {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00,
        };
        current_byte = (current_byte << 2) | bits;
        bits_used += 2;
        if bits_used == 8 {
            output.push(current_byte);
            current_byte = 0;
            bits_used = 0;
        }
    }
    if bits_used > 0 {
        current_byte <<= 8 - bits_used;
        output.push(current_byte);
    }
    output
}

/// Pack all quality strings through a given binning scheme.
fn pack_all_qualities(qualities: &[String], binning: QualityBinning) -> Vec<u8> {
    let mut packed = Vec::new();
    for qual in qualities {
        packed.extend_from_slice(&columnar::pack_qualities(qual, binning));
    }
    packed
}

// ---------------------------------------------------------------------------
// Whole-file benchmark (FQZ vs external tools)
// ---------------------------------------------------------------------------

/// Run whole-file compression benchmarks comparing FQZ vs other tools.
pub fn run_whole_file_benchmark(records: &[FastqRecord]) -> Result<WholeFileBenchmark> {
    let raw_fastq = reconstruct_fastq_bytes(records);
    let original_bytes = raw_fastq.len();
    let mut results = Vec::new();

    // FQZ (best defaults: template+zstd for headers, 2-bit+N-mask+zstd for sequences, BSC for qualities)
    match bench_fqz_whole(records, original_bytes) {
        Ok(r) => results.push(r),
        Err(e) => eprintln!("  [skip] FQZ whole-file benchmark failed: {}", e),
    }

    // gzip level 6 (default)
    match bench_gzip_whole("gzip -6", &raw_fastq, original_bytes, 6) {
        Ok(r) => results.push(r),
        Err(e) => eprintln!("  [skip] gzip -6 benchmark failed: {}", e),
    }

    // gzip level 9 (best)
    match bench_gzip_whole("gzip -9", &raw_fastq, original_bytes, 9) {
        Ok(r) => results.push(r),
        Err(e) => eprintln!("  [skip] gzip -9 benchmark failed: {}", e),
    }

    // bzip2 level 9
    match bench_bzip2_whole(&raw_fastq, original_bytes) {
        Ok(r) => results.push(r),
        Err(e) => eprintln!("  [skip] bzip2 benchmark failed: {}", e),
    }

    // bgzip (simulated with 64KB independent gzip blocks)
    match bench_bgzip_whole(&raw_fastq, original_bytes) {
        Ok(r) => results.push(r),
        Err(e) => eprintln!("  [skip] bgzip benchmark failed: {}", e),
    }

    // Sort by ratio descending
    results.sort_by(|a, b| b.ratio.partial_cmp(&a.ratio).unwrap_or(std::cmp::Ordering::Equal));

    Ok(WholeFileBenchmark {
        num_reads: records.len(),
        original_fastq_bytes: original_bytes,
        results,
    })
}

/// Reconstruct raw FASTQ text from records (for general-purpose compressor input).
fn reconstruct_fastq_bytes(records: &[FastqRecord]) -> Vec<u8> {
    let mut buf = Vec::new();
    for r in records {
        buf.extend_from_slice(r.id.as_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(r.sequence.as_bytes());
        buf.push(b'\n');
        buf.push(b'+');
        buf.push(b'\n');
        if let Some(ref q) = r.quality {
            buf.extend_from_slice(q.as_bytes());
        }
        buf.push(b'\n');
    }
    buf
}

/// Benchmark FQZ compression (best defaults).
fn bench_fqz_whole(records: &[FastqRecord], original_bytes: usize) -> Result<WholeFileResult> {
    let start = Instant::now();

    // Headers: template + zstd 9
    let (headers_compressed, _template) = super::compress_headers(records, 9)?;

    // Sequences: Raw + BSC adaptive (best from benchmark)
    let (sequences_compressed, _nmasks) = super::compress_sequences_raw_bsc(records)?;

    // Qualities: BSC (best lossless from benchmark)
    let qualities_compressed = super::compress_qualities(
        records,
        QualityBinning::None,
        3,
        QualityCompressor::Bsc,
    )?;

    let total_compressed = headers_compressed.len() + sequences_compressed.len()
        + qualities_compressed.len();
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;

    Ok(WholeFileResult {
        tool: "FQZ".to_string(),
        original_bytes,
        compressed_bytes: total_compressed,
        ratio: original_bytes as f64 / total_compressed as f64,
        compress_time_ms: elapsed,
        notes: "Template+zstd9 | Raw+BSC-adaptive | BSC".to_string(),
    })
}

/// Benchmark gzip compression at a given level.
fn bench_gzip_whole(
    name: &str,
    raw_data: &[u8],
    original_bytes: usize,
    level: u32,
) -> Result<WholeFileResult> {
    let start = Instant::now();
    let mut encoder = flate2::write::GzEncoder::new(
        Vec::new(),
        flate2::Compression::new(level),
    );
    encoder.write_all(raw_data)?;
    let compressed = encoder.finish()?;
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;

    Ok(WholeFileResult {
        tool: name.to_string(),
        original_bytes,
        compressed_bytes: compressed.len(),
        ratio: original_bytes as f64 / compressed.len() as f64,
        compress_time_ms: elapsed,
        notes: format!("flate2 level {}", level),
    })
}

/// Benchmark bzip2 compression (best level).
fn bench_bzip2_whole(raw_data: &[u8], original_bytes: usize) -> Result<WholeFileResult> {
    use bzip2::write::BzEncoder;
    use bzip2::Compression;

    let start = Instant::now();
    let mut encoder = BzEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(raw_data)?;
    let compressed = encoder.finish()?;
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;

    Ok(WholeFileResult {
        tool: "bzip2 -9".to_string(),
        original_bytes,
        compressed_bytes: compressed.len(),
        ratio: original_bytes as f64 / compressed.len() as f64,
        compress_time_ms: elapsed,
        notes: "bzip2 crate, best compression".to_string(),
    })
}

/// Benchmark bgzip (simulated: independent 64KB gzip blocks).
fn bench_bgzip_whole(raw_data: &[u8], original_bytes: usize) -> Result<WholeFileResult> {
    const BGZIP_BLOCK_SIZE: usize = 65536;

    let start = Instant::now();
    let mut compressed = Vec::new();

    for chunk in raw_data.chunks(BGZIP_BLOCK_SIZE) {
        let mut encoder = flate2::write::GzEncoder::new(
            Vec::new(),
            flate2::Compression::new(6),
        );
        encoder.write_all(chunk)?;
        let block = encoder.finish()?;
        compressed.extend_from_slice(&block);
    }

    let elapsed = start.elapsed().as_secs_f64() * 1000.0;

    Ok(WholeFileResult {
        tool: "bgzip".to_string(),
        original_bytes,
        compressed_bytes: compressed.len(),
        ratio: original_bytes as f64 / compressed.len() as f64,
        compress_time_ms: elapsed,
        notes: "Simulated BGZF (64KB independent blocks)".to_string(),
    })
}

// ---------------------------------------------------------------------------
// Report printing
// ---------------------------------------------------------------------------

fn format_size(bytes: usize) -> String {
    if bytes < 1024 {
        format!("{} B", bytes)
    } else if bytes < 1024 * 1024 {
        format!("{:.1} KB", bytes as f64 / 1024.0)
    } else {
        format!("{:.2} MB", bytes as f64 / (1024.0 * 1024.0))
    }
}

fn print_stream_table(title: &str, original_bytes: usize, num_reads: usize, results: &[StrategyResult]) {
    println!();
    println!("{} (original: {}, {} reads)", title, format_size(original_bytes), num_reads);
    println!("{}", "\u{2500}".repeat(105));
    println!(
        " {:<3} {:<42} {:>12} {:>8} {:>10} {:>9} {:>5}",
        "#", "Strategy", "Compressed", "Ratio", "bits/elem", "Time", "Lossy"
    );
    println!("{}", "\u{2500}".repeat(105));

    for (i, r) in results.iter().enumerate() {
        println!(
            " {:<3} {:<42} {:>12} {:>7.2}x {:>9.3} {:>7.1}ms {:>5}",
            i + 1,
            r.strategy,
            format_size(r.compressed_bytes),
            r.ratio,
            r.bits_per_element,
            r.compress_time_ms,
            if r.lossy { "yes" } else { "no" },
        );
    }
}

impl WholeFileBenchmark {
    pub fn print(&self) {
        println!();
        println!("====================================================================");
        println!("  WHOLE-FILE COMPRESSION COMPARISON");
        println!("====================================================================");
        println!();
        println!("  Reads: {}    Original FASTQ: {}", self.num_reads, format_size(self.original_fastq_bytes));
        println!();
        println!("{}", "\u{2500}".repeat(105));
        println!(
            " {:<3} {:<20} {:>12} {:>8} {:>10} {:>10}   {}",
            "#", "Tool", "Compressed", "Ratio", "bits/base", "Time", "Notes"
        );
        println!("{}", "\u{2500}".repeat(105));

        for (i, r) in self.results.iter().enumerate() {
            let bits_per_byte = (r.compressed_bytes as f64 * 8.0) / self.original_fastq_bytes as f64;
            println!(
                " {:<3} {:<20} {:>12} {:>7.2}x {:>9.3} {:>8.1}ms   {}",
                i + 1,
                r.tool,
                format_size(r.compressed_bytes),
                r.ratio,
                bits_per_byte,
                r.compress_time_ms,
                r.notes,
            );
        }
    }
}

impl BenchmarkReport {
    /// Print the full benchmark report with all strategies ranked per stream.
    pub fn print(&self) {
        println!();
        println!("====================================================================");
        println!("  FQZ COMPRESSION STRATEGY BENCHMARK");
        println!("====================================================================");
        println!();
        println!("  Reads: {}", self.num_reads);
        println!("  Total bases: {} ({})", self.total_bases, format_size(self.total_bases));
        println!("  Total quality scores: {} ({})", self.total_quality_scores, format_size(self.total_quality_scores));

        let header_orig = self.header_results.first().map(|r| r.original_bytes).unwrap_or(0);
        let seq_orig = self.sequence_results.first().map(|r| r.original_bytes).unwrap_or(0);
        let qual_orig = self.quality_results.first().map(|r| r.original_bytes).unwrap_or(0);

        print_stream_table("HEADERS", header_orig, self.num_reads, &self.header_results);
        print_stream_table("SEQUENCES", seq_orig, self.num_reads, &self.sequence_results);
        print_stream_table("QUALITIES", qual_orig, self.num_reads, &self.quality_results);

        // Print summary of best per stream
        println!();
        println!("====================================================================");
        println!("  BEST STRATEGIES");
        println!("====================================================================");

        if let Some(best) = self.header_results.first() {
            println!("  Headers:    {:<35} {:>7.2}x  {}", best.strategy, best.ratio, format_size(best.compressed_bytes));
        }
        if let Some(best) = self.sequence_results.first() {
            println!("  Sequences:  {:<35} {:>7.2}x  {}", best.strategy, best.ratio, format_size(best.compressed_bytes));
        }
        // Best lossless quality
        if let Some(best) = self.quality_results.iter().find(|r| !r.lossy) {
            println!("  Qualities (lossless):  {:<24} {:>7.2}x  {}", best.strategy, best.ratio, format_size(best.compressed_bytes));
        }
        // Best lossy quality
        if let Some(best) = self.quality_results.iter().find(|r| r.lossy) {
            println!("  Qualities (lossy):    {:<24} {:>7.2}x  {}", best.strategy, best.ratio, format_size(best.compressed_bytes));
        }

        // Combined estimate
        let best_header = self.header_results.first().map(|r| r.compressed_bytes).unwrap_or(header_orig);
        let best_seq = self.sequence_results.first().map(|r| r.compressed_bytes).unwrap_or(seq_orig);
        let best_qual_lossless = self.quality_results.iter().find(|r| !r.lossy).map(|r| r.compressed_bytes).unwrap_or(qual_orig);
        let total_orig = header_orig + seq_orig + qual_orig;
        let total_best = best_header + best_seq + best_qual_lossless;

        println!();
        println!("  Combined (best lossless): {} -> {} ({:.2}x)",
            format_size(total_orig), format_size(total_best),
            total_orig as f64 / total_best as f64);

        if let Some(ref wf) = self.whole_file {
            wf.print();
        }

        println!();
    }
}
