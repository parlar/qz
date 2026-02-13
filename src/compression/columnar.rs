/// Columnar compression for FASTQ files (patent-free alternative)
///
/// Strategy:
/// 1. Separate FASTQ into 3 streams (headers, sequences, quality)
/// 2. Encode each optimally (2-bit DNA, quality binning)
/// 3. Compress each with Zstd
/// 4. No reordering required!

use crate::io::FastqRecord;
use anyhow::{Context, Result};
use std::io::Write;
use tracing::info;

/// Quality binning schemes
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum QualityBinning {
    /// Illumina 8-level binning (3 bits per quality)
    Illumina8,
    /// Binary threshold (1 bit per quality)
    Binary { threshold: u8 },
    /// 4-level binning (2 bits per quality)
    FourLevel,
    /// No binning (keep original, 7 bits per quality)
    None,
}

impl QualityBinning {
    /// Get number of bits per quality score
    pub fn bits_per_quality(&self) -> usize {
        match self {
            QualityBinning::Illumina8 => 3,
            QualityBinning::Binary { .. } => 1,
            QualityBinning::FourLevel => 2,
            QualityBinning::None => 7,
        }
    }

    /// Encode a quality score
    pub fn encode(&self, phred: u8) -> u8 {
        match self {
            QualityBinning::Illumina8 => {
                // Illumina 8-level binning
                match phred {
                    0..=1 => 0,
                    2..=9 => 1,
                    10..=19 => 2,
                    20..=24 => 3,
                    25..=29 => 4,
                    30..=34 => 5,
                    35..=39 => 6,
                    _ => 7,
                }
            }
            QualityBinning::Binary { threshold } => {
                if phred >= *threshold { 1 } else { 0 }
            }
            QualityBinning::FourLevel => {
                match phred {
                    0..=9 => 0,
                    10..=19 => 1,
                    20..=29 => 2,
                    _ => 3,
                }
            }
            QualityBinning::None => phred.min(127),
        }
    }

    /// Decode a quality score
    pub fn decode(&self, encoded: u8) -> u8 {
        match self {
            QualityBinning::Illumina8 => {
                // Use representative values
                match encoded {
                    0 => 6,
                    1 => 6,
                    2 => 15,
                    3 => 22,
                    4 => 27,
                    5 => 33,
                    6 => 37,
                    7 => 40,
                    _ => 40,
                }
            }
            QualityBinning::Binary { threshold } => {
                if encoded == 1 { 40 } else { *threshold / 2 }
            }
            QualityBinning::FourLevel => {
                match encoded {
                    0 => 6,
                    1 => 15,
                    2 => 25,
                    3 => 37,
                    _ => 37,
                }
            }
            QualityBinning::None => encoded,
        }
    }
}

/// Pack quality scores with variable bit width
pub fn pack_qualities(qualities: &str, binning: QualityBinning) -> Vec<u8> {
    let bytes = qualities.as_bytes();
    let bits_per_qual = binning.bits_per_quality();

    let mut packed = Vec::new();
    let mut buffer = 0u64;
    let mut bits_in_buffer = 0;

    for &qual_ascii in bytes {
        // Convert ASCII to Phred (assuming Phred+33 encoding)
        let phred = qual_ascii.saturating_sub(33);
        let encoded = binning.encode(phred);

        // Add to buffer
        buffer |= (encoded as u64) << bits_in_buffer;
        bits_in_buffer += bits_per_qual;

        // Flush complete bytes
        while bits_in_buffer >= 8 {
            packed.push((buffer & 0xFF) as u8);
            buffer >>= 8;
            bits_in_buffer -= 8;
        }
    }

    // Flush remaining bits
    if bits_in_buffer > 0 {
        packed.push((buffer & 0xFF) as u8);
    }

    packed
}

/// Unpack quality scores
pub fn unpack_qualities(packed: &[u8], length: usize, binning: QualityBinning) -> String {
    let bits_per_qual = binning.bits_per_quality();
    let mask = (1u64 << bits_per_qual) - 1;

    let mut qualities = Vec::with_capacity(length);
    let mut buffer = 0u64;
    let mut bits_in_buffer = 0;
    let mut packed_idx = 0;

    for _ in 0..length {
        // Ensure we have enough bits
        while bits_in_buffer < bits_per_qual && packed_idx < packed.len() {
            buffer |= (packed[packed_idx] as u64) << bits_in_buffer;
            bits_in_buffer += 8;
            packed_idx += 1;
        }

        // Extract quality
        let encoded = (buffer & mask) as u8;
        let phred = binning.decode(encoded);
        qualities.push(phred + 33); // Convert to ASCII

        buffer >>= bits_per_qual;
        bits_in_buffer = bits_in_buffer.saturating_sub(bits_per_qual);
    }

    String::from_utf8_lossy(&qualities).to_string()
}

/// Write variable-length integer
fn write_varint<W: Write>(writer: &mut W, mut value: usize) -> std::io::Result<()> {
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

/// Columnar compression statistics
#[derive(Debug)]
pub struct ColumnarStats {
    #[allow(dead_code)]
    pub num_reads: usize,
    pub total_bases: usize,
    pub original_size: usize,
    pub headers_size: usize,
    pub sequences_size: usize,
    pub qualities_size: usize,
    pub total_compressed: usize,
}

impl ColumnarStats {
    pub fn compression_ratio(&self) -> f64 {
        if self.total_compressed > 0 {
            self.original_size as f64 / self.total_compressed as f64
        } else {
            0.0
        }
    }
}

/// Compress FASTQ records using columnar approach with N-mask encoding
pub fn compress_columnar(
    records: &[FastqRecord],
    binning: QualityBinning,
    compression_level: i32,
) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, ColumnarStats)> {
    use crate::compression::n_mask::encode_with_n_mask;

    info!("Compressing {} records using columnar method with N-mask encoding", records.len());

    let mut stats = ColumnarStats {
        num_reads: records.len(),
        total_bases: 0,
        original_size: 0,
        headers_size: 0,
        sequences_size: 0,
        qualities_size: 0,
        total_compressed: 0,
    };

    // Stream 1: Headers with length prefixes
    let mut header_stream = Vec::new();
    for record in records {
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.write_all(record.id.as_bytes())?;
        stats.original_size += record.id.len() + 1; // +1 for @ line
    }

    // Stream 2: Sequences (2-bit encoded with N-mask)
    let mut sequence_stream = Vec::new();
    let mut nmask_stream = Vec::new();
    for record in records {
        let encoding = encode_with_n_mask(&record.sequence);
        write_varint(&mut sequence_stream, encoding.length)?;
        sequence_stream.write_all(&encoding.sequence_2bit)?;
        nmask_stream.write_all(&encoding.n_mask)?;
        stats.total_bases += record.sequence.len();
        stats.original_size += record.sequence.len() + 1; // +1 for newline
    }

    // Stream 3: Qualities (binned) with length prefixes
    let mut quality_stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut quality_stream, qual.len())?;
            let packed = pack_qualities(qual, binning);
            quality_stream.write_all(&packed)?;
            stats.original_size += qual.len() + 2; // +2 for + line and qual line
        }
    }

    // Compress each stream with adaptive levels
    // Headers: highly redundant, benefit from higher levels (min level 15)
    let header_level = std::cmp::min(std::cmp::max(compression_level + 6, 15), 19);
    let headers_compressed = compress_stream(&header_stream, header_level)?;

    // Sequences: already 2-bit encoded, use user's level
    let sequences_compressed = compress_stream(&sequence_stream, compression_level)?;

    // N-masks: very small data, always use maximum compression
    let nmasks_compressed = compress_stream(&nmask_stream, 19)?;

    // Qualities: binned qualities compress better, add +3 for binned modes
    let quality_level = match binning {
        QualityBinning::None => compression_level,
        _ => std::cmp::min(compression_level + 3, 19), // Binned qualities are more compressible
    };
    let qualities_compressed = compress_stream(&quality_stream, quality_level)?;

    stats.headers_size = headers_compressed.len();
    stats.sequences_size = sequences_compressed.len() + nmasks_compressed.len();
    stats.qualities_size = qualities_compressed.len();
    stats.total_compressed = stats.headers_size + stats.sequences_size + stats.qualities_size;

    info!(
        "Columnar compression: {:.2}x ({} -> {} bytes)",
        stats.compression_ratio(),
        stats.original_size,
        stats.total_compressed
    );
    info!(
        "  Headers: {} bytes, Sequences: {} bytes, Qualities: {} bytes",
        stats.headers_size, stats.sequences_size, stats.qualities_size
    );

    Ok((headers_compressed, sequences_compressed, nmasks_compressed, qualities_compressed, stats))
}

/// Decompress columnar format back to FASTQ records
pub fn decompress_columnar(
    headers_compressed: &[u8],
    sequences_compressed: &[u8],
    nmasks_compressed: &[u8],
    qualities_compressed: &[u8],
    binning: QualityBinning,
) -> Result<Vec<FastqRecord>> {
    use crate::compression::n_mask::{decode_with_n_mask, NMaskEncoding};

    // Decompress streams
    let header_stream = decompress_stream(headers_compressed)?;
    let sequence_stream = decompress_stream(sequences_compressed)?;
    let nmask_stream = decompress_stream(nmasks_compressed)?;
    let quality_stream = decompress_stream(qualities_compressed)?;

    let mut records = Vec::new();
    let mut h_offset = 0;
    let mut s_offset = 0;
    let mut n_offset = 0;
    let mut q_offset = 0;

    // Reconstruct records
    while h_offset < header_stream.len() {
        // Read header
        let h_len = read_varint(&header_stream, &mut h_offset)
            .context("Failed to read header length")?;
        let id = String::from_utf8_lossy(&header_stream[h_offset..h_offset + h_len]).to_string();
        h_offset += h_len;

        // Read sequence with N-mask
        let s_len = read_varint(&sequence_stream, &mut s_offset)
            .context("Failed to read sequence length")?;
        let s_encoded_len = (s_len + 3) / 4; // 4 bases per byte (2-bit)
        let n_mask_len = (s_len + 7) / 8;     // 8 bases per byte (bitmap)

        let encoding = NMaskEncoding {
            sequence_2bit: sequence_stream[s_offset..s_offset + s_encoded_len].to_vec(),
            n_mask: nmask_stream[n_offset..n_offset + n_mask_len].to_vec(),
            length: s_len,
        };
        let sequence = decode_with_n_mask(&encoding);
        s_offset += s_encoded_len;
        n_offset += n_mask_len;

        // Read quality
        let quality = if q_offset < quality_stream.len() {
            let q_len = read_varint(&quality_stream, &mut q_offset)
                .context("Failed to read quality length")?;
            let bits_per_qual = binning.bits_per_quality();
            let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
            let quality_str = unpack_qualities(&quality_stream[q_offset..q_offset + q_encoded_len], q_len, binning);
            q_offset += q_encoded_len;
            Some(quality_str)
        } else {
            None
        };

        records.push(FastqRecord::new(id, sequence, quality));
    }

    Ok(records)
}

/// Compress a byte stream with zstd at specified level
fn compress_stream(data: &[u8], level: i32) -> Result<Vec<u8>> {
    zstd::bulk::compress(data, level)
        .map_err(|e| anyhow::anyhow!("Zstd compression failed: {}", e))
}

/// Decompress a zstd stream
fn decompress_stream(compressed: &[u8]) -> Result<Vec<u8>> {
    zstd::bulk::decompress(compressed, 100_000_000) // 100MB limit
        .map_err(|e| anyhow::anyhow!("Zstd decompression failed: {}", e))
}

