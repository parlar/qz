//! Streaming BSC decompression: parallel block decompression through bounded channels.

use anyhow::{Context, Result};
use std::time::Instant;
use tracing::info;
use crate::cli::DecompressConfig;
use super::*;
use super::codecs;

/// Block-by-block BSC stream decoder for memory-efficient decompression.
///
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
pub(super) fn can_stream_decompress(args: &DecompressConfig) -> Result<bool> {
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
    let nmasks_len = read_le_u64(&header, 36)?;

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
/// Peak memory: ~300 MB (3 streams x ~4 decompressed blocks x 25 MB + output buffer)
/// regardless of input size. Uses all available CPU cores via rayon.
pub(super) fn decompress_streaming_bsc(args: &DecompressConfig) -> Result<()> {
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

    let num_reads = read_le_u64(&header, 12)? as usize;
    let headers_len = read_le_u64(&header, 20)? as usize;
    let sequences_len = read_le_u64(&header, 28)? as usize;
    let nmasks_len = read_le_u64(&header, 36)? as usize;
    let qualities_len = read_le_u64(&header, 44)? as usize;

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
/// Parsed archive header metadata.
struct ArchiveHeader {
    encoding_type: u8,
    arithmetic_enabled: bool,
    read_lengths_opt: Option<Vec<usize>>,
    quality_binning: QualityBinning,
    quality_compressor: QualityCompressor,
    sequence_compressor: SequenceCompressor,
    header_compressor: HeaderCompressor,
    quality_model_opt: Option<quality_model::QualityModel>,
    quality_delta_enabled: bool,
    quality_dict_opt: Option<Vec<u8>>,
    read_id_template: read_id::ReadIdTemplate,
    num_reads: usize,
    /// Byte offset where stream data begins.
    data_offset: usize,
    headers_len: usize,
    sequences_len: usize,
    nmasks_len: usize,
    qualities_len: usize,
}

/// Parse the archive header from raw archive bytes.
///
/// Returns parsed metadata and the byte offset where stream data starts.
fn parse_archive_header(data: &[u8]) -> Result<ArchiveHeader> {
    if data.len() < 50 {
        anyhow::bail!("Invalid archive: too small");
    }

    let mut offset = 0;

    let encoding_type = data[offset];
    offset += 1;

    let arithmetic_enabled = data[offset] != 0;
    offset += 1;
    let read_lengths_opt = if arithmetic_enabled {
        let num_lengths = read_le_u32(data, offset)? as usize;
        offset += 4;
        let mut lengths = Vec::with_capacity(num_lengths);
        for _ in 0..num_lengths {
            let len = read_le_u32(data, offset)? as usize;
            offset += 4;
            lengths.push(len);
        }
        Some(lengths)
    } else {
        None
    };

    let quality_binning = code_to_binning(data[offset])?;
    offset += 1;

    let quality_compressor = if offset < data.len() && data[offset] <= 4 {
        let compressor = code_to_compressor(data[offset])?;
        offset += 1;
        compressor
    } else {
        QualityCompressor::Zstd
    };

    // Fqzcomp stores raw ASCII quality bytes (not bit-packed), so override binning
    let quality_binning = if quality_compressor == QualityCompressor::Fqzcomp {
        QualityBinning::None
    } else {
        quality_binning
    };

    let sequence_compressor = code_to_seq_compressor(data[offset])?;
    offset += 1;

    let header_compressor = code_to_header_compressor(data[offset])?;
    offset += 1;

    let quality_model_enabled = data[offset] != 0;
    offset += 1;
    let quality_model_opt = if quality_model_enabled {
        let model_size = read_le_u16(data, offset)? as usize;
        offset += 2;
        let model_bytes = &data[offset..offset + model_size];
        offset += model_size;
        Some(quality_model::deserialize_model(model_bytes)?)
    } else {
        None
    };

    let quality_delta_enabled = data[offset] != 0;
    offset += 1;

    let quality_dict_opt = if data[offset] != 0 {
        offset += 1;
        let dict_size = read_le_u32(data, offset)? as usize;
        offset += 4;
        let dict = data[offset..offset + dict_size].to_vec();
        offset += dict_size;
        Some(dict)
    } else {
        offset += 1;
        None
    };

    let template_prefix_len = read_le_u16(data, offset)? as usize;
    offset += 2;
    let template_prefix = String::from_utf8_lossy(&data[offset..offset + template_prefix_len]).to_string();
    offset += template_prefix_len;
    let template_has_comment = data[offset] != 0;
    offset += 1;

    let common_comment = if template_prefix_len > 0 && template_has_comment {
        let cc_len = read_le_u16(data, offset)? as usize;
        offset += 2;
        if cc_len > 0 {
            let cc = String::from_utf8_lossy(&data[offset..offset + cc_len]).to_string();
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

    let num_reads = read_le_u64(data, offset)? as usize;
    offset += 8;
    let headers_len = read_le_u64(data, offset)? as usize;
    offset += 8;
    let sequences_len = read_le_u64(data, offset)? as usize;
    offset += 8;
    let nmasks_len = read_le_u64(data, offset)? as usize;
    offset += 8;
    let qualities_len = read_le_u64(data, offset)? as usize;
    offset += 8;

    if data.len() < offset + headers_len + sequences_len + nmasks_len + qualities_len {
        anyhow::bail!("Invalid archive: data truncated");
    }

    Ok(ArchiveHeader {
        encoding_type,
        arithmetic_enabled,
        read_lengths_opt,
        quality_binning,
        quality_compressor,
        sequence_compressor,
        header_compressor,
        quality_model_opt,
        quality_delta_enabled,
        quality_dict_opt,
        read_id_template,
        num_reads,
        data_offset: offset,
        headers_len,
        sequences_len,
        nmasks_len,
        qualities_len,
    })
}

/// Decode packed quality data and attach quality strings to records.
///
/// Handles all quality encoding variants: quality_delta, fqzcomp, quality_model,
/// and standard binning mode. Modifies records in place.
fn decode_qualities_to_records(
    records: &mut [crate::io::FastqRecord],
    qualities_data: &[u8],
    hdr: &ArchiveHeader,
) -> Result<()> {
    if qualities_data.is_empty() {
        return Ok(());
    }

    if hdr.quality_delta_enabled {
        let mut encoded_deltas = Vec::new();
        let mut qual_offset = 0;
        while qual_offset < qualities_data.len() {
            let q_len = read_varint(qualities_data, &mut qual_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
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
        for record in records.iter_mut() {
            if qual_offset >= qualities_data.len() {
                break;
            }
            let q_len = read_varint(qualities_data, &mut qual_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;

            let quality_str = if hdr.quality_compressor == QualityCompressor::Fqzcomp {
                if qual_offset + q_len <= qualities_data.len() {
                    let raw = &qualities_data[qual_offset..qual_offset + q_len];
                    qual_offset += q_len;
                    unsafe { String::from_utf8_unchecked(raw.to_vec()) }
                } else {
                    break;
                }
            } else if let Some(ref model) = hdr.quality_model_opt {
                if qual_offset + q_len <= qualities_data.len() {
                    let deltas = quality_model::unpack_deltas(&qualities_data[qual_offset..qual_offset + q_len]);
                    qual_offset += q_len;
                    quality_model::decode_with_model(&deltas, model)
                } else {
                    break;
                }
            } else {
                let bits_per_qual = hdr.quality_binning.bits_per_quality();
                let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                if qual_offset + q_encoded_len <= qualities_data.len() {
                    let quality_str = columnar::unpack_qualities(
                        &qualities_data[qual_offset..qual_offset + q_encoded_len],
                        q_len, hdr.quality_binning,
                    );
                    qual_offset += q_encoded_len;
                    quality_str
                } else {
                    break;
                }
            };
            record.quality = Some(quality_str);
        }
    }

    Ok(())
}

pub(super) fn decompress(args: &DecompressConfig) -> Result<()> {
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

    // Parse archive header
    let hdr = parse_archive_header(&archive_data)?;
    let mut offset = hdr.data_offset;
    let headers = &archive_data[offset..offset + hdr.headers_len];
    offset += hdr.headers_len;
    let sequences = &archive_data[offset..offset + hdr.sequences_len];
    offset += hdr.sequences_len;
    let nmasks = &archive_data[offset..offset + hdr.nmasks_len];
    offset += hdr.nmasks_len;
    let qualities = &archive_data[offset..offset + hdr.qualities_len];
    let num_reads = hdr.num_reads;
    let template_prefix_len = hdr.read_id_template.prefix.len();
    let encoding_type = hdr.encoding_type;
    let arithmetic_enabled = hdr.arithmetic_enabled;
    let quality_binning = hdr.quality_binning;
    let quality_compressor = hdr.quality_compressor;
    let sequence_compressor = hdr.sequence_compressor;
    let header_compressor = hdr.header_compressor;
    let quality_delta_enabled = hdr.quality_delta_enabled;
    let qualities_len = hdr.qualities_len;
    let nmasks_len = hdr.nmasks_len;
    let read_id_template = &hdr.read_id_template;
    let read_lengths_opt = hdr.read_lengths_opt.clone();
    let quality_model_opt = &hdr.quality_model_opt;
    let quality_dict_opt = &hdr.quality_dict_opt;

    // Decompress based on encoding mode
    info!("Decompressing...");
    let records = if arithmetic_enabled {
        info!("Decompressing arithmetic-coded data...");

        // Decompress headers
        let read_ids = decompress_headers_dispatch(
            header_compressor, headers, template_prefix_len, &read_id_template, num_reads,
        )?;

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
        let read_ids = decompress_headers_dispatch(
            header_compressor, headers, template_prefix_len, &read_id_template, num_reads,
        )?;

        // Decompress sequences with de Bruijn graph
        let decoded_sequences = debruijn::decompress_sequences_debruijn(sequences, num_reads)
            .context("Failed to decompress de Bruijn-coded sequences")?;

        // Decompress qualities
        let qualities_data = if qualities_len == 0 {
            Vec::new()
        } else if let Some(dict) = quality_dict_opt {
            zstd_dict::decompress_with_dict(qualities, dict)
                .context("Failed to decompress quality scores (dictionary mode)")?
        } else {
            codecs::decompress_qualities_data(qualities, quality_compressor)
                .context("Failed to decompress quality scores")?
        };

        // Reconstruct records
        let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
            crate::io::FastqRecord::new(id, seq, None)
        }).collect();

        decode_qualities_to_records(&mut records, &qualities_data, &hdr)?;

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
                || codecs::decompress_headers_bsc(headers, num_reads)
                    .context("Failed to decompress headers"),
                || -> Result<Vec<u8>> {
                    if qualities_len == 0 {
                        Ok(Vec::new())
                    } else {
                        codecs::decompress_qualities_data(qualities, quality_compressor)
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

            decode_qualities_to_records(&mut records, &qualities_data, &hdr)?;

            records
        } else {
            // ── Factorize v1: 3 sub-streams (meta, raw_seq, delta_seq) ──
            let mut soff = 0usize;

            let mut read_sub = |label: &str| -> Result<&[u8]> {
                if soff + 8 > sequences.len() {
                    anyhow::bail!("Truncated factorized {} length", label);
                }
                let len = read_le_u64(sequences, soff)? as usize;
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
                    || codecs::decompress_headers_bsc(headers, num_reads)
                        .context("Failed to decompress headers"),
                    || -> Result<Vec<u8>> {
                        if qualities_len == 0 {
                            Ok(Vec::new())
                        } else {
                            codecs::decompress_qualities_data(qualities, quality_compressor)
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

            decode_qualities_to_records(&mut records, &qualities_data, &hdr)?;

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
            || codecs::decompress_headers_bsc(headers, num_reads)
                .context("Failed to decompress headers"),
            || -> Result<ultra::HarcDecoded> {
                if encoding_type == 8 {
                    ultra::decode_harc_sequences(sequences, num_reads)
                } else {
                    ultra::decode_reorder_local(sequences, num_reads)
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
                let qualities_data = codecs::decompress_qualities_data(qualities, quality_compressor)
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

        records.into_iter().enumerate()
            .map(|(i, r)| r.ok_or_else(|| anyhow::anyhow!("missing record at index {i} after un-permutation")))
            .collect::<Result<Vec<_>>>()?
    } else if encoding_type == 1 || encoding_type == 2 {
        // Delta or RLE mode: decompress all streams in parallel
        let ((header_result, seq_result), qual_result) = rayon::join(
            || rayon::join(
                || decompress_headers_dispatch(
                    header_compressor, headers, template_prefix_len, &read_id_template, num_reads,
                ),
                || match sequence_compressor {
                    SequenceCompressor::Bsc => bsc::decompress_parallel(sequences)
                        .context("Failed to decompress sequences (BSC)"),
                    SequenceCompressor::Zstd => decompress_zstd(sequences)
                        .context("Failed to decompress sequences (zstd)"),
                    SequenceCompressor::OpenZl => openzl::decompress_parallel(sequences)
                        .context("Failed to decompress sequences (OpenZL)"),
                },
            ),
            || if let Some(dict) = quality_dict_opt {
                zstd_dict::decompress_with_dict(qualities, dict)
                    .context("Failed to decompress quality scores (dictionary mode)")
            } else {
                codecs::decompress_qualities_data(qualities, quality_compressor)
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
            let len = read_le_u32(&seq_stream, s_offset)? as usize;
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

        decode_qualities_to_records(&mut records, &quality_stream, &hdr)?;

        records
    } else {
        // Normal mode: decompress all streams in parallel
        info!("Decompressing streams in parallel...");
        let use_quality_ctx = quality_compressor == QualityCompressor::QualityCtx;

        // Step 1: Parallel decompression of headers, sequences, qualities
        // For quality_ctx: skip quality decompression here (needs sequences first)
        let (header_result, (seq_result, qual_result)) = rayon::join(
            || decompress_headers_dispatch(
                header_compressor, headers, template_prefix_len, &read_id_template, num_reads,
            ),
            || rayon::join(
                || -> Result<Vec<String>> {
                    match sequence_compressor {
                        SequenceCompressor::Bsc => {
                            if nmasks_len > 0 {
                                codecs::decompress_sequences_2bit_bsc(sequences, nmasks, num_reads)
                                    .context("Failed to decompress 2-bit sequences (BSC)")
                            } else {
                                codecs::decompress_sequences_raw_bsc(sequences, num_reads, encoding_type)
                                    .context("Failed to decompress sequences (BSC)")
                            }
                        }
                        SequenceCompressor::OpenZl => {
                            codecs::decompress_sequences_raw_openzl(sequences, num_reads)
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
                    } else if let Some(dict) = quality_dict_opt {
                        zstd_dict::decompress_with_dict(qualities, dict)
                            .context("Failed to decompress quality scores (dictionary mode)")
                    } else {
                        codecs::decompress_qualities_data(qualities, quality_compressor)
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
                // Two-pass parallel quality decoding:
                // Pass 1: Sequential scan to locate record boundaries (fast varint parsing)
                use rayon::prelude::*;
                let bits_per_qual = quality_binning.bits_per_quality();
                let mut qual_entries: Vec<(usize, usize)> = Vec::with_capacity(records.len());
                let mut qual_offset = 0;
                for _ in 0..records.len() {
                    if qual_offset >= qualities_data.len() { break; }
                    let q_len = read_varint(&qualities_data, &mut qual_offset)
                        .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?;
                    let data_len = if quality_compressor == QualityCompressor::Fqzcomp || quality_model_opt.is_some() {
                        q_len
                    } else {
                        (q_len * bits_per_qual + 7) / 8
                    };
                    if qual_offset + data_len > qualities_data.len() { break; }
                    qual_entries.push((qual_offset, q_len));
                    qual_offset += data_len;
                }

                // Pass 2: Parallel decode (each record is independent)
                let decoded_quals: Vec<String> = qual_entries.par_iter().map(|&(data_off, q_len)| {
                    if quality_compressor == QualityCompressor::Fqzcomp {
                        let raw = &qualities_data[data_off..data_off + q_len];
                        unsafe { String::from_utf8_unchecked(raw.to_vec()) }
                    } else if let Some(model) = quality_model_opt {
                        let deltas = quality_model::unpack_deltas(&qualities_data[data_off..data_off + q_len]);
                        quality_model::decode_with_model(&deltas, model)
                    } else {
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        columnar::unpack_qualities(&qualities_data[data_off..data_off + q_encoded_len], q_len, quality_binning)
                    }
                }).collect();

                for (i, qual) in decoded_quals.into_iter().enumerate() {
                    if i < records.len() {
                        records[i].quality = Some(qual);
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

    // Batch FASTQ formatting into 2MB chunks to reduce write_all calls through dyn Write
    {
        const WRITE_BATCH: usize = 2 * 1024 * 1024;
        let mut buf = Vec::with_capacity(WRITE_BATCH + 1024);
        for record in &records {
            buf.extend_from_slice(record.id.as_bytes());
            buf.push(b'\n');
            buf.extend_from_slice(record.sequence.as_bytes());
            buf.extend_from_slice(b"\n+\n");
            if let Some(qual) = &record.quality {
                buf.extend_from_slice(qual.as_bytes());
            }
            buf.push(b'\n');
            if buf.len() >= WRITE_BATCH {
                output.write_all(&buf)?;
                buf.clear();
            }
        }
        if !buf.is_empty() {
            output.write_all(&buf)?;
        }
    }

    let elapsed = start_time.elapsed();
    info!("Decompression completed in {:.2}s", elapsed.as_secs_f64());

    Ok(())
}

