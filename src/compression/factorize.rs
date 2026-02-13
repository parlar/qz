/// Chunked read-to-read delta compression with inverted index and split streams.
///
/// For each 5M-record chunk:
///   Pass 1: Build inverted index mapping canonical k-mer hashes to read indices.
///   Pass 2: For each read (in order), find the best earlier match via the index,
///           encode as delta (0x00=match, actual base=mismatch) or raw.
///
/// Three separate sequence sub-streams for optimal BWT compression:
///   - meta: class byte (0=raw, 1=fwd delta, 2=rc delta) + varint(ref_index) for deltas
///   - raw_seq: raw reads only (pure DNA — optimal BWT context)
///   - delta_seq: delta reads only (mostly 0x00 — compresses extremely well)
///
/// Streams are written to temp files per chunk (bounded memory), then assembled into archive.

use anyhow::Result;
use tracing::info;
use std::time::Instant;

use super::*;

const CHUNK_SIZE: usize = 5_000_000;
const INDEX_K: usize = 21;
const INDEX_STRIDE: usize = 15;
const BUCKET_CAP: usize = 64;

pub(super) fn compress_factorize(args: &crate::cli::CompressArgs) -> Result<()> {
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

    info!("Factorized compression: chunked read-to-read delta with split streams ({} records/chunk)", CHUNK_SIZE);

    let mut reader = crate::io::FastqReader::from_path(&args.input[0], args.fasta)?;

    // Temp files for streaming compressed blocks
    let working_dir = &args.working_dir;
    let h_tmp_path = working_dir.join(".qz_factorize_h.tmp");
    let m_tmp_path = working_dir.join(".qz_factorize_m.tmp");
    let rs_tmp_path = working_dir.join(".qz_factorize_rs.tmp");
    let ds_tmp_path = working_dir.join(".qz_factorize_ds.tmp");
    let q_tmp_path = working_dir.join(".qz_factorize_q.tmp");

    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 { let _ = std::fs::remove_file(p); }
        }
    }
    let _cleanup = TmpCleanup(vec![
        h_tmp_path.clone(), m_tmp_path.clone(), rs_tmp_path.clone(),
        ds_tmp_path.clone(), q_tmp_path.clone(),
    ]);

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut m_tmp = BufWriter::new(std::fs::File::create(&m_tmp_path)?);
    let mut rs_tmp = BufWriter::new(std::fs::File::create(&rs_tmp_path)?);
    let mut ds_tmp = BufWriter::new(std::fs::File::create(&ds_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);

    let mut h_num_blocks: u32 = 0;
    let mut m_num_blocks: u32 = 0;
    let mut rs_num_blocks: u32 = 0;
    let mut ds_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;
    let mut num_reads: usize = 0;
    let mut original_size: usize = 0;
    let mut chunk_idx: usize = 0;
    let mut total_delta: usize = 0;
    let mut total_fwd: usize = 0;
    let mut total_rc: usize = 0;

    // Read first chunk
    let (mut cur_records, _, mut cur_orig) = read_chunk_records(&mut reader, CHUNK_SIZE)?;

    while !cur_records.is_empty() {
        let cur_reads = cur_records.len();
        info!("Chunk {}: {} reads", chunk_idx, cur_reads);

        // Pipeline: compress current chunk while reading next
        let (next_result, compress_result) = std::thread::scope(|scope| {
            let records = std::mem::take(&mut cur_records);

            let compress_handle = scope.spawn(move || -> Result<ChunkCompressed> {
                compress_chunk(records, quality_mode, quality_binning, no_quality, bsc_static)
            });

            // Read next chunk on main thread (overlapped with compression)
            let next = read_chunk_records(&mut reader, CHUNK_SIZE);

            let compressed = compress_handle.join().unwrap();
            (next, compressed)
        });

        let chunk = compress_result?;
        h_num_blocks += write_blocks_to_tmp(chunk.h_blocks, &mut h_tmp)?;
        m_num_blocks += write_blocks_to_tmp(chunk.m_blocks, &mut m_tmp)?;
        rs_num_blocks += write_blocks_to_tmp(chunk.rs_blocks, &mut rs_tmp)?;
        if !chunk.ds_blocks.is_empty() {
            ds_num_blocks += write_blocks_to_tmp(chunk.ds_blocks, &mut ds_tmp)?;
        }
        if !chunk.q_blocks.is_empty() {
            q_num_blocks += write_blocks_to_tmp(chunk.q_blocks, &mut q_tmp)?;
        }

        total_delta += chunk.delta_count;
        total_fwd += chunk.fwd_count;
        total_rc += chunk.rc_count;
        num_reads += cur_reads;
        original_size += cur_orig;
        chunk_idx += 1;

        let (nr, _, no) = next_result?;
        cur_records = nr;
        cur_orig = no;
    }

    // Flush temp files
    h_tmp.flush()?; m_tmp.flush()?; rs_tmp.flush()?; ds_tmp.flush()?; q_tmp.flush()?;
    drop(h_tmp); drop(m_tmp); drop(rs_tmp); drop(ds_tmp); drop(q_tmp);

    info!("Read {} records in {} chunks", num_reads, chunk_idx);
    info!("  {}/{} reads delta-encoded ({:.1}%), {} fwd + {} rc",
        total_delta, num_reads, 100.0 * total_delta as f64 / num_reads.max(1) as f64,
        total_fwd, total_rc);

    // Get temp file sizes
    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let m_tmp_size = std::fs::metadata(&m_tmp_path)?.len() as usize;
    let rs_tmp_size = std::fs::metadata(&rs_tmp_path)?.len() as usize;
    let ds_tmp_size = std::fs::metadata(&ds_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    let h_size = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let m_size = if m_num_blocks > 0 { 4 + m_tmp_size } else { 0 };
    let rs_size = if rs_num_blocks > 0 { 4 + rs_tmp_size } else { 0 };
    let ds_size = if ds_num_blocks > 0 { 4 + ds_tmp_size } else { 0 };
    let q_size = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    // Write final archive
    write_chunked_factorize_archive(
        &args.output, quality_binning, no_quality,
        num_reads, original_size, start_time,
        h_size, m_size, rs_size, ds_size, q_size,
        h_num_blocks, m_num_blocks, rs_num_blocks, ds_num_blocks, q_num_blocks,
        &h_tmp_path, &m_tmp_path, &rs_tmp_path, &ds_tmp_path, &q_tmp_path,
    )
}

// ── Per-chunk compression ─────────────────────────────────────────────────

struct ChunkCompressed {
    h_blocks: Vec<Vec<u8>>,
    m_blocks: Vec<Vec<u8>>,
    rs_blocks: Vec<Vec<u8>>,
    ds_blocks: Vec<Vec<u8>>,
    q_blocks: Vec<Vec<u8>>,
    delta_count: usize,
    fwd_count: usize,
    rc_count: usize,
}

fn compress_chunk(
    records: Vec<crate::io::FastqRecord>,
    quality_mode: crate::cli::QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    bsc_static: bool,
) -> Result<ChunkCompressed> {
    let n = records.len();

    // Pass 1: Build inverted index for this chunk
    let index = build_read_index(&records);

    // Pass 2: Encode reads (parallel)
    use rayon::prelude::*;

    let typical_len = records.iter().take(1000).map(|r| r.sequence.len()).max().unwrap_or(150);
    let max_dist = (typical_len as f64 * 0.30) as usize;

    struct ReadEncoded {
        header: Vec<u8>,
        meta: Vec<u8>,
        raw_seq: Vec<u8>,
        delta_seq: Vec<u8>,
        quality: Vec<u8>,
        is_delta: bool,
        orientation: u8,
    }

    let encoded: Vec<ReadEncoded> = records.par_iter().enumerate().map(|(i, record)| {
        let mut header = Vec::new();
        let mut meta = Vec::new();
        let mut raw_seq = Vec::new();
        let mut delta_seq = Vec::new();
        let mut quality = Vec::new();

        dna_utils::write_varint(&mut header, record.id.len() as u64);
        header.extend_from_slice(record.id.as_bytes());

        if !no_quality {
            if let Some(ref qual) = record.quality {
                let quantized = quality::quantize_quality(qual, quality_mode);
                dna_utils::write_varint(&mut quality, quantized.len() as u64);
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                quality.extend_from_slice(&packed);
            }
        }

        let read_seq = record.sequence.as_bytes();
        let best_ref = find_best_reference(i, read_seq, &records, &index, max_dist);

        match best_ref {
            Some((ref_idx, orientation, _hamming)) => {
                meta.push(orientation);
                dna_utils::write_varint(&mut meta, ref_idx as u64);

                let ref_seq = records[ref_idx].sequence.as_bytes();
                dna_utils::write_varint(&mut delta_seq, read_seq.len() as u64);

                if orientation == 1 {
                    for (&rb, &refb) in read_seq.iter().zip(ref_seq.iter()) {
                        delta_seq.push(if rb == refb { 0x00 } else { rb });
                    }
                    for &rb in &read_seq[ref_seq.len().min(read_seq.len())..] {
                        delta_seq.push(rb);
                    }
                } else {
                    let rn = ref_seq.len();
                    for (j, &rb) in read_seq.iter().enumerate() {
                        let rc_base = if j < rn { complement(ref_seq[rn - 1 - j]) } else { 0xFF };
                        delta_seq.push(if rb == rc_base { 0x00 } else { rb });
                    }
                }

                ReadEncoded { header, meta, raw_seq, delta_seq, quality,
                    is_delta: true, orientation }
            }
            None => {
                meta.push(0x00);
                dna_utils::write_varint(&mut raw_seq, read_seq.len() as u64);
                raw_seq.extend_from_slice(read_seq);
                ReadEncoded { header, meta, raw_seq, delta_seq, quality,
                    is_delta: false, orientation: 0 }
            }
        }
    }).collect();

    // Concatenate into streams
    let mut header_stream = Vec::with_capacity(n * 64);
    let mut meta_stream = Vec::with_capacity(n * 4);
    let mut raw_seq_stream = Vec::with_capacity(n * 152);
    let mut delta_seq_stream = Vec::new();
    let mut qual_stream = Vec::with_capacity(if no_quality { 0 } else { n * 136 });
    let mut delta_count = 0usize;
    let mut fwd_count = 0usize;
    let mut rc_count = 0usize;

    for enc in &encoded {
        header_stream.extend_from_slice(&enc.header);
        meta_stream.extend_from_slice(&enc.meta);
        raw_seq_stream.extend_from_slice(&enc.raw_seq);
        delta_seq_stream.extend_from_slice(&enc.delta_seq);
        qual_stream.extend_from_slice(&enc.quality);
        if enc.is_delta {
            delta_count += 1;
            if enc.orientation == 1 { fwd_count += 1; } else { rc_count += 1; }
        }
    }
    drop(encoded);

    info!("  {}/{} delta ({:.1}%), {} fwd + {} rc, meta: {} B, raw_seq: {} B, delta_seq: {} B",
        delta_count, n, 100.0 * delta_count as f64 / n as f64,
        fwd_count, rc_count,
        meta_stream.len(), raw_seq_stream.len(), delta_seq_stream.len());

    // BSC-compress streams (sequentially to limit memory, like chunked mode)
    let h_blocks = compress_stream_to_bsc_blocks(&header_stream, bsc_static)?;
    drop(header_stream);
    let m_blocks = compress_stream_to_bsc_blocks(&meta_stream, bsc_static)?;
    drop(meta_stream);
    let rs_blocks = compress_stream_to_bsc_blocks(&raw_seq_stream, bsc_static)?;
    drop(raw_seq_stream);
    let ds_blocks = if delta_seq_stream.is_empty() { Vec::new() }
                    else { compress_stream_to_bsc_blocks(&delta_seq_stream, bsc_static)? };
    drop(delta_seq_stream);
    let q_blocks = if no_quality || qual_stream.is_empty() { Vec::new() }
                   else { compress_stream_to_bsc_blocks(&qual_stream, bsc_static)? };
    drop(qual_stream);

    Ok(ChunkCompressed { h_blocks, m_blocks, rs_blocks, ds_blocks, q_blocks,
        delta_count, fwd_count, rc_count })
}

// ── Inverted index ────────────────────────────────────────────────────────

fn build_read_index(records: &[crate::io::FastqRecord]) -> rustc_hash::FxHashMap<u64, Vec<u32>> {
    let mut index: rustc_hash::FxHashMap<u64, Vec<u32>> =
        rustc_hash::FxHashMap::with_capacity_and_hasher(4_000_000, Default::default());

    for (i, record) in records.iter().enumerate() {
        let seq = record.sequence.as_bytes();
        if seq.len() < INDEX_K { continue; }
        let mut pos = 0;
        while pos + INDEX_K <= seq.len() {
            if let Some(fwd) = dna_utils::kmer_to_hash(&seq[pos..pos + INDEX_K]) {
                let rc = dna_utils::reverse_complement_hash(fwd, INDEX_K);
                let canonical = fwd.min(rc);
                let bucket = index.entry(canonical).or_default();
                if bucket.len() < BUCKET_CAP {
                    bucket.push(i as u32);
                }
            }
            pos += INDEX_STRIDE;
        }
    }
    index
}

fn find_best_reference(
    read_idx: usize,
    read: &[u8],
    records: &[crate::io::FastqRecord],
    index: &rustc_hash::FxHashMap<u64, Vec<u32>>,
    max_dist: usize,
) -> Option<(usize, u8, usize)> {
    if read.len() < INDEX_K { return None; }

    let mut candidate_counts: rustc_hash::FxHashMap<u32, u16> = rustc_hash::FxHashMap::default();

    let mut pos = 0;
    while pos + INDEX_K <= read.len() {
        if let Some(fwd) = dna_utils::kmer_to_hash(&read[pos..pos + INDEX_K]) {
            let rc = dna_utils::reverse_complement_hash(fwd, INDEX_K);
            let canonical = fwd.min(rc);
            if let Some(bucket) = index.get(&canonical) {
                for &j in bucket {
                    if (j as usize) < read_idx {
                        *candidate_counts.entry(j).or_insert(0) += 1;
                    }
                }
            }
        }
        pos += INDEX_STRIDE;
    }

    if candidate_counts.is_empty() { return None; }

    let mut candidates: Vec<(u32, u16)> = candidate_counts.into_iter().collect();
    candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    candidates.truncate(8);

    let mut best: Option<(usize, u8, usize)> = None;

    for (j, _score) in candidates {
        let ref_seq = records[j as usize].sequence.as_bytes();
        if ref_seq.len() != read.len() { continue; }

        let threshold = best.map(|b| b.2).unwrap_or(max_dist);
        let dist_fwd = hamming_capped(read, ref_seq, threshold);
        if dist_fwd <= threshold {
            best = Some((j as usize, 1, dist_fwd));
        }

        let threshold = best.map(|b| b.2).unwrap_or(max_dist);
        let dist_rc = hamming_rc_capped(read, ref_seq, threshold);
        if dist_rc < threshold {
            best = Some((j as usize, 2, dist_rc));
        }
    }

    if let Some((_, _, dist)) = best {
        let matching = read.len() - dist;
        if matching < 10 { return None; }
    }

    best
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        x => x,
    }
}

fn hamming_capped(a: &[u8], b: &[u8], cap: usize) -> usize {
    let mut dist = 0;
    for (&x, &y) in a.iter().zip(b.iter()) {
        if x != y { dist += 1; if dist > cap { return dist; } }
    }
    dist
}

fn hamming_rc_capped(a: &[u8], b: &[u8], cap: usize) -> usize {
    let n = a.len().min(b.len());
    let mut dist = 0;
    for i in 0..n {
        let rc_base = complement(b[n - 1 - i]);
        if a[i] != rc_base { dist += 1; if dist > cap { return dist; } }
    }
    dist
}

// ── Chunked archive writer ───────────────────────────────────────────────

/// Write chunked factorize archive (encoding_type=7).
///
/// sequences_region contains 3 sub-stream groups: meta, raw_seq, delta_seq.
/// Each group has: [u64 total_size] [u32 num_blocks] [blocks...]
fn write_chunked_factorize_archive(
    output_path: &std::path::Path,
    quality_binning: QualityBinning,
    _no_quality: bool,
    num_reads: usize,
    original_size: usize,
    start_time: Instant,
    h_size: usize, m_size: usize, rs_size: usize, ds_size: usize, q_size: usize,
    h_num_blocks: u32, m_num_blocks: u32, rs_num_blocks: u32, ds_num_blocks: u32, q_num_blocks: u32,
    h_tmp: &std::path::Path, m_tmp: &std::path::Path,
    rs_tmp: &std::path::Path, ds_tmp: &std::path::Path, q_tmp: &std::path::Path,
) -> Result<()> {
    use std::io::{Write, BufWriter};

    // sequences_len = 3 sub-stream headers (u64 each) + sub-stream data
    let factorized_seq_len = 3 * 8 + m_size + rs_size + ds_size;

    info!("Writing output file...");
    let mut out = BufWriter::new(std::fs::File::create(output_path)?);

    // Archive header (encoding_type=7)
    out.write_all(&[7u8])?;
    out.write_all(&[0u8])?;
    out.write_all(&[binning_to_code(quality_binning)])?;
    out.write_all(&[compressor_to_code(QualityCompressor::Bsc)])?;
    out.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?;
    out.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;
    out.write_all(&[0u8; 3])?;
    out.write_all(&0u16.to_le_bytes())?;
    out.write_all(&[0u8])?;

    out.write_all(&(num_reads as u64).to_le_bytes())?;
    out.write_all(&(h_size as u64).to_le_bytes())?;
    out.write_all(&(factorized_seq_len as u64).to_le_bytes())?;
    out.write_all(&0u64.to_le_bytes())?;
    out.write_all(&(q_size as u64).to_le_bytes())?;

    // Helper: copy num_blocks header + tmp file data to output
    let copy_stream = |num_blocks: u32, tmp_path: &std::path::Path, out: &mut BufWriter<std::fs::File>| -> Result<()> {
        if num_blocks == 0 { return Ok(()); }
        out.write_all(&num_blocks.to_le_bytes())?;
        let mut f = std::fs::File::open(tmp_path)?;
        std::io::copy(&mut f, out)?;
        Ok(())
    };

    // Headers
    copy_stream(h_num_blocks, h_tmp, &mut out)?;

    // Factorized sub-streams (each with u64 size prefix)
    out.write_all(&(m_size as u64).to_le_bytes())?;
    copy_stream(m_num_blocks, m_tmp, &mut out)?;
    out.write_all(&(rs_size as u64).to_le_bytes())?;
    copy_stream(rs_num_blocks, rs_tmp, &mut out)?;
    out.write_all(&(ds_size as u64).to_le_bytes())?;
    copy_stream(ds_num_blocks, ds_tmp, &mut out)?;

    // Qualities
    copy_stream(q_num_blocks, q_tmp, &mut out)?;

    out.flush()?;

    let metadata_size = 9 + 2 + 1 + 40;
    let total = h_size + factorized_seq_len + q_size + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total);
    info!("Compression ratio: {:.2}x", original_size as f64 / total as f64);
    info!("Stream breakdown:");
    info!("  Headers:    {} bytes ({:.1}%)", h_size, 100.0 * h_size as f64 / total as f64);
    info!("  Meta:       {} bytes ({:.1}%)", m_size, 100.0 * m_size as f64 / total as f64);
    info!("  Raw seqs:   {} bytes ({:.1}%)", rs_size, 100.0 * rs_size as f64 / total as f64);
    info!("  Delta seqs: {} bytes ({:.1}%)", ds_size, 100.0 * ds_size as f64 / total as f64);
    info!("  Qualities:  {} bytes ({:.1}%)", q_size, 100.0 * q_size as f64 / total as f64);
    info!("  Metadata:   {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total as f64);

    Ok(())
}

// ── Decoder ───────────────────────────────────────────────────────────────

/// Decode factorized sequences from 3 sub-streams (meta, raw_seq, delta_seq).
///
/// Must process sequentially since delta reads reference earlier decoded reads.
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
            let seq_len = read_varint(raw_seq_data, &mut raw_off)
                .ok_or_else(|| anyhow::anyhow!("Failed to read raw seq_len for read {}", i))?;
            if raw_off + seq_len > raw_seq_data.len() {
                anyhow::bail!("Truncated raw seq data for read {}", i);
            }
            let seq = raw_seq_data[raw_off..raw_off + seq_len].to_vec();
            raw_off += seq_len;
            sequences.push(seq);
        } else {
            let ref_idx = read_varint(meta_data, &mut meta_off)
                .ok_or_else(|| anyhow::anyhow!("Failed to read ref_index for read {}", i))?;
            if ref_idx >= sequences.len() {
                anyhow::bail!("Invalid ref_index {} for read {} (only {} decoded)",
                    ref_idx, i, sequences.len());
            }

            let seq_len = read_varint(delta_seq_data, &mut delta_off)
                .ok_or_else(|| anyhow::anyhow!("Failed to read delta seq_len for read {}", i))?;
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
