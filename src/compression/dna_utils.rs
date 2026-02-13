/// Shared DNA compression utilities used by spring_impl.rs and debruijn.rs

use anyhow::Result;

/// Compute Hamming distance between two sequences (must be same length)
pub fn hamming_distance(seq1: &[u8], seq2: &[u8]) -> usize {
    seq1.iter()
        .zip(seq2.iter())
        .filter(|(a, b)| a != b)
        .count()
}

/// Get reverse complement of a sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Hamming distance with early termination — returns None if distance exceeds threshold
pub fn hamming_distance_within(seq1: &[u8], seq2: &[u8], threshold: usize) -> Option<usize> {
    let mut dist = 0;
    for (&a, &b) in seq1.iter().zip(seq2.iter()) {
        if a != b {
            dist += 1;
            if dist > threshold {
                return None;
            }
        }
    }
    Some(dist)
}

/// Encode a k-mer as a 2-bit hash (4 bases → 8 bits, k=20 → 40 bits fits u64).
/// Returns None if the k-mer contains N or other ambiguous bases.
pub fn kmer_to_hash(seq: &[u8]) -> Option<u64> {
    let mut hash = 0u64;
    for &base in seq {
        let val = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        hash = (hash << 2) | val;
    }
    Some(hash)
}

/// Reverse complement a k-mer hash (A↔T, C↔G in 2-bit: 0↔3, 1↔2)
pub fn reverse_complement_hash(hash: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut h = hash;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (h & 3));
        h >>= 2;
    }
    rc
}

/// Convert a DNA base to a 0-3 index (A=0, C=1, G=2, T=3). N treated as A.
pub fn base_to_idx(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0, // N treated as A for voting
    }
}

/// Convert a 0-3 index back to a DNA base
pub fn idx_to_base(i: usize) -> u8 {
    match i { 0 => b'A', 1 => b'C', 2 => b'G', _ => b'T' }
}

/// Pack DNA sequences to 2-bit format (4 bases per byte)
/// - No N: 2-bit (4 bases/byte) A=0, C=1, G=2, T=3
/// - Has N: 4-bit (2 bases/byte) A=0, G=1, C=2, T=3, N=4 (C/G swapped!)
pub fn pack_dna_2bit(sequences: &[Vec<u8>]) -> Vec<u8> {
    let mut packed = Vec::new();

    // Write number of sequences
    write_varint(&mut packed, sequences.len() as u64);

    for seq in sequences {
        let len = seq.len();
        write_varint(&mut packed, len as u64);

        // Check if sequence contains N
        let has_n = seq.iter().any(|&b| !matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'));

        // Write encoding flag
        packed.push(if has_n { 1 } else { 0 });

        if has_n {
            // 4-bit encoding (2 bases per byte)
            // A=0, G=1, C=2, T=3, N=4
            for chunk in seq.chunks(2) {
                let mut byte = 0u8;
                for (i, &base) in chunk.iter().enumerate() {
                    let val = match base {
                        b'A' | b'a' => 0,
                        b'G' | b'g' => 1,
                        b'C' | b'c' => 2,
                        b'T' | b't' => 3,
                        _ => 4, // N
                    };
                    byte |= val << (4 * i); // 4 bits per base
                }
                packed.push(byte);
            }
        } else {
            // 2-bit encoding (4 bases per byte)
            // A=0, C=1, G=2, T=3
            for chunk in seq.chunks(4) {
                let mut byte = 0u8;
                for (i, &base) in chunk.iter().enumerate() {
                    let val = match base {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => 0, // Should never happen
                    };
                    byte |= val << (2 * i); // 2 bits per base
                }
                packed.push(byte);
            }
        }
    }

    packed
}

/// Unpack DNA sequences from 2-bit/4-bit packed format.
/// Returns (sequences, bytes_consumed)
pub fn unpack_dna_2bit(packed: &[u8], start_offset: usize) -> Result<(Vec<Vec<u8>>, usize)> {
    let mut offset = start_offset;
    let mut sequences = Vec::new();

    // Read number of sequences
    let num_seqs = read_varint(packed, &mut offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read sequence count"))? as usize;

    for _ in 0..num_seqs {
        let len = read_varint(packed, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length"))? as usize;

        if offset >= packed.len() {
            anyhow::bail!("Truncated: cannot read encoding flag");
        }

        // Read encoding flag
        let has_n = packed[offset] != 0;
        offset += 1;

        let mut seq = Vec::with_capacity(len);

        if has_n {
            // 4-bit decoding (2 bases per byte)
            // A=0, G=1, C=2, T=3, N=4
            let num_bytes = (len + 1) / 2;
            for byte_idx in 0..num_bytes {
                if offset >= packed.len() {
                    anyhow::bail!("Truncated 4-bit packed data");
                }

                let byte = packed[offset];
                offset += 1;

                let bases_in_byte = std::cmp::min(2, len - byte_idx * 2);
                for i in 0..bases_in_byte {
                    let val = (byte >> (4 * i)) & 0xF;
                    let base = match val {
                        0 => b'A',
                        1 => b'G',
                        2 => b'C',
                        3 => b'T',
                        4 => b'N',
                        _ => b'N',
                    };
                    seq.push(base);
                }
            }
        } else {
            // 2-bit decoding (4 bases per byte)
            // A=0, C=1, G=2, T=3
            let num_bytes = (len + 3) / 4;
            for byte_idx in 0..num_bytes {
                if offset >= packed.len() {
                    anyhow::bail!("Truncated 2-bit packed data");
                }

                let byte = packed[offset];
                offset += 1;

                let bases_in_byte = std::cmp::min(4, len - byte_idx * 4);
                for i in 0..bases_in_byte {
                    let val = (byte >> (2 * i)) & 0b11;
                    let base = match val {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        3 => b'T',
                        _ => b'N',
                    };
                    seq.push(base);
                }
            }
        }

        sequences.push(seq);
    }

    let bytes_consumed = offset - start_offset;
    Ok((sequences, bytes_consumed))
}

/// Write a variable-length integer (LEB128 encoding)
pub fn write_varint(output: &mut Vec<u8>, mut value: u64) {
    while value >= 128 {
        output.push((value & 0x7F) as u8 | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
}

/// Read a variable-length integer (LEB128 encoding)
pub fn read_varint(data: &[u8], offset: &mut usize) -> Option<u64> {
    let mut value = 0u64;
    let mut shift = 0;

    loop {
        if *offset >= data.len() {
            return None;
        }

        let byte = data[*offset];
        *offset += 1;

        value |= ((byte & 0x7F) as u64) << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
    }
}

/// Write a little-endian u64
pub fn write_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

/// Read a little-endian u64
pub fn read_u64(data: &[u8], offset: &mut usize) -> Option<u64> {
    if *offset + 8 > data.len() {
        return None;
    }
    let bytes: [u8; 8] = data[*offset..*offset + 8].try_into().ok()?;
    *offset += 8;
    Some(u64::from_le_bytes(bytes))
}
