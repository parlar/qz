//! Delta encoding for DNA sequences
//!
//! Encodes sequences as differences from a reference, which can improve
//! compression for highly similar reads (e.g., PCR duplicates, technical replicates)

use anyhow::Result;

/// Write a variable-length integer to a byte buffer.
fn push_varint(buf: &mut Vec<u8>, mut value: usize) {
    while value >= 0x80 {
        buf.push(((value & 0x7F) | 0x80) as u8);
        value >>= 7;
    }
    buf.push(value as u8);
}

/// Read a variable-length integer from a byte slice at the given offset.
/// Returns None if the data is truncated or the varint exceeds 10 bytes.
fn read_varint(data: &[u8], offset: &mut usize) -> Option<usize> {
    let mut value = 0usize;
    let mut shift = 0;
    loop {
        if *offset >= data.len() || shift >= 70 {
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

/// Apply simple delta encoding to sequences before compression
/// Finds common prefixes and suffixes to reduce redundancy
pub fn apply_delta_encoding(sequences: &[String]) -> Vec<Vec<u8>> {
    if sequences.is_empty() {
        return Vec::new();
    }

    let mut encoded = Vec::with_capacity(sequences.len());

    // First sequence is the reference
    encoded.push(sequences[0].as_bytes().to_vec());

    // Encode subsequent sequences as deltas from previous
    for i in 1..sequences.len() {
        let prev = sequences[i - 1].as_bytes();
        let curr = sequences[i].as_bytes();

        // Find common prefix length
        let prefix_len = prev.iter()
            .zip(curr.iter())
            .take_while(|(a, b)| a == b)
            .count();

        // Find common suffix length (excluding prefix)
        let suffix_len = if prefix_len < prev.len().min(curr.len()) {
            prev[prefix_len..].iter().rev()
                .zip(curr[prefix_len..].iter().rev())
                .take_while(|(a, b)| a == b)
                .count()
        } else {
            0
        };

        // Encode as: [prefix_len:varint][suffix_len:varint][middle_bytes]
        let middle_start = prefix_len;
        let middle_end = curr.len().saturating_sub(suffix_len);

        let mut delta = Vec::new();
        push_varint(&mut delta, prefix_len);
        push_varint(&mut delta, suffix_len);
        if middle_end > middle_start {
            delta.extend_from_slice(&curr[middle_start..middle_end]);
        }

        encoded.push(delta);
    }

    encoded
}

/// Decode delta-encoded sequences
pub fn decode_delta_encoding(encoded: &[Vec<u8>]) -> Result<Vec<String>> {
    if encoded.is_empty() {
        return Ok(Vec::new());
    }

    let mut decoded = Vec::with_capacity(encoded.len());

    // First sequence is the reference
    decoded.push(String::from_utf8_lossy(&encoded[0]).to_string());

    // Decode subsequent sequences
    for i in 1..encoded.len() {
        let prev = decoded[i - 1].as_bytes();
        let delta = &encoded[i];

        let mut offset = 0;
        let prefix_len = read_varint(delta, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Invalid delta encoding: truncated prefix length"))?;
        let suffix_len = read_varint(delta, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Invalid delta encoding: truncated suffix length"))?;
        let middle = &delta[offset..];

        // Validate lengths against previous sequence
        if prefix_len > prev.len() {
            anyhow::bail!(
                "Invalid delta encoding: prefix_len {} exceeds previous sequence length {}",
                prefix_len, prev.len()
            );
        }
        if suffix_len > prev.len().saturating_sub(prefix_len) {
            anyhow::bail!(
                "Invalid delta encoding: suffix_len {} exceeds available bytes in previous sequence",
                suffix_len
            );
        }

        // Reconstruct: prefix + middle + suffix
        let mut reconstructed = Vec::with_capacity(prefix_len + middle.len() + suffix_len);
        reconstructed.extend_from_slice(&prev[..prefix_len]);
        reconstructed.extend_from_slice(middle);
        if suffix_len > 0 {
            reconstructed.extend_from_slice(&prev[prev.len() - suffix_len..]);
        }

        decoded.push(String::from_utf8_lossy(&reconstructed).to_string());
    }

    Ok(decoded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_identical() {
        let sequences = vec![
            "ACGTACGT".to_string(),
            "ACGTACGT".to_string(),
        ];
        let encoded = apply_delta_encoding(&sequences);
        let decoded = decode_delta_encoding(&encoded).unwrap();
        assert_eq!(sequences, decoded);
    }

    #[test]
    fn test_delta_similar() {
        let sequences = vec![
            "ACGTACGTACGT".to_string(),
            "ACGTNNNNACGT".to_string(), // Same prefix/suffix, different middle
        ];
        let encoded = apply_delta_encoding(&sequences);
        let decoded = decode_delta_encoding(&encoded).unwrap();
        assert_eq!(sequences, decoded);
    }

    #[test]
    fn test_delta_different() {
        let sequences = vec![
            "AAAA".to_string(),
            "TTTT".to_string(),
        ];
        let encoded = apply_delta_encoding(&sequences);
        let decoded = decode_delta_encoding(&encoded).unwrap();
        assert_eq!(sequences, decoded);
    }

    #[test]
    fn test_delta_long_reads() {
        // Reads >255 bp with shared prefix/suffix that exceed u8 range
        let base = "ACGT".repeat(200); // 800 bp
        let mut variant = base.clone();
        // Change one base in the middle
        unsafe { variant.as_bytes_mut()[400] = b'N'; }
        let sequences = vec![base, variant];
        let encoded = apply_delta_encoding(&sequences);
        let decoded = decode_delta_encoding(&encoded).unwrap();
        assert_eq!(sequences, decoded);
    }

    #[test]
    fn test_delta_empty_middle() {
        let sequences = vec![
            "ACGTACGT".to_string(),
            "ACGTACGT".to_string(), // identical
        ];
        let encoded = apply_delta_encoding(&sequences);
        let decoded = decode_delta_encoding(&encoded).unwrap();
        assert_eq!(sequences, decoded);
        // Delta should just be two varints (prefix=8, suffix=0) and no middle
        assert!(encoded[1].len() <= 3);
    }
}
