//! Delta encoding for DNA sequences
//!
//! Encodes sequences as differences from a reference, which can improve
//! compression for highly similar reads (e.g., PCR duplicates, technical replicates)

use anyhow::Result;

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

        // Encode as: [prefix_len:u8][suffix_len:u8][middle_bytes]
        let middle_start = prefix_len;
        let middle_end = curr.len().saturating_sub(suffix_len);

        let mut delta = Vec::new();
        delta.push(prefix_len.min(255) as u8);
        delta.push(suffix_len.min(255) as u8);
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

        if delta.len() < 2 {
            anyhow::bail!("Invalid delta encoding");
        }

        let prefix_len = delta[0] as usize;
        let suffix_len = delta[1] as usize;
        let middle = &delta[2..];

        // Reconstruct: prefix + middle + suffix
        let mut reconstructed = Vec::new();
        reconstructed.extend_from_slice(&prev[..prefix_len.min(prev.len())]);
        reconstructed.extend_from_slice(middle);
        if suffix_len > 0 && prev.len() >= suffix_len {
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
}
