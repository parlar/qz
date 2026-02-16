//! Run-Length Encoding (RLE) for DNA sequences
//!
//! Compresses homopolymer runs (e.g., AAAAA → A5)
//! Useful for sequences with long repeats like PacBio/Nanopore data

use anyhow::Result;

/// Apply run-length encoding to a DNA sequence
/// Format: base byte, then if run > 1: (0x80 | min(run, 127)), and if run > 127
/// additional continuation bytes.
pub fn apply_rle_encoding(sequence: &str) -> Vec<u8> {
    if sequence.is_empty() {
        return Vec::new();
    }

    let bytes = sequence.as_bytes();
    let mut encoded = Vec::new();
    let mut i = 0;

    while i < bytes.len() {
        let base = bytes[i];
        let mut run_length = 1usize;

        // Count consecutive identical bases (no cap)
        while i + run_length < bytes.len() && bytes[i + run_length] == base {
            run_length += 1;
        }

        // Store base
        encoded.push(base);

        // Store run length if > 1 using a variable-length scheme:
        // First byte: 0x80 | (min(remaining, 127))
        // If remaining > 127, additional bytes follow with the same pattern
        // A byte without 0x80 set (or no byte) means no more run length data
        if run_length > 1 {
            let mut remaining = run_length;
            loop {
                if remaining <= 127 {
                    encoded.push(0x80 | (remaining as u8));
                    break;
                } else {
                    // Signal "127 more and keep going" with 0xFF
                    encoded.push(0xFF);
                    remaining -= 127;
                }
            }
        }

        i += run_length;
    }

    encoded
}

/// Decode run-length encoded sequence
pub fn decode_rle_encoding(encoded: &[u8]) -> Result<String> {
    let mut decoded = Vec::new();
    let mut i = 0;

    while i < encoded.len() {
        let base = encoded[i];
        i += 1;

        // Accumulate run length from continuation bytes
        let mut run_length = 1usize;
        if i < encoded.len() && (encoded[i] & 0x80) != 0 {
            run_length = 0;
            while i < encoded.len() && (encoded[i] & 0x80) != 0 {
                let chunk = (encoded[i] & 0x7F) as usize;
                run_length += chunk;
                i += 1;
                if chunk < 127 {
                    break; // Last run-length byte
                }
            }
        }

        // Expand run
        decoded.resize(decoded.len() + run_length, base);
    }

    Ok(String::from_utf8_lossy(&decoded).to_string())
}

/// Apply RLE to multiple sequences
pub fn apply_rle_to_sequences(sequences: &[String]) -> Vec<Vec<u8>> {
    sequences.iter()
        .map(|seq| apply_rle_encoding(seq))
        .collect()
}

/// Decode RLE from multiple sequences
pub fn decode_rle_from_sequences(encoded: &[Vec<u8>]) -> Result<Vec<String>> {
    encoded.iter()
        .map(|enc| decode_rle_encoding(enc))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rle_no_runs() {
        let seq = "ACGTACGT";
        let encoded = apply_rle_encoding(seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }

    #[test]
    fn test_rle_with_runs() {
        let seq = "AAAAACCCGGGTTT";
        let encoded = apply_rle_encoding(seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
        assert!(encoded.len() < seq.len());
    }

    #[test]
    fn test_rle_long_run() {
        let seq = "A".repeat(100);
        let encoded = apply_rle_encoding(&seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }

    #[test]
    fn test_rle_very_long_run() {
        // Run exceeding old 255 limit
        let seq = "A".repeat(500);
        let encoded = apply_rle_encoding(&seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }

    #[test]
    fn test_rle_extreme_run() {
        // Very long homopolymer (e.g. PacBio/Nanopore)
        let seq = "T".repeat(10_000);
        let encoded = apply_rle_encoding(&seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
        // Should be much smaller than original
        assert!(encoded.len() < 200);
    }

    #[test]
    fn test_rle_mixed() {
        let seq = "ACGTAAAACCCCNNNN";
        let encoded = apply_rle_encoding(seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }

    #[test]
    fn test_rle_single_base() {
        let seq = "A";
        let encoded = apply_rle_encoding(seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }
}
