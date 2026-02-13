//! Run-Length Encoding (RLE) for DNA sequences
//!
//! Compresses homopolymer runs (e.g., AAAAA → A5)
//! Useful for sequences with long repeats like PacBio/Nanopore data

use anyhow::Result;

/// Apply run-length encoding to a DNA sequence
/// Format: base byte + run length (if > 1)
pub fn apply_rle_encoding(sequence: &str) -> Vec<u8> {
    if sequence.is_empty() {
        return Vec::new();
    }

    let bytes = sequence.as_bytes();
    let mut encoded = Vec::new();
    let mut i = 0;

    while i < bytes.len() {
        let base = bytes[i];
        let mut run_length = 1;

        // Count consecutive identical bases
        while i + run_length < bytes.len() && bytes[i + run_length] == base {
            run_length += 1;
            if run_length >= 255 {
                break; // Max run length is 255
            }
        }

        // Store base
        encoded.push(base);

        // Store run length if > 1
        if run_length > 1 {
            encoded.push(0x80 | (run_length as u8)); // High bit set = run length follows
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
        if i >= encoded.len() {
            break;
        }

        let base = encoded[i];
        i += 1;

        // Check if next byte is a run length (high bit set)
        let run_length = if i < encoded.len() && (encoded[i] & 0x80) != 0 {
            let length = (encoded[i] & 0x7F) as usize;
            i += 1;
            length
        } else {
            1
        };

        // Expand run
        for _ in 0..run_length {
            decoded.push(base);
        }
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
        // Encoding should be smaller: A + 5, C + 3, G + 3, T + 3 = 8 bytes vs 14 original
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
    fn test_rle_mixed() {
        let seq = "ACGTAAAACCCCNNNN";
        let encoded = apply_rle_encoding(seq);
        let decoded = decode_rle_encoding(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }
}
