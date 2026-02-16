use crate::cli::QualityMode;
use std::borrow::Cow;
use tracing::info;

/// Illumina 8-level binning table
const ILLUMINA_BINS: [u8; 94] = [
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 0-9
    39, 39, 39, 39, 39, 39, 39, 39, 39, 39, // 10-19
    48, 48, 48, 48, 48, 48, 48, 48, 48, 48, // 20-29
    55, 55, 55, 55, 55, 55, 55, 55, 55, 55, // 30-39
    60, 60, 60, 60, 60, 60, 60, 60, 60, 60, // 40-49
    66, 66, 66, 66, 66, 66, 66, 66, 66, 66, // 50-59
    70, 70, 70, 70, 70, 70, 70, 70, 70, 70, // 60-69
    73, 73, 73, 73, 73, 73, 73, 73, 73, 73, // 70-79
    73, 73, 73, 73, 73, 73, 73, 73, 73, 73, // 80-89
    73, 73, 73, 73, // 90-93
];

/// Quantize quality string based on mode
pub fn quantize_quality<'a>(quality: &'a str, mode: QualityMode) -> Cow<'a, str> {
    match mode {
        QualityMode::Lossless => Cow::Borrowed(quality),
        QualityMode::IlluminaBin => Cow::Owned(quantize_illumina(quality)),
        QualityMode::Illumina4 => Cow::Owned(quantize_four_level(quality)),
        QualityMode::Binary => Cow::Owned(quantize_binary(quality, 20, 40, 6)),
        QualityMode::Qvz => {
            // TODO: Implement QVZ compression
            info!("QVZ mode not yet implemented, using lossless");
            Cow::Borrowed(quality)
        }
        QualityMode::Discard => {
            // Discard quality: replace all with minimum quality '!'
            Cow::Owned("!".repeat(quality.len()))
        }
    }
}

/// Apply 4-level binning (more aggressive than Illumina 8-level)
fn quantize_four_level(quality: &str) -> String {
    quality
        .bytes()
        .map(|qv| {
            let phred = (qv as i32).saturating_sub(33).max(0).min(93);
            let binned = if phred < 10 {
                6  // Low quality
            } else if phred < 20 {
                15 // Medium quality
            } else if phred < 30 {
                25 // Good quality
            } else {
                37 // High quality
            };
            (binned + 33) as u8 as char
        })
        .collect()
}

/// Apply Illumina 8-level binning
fn quantize_illumina(quality: &str) -> String {
    quality
        .bytes()
        .map(|qv| {
            let idx = (qv as usize).saturating_sub(33).min(93);
            ILLUMINA_BINS[idx] as char
        })
        .collect()
}

/// Binary thresholding: high if >= threshold, low otherwise
fn quantize_binary(quality: &str, threshold: u8, high: u8, low: u8) -> String {
    quality
        .bytes()
        .map(|qv| {
            let phred = qv.saturating_sub(33);
            if phred >= threshold {
                (33 + high) as char
            } else {
                (33 + low) as char
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lossless() {
        let quality = "IIIIIIII";
        assert_eq!(quantize_quality(quality, QualityMode::Lossless), quality);
    }

    #[test]
    fn test_illumina_binning() {
        // Quality 40 (Phred 7) should map to bin
        let quality = "!!!!!!!!"; // Phred 0
        let binned = quantize_illumina(quality);
        assert_eq!(binned.chars().next().unwrap(), '!'); // Should map to lowest bin
    }

    #[test]
    fn test_binary_threshold() {
        // Mix of high and low quality
        let quality = "!!!!IIII"; // Phred 0 and 40
        let binned = quantize_binary(quality, 20, 40, 6);
        // First 4 should be low (6), last 4 should be high (40)
        assert!(binned.starts_with("''''"));
        assert!(binned.ends_with("IIII"));
    }
}
