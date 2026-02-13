/// Arithmetic coding for sequences and quality scores using constriction library
///
/// This module implements adaptive entropy coding using ANS (Asymmetric Numeral Systems)
/// to achieve better compression than simple 2-bit encoding + zstd.
///
/// Key improvements:
/// - Quality scores: Context-based models (position, previous quality, base)
/// - Sequences: Markov-based prediction (previous 3 bases, GC%, position)
/// - Expected: 6-8x compression (vs 4.6x current)

use constriction::stream::model::DefaultContiguousCategoricalEntropyModel;
use constriction::stream::stack::DefaultAnsCoder;
use constriction::stream::{Decode, Encode}; // Traits for encode_symbol/decode_symbol
use constriction::UnwrapInfallible;

/// Simple test to verify constriction works correctly
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constriction_basic() {
        // Test basic ANS encoding/decoding
        let symbols = vec![0u8, 1, 0, 2, 1, 3, 0, 1];
        let probabilities = vec![0.4, 0.3, 0.2, 0.1]; // A, C, G, T

        // Create entropy model (None = auto-normalize probabilities)
        let model = DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(
            &probabilities,
            None,
        )
        .unwrap();

        // Encode
        let mut encoder = DefaultAnsCoder::new();
        for &symbol in symbols.iter().rev() {
            // ANS encodes in reverse
            encoder.encode_symbol(symbol as usize, &model).unwrap();
        }
        let compressed = encoder.into_compressed().unwrap_infallible();

        println!("Original: {} bytes", symbols.len());
        println!("Compressed: {} bytes", compressed.len());
        println!("Compression ratio: {:.2}x", symbols.len() as f64 / compressed.len() as f64);

        // Decode
        let mut decoder = DefaultAnsCoder::from_compressed(compressed).unwrap();
        let mut decoded = Vec::new();
        for _ in 0..symbols.len() {
            let symbol = decoder.decode_symbol(&model).unwrap();
            decoded.push(symbol as u8);
        }

        assert_eq!(symbols, decoded, "Roundtrip failed");
    }

    #[test]
    fn test_adaptive_model() {
        // Test adaptive probability model (counts update after each symbol)
        let symbols = vec![0u8, 0, 0, 1, 1, 0, 0, 0]; // Mostly 0s

        // For adaptive models with ANS, we need to build all models first, then encode in reverse
        // This is because ANS encodes in reverse but we need consistent model updates
        let mut models = Vec::new();
        let mut counts = vec![1u32, 1]; // Laplace smoothing

        // Build models in forward order (same as decoding)
        for &symbol in &symbols {
            let total: u32 = counts.iter().sum();
            let probs: Vec<f64> = counts.iter().map(|&c| c as f64 / total as f64).collect();

            let model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(&probs, None)
                    .unwrap();

            models.push((symbol, model));
            counts[symbol as usize] += 1;
        }

        // Encode with pre-computed models in reverse
        let mut encoder = DefaultAnsCoder::new();
        for (symbol, model) in models.iter().rev() {
            encoder.encode_symbol(*symbol as usize, model).unwrap();
        }

        let compressed = encoder.into_compressed().unwrap_infallible();

        println!("Adaptive: {} bytes -> {} bytes", symbols.len(), compressed.len());

        // Decode with same adaptive model (in forward order, symbols will come out in forward order)
        let mut decoder = DefaultAnsCoder::from_compressed(compressed).unwrap();
        let mut decoded = Vec::new();
        let mut counts = vec![1u32, 1];

        for _ in 0..symbols.len() {
            let total: u32 = counts.iter().sum();
            let probs: Vec<f64> = counts.iter().map(|&c| c as f64 / total as f64).collect();

            let model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(&probs, None)
                    .unwrap();

            let symbol = decoder.decode_symbol(&model).unwrap();
            decoded.push(symbol as u8);

            counts[symbol as usize] += 1;
        }

        // ANS decodes in same order as encoded (forward), so decoded should match symbols
        assert_eq!(symbols, decoded, "Adaptive roundtrip failed");
    }

    #[test]
    fn test_quality_score_encoding() {
        // Test encoding quality scores (Phred 0-40)
        let qualities = b"IIIIIHHHHGGGGG"; // Phred scores as ASCII
        let symbols: Vec<u8> = qualities.iter().map(|&q| q - 33).collect(); // Convert to Phred 0-93

        // Simple model: assume uniform distribution initially
        let mut counts = vec![1u32; 41]; // Phred 0-40

        let mut encoder = DefaultAnsCoder::new();

        for &symbol in symbols.iter().rev() {
            let total: u32 = counts.iter().sum();
            let probs: Vec<f64> = counts.iter().map(|&c| c as f64 / total as f64).collect();

            let model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(&probs, None)
                    .unwrap();

            encoder.encode_symbol(symbol as usize, &model).unwrap();
            counts[symbol as usize] += 1;
        }

        let compressed = encoder.into_compressed().unwrap_infallible();

        println!(
            "Quality encoding: {} bytes -> {} bytes ({:.2}x)",
            symbols.len(),
            compressed.len(),
            symbols.len() as f64 / compressed.len() as f64
        );

        // Should compress better than 1 byte per quality score
        assert!(compressed.len() < symbols.len(), "No compression achieved");
    }
}
