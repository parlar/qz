use crate::io::FastqRecord;
use anyhow::Result;
use rayon::prelude::*;
use std::cmp::Ordering;
use tracing::info;

/// Compute Hamming distance between two sequences (for equal-length sequences)
fn hamming_distance(s1: &str, s2: &str) -> usize {
    s1.chars()
        .zip(s2.chars())
        .filter(|(a, b)| a != b)
        .count()
}

/// Compute a hash for k-mer based clustering
fn compute_kmer_hash(seq: &str, k: usize) -> u64 {
    if seq.len() < k {
        return 0;
    }

    // Use first k-mer as a simple hash
    let kmer = &seq[..k];
    let mut hash = 0u64;

    for (i, c) in kmer.chars().enumerate() {
        let val = match c {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => 0,
        };
        hash |= val << (i * 2);
    }

    hash
}

/// Order records by k-mer similarity for better compression
pub fn reorder_by_similarity(mut records: Vec<FastqRecord>, k: usize) -> Result<Vec<FastqRecord>> {
    if records.is_empty() {
        return Ok(records);
    }

    info!(
        "Reordering {} records by {}-mer similarity...",
        records.len(),
        k
    );

    // Compute k-mer hashes in parallel
    let mut indexed_records: Vec<_> = records
        .par_iter()
        .enumerate()
        .map(|(idx, record)| {
            let hash = compute_kmer_hash(&record.sequence, k);
            (idx, hash)
        })
        .collect();

    // Sort by k-mer hash to group similar sequences
    indexed_records.sort_by(|a, b| {
        a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)) // Secondary sort by original index for stability
    });

    // Reorder records based on sorted indices
    let original_records = records.clone();
    for (new_pos, (old_idx, _)) in indexed_records.iter().enumerate() {
        records[new_pos] = original_records[*old_idx].clone();
    }

    info!("Reordering complete");
    Ok(records)
}

/// Advanced reordering using greedy nearest-neighbor approach
pub fn reorder_greedy(records: Vec<FastqRecord>) -> Result<Vec<FastqRecord>> {
    if records.len() <= 1 {
        return Ok(records);
    }

    info!("Reordering {} records using greedy algorithm...", records.len());

    let n = records.len();
    let mut ordered = Vec::with_capacity(n);
    let mut used = vec![false; n];

    // Start with first record
    ordered.push(records[0].clone());
    used[0] = true;

    // Greedy: always pick the most similar unused record to the last added
    for _ in 1..n {
        let last_seq = &ordered.last().unwrap().sequence;
        let last_len = last_seq.len();

        // Find the most similar unused record
        let (best_idx, _best_dist) = records
            .iter()
            .enumerate()
            .filter(|(idx, _)| !used[*idx])
            .map(|(idx, record)| {
                // Use prefix similarity for different length reads
                let min_len = last_len.min(record.sequence.len());
                let dist = hamming_distance(&last_seq[..min_len], &record.sequence[..min_len]);
                (idx, dist)
            })
            .min_by_key(|(_, dist)| *dist)
            .unwrap_or((0, usize::MAX));

        if best_idx < n {
            ordered.push(records[best_idx].clone());
            used[best_idx] = true;
        }
    }

    info!("Greedy reordering complete");
    Ok(ordered)
}

/// Reorder by GC content (simpler, faster alternative)
pub fn reorder_by_gc_content(mut records: Vec<FastqRecord>) -> Result<Vec<FastqRecord>> {
    info!("Reordering {} records by GC content...", records.len());

    records.par_sort_by(|a, b| {
        let gc_a = compute_gc_content(&a.sequence);
        let gc_b = compute_gc_content(&b.sequence);

        gc_a.partial_cmp(&gc_b)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.sequence.len().cmp(&b.sequence.len()))
    });

    info!("GC content reordering complete");
    Ok(records)
}

/// Compute GC content percentage
fn compute_gc_content(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }

    let gc_count = seq
        .chars()
        .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
        .count();

    (gc_count as f64) / (seq.len() as f64)
}

/// Choose best reordering strategy based on dataset characteristics
pub fn reorder_auto(records: Vec<FastqRecord>) -> Result<Vec<FastqRecord>> {
    let n = records.len();

    if n <= 100 {
        // For small datasets, use greedy (more accurate)
        info!("Small dataset detected, using greedy reordering");
        reorder_greedy(records)
    } else if n <= 100000 {
        // For medium datasets, use k-mer based (good balance)
        info!("Medium dataset detected, using k-mer reordering");
        reorder_by_similarity(records, 12)
    } else {
        // For large datasets, use GC content (fastest)
        info!("Large dataset detected, using GC content reordering");
        reorder_by_gc_content(records)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_distance() {
        assert_eq!(hamming_distance("ACGT", "ACGT"), 0);
        assert_eq!(hamming_distance("ACGT", "ACGA"), 1);
        assert_eq!(hamming_distance("ACGT", "TGCA"), 4);
    }

    #[test]
    fn test_gc_content() {
        assert_eq!(compute_gc_content("AAAA"), 0.0);
        assert_eq!(compute_gc_content("GGGG"), 1.0);
        assert_eq!(compute_gc_content("ACGT"), 0.5);
    }

    #[test]
    fn test_kmer_hash() {
        let hash1 = compute_kmer_hash("ACGTACGT", 4);
        let hash2 = compute_kmer_hash("ACGTTTTT", 4);
        assert_eq!(hash1, hash2); // Same first k-mer

        let hash3 = compute_kmer_hash("TGCA", 4);
        assert_ne!(hash1, hash3); // Different k-mer
    }

    #[test]
    fn test_reorder_by_gc() {
        let records = vec![
            FastqRecord::new("@r1".to_string(), "GGGG".to_string(), None), // GC=100%
            FastqRecord::new("@r2".to_string(), "AAAA".to_string(), None), // GC=0%
            FastqRecord::new("@r3".to_string(), "ACGT".to_string(), None), // GC=50%
        ];

        let reordered = reorder_by_gc_content(records).unwrap();

        // Should be sorted by GC content
        assert_eq!(reordered[0].sequence, "AAAA");
        assert_eq!(reordered[1].sequence, "ACGT");
        assert_eq!(reordered[2].sequence, "GGGG");
    }
}
