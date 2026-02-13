/// Loose-template sequence compression
///
/// Key insight: references don't need to be true genomic assemblies.
/// Even "imperfect" templates (up to ~20% mismatch) reduce encoding cost
/// vs raw bases, because encoding mismatches is cheaper than encoding
/// every base from scratch.
///
/// Approach:
/// 1. Build de Bruijn graph with SMALL k (11-15) → dense graph even at low coverage
/// 2. Extract unitigs as "loose templates" (not true contigs)
/// 3. Map ALL reads with HIGH mismatch tolerance (≤30 mismatches out of 150bp = 20%)
/// 4. Delta-encode each read against its best-matching template
/// 5. BSC-compress all streams
///
/// This works even at 0.25x coverage because small k-mers appear in many
/// reads regardless of true genomic overlap.

use super::bsc;
use super::dna_utils::*;
use anyhow::Result;
use minimizer_iter::MinimizerBuilder;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

// ── Configurable parameters ─────────────────────────────────────────────────

/// Parameters for template-based compression
pub struct TemplateParams {
    pub k: usize,
    pub mismatch_threshold: usize,
    pub tip_trim_threshold: u16,
    /// Use parallel template extraction + reduced mapping threshold for speed
    pub fast: bool,
    /// Use syncmer-based k-mer selection instead of fixed-stride for template indexing
    pub use_syncmers: bool,
}

impl Default for TemplateParams {
    fn default() -> Self {
        Self {
            k: 13,
            mismatch_threshold: 30,    // ~20% of 150bp
            tip_trim_threshold: 1,
            fast: false,
            use_syncmers: false,
        }
    }
}

// ── Data structures ─────────────────────────────────────────────────────────

struct TemplateGraph {
    kmers: FxHashMap<u64, u16>,
    k: usize,
}

struct Template {
    sequence: Vec<u8>,
}

#[derive(Clone)]
struct ReadMapping {
    template_id: usize,
    position: u32,
    is_rc: bool,
}

struct TemplateKmerHit {
    template_id: u32,
    position: u32,
    is_forward: bool,
}

// ── Graph construction ──────────────────────────────────────────────────────

fn canonical_hash(seq: &[u8], k: usize) -> Option<u64> {
    let hash = kmer_to_hash(seq)?;
    let rc = reverse_complement_hash(hash, k);
    Some(hash.min(rc))
}

impl TemplateGraph {
    fn from_sequences(sequences: &[Vec<u8>], k: usize) -> Self {
        let chunk_size = (sequences.len() / rayon::current_num_threads().max(1)).max(1000);

        let local_maps: Vec<FxHashMap<u64, u16>> = sequences
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut local: FxHashMap<u64, u16> = FxHashMap::default();
                local.reserve(chunk.len() * (150 - k + 1));
                for seq in chunk {
                    if seq.len() < k { continue; }
                    for i in 0..=seq.len() - k {
                        if let Some(canon) = canonical_hash(&seq[i..i + k], k) {
                            let count = local.entry(canon).or_insert(0);
                            *count = count.saturating_add(1);
                        }
                    }
                }
                local
            })
            .collect();

        let estimated_size = local_maps.iter().map(|m| m.len()).max().unwrap_or(0);
        let mut kmers: FxHashMap<u64, u16> = FxHashMap::default();
        kmers.reserve(estimated_size);
        for local in local_maps {
            for (hash, count) in local {
                let entry = kmers.entry(hash).or_insert(0);
                *entry = entry.saturating_add(count);
            }
        }

        eprintln!("  Template graph: {} unique canonical {}-mers from {} reads",
                 kmers.len(), k, sequences.len());

        TemplateGraph { kmers, k }
    }

    fn tip_trim(&mut self, min_count: u16) {
        let before = self.kmers.len();
        self.kmers.retain(|_, count| *count > min_count);
        eprintln!("  Tip trimming (count <= {}): {} → {} k-mers",
                 min_count, before, self.kmers.len());
    }

    fn contains(&self, canon_hash: u64) -> bool {
        self.kmers.contains_key(&canon_hash)
    }

    fn hash_to_sequence(hash: u64, k: usize) -> Vec<u8> {
        let mut seq = vec![0u8; k];
        let mut h = hash;
        for i in (0..k).rev() {
            seq[i] = idx_to_base(h as usize & 3);
            h >>= 2;
        }
        seq
    }

    fn successors(&self, kmer_seq: &[u8]) -> Vec<(u64, u8, u16)> {
        let k = self.k;
        let suffix = &kmer_seq[1..];
        let mut result = Vec::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            let mut candidate = Vec::with_capacity(k);
            candidate.extend_from_slice(suffix);
            candidate.push(base);
            if let Some(canon) = canonical_hash(&candidate, k) {
                if let Some(&count) = self.kmers.get(&canon) {
                    result.push((canon, base, count));
                }
            }
        }
        result
    }

    fn predecessors(&self, kmer_seq: &[u8]) -> Vec<(u64, u8, u16)> {
        let k = self.k;
        let prefix = &kmer_seq[..k - 1];
        let mut result = Vec::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            let mut candidate = Vec::with_capacity(k);
            candidate.push(base);
            candidate.extend_from_slice(prefix);
            if let Some(canon) = canonical_hash(&candidate, k) {
                if let Some(&count) = self.kmers.get(&canon) {
                    result.push((canon, base, count));
                }
            }
        }
        result
    }

    /// Get the highest-coverage successor (greedy walk)
    fn best_successor(&self, kmer_seq: &[u8]) -> Option<(u64, u8)> {
        self.successors(kmer_seq)
            .into_iter()
            .max_by_key(|&(_, _, count)| count)
            .map(|(hash, base, _)| (hash, base))
    }

    /// Get the highest-coverage predecessor (greedy walk)
    fn best_predecessor(&self, kmer_seq: &[u8]) -> Option<(u64, u8)> {
        self.predecessors(kmer_seq)
            .into_iter()
            .max_by_key(|&(_, _, count)| count)
            .map(|(hash, base, _)| (hash, base))
    }
}

// ── Template (greedy walk) extraction ───────────────────────────────────────

#[inline(always)]
fn base_to_2bit(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Single greedy walk using rolling fwd+rc hash: O(1) per step instead of O(k).
fn walk_from_seed(
    graph: &TemplateGraph,
    start_hash: u64,
    k: usize,
    mask: u64,
    rc_shift: usize, // = 2 * (k - 1)
    max_walk_length: usize,
    min_length: usize,
) -> Option<Template> {
    let start_seq = TemplateGraph::hash_to_sequence(start_hash, k);
    let start_seq = if kmer_to_hash(&start_seq) == Some(start_hash) {
        start_seq
    } else {
        reverse_complement(&TemplateGraph::hash_to_sequence(start_hash, k))
    };
    let start_fwd = kmer_to_hash(&start_seq).unwrap();
    let start_rc = reverse_complement_hash(start_fwd, k);

    let mut local_visited: FxHashSet<u64> = FxHashSet::default();
    local_visited.insert(start_fwd.min(start_rc));

    // Walk RIGHT: rolling fwd+rc hash
    let mut right_bases: Vec<u8> = Vec::with_capacity(1000);
    let mut cur_fwd = start_fwd;
    let mut cur_rc = start_rc;
    while right_bases.len() < max_walk_length {
        let mut best_canon = 0u64;
        let mut best_base = 0u8;
        let mut best_count = 0u16;
        let mut best_fwd = 0u64;
        let mut best_rc = 0u64;

        for bval in 0u64..4 {
            let cand_fwd = ((cur_fwd << 2) | bval) & mask;
            let cand_rc = (cur_rc >> 2) | ((3 - bval) << rc_shift);
            let canon = cand_fwd.min(cand_rc);

            if let Some(&count) = graph.kmers.get(&canon) {
                if count > best_count {
                    best_count = count;
                    best_canon = canon;
                    best_base = idx_to_base(bval as usize);
                    best_fwd = cand_fwd;
                    best_rc = cand_rc;
                }
            }
        }

        if best_count == 0 || local_visited.contains(&best_canon) { break; }
        local_visited.insert(best_canon);
        right_bases.push(best_base);
        cur_fwd = best_fwd;
        cur_rc = best_rc;
    }

    // Walk LEFT: rolling fwd+rc hash
    let mut left_bases: Vec<u8> = Vec::with_capacity(1000);
    cur_fwd = start_fwd;
    cur_rc = start_rc;
    while left_bases.len() < max_walk_length {
        let mut best_canon = 0u64;
        let mut best_base = 0u8;
        let mut best_count = 0u16;
        let mut best_fwd = 0u64;
        let mut best_rc = 0u64;

        for bval in 0u64..4 {
            let cand_fwd = (bval << rc_shift) | (cur_fwd >> 2);
            let cand_rc = ((cur_rc << 2) | (3 - bval)) & mask;
            let canon = cand_fwd.min(cand_rc);

            if let Some(&count) = graph.kmers.get(&canon) {
                if count > best_count {
                    best_count = count;
                    best_canon = canon;
                    best_base = idx_to_base(bval as usize);
                    best_fwd = cand_fwd;
                    best_rc = cand_rc;
                }
            }
        }

        if best_count == 0 || local_visited.contains(&best_canon) { break; }
        local_visited.insert(best_canon);
        left_bases.push(best_base);
        cur_fwd = best_fwd;
        cur_rc = best_rc;
    }

    left_bases.reverse();
    let total_len = left_bases.len() + k + right_bases.len();

    if total_len >= min_length {
        let mut seq = Vec::with_capacity(total_len);
        seq.extend_from_slice(&left_bases);
        seq.extend_from_slice(&start_seq);
        seq.extend_from_slice(&right_bases);
        Some(Template { sequence: seq })
    } else {
        None
    }
}

/// Extract templates via sequential greedy walks with diversity tracking.
/// Uses rolling hash walks for speed + used_seeds set to ensure templates
/// cover distinct regions of k-mer space.
fn extract_templates(graph: &TemplateGraph, min_length: usize) -> Vec<Template> {
    let k = graph.k;
    let max_walk_length = 5_000;
    let max_templates = 50_000;
    let min_seed_count = 3u16;
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    // Pre-filter: only collect seeds with count ≥ 3 before sorting (was 52M → 985K)
    let mut seeds: Vec<(u64, u16)> = graph.kmers.iter()
        .filter(|&(_, c)| *c >= min_seed_count)
        .map(|(&h, &c)| (h, c))
        .collect();
    seeds.sort_unstable_by(|a, b| b.1.cmp(&a.1));

    eprintln!("  {} seeds with count ≥ {} (out of {} total k-mers)",
             seeds.len(), min_seed_count, graph.kmers.len());

    let rc_shift = 2 * (k - 1);
    let mut used_seeds: FxHashSet<u64> = FxHashSet::default();
    let mut templates: Vec<Template> = Vec::new();

    for &(start_hash, start_count) in &seeds {
        if templates.len() >= max_templates { break; }
        if start_count < min_seed_count { break; }
        if used_seeds.contains(&start_hash) { continue; }

        if let Some(template) = walk_from_seed(graph, start_hash, k, mask, rc_shift, max_walk_length, min_length) {
            // Mark template k-mers as used seeds using rolling fwd+rc hash
            let seq = &template.sequence;
            if seq.len() >= k {
                let mut fwd = 0u64;
                for i in 0..k {
                    fwd = (fwd << 2) | base_to_2bit(seq[i]);
                }
                fwd &= mask;
                let mut rc = reverse_complement_hash(fwd, k); // O(k) once
                used_seeds.insert(fwd.min(rc));

                for i in k..seq.len() {
                    let bval = base_to_2bit(seq[i]);
                    fwd = ((fwd << 2) | bval) & mask;
                    rc = (rc >> 2) | ((3 - bval) << rc_shift);
                    used_seeds.insert(fwd.min(rc));
                }
            }
            templates.push(template);
        }
    }

    templates.sort_unstable_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

    let total_bases: usize = templates.iter().map(|t| t.sequence.len()).sum();
    let max_len = templates.first().map(|t| t.sequence.len()).unwrap_or(0);
    let long_templates = templates.iter().filter(|t| t.sequence.len() >= 150).count();
    eprintln!("  Extracted {} templates ({} total bases, longest={}, {} ≥ 150bp)",
             templates.len(), total_bases, max_len, long_templates);

    templates
}

/// Extend template tips through count=1 k-mers in the full graph.
/// After extracting from count≥2 filtered graph, templates are truncated at
/// singleton boundaries. This function extends them using the original graph.
fn extend_tips(templates: &mut [Template], graph: &TemplateGraph) -> (usize, usize) {
    let k = graph.k;
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let rc_shift = 2 * (k - 1);
    let max_extension = 500;

    let mut extended_count = 0usize;
    let mut total_bases_added = 0usize;

    for template in templates.iter_mut() {
        let seq = &template.sequence;
        if seq.len() < k { continue; }

        let orig_len = seq.len();

        // Build visited set from existing template k-mers (avoid re-entering)
        let mut visited: FxHashSet<u64> = FxHashSet::default();
        let mut fwd = 0u64;
        for i in 0..k {
            fwd = (fwd << 2) | base_to_2bit(seq[i]);
        }
        fwd &= mask;
        let mut rc = reverse_complement_hash(fwd, k);
        visited.insert(fwd.min(rc));
        for i in k..seq.len() {
            let bval = base_to_2bit(seq[i]);
            fwd = ((fwd << 2) | bval) & mask;
            rc = (rc >> 2) | ((3 - bval) << rc_shift);
            visited.insert(fwd.min(rc));
        }

        // Extend RIGHT from last k-mer
        // fwd/rc are already set to the last k-mer from the loop above
        let mut right_ext: Vec<u8> = Vec::new();
        while right_ext.len() < max_extension {
            let mut best_count = 0u16;
            let mut best_base = 0u8;
            let mut best_fwd = 0u64;
            let mut best_rc = 0u64;
            let mut best_canon = 0u64;

            for bval in 0u64..4 {
                let cand_fwd = ((fwd << 2) | bval) & mask;
                let cand_rc = (rc >> 2) | ((3 - bval) << rc_shift);
                let canon = cand_fwd.min(cand_rc);
                if let Some(&count) = graph.kmers.get(&canon) {
                    if count > best_count {
                        best_count = count;
                        best_base = idx_to_base(bval as usize);
                        best_fwd = cand_fwd;
                        best_rc = cand_rc;
                        best_canon = canon;
                    }
                }
            }

            if best_count == 0 || visited.contains(&best_canon) { break; }
            visited.insert(best_canon);
            right_ext.push(best_base);
            fwd = best_fwd;
            rc = best_rc;
        }

        // Extend LEFT from first k-mer
        fwd = 0u64;
        for i in 0..k {
            fwd = (fwd << 2) | base_to_2bit(seq[i]);
        }
        fwd &= mask;
        rc = reverse_complement_hash(fwd, k);

        let mut left_ext: Vec<u8> = Vec::new();
        while left_ext.len() < max_extension {
            let mut best_count = 0u16;
            let mut best_base = 0u8;
            let mut best_fwd = 0u64;
            let mut best_rc = 0u64;
            let mut best_canon = 0u64;

            for bval in 0u64..4 {
                let cand_fwd = (bval << rc_shift) | (fwd >> 2);
                let cand_rc = ((rc << 2) | (3 - bval)) & mask;
                let canon = cand_fwd.min(cand_rc);
                if let Some(&count) = graph.kmers.get(&canon) {
                    if count > best_count {
                        best_count = count;
                        best_base = idx_to_base(bval as usize);
                        best_fwd = cand_fwd;
                        best_rc = cand_rc;
                        best_canon = canon;
                    }
                }
            }

            if best_count == 0 || visited.contains(&best_canon) { break; }
            visited.insert(best_canon);
            left_ext.push(best_base);
            fwd = best_fwd;
            rc = best_rc;
        }
        left_ext.reverse();

        if !left_ext.is_empty() || !right_ext.is_empty() {
            let mut new_seq = Vec::with_capacity(left_ext.len() + seq.len() + right_ext.len());
            new_seq.extend_from_slice(&left_ext);
            new_seq.extend_from_slice(seq);
            new_seq.extend_from_slice(&right_ext);
            total_bases_added += new_seq.len() - orig_len;
            template.sequence = new_seq;
            extended_count += 1;
        }
    }

    (extended_count, total_bases_added)
}

/// Error-correct templates using k-mer voting from the graph.
/// Finds positions where a k-mer has count=1 but a neighbor has count≥3,
/// tries all 3 alternative bases and picks the one with highest k-mer count.
fn error_correct_templates(templates: &mut [Template], graph: &TemplateGraph) -> usize {
    let k = graph.k;
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let rc_shift = 2 * (k - 1);
    let mut total_corrections = 0usize;

    for template in templates.iter_mut() {
        if template.sequence.len() < k + 2 { continue; }

        // Phase 1: compute rolling hashes and counts (immutable borrow of sequence)
        let n_kmers = template.sequence.len() - k + 1;
        let mut hashes_fwd: Vec<u64> = Vec::with_capacity(n_kmers);
        let mut hashes_rc: Vec<u64> = Vec::with_capacity(n_kmers);
        let mut counts: Vec<u16> = Vec::with_capacity(n_kmers);

        let mut fwd = 0u64;
        for i in 0..k {
            fwd = (fwd << 2) | base_to_2bit(template.sequence[i]);
        }
        fwd &= mask;
        let mut rc = reverse_complement_hash(fwd, k);
        let canon = fwd.min(rc);
        hashes_fwd.push(fwd);
        hashes_rc.push(rc);
        counts.push(*graph.kmers.get(&canon).unwrap_or(&0));

        for i in k..template.sequence.len() {
            let bval = base_to_2bit(template.sequence[i]);
            fwd = ((fwd << 2) | bval) & mask;
            rc = (rc >> 2) | ((3 - bval) << rc_shift);
            let canon = fwd.min(rc);
            hashes_fwd.push(fwd);
            hashes_rc.push(rc);
            counts.push(*graph.kmers.get(&canon).unwrap_or(&0));
        }

        // Phase 2: collect corrections (no borrow of sequence)
        let mut corrections: Vec<(usize, u8)> = Vec::new();
        let mut corrected_pos: FxHashSet<usize> = FxHashSet::default();

        for pos in 0..n_kmers {
            if counts[pos] >= 2 { continue; }

            // Check if left neighbor has high count → error at pos+k-1 (trailing edge)
            if pos > 0 && counts[pos - 1] >= 3 {
                let error_idx = pos + k - 1;
                if error_idx >= template.sequence.len() || corrected_pos.contains(&error_idx) { continue; }

                let cur_fwd = hashes_fwd[pos];
                let cur_rc = hashes_rc[pos];
                let orig_bval = base_to_2bit(template.sequence[error_idx]);
                let mut best_count = counts[pos];
                let mut best_bval = orig_bval;

                for bval in 0u64..4 {
                    if bval == orig_bval { continue; }
                    let new_fwd = (cur_fwd & !3u64) | bval;
                    let new_rc = (cur_rc & !(3u64 << rc_shift)) | ((3 - bval) << rc_shift);
                    let new_canon = new_fwd.min(new_rc);
                    if let Some(&count) = graph.kmers.get(&new_canon) {
                        if count > best_count {
                            best_count = count;
                            best_bval = bval;
                        }
                    }
                }

                if best_bval != orig_bval && best_count >= 3 {
                    corrections.push((error_idx, idx_to_base(best_bval as usize)));
                    corrected_pos.insert(error_idx);
                }
            }

            // Check if right neighbor has high count → error at pos (leading edge)
            if pos + 1 < n_kmers && counts[pos + 1] >= 3 {
                let error_idx = pos;
                if corrected_pos.contains(&error_idx) { continue; }

                let cur_fwd = hashes_fwd[pos];
                let cur_rc = hashes_rc[pos];
                let orig_bval = base_to_2bit(template.sequence[error_idx]);
                let mut best_count = counts[pos];
                let mut best_bval = orig_bval;

                for bval in 0u64..4 {
                    if bval == orig_bval { continue; }
                    let new_fwd = (cur_fwd & !(3u64 << rc_shift)) | (bval << rc_shift);
                    let new_rc = (cur_rc & !3u64) | (3 - bval);
                    let new_canon = new_fwd.min(new_rc);
                    if let Some(&count) = graph.kmers.get(&new_canon) {
                        if count > best_count {
                            best_count = count;
                            best_bval = bval;
                        }
                    }
                }

                if best_bval != orig_bval && best_count >= 3 {
                    corrections.push((error_idx, idx_to_base(best_bval as usize)));
                    corrected_pos.insert(error_idx);
                }
            }
        }

        // Phase 3: apply corrections (mutable borrow of sequence)
        for (idx, base) in &corrections {
            template.sequence[*idx] = *base;
        }
        total_corrections += corrections.len();
    }

    total_corrections
}

/// Fast template extraction: filter graph to count≥2 (fits L3 cache), then
/// walk sequentially with used_seeds tracking for diversity.
/// Much faster than full extraction because 52M→3M k-mers fits in cache.
/// After extraction, extends tips through count=1 k-mers and error-corrects.
fn extract_templates_fast(graph: &mut TemplateGraph, min_length: usize) -> Vec<Template> {
    let before = graph.kmers.len();

    // Save full graph, replace with filtered for cache-friendly extraction
    let full_kmers = std::mem::take(&mut graph.kmers);
    graph.kmers = full_kmers.iter()
        .filter(|&(_, c)| *c >= 2)
        .map(|(&h, &c)| (h, c))
        .collect();

    let count2 = graph.kmers.len();
    eprintln!("  Fast filter: {} → {} k-mers (count ≥ 2)", before, count2);

    let mut templates = extract_templates(graph, min_length);

    // Only do tip extension + error correction if full graph fits in L3 cache (~30MB)
    // Each k-mer entry: ~20 bytes in FxHashMap → 50M entries ≈ 1GB
    // At >60M k-mers, the graph is too large for cache-friendly walking
    let max_graph_for_tips = 60_000_000;

    if before <= max_graph_for_tips {
        // Restore full graph for tip extension + error correction
        graph.kmers = full_kmers;
        let t = std::time::Instant::now();

        let (ext_count, ext_bases) = extend_tips(&mut templates, graph);
        eprintln!("  Tip extension: {} templates extended, {} bases added ({:.1}s)",
                 ext_count, ext_bases, t.elapsed().as_secs_f64());

        let corrections = error_correct_templates(&mut templates, graph);
        eprintln!("  Error correction: {} bases corrected", corrections);
    } else {
        eprintln!("  Skipping tip extension (graph too large: {}M k-mers > {}M threshold)",
                 before / 1_000_000, max_graph_for_tips / 1_000_000);
        // Error-correct using only count≥2 graph (still in memory, much smaller)
        let corrections = error_correct_templates(&mut templates, graph);
        eprintln!("  Error correction (count≥2 only): {} bases corrected", corrections);
        // Restore full kmers for potential future use (callers may need it)
        graph.kmers = full_kmers;
    }

    templates
}

// ── Read mapping with high mismatch tolerance ───────────────────────────────

fn build_template_index(
    templates: &[Template],
    k: usize,
) -> FxHashMap<u64, Vec<TemplateKmerHit>> {
    let total_bases: usize = templates.iter().map(|t| t.sequence.len()).sum();
    // Sample every `step` positions — for small k, index is very large
    let step = if total_bases > 100_000_000 { 8 }
               else if total_bases > 10_000_000 { 4 }
               else { 2 };

    let mut index: FxHashMap<u64, Vec<TemplateKmerHit>> = FxHashMap::default();
    index.reserve(total_bases / step);

    for (tid, template) in templates.iter().enumerate() {
        if template.sequence.len() < k { continue; }
        // Only index templates ≥ read length (short ones aren't useful)
        // Actually, index all — even shorter templates can match partial reads
        let mut pos = 0;
        while pos + k <= template.sequence.len() {
            let kmer = &template.sequence[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let is_forward = fwd <= rc;
                let bucket = index.entry(canon).or_default();
                // Cap entries per k-mer to avoid slow lookups on repetitive k-mers
                if bucket.len() < 100 {
                    bucket.push(TemplateKmerHit {
                        template_id: tid as u32,
                        position: pos as u32,
                        is_forward,
                    });
                }
            }
            pos += step;
        }
    }

    eprintln!("  Template index: {} entries (step={})",
             index.values().map(|v| v.len()).sum::<usize>(), step);

    index
}

fn map_reads_to_templates(
    sequences: &[Vec<u8>],
    templates: &[Template],
    index: &FxHashMap<u64, Vec<TemplateKmerHit>>,
    k: usize,
    mismatch_threshold: usize,
) -> Vec<Option<ReadMapping>> {
    sequences
        .par_iter()
        .map(|read| {
            if read.len() < k { return None; }

            // Use more anchor positions for small k (more hits per k-mer)
            let num_anchors = 8;
            let step = read.len().saturating_sub(k) / (num_anchors - 1).max(1);

            let mut best: Option<ReadMapping> = None;
            let mut best_dist = usize::MAX;

            for anchor_idx in 0..num_anchors {
                let anchor_pos = (anchor_idx * step).min(read.len().saturating_sub(k));
                if anchor_pos + k > read.len() { continue; }
                let kmer = &read[anchor_pos..anchor_pos + k];

                let fwd = match kmer_to_hash(kmer) {
                    Some(h) => h,
                    None => continue,
                };
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let read_is_fwd = fwd <= rc;

                let hits = match index.get(&canon) {
                    Some(h) => h,
                    None => continue,
                };

                // Limit per-anchor work
                let max_hits = 50;
                for (hit_idx, hit) in hits.iter().enumerate() {
                    if hit_idx >= max_hits { break; }

                    let template = &templates[hit.template_id as usize];
                    let maps_fwd = read_is_fwd == hit.is_forward;

                    let read_start: i64 = if maps_fwd {
                        hit.position as i64 - anchor_pos as i64
                    } else {
                        hit.position as i64 - (read.len() as i64 - anchor_pos as i64 - k as i64)
                    };

                    if read_start < 0 { continue; }
                    let read_start = read_start as usize;
                    if read_start + read.len() > template.sequence.len() { continue; }

                    let region = &template.sequence[read_start..read_start + read.len()];

                    let dist = if maps_fwd {
                        hamming_distance_within(region, read, mismatch_threshold)
                    } else {
                        let rc_region = reverse_complement(region);
                        hamming_distance_within(&rc_region, read, mismatch_threshold)
                    };

                    if let Some(d) = dist {
                        if d < best_dist {
                            best_dist = d;
                            best = Some(ReadMapping {
                                template_id: hit.template_id as usize,
                                position: read_start as u32,
                                is_rc: !maps_fwd,
                            });
                            if d == 0 { return best; }
                        }
                    }
                }

                if best_dist == 0 { break; }
            }

            best
        })
        .collect()
}

// ── Syncmer-based read mapping ──────────────────────────────────────────────

fn build_template_index_syncmer(
    templates: &[Template],
    k: usize,
) -> FxHashMap<u64, Vec<TemplateKmerHit>> {
    // Closed syncmers with s close to k for high density
    // density ≈ 2/(k-s+1) for closed syncmers
    // s = k-3 → density ≈ 2/4 = 50%, matching stride step=2
    let s = k.saturating_sub(3).max(5);
    let t_end = k - s; // target positions: [0, k-s] for closed syncmers

    let mut index: FxHashMap<u64, Vec<TemplateKmerHit>> = FxHashMap::default();

    let mut total_indexed = 0usize;
    for (tid, template) in templates.iter().enumerate() {
        if template.sequence.len() < k { continue; }

        // Closed syncmers: smallest s-mer at position 0 or k-s
        let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], &template.sequence);

        for pos in positions {
            if pos + k > template.sequence.len() { continue; }
            let kmer = &template.sequence[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let is_forward = fwd <= rc;
                let bucket = index.entry(canon).or_default();
                if bucket.len() < 100 {
                    bucket.push(TemplateKmerHit {
                        template_id: tid as u32,
                        position: pos as u32,
                        is_forward,
                    });
                    total_indexed += 1;
                }
            }
        }
    }

    eprintln!("  Template index (closed syncmer): {} entries, {} unique k-mers (s={}, targets=[0,{}])",
             total_indexed, index.len(), s, t_end);

    index
}

fn map_reads_to_templates_syncmer(
    sequences: &[Vec<u8>],
    templates: &[Template],
    index: &FxHashMap<u64, Vec<TemplateKmerHit>>,
    k: usize,
    mismatch_threshold: usize,
) -> Vec<Option<ReadMapping>> {
    let s = k.saturating_sub(3).max(5);
    let t_end = k - s;

    sequences
        .par_iter()
        .map(|read| {
            if read.len() < k { return None; }

            // Get syncmer positions from the read (content-addressed selection)
            let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], read);

            if positions.is_empty() { return None; }

            let mut best: Option<ReadMapping> = None;
            let mut best_dist = usize::MAX;

            // Evenly space anchors across all syncmer positions
            let max_anchors = 16;
            let n_pos = positions.len();
            let step_pos = if n_pos > max_anchors { n_pos / max_anchors } else { 1 };
            for idx in (0..n_pos).step_by(step_pos).take(max_anchors) {
                let anchor_pos = positions[idx];
                if anchor_pos + k > read.len() { continue; }
                let kmer = &read[anchor_pos..anchor_pos + k];

                let fwd = match kmer_to_hash(kmer) {
                    Some(h) => h,
                    None => continue,
                };
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let read_is_fwd = fwd <= rc;

                let hits = match index.get(&canon) {
                    Some(h) => h,
                    None => continue,
                };

                let max_hits = 50;
                for (hit_idx, hit) in hits.iter().enumerate() {
                    if hit_idx >= max_hits { break; }

                    let template = &templates[hit.template_id as usize];
                    let maps_fwd = read_is_fwd == hit.is_forward;

                    let read_start: i64 = if maps_fwd {
                        hit.position as i64 - anchor_pos as i64
                    } else {
                        hit.position as i64 - (read.len() as i64 - anchor_pos as i64 - k as i64)
                    };

                    if read_start < 0 { continue; }
                    let read_start = read_start as usize;
                    if read_start + read.len() > template.sequence.len() { continue; }

                    let region = &template.sequence[read_start..read_start + read.len()];

                    let dist = if maps_fwd {
                        hamming_distance_within(region, read, mismatch_threshold)
                    } else {
                        let rc_region = reverse_complement(region);
                        hamming_distance_within(&rc_region, read, mismatch_threshold)
                    };

                    if let Some(d) = dist {
                        if d < best_dist {
                            best_dist = d;
                            best = Some(ReadMapping {
                                template_id: hit.template_id as usize,
                                position: read_start as u32,
                                is_rc: !maps_fwd,
                            });
                            if d == 0 { return best; }
                        }
                    }
                }

                if best_dist == 0 { break; }
            }

            best
        })
        .collect()
}

// ── Encoding ────────────────────────────────────────────────────────────────

fn encode_reads(
    sequences: &[Vec<u8>],
    templates: &[Template],
    mappings: &[Option<ReadMapping>],
) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    // Determine which templates are actually used and remap IDs
    let mut used_ids: FxHashSet<usize> = FxHashSet::default();
    for m in mappings.iter().flatten() {
        used_ids.insert(m.template_id);
    }

    let mut id_map: FxHashMap<usize, usize> = FxHashMap::default();
    let mut used_templates: Vec<&Template> = Vec::new();
    // Sort to ensure deterministic ordering
    let mut sorted_ids: Vec<usize> = used_ids.into_iter().collect();
    sorted_ids.sort();
    for &uid in &sorted_ids {
        id_map.insert(uid, used_templates.len());
        used_templates.push(&templates[uid]);
    }

    // Stream 0: template consensus sequences (2-bit packed)
    let consensus_seqs: Vec<Vec<u8>> = used_templates.iter().map(|t| t.sequence.clone()).collect();
    let consensus_bytes = pack_dna_2bit(&consensus_seqs);

    // Stream 1: metadata (mapped/singleton flag + contig_id + position + orientation)
    let mut meta_bytes: Vec<u8> = Vec::new();
    // Stream 2: mismatches
    let mut mismatch_bytes: Vec<u8> = Vec::new();
    // Stream 3: singletons
    let mut singleton_seqs: Vec<Vec<u8>> = Vec::new();

    write_varint(&mut meta_bytes, sequences.len() as u64);

    let mut total_mismatches = 0u64;

    for (read_idx, mapping) in mappings.iter().enumerate() {
        match mapping {
            Some(m) => {
                let mapped_id = id_map[&m.template_id];
                meta_bytes.push(1); // mapped
                write_varint(&mut meta_bytes, mapped_id as u64);
                write_varint(&mut meta_bytes, m.position as u64);
                meta_bytes.push(if m.is_rc { 1 } else { 0 });
                write_varint(&mut meta_bytes, sequences[read_idx].len() as u64);

                let read = &sequences[read_idx];
                let template = &templates[m.template_id];
                let region = &template.sequence[m.position as usize..m.position as usize + read.len()];

                let (compare_read, compare_ref) = if m.is_rc {
                    (reverse_complement(read), region.to_vec())
                } else {
                    (read.clone(), region.to_vec())
                };

                let mut mismatches: Vec<(u16, u8)> = Vec::new();
                for (i, (&r, &c)) in compare_read.iter().zip(compare_ref.iter()).enumerate() {
                    if r != c {
                        mismatches.push((i as u16, r));
                    }
                }

                total_mismatches += mismatches.len() as u64;

                write_varint(&mut mismatch_bytes, mismatches.len() as u64);
                let mut prev_pos: u16 = 0;
                for (pos, base) in &mismatches {
                    write_varint(&mut mismatch_bytes, (*pos - prev_pos) as u64);
                    mismatch_bytes.push(*base);
                    prev_pos = *pos;
                }
            }
            None => {
                meta_bytes.push(0); // singleton
                singleton_seqs.push(sequences[read_idx].clone());
            }
        }
    }

    let singleton_bytes = pack_dna_2bit(&singleton_seqs);

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    let avg_mm = if mapped > 0 { total_mismatches as f64 / mapped as f64 } else { 0.0 };
    eprintln!("  Encoded: {} mapped (avg {:.1} mismatches), {} singletons",
             mapped, avg_mm, singleton_seqs.len());

    (consensus_bytes, meta_bytes, mismatch_bytes, singleton_bytes)
}

// ── Public API ──────────────────────────────────────────────────────────────

/// Compress sequences using loose-template approach.
/// Returns compressed bytes in DBG1 format (compatible with debruijn.rs decompressor).
pub fn compress_sequences_template(
    sequences: &[String],
    params: &TemplateParams,
) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    if seqs.is_empty() { return Ok(Vec::new()); }

    let k = params.k;
    let mismatch_threshold = params.mismatch_threshold;

    eprintln!("Template compression: {} reads, k={}, max_mismatches={}",
             seqs.len(), k, mismatch_threshold);

    // Step 1: Build de Bruijn graph
    eprintln!("  Step 1: Building graph (k={})...", k);
    let t = std::time::Instant::now();
    let mut graph = TemplateGraph::from_sequences(&seqs, k);

    // Adaptive tip trimming
    let total_kmers = graph.kmers.len();
    let singletons = graph.kmers.values().filter(|&&c| c == 1).count();
    let singleton_frac = singletons as f64 / total_kmers.max(1) as f64;
    eprintln!("  Singleton fraction: {:.1}%", singleton_frac * 100.0);

    if params.tip_trim_threshold > 0 && singleton_frac < 0.95 {
        graph.tip_trim(params.tip_trim_threshold);
    } else {
        eprintln!("  Skipping tip trimming (shallow coverage or disabled)");
    }
    eprintln!("  Graph built: {:.1}s", t.elapsed().as_secs_f64());

    // Step 2: Extract templates via greedy walks
    let min_length = seqs.iter().map(|s| s.len()).max().unwrap_or(150);
    eprintln!("  Step 2: Extracting templates (min_length={})...", min_length);
    let t = std::time::Instant::now();
    let templates = extract_templates(&graph, min_length);
    // Drop graph to free memory
    drop(graph);
    eprintln!("  Templates extracted: {:.1}s", t.elapsed().as_secs_f64());

    // Step 3: Build index and map reads
    eprintln!("  Step 3: Mapping reads (mismatch_threshold={})...", mismatch_threshold);
    let t = std::time::Instant::now();
    let index = build_template_index(&templates, k);
    let mappings = map_reads_to_templates(&seqs, &templates, &index, k, mismatch_threshold);
    // Drop index to free memory
    drop(index);

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    eprintln!("  Mapping: {}/{} reads mapped ({:.1}%) ({:.1}s)",
             mapped, seqs.len(), mapped as f64 / seqs.len() as f64 * 100.0,
             t.elapsed().as_secs_f64());

    // Step 4: Encode
    eprintln!("  Step 4: Encoding...");
    let (consensus_raw, meta_raw, mismatch_raw, singleton_raw) =
        encode_reads(&seqs, &templates, &mappings);

    // Step 5: BSC-compress
    eprintln!("  Step 5: BSC-compressing...");
    let bsc_consensus = bsc::compress_parallel(&consensus_raw)?;
    let bsc_meta = bsc::compress_parallel(&meta_raw)?;
    let bsc_mismatches = bsc::compress_parallel(&mismatch_raw)?;
    let bsc_singletons = bsc::compress_parallel(&singleton_raw)?;

    eprintln!("  Raw: consensus={}, meta={}, mismatches={}, singletons={}",
             consensus_raw.len(), meta_raw.len(), mismatch_raw.len(), singleton_raw.len());
    eprintln!("  BSC: consensus={}, meta={}, mismatches={}, singletons={}",
             bsc_consensus.len(), bsc_meta.len(), bsc_mismatches.len(), bsc_singletons.len());

    // Serialize (DBG1 format — compatible with debruijn.rs decompressor)
    let mut output = Vec::new();
    output.extend_from_slice(b"DBG1");
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_consensus.len() as u64);
    write_u64(&mut output, bsc_meta.len() as u64);
    write_u64(&mut output, bsc_mismatches.len() as u64);
    write_u64(&mut output, bsc_singletons.len() as u64);
    output.extend_from_slice(&bsc_consensus);
    output.extend_from_slice(&bsc_meta);
    output.extend_from_slice(&bsc_mismatches);
    output.extend_from_slice(&bsc_singletons);

    let total = output.len();
    let original: usize = sequences.iter().map(|s| s.len()).sum();
    eprintln!("  Template: {} → {} bytes ({:.2}x)",
             original, total, original as f64 / total as f64);

    Ok(output)
}

/// Compress sequences using template substitution approach.
/// Instead of splitting into 4 streams, replaces matching bases with \0
/// and feeds the ENTIRE sequence block to BSC as one stream.
/// This preserves BWT context across all reads.
pub fn compress_sequences_template_subst(
    sequences: &[String],
    params: &TemplateParams,
) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    if seqs.is_empty() { return Ok(Vec::new()); }

    let k = params.k;
    let mismatch_threshold = params.mismatch_threshold;

    eprintln!("Template-subst: {} reads, k={}, max_mismatches={}",
             seqs.len(), k, mismatch_threshold);

    // Steps 1-3: identical to compress_sequences_template
    eprintln!("  Step 1: Building graph (k={})...", k);
    let t = std::time::Instant::now();
    let mut graph = TemplateGraph::from_sequences(&seqs, k);

    let total_kmers = graph.kmers.len();
    let singletons = graph.kmers.values().filter(|&&c| c == 1).count();
    let singleton_frac = singletons as f64 / total_kmers.max(1) as f64;

    if params.tip_trim_threshold > 0 && singleton_frac < 0.95 {
        graph.tip_trim(params.tip_trim_threshold);
    }
    eprintln!("  Graph built: {:.1}s", t.elapsed().as_secs_f64());

    let min_length = seqs.iter().map(|s| s.len()).max().unwrap_or(150);
    eprintln!("  Step 2: Extracting templates (min_length={})...", min_length);
    let t = std::time::Instant::now();
    let templates = extract_templates(&graph, min_length);
    drop(graph);
    eprintln!("  Templates extracted: {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("  Step 3: Mapping reads (mismatch_threshold={})...", mismatch_threshold);
    let t = std::time::Instant::now();
    let index = build_template_index(&templates, k);
    let mappings = map_reads_to_templates(&seqs, &templates, &index, k, mismatch_threshold);
    drop(index);

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    eprintln!("  Mapping: {}/{} mapped ({:.1}%) ({:.1}s)",
             mapped, seqs.len(), mapped as f64 / seqs.len() as f64 * 100.0,
             t.elapsed().as_secs_f64());

    // Step 4: Build substituted stream
    // For mapped reads: replace matching bases with \0, keep mismatches as original
    // For singletons: keep original bases unchanged
    eprintln!("  Step 4: Building substituted stream...");
    let total_bases: usize = seqs.iter().map(|s| s.len()).sum();
    let mut subst_stream = Vec::with_capacity(total_bases);
    let mut meta_bytes: Vec<u8> = Vec::new();
    let mut total_subst = 0u64;

    // Determine which templates are used and remap IDs
    let mut used_ids: FxHashSet<usize> = FxHashSet::default();
    for m in mappings.iter().flatten() {
        used_ids.insert(m.template_id);
    }
    let mut id_map: FxHashMap<usize, usize> = FxHashMap::default();
    let mut used_templates: Vec<&Template> = Vec::new();
    let mut sorted_ids: Vec<usize> = used_ids.into_iter().collect();
    sorted_ids.sort();
    for &uid in &sorted_ids {
        id_map.insert(uid, used_templates.len());
        used_templates.push(&templates[uid]);
    }

    write_varint(&mut meta_bytes, seqs.len() as u64);

    for (read_idx, mapping) in mappings.iter().enumerate() {
        let read = &seqs[read_idx];
        match mapping {
            Some(m) => {
                let mapped_id = id_map[&m.template_id];
                meta_bytes.push(1); // mapped
                write_varint(&mut meta_bytes, mapped_id as u64);
                write_varint(&mut meta_bytes, m.position as u64);
                meta_bytes.push(if m.is_rc { 1 } else { 0 });

                let template = &templates[m.template_id];
                let region = &template.sequence[m.position as usize..m.position as usize + read.len()];

                let compare_ref = if m.is_rc {
                    reverse_complement(region)
                } else {
                    region.to_vec()
                };

                for (i, &read_base) in read.iter().enumerate() {
                    if read_base == compare_ref[i] {
                        subst_stream.push(0u8); // match → \0
                        total_subst += 1;
                    } else {
                        subst_stream.push(read_base); // mismatch → keep original
                    }
                }
            }
            None => {
                meta_bytes.push(0); // singleton
                subst_stream.extend_from_slice(read);
            }
        }
    }

    let avg_subst = if mapped > 0 { total_subst as f64 / mapped as f64 } else { 0.0 };
    eprintln!("  Substituted: {} bases → \\0 (avg {:.1}/read across {} mapped reads)",
             total_subst, avg_subst, mapped);

    // Template consensus for decompression
    let consensus_seqs: Vec<Vec<u8>> = used_templates.iter().map(|t| t.sequence.clone()).collect();
    let consensus_bytes = pack_dna_2bit(&consensus_seqs);

    // Step 5: BSC-compress all streams
    eprintln!("  Step 5: BSC-compressing...");
    let bsc_subst = bsc::compress_parallel_adaptive(&subst_stream)?;
    let bsc_meta = bsc::compress_parallel(&meta_bytes)?;
    let bsc_consensus = bsc::compress_parallel(&consensus_bytes)?;

    eprintln!("  Raw: subst_stream={}, meta={}, consensus={}",
             subst_stream.len(), meta_bytes.len(), consensus_bytes.len());
    eprintln!("  BSC: subst_stream={}, meta={}, consensus={}",
             bsc_subst.len(), bsc_meta.len(), bsc_consensus.len());

    // Serialize
    let mut output = Vec::new();
    output.extend_from_slice(b"TSB1"); // Template SuBstitution v1
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_subst.len() as u64);
    write_u64(&mut output, bsc_meta.len() as u64);
    write_u64(&mut output, bsc_consensus.len() as u64);
    output.extend_from_slice(&bsc_subst);
    output.extend_from_slice(&bsc_meta);
    output.extend_from_slice(&bsc_consensus);

    let total = output.len();
    eprintln!("  Template-subst: {} → {} bytes ({:.2}x)",
             total_bases, total, total_bases as f64 / total as f64);

    Ok(output)
}

// ── Minimizer-based read-as-template substitution ───────────────────────────

/// Compute canonical minimizer hashes for a read using minimizer-iter.
/// Returns deduplicated sorted list of canonical (min(fwd, rc)) k-mer hashes.
fn compute_read_minimizers(seq: &[u8], k: usize, w: u16) -> Vec<u64> {
    if seq.len() < k + w as usize - 1 {
        return Vec::new();
    }

    let positions: Vec<usize> = MinimizerBuilder::<u64>::new()
        .minimizer_size(k)
        .width(w)
        .canonical()
        .iter_pos(seq)
        .map(|(pos, _is_rc)| pos)
        .collect();

    let mut hashes: Vec<u64> = Vec::with_capacity(positions.len());
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            hashes.push(fwd.min(rc));
        }
    }

    hashes.sort_unstable();
    hashes.dedup();
    hashes
}

/// Hamming distance against reverse complement without allocating.
/// Compares seq1 reversed-complemented against seq2.
#[inline]
fn hamming_distance_rc(seq1: &[u8], seq2: &[u8], threshold: usize) -> Option<usize> {
    if seq1.len() != seq2.len() { return None; }
    let n = seq1.len();
    let mut dist = 0;
    for i in 0..n {
        let rc_base = match seq1[n - 1 - i] {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            x => x,
        };
        if rc_base != seq2[i] {
            dist += 1;
            if dist > threshold { return None; }
        }
    }
    Some(dist)
}

/// Compress sequences using minimizer-based read-as-template substitution.
///
/// Uses reads themselves as templates (no de Bruijn graph):
/// 1. Compute minimizers per read (O(n) via minimizer-iter)
/// 2. Greedy clustering: each read maps to best-matching earlier read
/// 3. Substitute matching bases with \0
/// 4. BSC compress the single stream
pub fn compress_sequences_minimizer_subst(
    sequences: &[String],
    mini_k: usize,
    mini_w: u16,
    mismatch_threshold: usize,
) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    if seqs.is_empty() { return Ok(Vec::new()); }

    let max_candidates = 10usize;
    let max_entries_per_minimizer = 50usize;

    let total_bases: usize = seqs.iter().map(|s| s.len()).sum();
    eprintln!("Minimizer-subst: {} reads, mini_k={}, mini_w={}, max_mm={}",
             seqs.len(), mini_k, mini_w, mismatch_threshold);

    // Step 1: Compute minimizers for all reads in parallel
    eprintln!("  Step 1: Computing minimizers...");
    let t = std::time::Instant::now();
    let all_minimizers: Vec<Vec<u64>> = seqs
        .par_iter()
        .map(|seq| compute_read_minimizers(seq, mini_k, mini_w))
        .collect();
    let total_mini: usize = all_minimizers.iter().map(|m| m.len()).sum();
    let avg_minimizers = total_mini as f64 / all_minimizers.len() as f64;
    eprintln!("  Minimizers computed: {:.1}s (avg {:.1} unique per read)",
             t.elapsed().as_secs_f64(), avg_minimizers);

    // Step 2: Greedy sequential clustering
    eprintln!("  Step 2: Greedy clustering...");
    let t = std::time::Instant::now();

    let mut index: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    index.reserve(total_mini);

    // mapping[i] = Some((ref_read_id, is_rc)) if mapped, None if representative
    let mut mappings: Vec<Option<(u32, bool)>> = vec![None; seqs.len()];
    let mut num_mapped = 0usize;

    for i in 0..seqs.len() {
        let minimizers = &all_minimizers[i];

        if !minimizers.is_empty() {
            // Count shared minimizers per candidate read
            let mut candidate_counts: FxHashMap<u32, u16> = FxHashMap::default();
            for &canon in minimizers {
                if let Some(read_ids) = index.get(&canon) {
                    for &rid in read_ids {
                        *candidate_counts.entry(rid).or_insert(0) += 1;
                    }
                }
            }

            if !candidate_counts.is_empty() {
                // Sort by shared count descending, take top candidates
                let mut candidates: Vec<(u32, u16)> = candidate_counts.into_iter().collect();
                candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));
                candidates.truncate(max_candidates);

                let mut best_ref: Option<(u32, bool)> = None;
                let mut best_dist = usize::MAX;

                for &(rid, _shared) in &candidates {
                    let ref_read = &seqs[rid as usize];
                    let query_read = &seqs[i];

                    if ref_read.len() != query_read.len() { continue; }

                    let threshold = best_dist.min(mismatch_threshold);

                    // Forward match
                    if let Some(d) = hamming_distance_within(ref_read, query_read, threshold) {
                        if d < best_dist {
                            best_dist = d;
                            best_ref = Some((rid, false));
                            if d == 0 { break; }
                        }
                    }

                    // Reverse complement match (allocation-free)
                    if best_dist > 0 {
                        let threshold = best_dist.min(mismatch_threshold);
                        if let Some(d) = hamming_distance_rc(ref_read, query_read, threshold) {
                            if d < best_dist {
                                best_dist = d;
                                best_ref = Some((rid, true));
                                if d == 0 { break; }
                            }
                        }
                    }
                }

                if let Some(best) = best_ref {
                    mappings[i] = Some(best);
                    num_mapped += 1;
                }
            }
        }

        // Add ALL reads' minimizers to the index (allows chains, better clustering)
        for &canon in minimizers {
            let bucket = index.entry(canon).or_default();
            if bucket.len() < max_entries_per_minimizer {
                bucket.push(i as u32);
            }
        }

        if (i + 1) % 100_000 == 0 {
            eprintln!("    ... {}/{} reads, {} mapped ({:.1}%)",
                     i + 1, seqs.len(), num_mapped,
                     num_mapped as f64 / (i + 1) as f64 * 100.0);
        }
    }

    let num_reps = seqs.len() - num_mapped;
    eprintln!("  Clustering: {}/{} mapped ({:.1}%), {} reps ({:.1}s)",
             num_mapped, seqs.len(),
             num_mapped as f64 / seqs.len() as f64 * 100.0,
             num_reps, t.elapsed().as_secs_f64());

    // Step 3: Build substituted stream
    eprintln!("  Step 3: Building substituted stream...");
    let mut subst_stream = Vec::with_capacity(total_bases);
    let mut meta_bytes: Vec<u8> = Vec::new();
    let mut total_subst = 0u64;

    write_varint(&mut meta_bytes, seqs.len() as u64);

    for (i, mapping) in mappings.iter().enumerate() {
        let read = &seqs[i];
        match mapping {
            Some(m) => {
                let (ref_id, is_rc) = *m;
                meta_bytes.push(1); // mapped
                write_varint(&mut meta_bytes, ref_id as u64);
                meta_bytes.push(if is_rc { 1 } else { 0 });

                let ref_read = &seqs[ref_id as usize];

                if is_rc {
                    // Compare against RC of reference, allocation-free
                    let n = ref_read.len().min(read.len());
                    for j in 0..n {
                        let rc_base = match ref_read[n - 1 - j] {
                            b'A' | b'a' => b'T',
                            b'C' | b'c' => b'G',
                            b'G' | b'g' => b'C',
                            b'T' | b't' => b'A',
                            x => x,
                        };
                        if read[j] == rc_base {
                            subst_stream.push(0u8);
                            total_subst += 1;
                        } else {
                            subst_stream.push(read[j]);
                        }
                    }
                } else {
                    for (j, &read_base) in read.iter().enumerate() {
                        if j < ref_read.len() && read_base == ref_read[j] {
                            subst_stream.push(0u8);
                            total_subst += 1;
                        } else {
                            subst_stream.push(read_base);
                        }
                    }
                }
            }
            None => {
                meta_bytes.push(0); // representative
                subst_stream.extend_from_slice(read);
            }
        }
    }

    let avg_subst = if num_mapped > 0 { total_subst as f64 / num_mapped as f64 } else { 0.0 };
    eprintln!("  Substituted: {} bases → \\0 (avg {:.1}/read across {} mapped reads)",
             total_subst, avg_subst, num_mapped);

    // Step 4: BSC compress
    eprintln!("  Step 4: BSC-compressing...");
    let t = std::time::Instant::now();
    let bsc_subst = bsc::compress_parallel_adaptive(&subst_stream)?;
    let bsc_meta = bsc::compress_parallel(&meta_bytes)?;
    eprintln!("  BSC done: {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("  Raw: subst_stream={}, meta={}",
             subst_stream.len(), meta_bytes.len());
    eprintln!("  BSC: subst_stream={}, meta={}",
             bsc_subst.len(), bsc_meta.len());

    // Serialize (MSB1 = Minimizer SuBstitution v1)
    let mut output = Vec::new();
    output.extend_from_slice(b"MSB1");
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_subst.len() as u64);
    write_u64(&mut output, bsc_meta.len() as u64);
    output.extend_from_slice(&bsc_subst);
    output.extend_from_slice(&bsc_meta);

    let total = output.len();
    eprintln!("  Minimizer-subst: {} → {} bytes ({:.2}x)",
             total_bases, total, total_bases as f64 / total as f64);

    Ok(output)
}

// ── Hybrid encoding (SPRING-style) ──────────────────────────────────────────

/// Compress sequences using hybrid encoding (SPRING-style).
///
/// Like compress_sequences_template_subst, uses the de Bruijn graph pipeline
/// (steps 1-3) to build templates and map reads. But instead of substituting
/// matching bases with \0 in a single stream, it:
///
/// - **Singletons**: kept as raw ASCII in one large BSC stream (preserves BWT context)
/// - **Mapped reads**: REMOVED from BSC stream entirely, encoded as compact metadata:
///   - Flag stream: 1 byte per read (0=singleton, 1=mapped)
///   - Mapping metadata: (template_id, position, is_rc, mismatch_count) per mapped read
///   - Mismatch positions: delta-encoded per mapped read
///   - Mismatch bases: one base per mismatch
///   - Consensus: 2-bit packed templates
///
/// This avoids inflating the BSC stream with \0 bytes for mapped reads.
pub fn compress_sequences_template_hybrid(
    sequences: &[String],
    params: &TemplateParams,
) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    if seqs.is_empty() { return Ok(Vec::new()); }

    let k = params.k;
    let mismatch_threshold = params.mismatch_threshold;

    eprintln!("Template-hybrid: {} reads, k={}, max_mismatches={}",
             seqs.len(), k, mismatch_threshold);

    // Steps 1-3: Build graph, extract templates, map reads
    eprintln!("  Step 1: Building graph (k={})...", k);
    let t = std::time::Instant::now();
    let mut graph = TemplateGraph::from_sequences(&seqs, k);

    if !params.fast {
        let total_kmers = graph.kmers.len();
        let singletons_count = graph.kmers.values().filter(|&&c| c == 1).count();
        let singleton_frac = singletons_count as f64 / total_kmers.max(1) as f64;
        if params.tip_trim_threshold > 0 && singleton_frac < 0.95 {
            graph.tip_trim(params.tip_trim_threshold);
        }
    }
    eprintln!("  Graph built: {:.1}s", t.elapsed().as_secs_f64());

    let min_length = seqs.iter().map(|s| s.len()).max().unwrap_or(150);
    eprintln!("  Step 2: Extracting templates (min_length={})...", min_length);
    let t = std::time::Instant::now();
    let templates = if params.fast {
        extract_templates_fast(&mut graph, min_length)
    } else {
        extract_templates(&graph, min_length)
    };
    drop(graph);
    eprintln!("  Templates extracted: {:.1}s", t.elapsed().as_secs_f64());

    let map_threshold = mismatch_threshold;
    eprintln!("  Step 3: Mapping reads (mismatch_threshold={}, syncmers={})...",
             map_threshold, params.use_syncmers);
    let t = std::time::Instant::now();
    let mappings = if params.use_syncmers {
        let index = build_template_index_syncmer(&templates, k);
        let m = map_reads_to_templates_syncmer(&seqs, &templates, &index, k, map_threshold);
        drop(index);
        m
    } else {
        let index = build_template_index(&templates, k);
        let m = map_reads_to_templates(&seqs, &templates, &index, k, map_threshold);
        drop(index);
        m
    };

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    eprintln!("  Mapping: {}/{} mapped ({:.1}%) ({:.1}s)",
             mapped, seqs.len(), mapped as f64 / seqs.len() as f64 * 100.0,
             t.elapsed().as_secs_f64());

    // Step 4: Hybrid encoding — separate singletons from mapped reads
    eprintln!("  Step 4: Hybrid encoding...");

    // Determine which templates are used and remap IDs
    let mut used_ids: FxHashSet<usize> = FxHashSet::default();
    for m in mappings.iter().flatten() {
        used_ids.insert(m.template_id);
    }
    let mut id_map: FxHashMap<usize, usize> = FxHashMap::default();
    let mut used_templates: Vec<&Template> = Vec::new();
    let mut sorted_ids: Vec<usize> = used_ids.into_iter().collect();
    sorted_ids.sort();
    for &uid in &sorted_ids {
        id_map.insert(uid, used_templates.len());
        used_templates.push(&templates[uid]);
    }

    // Stream 1: Flag stream — 1 byte per read (0=singleton, 1=mapped)
    let mut flag_stream: Vec<u8> = Vec::with_capacity(seqs.len());

    // Stream 2: Singleton bases — raw ASCII, one big BSC block
    let mut singleton_stream: Vec<u8> = Vec::new();

    // Stream 3: Mapping metadata — (template_id, position, is_rc) per mapped read
    let mut map_meta_stream: Vec<u8> = Vec::new();

    // Stream 4: Mismatch positions — delta-encoded per mapped read
    let mut mismatch_pos_stream: Vec<u8> = Vec::new();

    // Stream 5: Mismatch bases — one base per mismatch
    let mut mismatch_base_stream: Vec<u8> = Vec::new();

    let mut total_mismatches = 0u64;

    for (read_idx, mapping) in mappings.iter().enumerate() {
        let read = &seqs[read_idx];
        match mapping {
            Some(m) => {
                flag_stream.push(1);

                let mapped_id = id_map[&m.template_id];
                write_varint(&mut map_meta_stream, mapped_id as u64);
                write_varint(&mut map_meta_stream, m.position as u64);
                map_meta_stream.push(if m.is_rc { 1 } else { 0 });

                // Compute mismatches
                let template = &templates[m.template_id];
                let region = &template.sequence[m.position as usize..m.position as usize + read.len()];

                let compare_ref = if m.is_rc {
                    reverse_complement(region)
                } else {
                    region.to_vec()
                };

                let mut mm_count = 0u32;
                let mut prev_pos = 0u16;
                // Collect mismatches first to write count before positions
                let mut mm_positions: Vec<u16> = Vec::new();
                let mut mm_bases: Vec<u8> = Vec::new();

                for (i, (&r, &c)) in read.iter().zip(compare_ref.iter()).enumerate() {
                    if r != c {
                        mm_positions.push(i as u16);
                        mm_bases.push(r);
                        mm_count += 1;
                    }
                }

                total_mismatches += mm_count as u64;

                // Write mismatch count
                write_varint(&mut mismatch_pos_stream, mm_count as u64);

                // Delta-encoded positions
                for &pos in &mm_positions {
                    write_varint(&mut mismatch_pos_stream, (pos - prev_pos) as u64);
                    prev_pos = pos;
                }

                // Mismatch bases (raw, not delta-encoded — only 4 possible values)
                mismatch_base_stream.extend_from_slice(&mm_bases);
            }
            None => {
                flag_stream.push(0);
                singleton_stream.extend_from_slice(read);
            }
        }
    }

    let avg_mm = if mapped > 0 { total_mismatches as f64 / mapped as f64 } else { 0.0 };
    let num_singletons = seqs.len() - mapped;
    eprintln!("  Hybrid: {} mapped (avg {:.1} mm), {} singletons",
             mapped, avg_mm, num_singletons);
    eprintln!("  Raw streams: flags={}, singletons={}, map_meta={}, mm_pos={}, mm_bases={}",
             flag_stream.len(), singleton_stream.len(), map_meta_stream.len(),
             mismatch_pos_stream.len(), mismatch_base_stream.len());

    // Template consensus for decompression
    let consensus_seqs: Vec<Vec<u8>> = used_templates.iter().map(|t| t.sequence.clone()).collect();
    let consensus_bytes = pack_dna_2bit(&consensus_seqs);

    // Step 5: BSC-compress all streams
    eprintln!("  Step 5: BSC-compressing...");
    let t = std::time::Instant::now();

    // Use adaptive BSC for the big singleton stream (like baseline)
    let bsc_singletons = bsc::compress_parallel_adaptive(&singleton_stream)?;
    let bsc_flags = bsc::compress_parallel(&flag_stream)?;
    let bsc_map_meta = bsc::compress_parallel(&map_meta_stream)?;
    let bsc_mm_pos = bsc::compress_parallel(&mismatch_pos_stream)?;
    let bsc_mm_bases = bsc::compress_parallel(&mismatch_base_stream)?;
    let bsc_consensus = bsc::compress_parallel(&consensus_bytes)?;

    eprintln!("  BSC done: {:.1}s", t.elapsed().as_secs_f64());
    eprintln!("  BSC: singletons={}, flags={}, map_meta={}, mm_pos={}, mm_bases={}, consensus={}",
             bsc_singletons.len(), bsc_flags.len(), bsc_map_meta.len(),
             bsc_mm_pos.len(), bsc_mm_bases.len(), bsc_consensus.len());

    let total_compressed = bsc_singletons.len() + bsc_flags.len() + bsc_map_meta.len()
        + bsc_mm_pos.len() + bsc_mm_bases.len() + bsc_consensus.len();
    let total_bases: usize = seqs.iter().map(|s| s.len()).sum();
    eprintln!("  Total: {} bytes ({:.2}x)", total_compressed,
             total_bases as f64 / total_compressed as f64);

    // Serialize (THB1 = Template HyBrid v1)
    let mut output = Vec::new();
    output.extend_from_slice(b"THB1");
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_singletons.len() as u64);
    write_u64(&mut output, bsc_flags.len() as u64);
    write_u64(&mut output, bsc_map_meta.len() as u64);
    write_u64(&mut output, bsc_mm_pos.len() as u64);
    write_u64(&mut output, bsc_mm_bases.len() as u64);
    write_u64(&mut output, bsc_consensus.len() as u64);
    output.extend_from_slice(&bsc_singletons);
    output.extend_from_slice(&bsc_flags);
    output.extend_from_slice(&bsc_map_meta);
    output.extend_from_slice(&bsc_mm_pos);
    output.extend_from_slice(&bsc_mm_bases);
    output.extend_from_slice(&bsc_consensus);

    let total = output.len();
    eprintln!("  Template-hybrid: {} → {} bytes ({:.2}x)",
             total_bases, total, total_bases as f64 / total as f64);

    Ok(output)
}

// ── Soft reorder experiments ────────────────────────────────────────────────

/// Soft reorder: sort reads by minimizer hash, BSC compress, store permutation.
/// Returns (bsc_reordered_stream, bsc_permutation, stats_string).
pub fn compress_sequences_reordered(
    sequences: &[String],
) -> Result<(Vec<u8>, Vec<u8>)> {
    let mini_k = 15usize;
    let mini_w = 11u16;

    eprintln!("  Computing sort keys (mini_k={}, mini_w={})...", mini_k, mini_w);
    let t = std::time::Instant::now();

    // Sort key = smallest canonical minimizer hash
    let sort_keys: Vec<u64> = sequences
        .par_iter()
        .map(|seq| {
            let mins = compute_read_minimizers(seq.as_bytes(), mini_k, mini_w);
            mins.into_iter().min().unwrap_or(u64::MAX)
        })
        .collect();
    eprintln!("  Sort keys: {:.2}s", t.elapsed().as_secs_f64());

    // Stable argsort
    let t = std::time::Instant::now();
    let mut perm: Vec<u32> = (0..sequences.len() as u32).collect();
    perm.sort_by_key(|&i| sort_keys[i as usize]);
    eprintln!("  Argsort: {:.2}s", t.elapsed().as_secs_f64());

    // Write reordered reads
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let mut reordered = Vec::with_capacity(total_bases);
    for &idx in &perm {
        reordered.extend_from_slice(sequences[idx as usize].as_bytes());
    }

    // BSC compress reordered stream
    let t = std::time::Instant::now();
    let bsc_data = bsc::compress_parallel_adaptive(&reordered)?;
    eprintln!("  BSC reordered: {:.2}s ({} → {} bytes)",
             t.elapsed().as_secs_f64(), total_bases, bsc_data.len());

    // Permutation: store as u32 LE bytes, BSC'd
    let perm_bytes: Vec<u8> = perm.iter()
        .flat_map(|&idx| idx.to_le_bytes())
        .collect();
    let bsc_perm = bsc::compress_parallel(&perm_bytes)?;
    eprintln!("  Permutation: {} raw → {} BSC'd",
             perm_bytes.len(), bsc_perm.len());

    Ok((bsc_data, bsc_perm))
}

/// Columnar layout: transpose reads (column-major) then BSC.
/// No extra metadata needed — decompressor just transposes back.
pub fn compress_sequences_columnar(
    sequences: &[String],
) -> Result<Vec<u8>> {
    let n = sequences.len();
    if n == 0 { return Ok(Vec::new()); }
    let read_len = sequences[0].len();
    let total_bases = n * read_len;

    eprintln!("  Transposing {} reads x {} bp...", n, read_len);
    let t = std::time::Instant::now();

    let mut columnar = Vec::with_capacity(total_bases);
    // Column-major: all bases at position 0, then position 1, etc.
    let seqs_bytes: Vec<&[u8]> = sequences.iter().map(|s| s.as_bytes()).collect();
    for col in 0..read_len {
        for seq in &seqs_bytes {
            columnar.push(seq[col]);
        }
    }
    eprintln!("  Transpose: {:.3}s", t.elapsed().as_secs_f64());

    let t = std::time::Instant::now();
    let bsc_data = bsc::compress_parallel_adaptive(&columnar)?;
    eprintln!("  BSC columnar: {:.2}s ({} → {} bytes)",
             t.elapsed().as_secs_f64(), total_bases, bsc_data.len());

    Ok(bsc_data)
}

/// Anchor-dictionary compression: encode reads against previously-seen similar reads.
///
/// Uses syncmer anchors (k=31) to build an online dictionary of seen reads.
/// For each new read:
/// 1. Extract syncmer anchor hashes
/// 2. Look up anchors in dictionary → candidate previous reads with shared anchors
/// 3. Rank candidates by shared anchor count, verify top candidates via Hamming distance
/// 4. If best match is within threshold, encode as \0-substitution against match
/// 5. Add current read's anchors to dictionary
///
/// No de Bruijn graph, no template extraction — just read-to-read matching.
/// This is fast and scales well to large datasets.
pub fn compress_sequences_anchor_dict(
    sequences: &[String],
    anchor_k: usize,
    mismatch_threshold: usize,
) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    if seqs.is_empty() { return Ok(Vec::new()); }

    let max_candidates = 10usize;
    let max_entries_per_anchor = 50usize;

    let total_bases: usize = seqs.iter().map(|s| s.len()).sum();
    let s = anchor_k.saturating_sub(3).max(5);
    let t_end = anchor_k - s;

    eprintln!("Anchor-dict: {} reads, k={}, s={}, max_mm={}",
             seqs.len(), anchor_k, s, mismatch_threshold);

    // Step 1: Compute syncmer anchors for all reads in parallel
    eprintln!("  Step 1: Computing syncmer anchors...");
    let t = std::time::Instant::now();
    let all_anchors: Vec<Vec<u64>> = seqs
        .par_iter()
        .map(|seq| {
            if seq.len() < anchor_k { return Vec::new(); }
            let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], seq);
            let mut hashes = Vec::with_capacity(positions.len());
            for pos in positions {
                if pos + anchor_k > seq.len() { continue; }
                let kmer = &seq[pos..pos + anchor_k];
                if let Some(fwd) = kmer_to_hash(kmer) {
                    let rc = reverse_complement_hash(fwd, anchor_k);
                    hashes.push(fwd.min(rc));
                }
            }
            hashes.sort_unstable();
            hashes.dedup();
            hashes
        })
        .collect();
    let total_anchors: usize = all_anchors.iter().map(|a| a.len()).sum();
    let avg_anchors = total_anchors as f64 / all_anchors.len().max(1) as f64;
    eprintln!("  Anchors computed: {:.1}s (avg {:.1} unique per read)",
             t.elapsed().as_secs_f64(), avg_anchors);

    // Step 2: Sequential dictionary matching
    eprintln!("  Step 2: Dictionary matching...");
    let t = std::time::Instant::now();

    let mut index: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    index.reserve(total_anchors);

    // mapping[i] = Some((ref_read_id, is_rc)) if mapped, None if representative
    let mut mappings: Vec<Option<(u32, bool)>> = vec![None; seqs.len()];
    let mut num_mapped = 0usize;
    let mut total_dist = 0u64;

    for i in 0..seqs.len() {
        let anchors = &all_anchors[i];

        if !anchors.is_empty() {
            // Count shared anchors per candidate read
            let mut candidate_counts: FxHashMap<u32, u16> = FxHashMap::default();
            for &canon in anchors {
                if let Some(read_ids) = index.get(&canon) {
                    for &rid in read_ids {
                        *candidate_counts.entry(rid).or_insert(0) += 1;
                    }
                }
            }

            if !candidate_counts.is_empty() {
                // Sort by shared count descending, take top candidates
                let mut candidates: Vec<(u32, u16)> = candidate_counts.into_iter().collect();
                candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));
                candidates.truncate(max_candidates);

                let mut best_ref: Option<(u32, bool)> = None;
                let mut best_dist = usize::MAX;

                for &(rid, _shared) in &candidates {
                    let ref_read = &seqs[rid as usize];
                    let query_read = &seqs[i];

                    if ref_read.len() != query_read.len() { continue; }

                    let threshold = best_dist.min(mismatch_threshold);

                    // Forward match
                    if let Some(d) = hamming_distance_within(ref_read, query_read, threshold) {
                        if d < best_dist {
                            best_dist = d;
                            best_ref = Some((rid, false));
                            if d == 0 { break; }
                        }
                    }

                    // RC match (allocation-free)
                    if best_dist > 0 {
                        let threshold = best_dist.min(mismatch_threshold);
                        if let Some(d) = hamming_distance_rc(ref_read, query_read, threshold) {
                            if d < best_dist {
                                best_dist = d;
                                best_ref = Some((rid, true));
                                if d == 0 { break; }
                            }
                        }
                    }
                }

                if let Some(best) = best_ref {
                    mappings[i] = Some(best);
                    num_mapped += 1;
                    total_dist += best_dist as u64;
                }
            }
        }

        // Add ALL reads' anchors to the index
        for &canon in anchors {
            let bucket = index.entry(canon).or_default();
            if bucket.len() < max_entries_per_anchor {
                bucket.push(i as u32);
            }
        }

        if (i + 1) % 100_000 == 0 {
            eprintln!("    ... {}/{} reads, {} matched ({:.1}%)",
                     i + 1, seqs.len(), num_mapped,
                     num_mapped as f64 / (i + 1) as f64 * 100.0);
        }
    }

    let avg_dist = if num_mapped > 0 { total_dist as f64 / num_mapped as f64 } else { 0.0 };
    eprintln!("  Matching: {}/{} matched ({:.1}%), avg dist={:.1} ({:.1}s)",
             num_mapped, seqs.len(),
             num_mapped as f64 / seqs.len() as f64 * 100.0,
             avg_dist, t.elapsed().as_secs_f64());

    // Step 3: Build substituted stream (\0 for matching bases)
    eprintln!("  Step 3: Building substituted stream...");
    let mut subst_stream = Vec::with_capacity(total_bases);
    let mut meta_bytes: Vec<u8> = Vec::new();
    let mut total_subst = 0u64;

    write_varint(&mut meta_bytes, seqs.len() as u64);

    for (i, mapping) in mappings.iter().enumerate() {
        let read = &seqs[i];
        if let Some(&(ref_id, is_rc)) = mapping.as_ref() {
            meta_bytes.push(1);
            write_varint(&mut meta_bytes, ref_id as u64);
            meta_bytes.push(if is_rc { 1 } else { 0 });

            let ref_read = &seqs[ref_id as usize];

            if is_rc {
                let n = ref_read.len().min(read.len());
                for j in 0..n {
                    let rc_base = match ref_read[n - 1 - j] {
                        b'A' | b'a' => b'T',
                        b'C' | b'c' => b'G',
                        b'G' | b'g' => b'C',
                        b'T' | b't' => b'A',
                        x => x,
                    };
                    if read[j] == rc_base {
                        subst_stream.push(0u8);
                        total_subst += 1;
                    } else {
                        subst_stream.push(read[j]);
                    }
                }
            } else {
                for (j, &read_base) in read.iter().enumerate() {
                    if j < ref_read.len() && read_base == ref_read[j] {
                        subst_stream.push(0u8);
                        total_subst += 1;
                    } else {
                        subst_stream.push(read_base);
                    }
                }
            }
        } else {
            meta_bytes.push(0);
            subst_stream.extend_from_slice(read);
        }
    }

    let avg_subst = if num_mapped > 0 { total_subst as f64 / num_mapped as f64 } else { 0.0 };
    eprintln!("  Substituted: {} bases → \\0 (avg {:.1}/read across {} matched reads)",
             total_subst, avg_subst, num_mapped);

    // Step 4: BSC compress
    eprintln!("  Step 4: BSC-compressing...");
    let t = std::time::Instant::now();
    let bsc_subst = bsc::compress_parallel_adaptive(&subst_stream)?;
    let bsc_meta = bsc::compress_parallel(&meta_bytes)?;
    eprintln!("  BSC done: {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("  Raw: subst_stream={}, meta={}",
             subst_stream.len(), meta_bytes.len());
    eprintln!("  BSC: subst_stream={}, meta={}",
             bsc_subst.len(), bsc_meta.len());

    // Serialize (ADB1 = Anchor Dictionary v1)
    let mut output = Vec::new();
    output.extend_from_slice(b"ADB1");
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_subst.len() as u64);
    write_u64(&mut output, bsc_meta.len() as u64);
    output.extend_from_slice(&bsc_subst);
    output.extend_from_slice(&bsc_meta);

    let total = output.len();
    eprintln!("  Anchor-dict: {} → {} bytes ({:.2}x)",
             total_bases, total, total_bases as f64 / total as f64);

    Ok(output)
}

/// Sort reads by syncmer-based sort key, then BSC.
/// Uses primary syncmer hash as sort key to group reads from same genomic region.
/// Stores permutation for decompression (order reconstruction).
pub fn compress_sequences_syncmer_reorder(
    sequences: &[String],
    anchor_k: usize,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let n = sequences.len();
    if n == 0 { return Ok((Vec::new(), Vec::new())); }

    let s = anchor_k.saturating_sub(3).max(5);
    let t_end = anchor_k - s;

    eprintln!("  Computing syncmer sort keys (k={}, s={})...", anchor_k, s);
    let t = std::time::Instant::now();

    // Sort key = smallest canonical syncmer hash in each read
    let sort_keys: Vec<u64> = sequences
        .par_iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            if seq_bytes.len() < anchor_k { return u64::MAX; }
            let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], seq_bytes);
            let mut min_hash = u64::MAX;
            for pos in positions {
                if pos + anchor_k > seq_bytes.len() { continue; }
                let kmer = &seq_bytes[pos..pos + anchor_k];
                if let Some(fwd) = kmer_to_hash(kmer) {
                    let rc = reverse_complement_hash(fwd, anchor_k);
                    let canon = fwd.min(rc);
                    if canon < min_hash { min_hash = canon; }
                }
            }
            min_hash
        })
        .collect();
    eprintln!("  Sort keys: {:.2}s", t.elapsed().as_secs_f64());

    let t = std::time::Instant::now();
    let mut perm: Vec<u32> = (0..n as u32).collect();
    perm.sort_by_key(|&i| sort_keys[i as usize]);
    eprintln!("  Argsort: {:.2}s", t.elapsed().as_secs_f64());

    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let mut reordered = Vec::with_capacity(total_bases);
    for &idx in &perm {
        reordered.extend_from_slice(sequences[idx as usize].as_bytes());
    }

    let t = std::time::Instant::now();
    let bsc_data = bsc::compress_parallel_adaptive(&reordered)?;
    eprintln!("  BSC syncmer-reorder: {:.2}s ({} → {} bytes, {:.2}x)",
             t.elapsed().as_secs_f64(), total_bases, bsc_data.len(),
             total_bases as f64 / bsc_data.len() as f64);

    let perm_bytes: Vec<u8> = perm.iter()
        .flat_map(|&idx| idx.to_le_bytes())
        .collect();
    let bsc_perm = bsc::compress_parallel(&perm_bytes)?;
    eprintln!("  Permutation: {} raw → {} BSC'd", perm_bytes.len(), bsc_perm.len());

    Ok((bsc_data, bsc_perm))
}

/// Sort reads by syncmer + \0 substitution: reorder AND encode against neighbors.
/// After reordering by syncmer hash, encode each read as delta from its predecessor
/// (since nearby reads in sorted order likely overlap on the genome).
pub fn compress_sequences_syncmer_reorder_delta(
    sequences: &[String],
    anchor_k: usize,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let n = sequences.len();
    if n == 0 { return Ok((Vec::new(), Vec::new())); }

    let s = anchor_k.saturating_sub(3).max(5);
    let t_end = anchor_k - s;

    eprintln!("  Computing syncmer sort keys (k={}, s={})...", anchor_k, s);
    let t = std::time::Instant::now();

    let sort_keys: Vec<u64> = sequences
        .par_iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            if seq_bytes.len() < anchor_k { return u64::MAX; }
            let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], seq_bytes);
            let mut min_hash = u64::MAX;
            for pos in positions {
                if pos + anchor_k > seq_bytes.len() { continue; }
                let kmer = &seq_bytes[pos..pos + anchor_k];
                if let Some(fwd) = kmer_to_hash(kmer) {
                    let rc = reverse_complement_hash(fwd, anchor_k);
                    let canon = fwd.min(rc);
                    if canon < min_hash { min_hash = canon; }
                }
            }
            min_hash
        })
        .collect();
    eprintln!("  Sort keys: {:.2}s", t.elapsed().as_secs_f64());

    let t = std::time::Instant::now();
    let mut perm: Vec<u32> = (0..n as u32).collect();
    perm.sort_by_key(|&i| sort_keys[i as usize]);
    eprintln!("  Argsort: {:.2}s", t.elapsed().as_secs_f64());

    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();

    // Build reordered stream with \0 substitution against predecessor
    let t = std::time::Instant::now();
    let read_len = sequences[0].len();
    let mut subst_stream = Vec::with_capacity(total_bases);
    let mut total_subst = 0u64;
    let mut num_deltas = 0u64;

    // First read in sorted order goes as-is
    subst_stream.extend_from_slice(sequences[perm[0] as usize].as_bytes());

    for w in 1..n {
        let cur = sequences[perm[w] as usize].as_bytes();
        let prev = sequences[perm[w - 1] as usize].as_bytes();

        // Only delta-encode if they share the same sort key (same primary syncmer)
        if sort_keys[perm[w] as usize] == sort_keys[perm[w - 1] as usize]
            && cur.len() == prev.len()
        {
            // Try fwd match
            let mut fwd_dist = 0usize;
            for i in 0..cur.len() {
                if cur[i] != prev[i] { fwd_dist += 1; }
            }

            // Try RC match
            let mut rc_dist = 0usize;
            let prev_len = prev.len();
            for i in 0..cur.len() {
                let rc_base = match prev[prev_len - 1 - i] {
                    b'A' | b'a' => b'T',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    b'T' | b't' => b'A',
                    x => x,
                };
                if cur[i] != rc_base { rc_dist += 1; }
            }

            let threshold = read_len / 4; // 25% mismatch threshold
            let (use_rc, dist) = if rc_dist < fwd_dist { (true, rc_dist) } else { (false, fwd_dist) };

            if dist <= threshold {
                num_deltas += 1;
                if use_rc {
                    for i in 0..cur.len() {
                        let rc_base = match prev[prev_len - 1 - i] {
                            b'A' | b'a' => b'T',
                            b'C' | b'c' => b'G',
                            b'G' | b'g' => b'C',
                            b'T' | b't' => b'A',
                            x => x,
                        };
                        if cur[i] == rc_base {
                            subst_stream.push(0u8);
                            total_subst += 1;
                        } else {
                            subst_stream.push(cur[i]);
                        }
                    }
                } else {
                    for i in 0..cur.len() {
                        if cur[i] == prev[i] {
                            subst_stream.push(0u8);
                            total_subst += 1;
                        } else {
                            subst_stream.push(cur[i]);
                        }
                    }
                }
                continue;
            }
        }

        // No delta: emit raw
        subst_stream.extend_from_slice(cur);
    }

    let avg_subst = if num_deltas > 0 { total_subst as f64 / num_deltas as f64 } else { 0.0 };
    eprintln!("  Delta encoding: {} deltas out of {} ({:.1}%), avg {:.1} subst/read ({:.2}s)",
             num_deltas, n, num_deltas as f64 / n as f64 * 100.0,
             avg_subst, t.elapsed().as_secs_f64());

    let t = std::time::Instant::now();
    let bsc_data = bsc::compress_parallel_adaptive(&subst_stream)?;
    eprintln!("  BSC reorder+delta: {:.2}s ({} → {} bytes, {:.2}x)",
             t.elapsed().as_secs_f64(), total_bases, bsc_data.len(),
             total_bases as f64 / bsc_data.len() as f64);

    let perm_bytes: Vec<u8> = perm.iter()
        .flat_map(|&idx| idx.to_le_bytes())
        .collect();
    let bsc_perm = bsc::compress_parallel(&perm_bytes)?;

    Ok((bsc_data, bsc_perm))
}

/// Reordered + columnar: sort reads by minimizer, then transpose, then BSC.
pub fn compress_sequences_reordered_columnar(
    sequences: &[String],
) -> Result<(Vec<u8>, Vec<u8>)> {
    let mini_k = 15usize;
    let mini_w = 11u16;

    let sort_keys: Vec<u64> = sequences
        .par_iter()
        .map(|seq| {
            let mins = compute_read_minimizers(seq.as_bytes(), mini_k, mini_w);
            mins.into_iter().min().unwrap_or(u64::MAX)
        })
        .collect();

    let mut perm: Vec<u32> = (0..sequences.len() as u32).collect();
    perm.sort_by_key(|&i| sort_keys[i as usize]);

    let n = sequences.len();
    let read_len = sequences[0].len();
    let total_bases = n * read_len;

    // Reorder then transpose
    let mut columnar = Vec::with_capacity(total_bases);
    for col in 0..read_len {
        for &idx in &perm {
            columnar.push(sequences[idx as usize].as_bytes()[col]);
        }
    }

    let bsc_data = bsc::compress_parallel_adaptive(&columnar)?;

    let perm_bytes: Vec<u8> = perm.iter()
        .flat_map(|&idx| idx.to_le_bytes())
        .collect();
    let bsc_perm = bsc::compress_parallel(&perm_bytes)?;

    eprintln!("  Reorder+columnar: {} → {} + {} perm",
             total_bases, bsc_data.len(), bsc_perm.len());

    Ok((bsc_data, bsc_perm))
}
