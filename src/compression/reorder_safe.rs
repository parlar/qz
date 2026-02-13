//! Patent-safe read reordering strategies
//!
//! This module implements read reordering based on METADATA and SIMPLE STATISTICS,
//! NOT sequence similarity or homology detection.
//!
//! All methods here are patent-safe because they:
//! - Sort by metadata (flowcell coordinates, read length)
//! - Use simple aggregate statistics (GC content)
//! - Apply deterministic operations (lexicographic sorting)
//!
//! These do NOT detect sequence similarity or cluster by homology,
//! thus avoiding US20200058379A1 (SPRING patent).

use crate::io::fastq::FastqRecord;
use anyhow::Result;

/// Patent-safe reordering modes
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReorderMode {
    /// No reordering (preserve input order)
    None,
    /// Sort by flowcell tile and X/Y coordinates (metadata-based)
    Flowcell,
    /// Sort by GC content percentage (simple statistic)
    GcContent,
    /// Sort by read length (trivial sorting)
    Length,
    /// Sort lexicographically by sequence (deterministic)
    Lexicographic,
    /// Smart: flowcell + GC + length (multi-level metadata sort)
    Smart,
}

/// Flowcell coordinates extracted from Illumina read ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct FlowcellCoords {
    pub lane: u32,
    pub tile: u32,
    pub x: u32,
    pub y: u32,
}

/// Parse flowcell coordinates from Illumina read ID
/// Format: @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y [COMMENT]
pub fn parse_flowcell_coords(read_id: &str) -> Option<FlowcellCoords> {
    let id_without_at = read_id.strip_prefix('@').unwrap_or(read_id);

    // Split off comment if present
    let base_id = if let Some(space_pos) = id_without_at.find(' ') {
        &id_without_at[..space_pos]
    } else {
        id_without_at
    };

    // Parse colon-separated fields
    let fields: Vec<&str> = base_id.split(':').collect();

    if fields.len() < 7 {
        return None; // Not standard Illumina format
    }

    // INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y
    let lane = fields[3].parse().ok()?;
    let tile = fields[4].parse().ok()?;
    let x = fields[5].parse().ok()?;
    let y = fields[6].parse().ok()?;

    Some(FlowcellCoords { lane, tile, x, y })
}

/// Calculate GC content percentage (0-100)
pub fn gc_content(sequence: &str) -> u8 {
    if sequence.is_empty() {
        return 0;
    }

    let gc_count = sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
        .count();

    ((gc_count * 100) / sequence.len()).min(100) as u8
}

/// Reorder reads using patent-safe strategy
pub fn reorder_reads(reads: &mut [FastqRecord], mode: ReorderMode) -> Result<()> {
    match mode {
        ReorderMode::None => {
            // No reordering
            Ok(())
        }
        ReorderMode::Flowcell => {
            // Sort by flowcell coordinates (LANE, TILE, X, Y)
            reads.sort_by_key(|r| parse_flowcell_coords(&r.id));
            Ok(())
        }
        ReorderMode::GcContent => {
            // Sort by GC content percentage
            reads.sort_by_key(|r| gc_content(&r.sequence));
            Ok(())
        }
        ReorderMode::Length => {
            // Sort by read length
            reads.sort_by_key(|r| r.sequence.len());
            Ok(())
        }
        ReorderMode::Lexicographic => {
            // Sort alphabetically by sequence
            reads.sort_by(|a, b| a.sequence.cmp(&b.sequence));
            Ok(())
        }
        ReorderMode::Smart => {
            // Multi-level sort: flowcell coords → GC bin → length → X/Y
            reads.sort_by_key(|r| {
                let coords = parse_flowcell_coords(&r.id);
                let gc_bin = gc_content(&r.sequence) / 10; // 10% bins
                let length = r.sequence.len();

                if let Some(c) = coords {
                    (Some(c.lane), Some(c.tile), gc_bin, length, Some(c.x), Some(c.y))
                } else {
                    (None, None, gc_bin, length, None, None)
                }
            });
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_flowcell_coords() {
        let id = "@INSTRUMENT:123:FLOWCELL:1:2101:10000:1000";
        let coords = parse_flowcell_coords(id).unwrap();

        assert_eq!(coords.lane, 1);
        assert_eq!(coords.tile, 2101);
        assert_eq!(coords.x, 10000);
        assert_eq!(coords.y, 1000);
    }

    #[test]
    fn test_parse_flowcell_coords_with_comment() {
        let id = "@INSTRUMENT:123:FLOWCELL:1:2101:10000:1000 1:N:0:ATCG";
        let coords = parse_flowcell_coords(id).unwrap();

        assert_eq!(coords.lane, 1);
        assert_eq!(coords.tile, 2101);
    }

    #[test]
    fn test_parse_non_illumina_format() {
        let id = "READ_1";
        let coords = parse_flowcell_coords(id);
        assert!(coords.is_none());
    }

    #[test]
    fn test_gc_content() {
        assert_eq!(gc_content("AAAA"), 0);
        assert_eq!(gc_content("GGGG"), 100);
        assert_eq!(gc_content("CCCC"), 100);
        assert_eq!(gc_content("ATCG"), 50);
        assert_eq!(gc_content("ACGT"), 50);
        assert_eq!(gc_content(""), 0);
    }

    #[test]
    fn test_flowcell_reordering() {
        let mut reads = vec![
            FastqRecord {
                id: "@INST:1:FC:1:2101:10005:1000".to_string(),
                sequence: "ACGT".to_string(),
                quality: Some("IIII".to_string()),
            },
            FastqRecord {
                id: "@INST:1:FC:1:2101:10001:1000".to_string(),
                sequence: "GGGG".to_string(),
                quality: Some("IIII".to_string()),
            },
            FastqRecord {
                id: "@INST:1:FC:1:2101:10003:1000".to_string(),
                sequence: "TTTT".to_string(),
                quality: Some("IIII".to_string()),
            },
        ];

        reorder_reads(&mut reads, ReorderMode::Flowcell).unwrap();

        // Should be sorted by X coordinate
        assert_eq!(reads[0].id, "@INST:1:FC:1:2101:10001:1000");
        assert_eq!(reads[1].id, "@INST:1:FC:1:2101:10003:1000");
        assert_eq!(reads[2].id, "@INST:1:FC:1:2101:10005:1000");
    }

    #[test]
    fn test_gc_reordering() {
        let mut reads = vec![
            FastqRecord {
                id: "READ1".to_string(),
                sequence: "GGGG".to_string(), // GC = 100%
                quality: Some("IIII".to_string()),
            },
            FastqRecord {
                id: "READ2".to_string(),
                sequence: "AAAA".to_string(), // GC = 0%
                quality: Some("IIII".to_string()),
            },
            FastqRecord {
                id: "READ3".to_string(),
                sequence: "ACGT".to_string(), // GC = 50%
                quality: Some("IIII".to_string()),
            },
        ];

        reorder_reads(&mut reads, ReorderMode::GcContent).unwrap();

        // Should be sorted by GC content (ascending)
        assert_eq!(reads[0].sequence, "AAAA"); // 0%
        assert_eq!(reads[1].sequence, "ACGT"); // 50%
        assert_eq!(reads[2].sequence, "GGGG"); // 100%
    }

    #[test]
    fn test_length_reordering() {
        let mut reads = vec![
            FastqRecord {
                id: "READ1".to_string(),
                sequence: "ACGTACGT".to_string(), // Length 8
                quality: Some("IIIIIIII".to_string()),
            },
            FastqRecord {
                id: "READ2".to_string(),
                sequence: "ACGT".to_string(), // Length 4
                quality: Some("IIII".to_string()),
            },
            FastqRecord {
                id: "READ3".to_string(),
                sequence: "ACGTACGTAC".to_string(), // Length 10
                quality: Some("IIIIIIIIII".to_string()),
            },
        ];

        reorder_reads(&mut reads, ReorderMode::Length).unwrap();

        // Should be sorted by length (ascending)
        assert_eq!(reads[0].sequence.len(), 4);
        assert_eq!(reads[1].sequence.len(), 8);
        assert_eq!(reads[2].sequence.len(), 10);
    }
}
