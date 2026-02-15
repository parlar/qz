use qz::cli::{CompressArgs, DecompressArgs, HeaderCompressor, QualityCompressor, QualityMode, ReorderMode, SequenceCompressor};
use std::fs;
use tempfile::TempDir;

/// Helper to create default DecompressArgs
fn decompress_args(input: std::path::PathBuf, output: std::path::PathBuf, working_dir: std::path::PathBuf) -> DecompressArgs {
    DecompressArgs {
        input,
        output: vec![output],
        working_dir,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
    }
}

#[test]
fn test_compress_decompress_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len());

    for (orig, decomp) in original_lines.iter().zip(decompressed_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_decompress_gzipped_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };
    qz::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq.gz");
    let d_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: true,
        gzip_level: 6,
    };
    qz::compression::decompress(&d_args).unwrap();

    assert!(output_fastq.exists());

    let data = fs::read(&output_fastq).unwrap();
    assert!(data.len() >= 2);
    assert_eq!(data[0], 0x1f);
    assert_eq!(data[1], 0x8b);
}

#[test]
fn test_decompress_no_quality() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        no_quality: true,
        quality_compressor: QualityCompressor::Zstd,
        sequence_compressor: SequenceCompressor::Zstd,
        header_compressor: HeaderCompressor::Zstd,
        ..CompressArgs::default()
    };
    qz::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("IIIIIIII")); // Dummy quality
}

// ========================================
// ERROR HANDLING TESTS
// ========================================

#[test]
fn test_empty_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("empty.fastq");
    fs::write(&input_fastq, "").unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };

    // Should handle empty file gracefully (either error or create empty archive)
    let result = qz::compression::compress(&compress_args);
    if result.is_ok() {
        assert!(archive_path.exists());
    }
}

#[test]
fn test_single_read() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("single.fastq");
    fs::write(&input_fastq, "@read1\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("read1"));
    assert!(content.contains("ACGT"));
}

#[test]
fn test_long_read() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let long_seq = "ACGT".repeat(250);
    let long_qual = "I".repeat(1000);
    let fastq_data = format!("@long_read\n{}\n+\n{}\n", long_seq, long_qual);

    let input_fastq = temp_path.join("long.fastq");
    fs::write(&input_fastq, fastq_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains(&long_seq));
}

#[test]
fn test_reads_with_n_bases() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTNNNNACGT\n+\nIIII####IIII\n@read2\nNNNNNNNN\n+\n########\n";
    let input_fastq = temp_path.join("with_n.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("ACGTNNNNA")); // N bases should be preserved
}

// ========================================
// QUALITY MODE TESTS
// ========================================

#[test]
fn test_illumina_binning_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIHHGGFF\n";
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_mode: QualityMode::IlluminaBin,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());
    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("ACGTACGT"));
}

#[test]
fn test_binary_quality_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIHHGGFF\n";
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_mode: QualityMode::Binary,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());
}

// ========================================
// COMPRESSION LEVEL TESTS
// ========================================

#[test]
fn test_different_compression_levels() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n".repeat(100);
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, &test_data).unwrap();

    let levels = vec![1, 3, 9, 19];
    let mut sizes = Vec::new();

    for level in levels {
        let archive_path = temp_path.join(format!("test_level_{}.qz", level));
        let compress_args = CompressArgs {
            input: vec![input_fastq.clone()],
            output: archive_path.clone(),
            working_dir: temp_path.clone(),
            threads: 1,
            compression_level: level,
            ..CompressArgs::default()
        };

        qz::compression::compress(&compress_args).unwrap();
        let size = fs::metadata(&archive_path).unwrap().len();
        sizes.push((level, size));
    }

    assert_eq!(sizes.len(), 4);

    for (level, size) in &sizes {
        println!("Level {}: {} bytes", level, size);
    }
}

// ========================================
// MULTI-THREADING TEST
// ========================================

#[test]
fn test_multithreading_deterministic() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n".repeat(100);
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, &test_data).unwrap();

    // Compress with 1 thread
    let archive1 = temp_path.join("test1.qz");
    let compress_args1 = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive1.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressArgs::default()
    };
    qz::compression::compress(&compress_args1).unwrap();

    // Compress with 4 threads
    let archive2 = temp_path.join("test2.qz");
    let compress_args2 = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive2.clone(),
        working_dir: temp_path.clone(),
        threads: 4,
        ..CompressArgs::default()
    };
    qz::compression::compress(&compress_args2).unwrap();

    // Decompress both
    let output1 = temp_path.join("output1.fastq");
    let output2 = temp_path.join("output2.fastq");

    qz::compression::decompress(&decompress_args(archive1, output1.clone(), temp_path.clone())).unwrap();
    qz::compression::decompress(&decompress_args(archive2, output2.clone(), temp_path)).unwrap();

    // Both should produce identical output
    let content1 = fs::read_to_string(&output1).unwrap();
    let content2 = fs::read_to_string(&output2).unwrap();
    assert_eq!(content1, content2, "Multi-threaded compression should produce same output as single-threaded");
}

// ========== ERROR HANDLING TESTS ==========

#[test]
fn test_decompress_truncated_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create a valid archive
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        ..CompressArgs::default()
    };
    qz::compression::compress(&compress_args).unwrap();

    // Truncate to half size
    let archive_data = fs::read(&archive_path).unwrap();
    let truncated_path = temp_path.join("truncated.qz");
    fs::write(&truncated_path, &archive_data[..archive_data.len() / 2]).unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result = qz::compression::decompress(&decompress_args(
        truncated_path, output_fastq, temp_path,
    ));
    assert!(result.is_err(), "Decompressing a truncated archive should fail");
}

#[test]
fn test_decompress_corrupted_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create a valid archive
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        ..CompressArgs::default()
    };
    qz::compression::compress(&compress_args).unwrap();

    // Corrupt bytes in the data section
    let mut archive_data = fs::read(&archive_path).unwrap();
    let mid = archive_data.len() / 2;
    for i in mid..std::cmp::min(mid + 32, archive_data.len()) {
        archive_data[i] ^= 0xFF;
    }
    let corrupted_path = temp_path.join("corrupted.qz");
    fs::write(&corrupted_path, &archive_data).unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result = qz::compression::decompress(&decompress_args(
        corrupted_path, output_fastq, temp_path,
    ));
    assert!(result.is_err(), "Decompressing a corrupted archive should fail");
}

#[test]
fn test_decompress_zero_byte_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let empty_archive = temp_path.join("empty.qz");
    fs::write(&empty_archive, b"").unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result = qz::compression::decompress(&decompress_args(
        empty_archive, output_fastq, temp_path,
    ));
    assert!(result.is_err(), "Decompressing a zero-byte file should fail");
}

#[test]
fn test_twobit_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGTNNACGT\n+\nIIIIIIIIIIIIII\n@read2\nTGCATGCAATGCAT\n+\nHHHHHHHHHHHHHH\n@read3\nNNNACGTACGTNNN\n+\nAAAAAAAAAAAABB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        twobit: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    assert_eq!(original.trim(), decompressed.trim());
}

#[test]
fn test_header_template_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Illumina-style headers with common prefix and varying coordinates
    let test_data = "@INSTRUMENT:123:FLOWCELL:1:2101:1000:2000 1:N:0:ATCG\nACGTACGT\n+\nIIIIIIII\n\
                     @INSTRUMENT:123:FLOWCELL:1:2101:1001:2001 1:N:0:ATCG\nTGCATGCA\n+\nHHHHHHHH\n\
                     @INSTRUMENT:123:FLOWCELL:1:2101:1002:2002 1:N:0:ATCG\nAAAACCCC\n+\nBBBBBBBB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        header_template: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    assert_eq!(original.trim(), decompressed.trim());
}

#[test]
fn test_quality_modeling_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n@read3\nAAAACCCC\n+\nBBBBAAAA\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_modeling: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    // Quality modeling is lossy (i8 delta clamping), so check sequences match and qualities are close
    let original_content = fs::read_to_string(&input_fastq).unwrap();
    let original_lines: Vec<&str> = original_content.lines().collect();
    let decompressed_content = fs::read_to_string(&output_fastq).unwrap();
    let decompressed_lines: Vec<&str> = decompressed_content.lines().collect();

    // Check same number of lines
    assert_eq!(original_lines.len(), decompressed_lines.len());

    // Check headers and sequences match exactly
    for (i, (orig, decomp)) in original_lines.iter().zip(decompressed_lines.iter()).enumerate() {
        let line_type = i % 4;
        if line_type == 0 || line_type == 1 || line_type == 2 {
            // Header, sequence, separator should match exactly
            assert_eq!(orig.trim(), decomp.trim(), "Line {} mismatch", i);
        }
        // Quality lines may differ slightly due to modeling
    }
}

#[test]
fn test_openzl_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@read2\nTGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHH\n@read3\nAAAACCCCGGGGTTTT\n+\nFFFFFFFFFFFFFFFF\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_compressor: QualityCompressor::OpenZl,
        sequence_compressor: SequenceCompressor::OpenZl,
        header_compressor: HeaderCompressor::OpenZl,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len());

    for (orig, decomp) in original_lines.iter().zip(decompressed_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_fqzcomp_quality_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Multiple reads with varying quality profiles to test mean-quality sorting
    let test_data = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
                     @read2\nTGCATGCATGCATGCA\n+\n!!!!!!!!!!!!!!!!\n\
                     @read3\nAAAACCCCGGGGTTTT\n+\nFFFFFFFFFFFFFFFF\n\
                     @read4\nACGTTGCAACGTTGCA\n+\nHHHHAAAAHHHHAAAA\n\
                     @read5\nGGGGCCCCAAAATTTT\n+\nBBBBBBBBBBBBBBBB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_compressor: QualityCompressor::Fqzcomp,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len(),
        "Line count mismatch: expected {}, got {}", original_lines.len(), decompressed_lines.len());

    for (i, (orig, decomp)) in original_lines.iter().zip(decompressed_lines.iter()).enumerate() {
        assert_eq!(orig.trim(), decomp.trim(),
            "Line {} mismatch:\n  original:     {:?}\n  decompressed: {:?}", i, orig, decomp);
    }
}

#[test]
fn test_openzl_mixed_compressors() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_compressor: QualityCompressor::OpenZl,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::OpenZl,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len());

    for (orig, decomp) in original_lines.iter().zip(decompressed_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_reorder_local_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
@read4\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        reorder: Some(ReorderMode::Local),
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    // With reorder, reads may be in a different order but all reads must be present
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());

    let mut orig_records: Vec<String> = orig_lines.chunks(4).map(|c| c.join("\n")).collect();
    let mut decomp_records: Vec<String> = decomp_lines.chunks(4).map(|c| c.join("\n")).collect();
    orig_records.sort();
    decomp_records.sort();
    assert_eq!(orig_records, decomp_records);
}

#[test]
fn test_reorder_global_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
@read4\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        reorder: Some(ReorderMode::Global),
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path.clone())).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());

    let mut orig_records: Vec<String> = orig_lines.chunks(4).map(|c| c.join("\n")).collect();
    let mut decomp_records: Vec<String> = decomp_lines.chunks(4).map(|c| c.join("\n")).collect();
    orig_records.sort();
    decomp_records.sort();
    assert_eq!(orig_records, decomp_records);
}

#[test]
fn test_sequence_hints_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Reads >= 21bp for valid syncmers, plus one short read to test 0xFF fallback
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n\
@read4\nACGT\n+\nIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        sequence_hints: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 4
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[0], 4, "encoding_type should be 4 for sequence hints");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());
    for (orig, decomp) in orig_lines.iter().zip(decomp_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_sequence_delta_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Include near-duplicate reads to exercise delta path, plus a unique read
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nACGTACGTACGTACGTACGTACGTATGTACGT\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n\
@read4\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read5\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n+\nDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n\
@read6\nACGT\n+\nIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        sequence_delta: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 5
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[0], 5, "encoding_type should be 5 for inline delta");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());
    for (orig, decomp) in orig_lines.iter().zip(decomp_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_rc_canon_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("test.fastq");

    // Create test data with reads that benefit from RC canonicalization
    fs::write(&input_fastq,
        "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
         @read2\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n\
         @read3\nTTTTAAAACCCCGGGG\n+\nAAAAAAAAAAAAAAAA\n\
         @read4\nCCCCGGGGTTTTAAAA\n+\nBBBBBBBBBBBBBBBB\n"
    ).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        rc_canon: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 6
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[0], 6, "encoding_type should be 6 for rc_canon");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (orig, decomp)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(orig.trim(), decomp.trim(), "mismatch at line {}", i);
    }
}

#[test]
fn test_factorize_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Include near-duplicate reads to exercise the factorize match path
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nACGTACGTACGTACGTACGTACGTATGTACGT\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n\
@read4\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read5\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n+\nDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n\
@read6\nACGT\n+\nIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        factorize: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 7
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[0], 7, "encoding_type should be 7 for factorize");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (orig, decomp)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(orig.trim(), decomp.trim(), "mismatch at line {}", i);
    }
}

#[test]
fn test_local_reorder_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read3\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read4\nGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read5\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read6\nTGCATGCATGCATGCATGCATGCATGCATGCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        local_reorder: true,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[0], 8, "encoding_type should be 8 for local-reorder");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();
    assert!(output_fastq.exists());

    // Local reorder restores original order — verify exact match
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}

#[test]
fn test_ultra_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read4\nGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n\
@read5\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read6\nTGCATGCATGCATGCATGCATGCATGCATGCT\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ultra: Some(3),
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[0], 9, "encoding_type should be 9 for ultra");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();
    assert!(output_fastq.exists());

    // Ultra preserves original order — verify exact match
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}

#[test]
fn test_quality_ctx_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Generate 20 fixed-length reads (quality_ctx requires fixed length)
    let mut test_data = String::new();
    let seq_base = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp
    let qual_base = "IIHHJJBBCCDDEEFFIIHHJJBBCCDDEEFF"; // 32 chars
    assert_eq!(seq_base.len(), 32);
    assert_eq!(qual_base.len(), 32);
    for i in 0..20 {
        // Rotate sequence to get variety
        let rot = i % seq_base.len();
        let seq: String = seq_base.bytes().cycle().skip(rot).take(32).map(|b| b as char).collect();
        test_data.push_str(&format!("@read{}\n{}\n+\n{}\n", i, seq, qual_base));
    }
    fs::write(&input_fastq, &test_data).unwrap();

    // Force quality_ctx (bypasses the 100K-read auto-selection threshold)
    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_compressor: QualityCompressor::QualityCtx,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify archive header has quality_compressor = 4 (QualityCtx)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[3], 4, "quality_compressor code should be 4 (QualityCtx)");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}

#[test]
fn test_quality_ctx_variable_length_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Variable-length reads: 20bp, 32bp, 15bp, 40bp, 32bp
    let test_data = "\
@read0\nACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFIIHH\n\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFIIHHJJBBCCDDEEFF\n\
@read2\nACGTACGTACGTACG\n+\nIIHHJJBBCCDDEEF\n\
@read3\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFIIHHJJBBCCDDEEFFIIHHJJBB\n\
@read4\nTTTGGGCCCAAAGGGTTTGGGCCCAAAGGGTT\n+\nBBBBAAAABBBBAAAABBBBAAAABBBBAAAABB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_compressor: QualityCompressor::QualityCtx,
        ..CompressArgs::default()
    };

    qz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify archive uses QualityCtx (not fallen back to BSC)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[3], 4, "quality_compressor should be 4 (QualityCtx), not BSC fallback");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();
    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}
