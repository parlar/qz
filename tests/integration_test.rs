use fqz::cli::{CompressArgs, DecompressArgs, HeaderCompressor, QualityCompressor, QualityMode, ReorderMode, SequenceCompressor};
use std::fs;
use tempfile::TempDir;

#[test]
fn test_compress_decompress_roundtrip() {
    // Create temporary directory for test files
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create a small test FASTQ file
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    // Compress the file
    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    fqz::compression::compress(&compress_args).unwrap();

    // Verify archive was created
    assert!(archive_path.exists());

    // Decompress the archive
    let output_fastq = temp_path.join("decompressed.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path.clone(),
        output: vec![output_fastq.clone()],
        working_dir: temp_path.clone(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };

    fqz::compression::decompress(&decompress_args).unwrap();

    // Verify decompressed file exists
    assert!(output_fastq.exists());

    // Compare original and decompressed content
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    // Parse both files into records for comparison
    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len());

    // Compare each record (allowing for whitespace differences)
    for (orig, decomp) in original_lines.iter().zip(decompressed_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_decompress_gzipped_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create test FASTQ
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    // Compress
    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };
    fqz::compression::compress(&compress_args).unwrap();

    // Decompress to gzipped output
    let output_fastq = temp_path.join("decompressed.fastq.gz");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: true,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();

    // Verify gzipped file was created
    assert!(output_fastq.exists());

    // Verify it's a valid gzip file by checking magic bytes
    let data = fs::read(&output_fastq).unwrap();
    assert!(data.len() >= 2);
    assert_eq!(data[0], 0x1f);
    assert_eq!(data[1], 0x8b);
}

#[test]
fn test_decompress_no_quality() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create test FASTQ
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    // Compress without quality — use Zstd since BSC needs minimum data size
    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: true,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Zstd,
        sequence_compressor: SequenceCompressor::Zstd,
        header_compressor: HeaderCompressor::Zstd,
        bsc_static: false,
        chunked: false,
    };
    fqz::compression::compress(&compress_args).unwrap();

    // Decompress
    let output_fastq = temp_path.join("decompressed.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();

    // Verify decompressed file has dummy quality values
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

    // Create empty FASTQ file
    let input_fastq = temp_path.join("empty.fastq");
    fs::write(&input_fastq, "").unwrap();

    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    // Should handle empty file gracefully (either error or create empty archive)
    let result = fqz::compression::compress(&compress_args);
    // Accept either success with empty archive or specific error
    if result.is_ok() {
        assert!(archive_path.exists());
    }
}

#[test]
fn test_single_read() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Single read
    let input_fastq = temp_path.join("single.fastq");
    fs::write(&input_fastq, "@read1\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    fqz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Decompress and verify
    let output_fastq = temp_path.join("output.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("read1"));
    assert!(content.contains("ACGT"));
}

#[test]
fn test_long_read() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create a very long read (1000bp)
    let long_seq = "ACGT".repeat(250);
    let long_qual = "I".repeat(1000);
    let fastq_data = format!("@long_read\n{}\n+\n{}\n", long_seq, long_qual);
    
    let input_fastq = temp_path.join("long.fastq");
    fs::write(&input_fastq, fastq_data).unwrap();

    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: true,  // Enable long read mode
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    fqz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Decompress and verify
    let output_fastq = temp_path.join("output.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains(&long_seq));
}

#[test]
fn test_reads_with_n_bases() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Reads with N bases
    let test_data = "@read1\nACGTNNNNACGT\n+\nIIII####IIII\n@read2\nNNNNNNNN\n+\n########\n";
    let input_fastq = temp_path.join("with_n.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    fqz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Decompress and verify N bases preserved
    let output_fastq = temp_path.join("output.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();
    
    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("ACGTNNNNA"));  // N bases should be preserved
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

    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::IlluminaBin,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    fqz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Decompress
    let output_fastq = temp_path.join("output.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();
    
    // Verify file exists (quality may differ due to binning)
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

    let archive_path = temp_path.join("test.fqz");
    let compress_args = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Binary,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };

    fqz::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Decompress
    let output_fastq = temp_path.join("output.fastq");
    let decompress_args = DecompressArgs {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args).unwrap();
    
    // Verify file exists
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
        let archive_path = temp_path.join(format!("test_level_{}.fqz", level));
        let compress_args = CompressArgs {
            input: vec![input_fastq.clone()],
            output: archive_path.clone(),
            working_dir: temp_path.clone(),
            num_threads: 1,
            allow_reordering: false,
            no_quality: false,
            no_ids: false,
            long_mode: false,
            gzipped: false,
            fasta: false,
            quality_mode: QualityMode::Lossless,
            delta_encoding: false,
            rle_encoding: false,
            quality_modeling: false,
            quality_delta: false,
            dict_training: false,
            dict_size: 64,
            compression_level: level,
            threads: 1,
            use_zstd: false,
            reorder_by: ReorderMode::None,
            arithmetic: false,
            kmer_reference: false,
            debruijn: false,
            kmer_size: 20,
            quality_compressor: QualityCompressor::Bsc,
            sequence_compressor: SequenceCompressor::Bsc,
            header_compressor: HeaderCompressor::Bsc,
            bsc_static: false,
            chunked: false,
        };

        fqz::compression::compress(&compress_args).unwrap();
        let size = fs::metadata(&archive_path).unwrap().len();
        sizes.push((level, size));
    }

    // Verify all succeeded
    assert_eq!(sizes.len(), 4);
    
    // Generally higher levels should produce smaller or equal files
    // (though not guaranteed on small files)
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
    let archive1 = temp_path.join("test1.fqz");
    let compress_args1 = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive1.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 1,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };
    fqz::compression::compress(&compress_args1).unwrap();

    // Compress with 4 threads
    let archive2 = temp_path.join("test2.fqz");
    let compress_args2 = CompressArgs {
        input: vec![input_fastq.clone()],
        output: archive2.clone(),
        working_dir: temp_path.clone(),
        num_threads: 1,
        allow_reordering: false,
        no_quality: false,
        no_ids: false,
        long_mode: false,
        gzipped: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        delta_encoding: false,
        rle_encoding: false,
        quality_modeling: false,
        quality_delta: false,
        dict_training: false,
        dict_size: 64,
        compression_level: 3,
        threads: 4,
        use_zstd: false,
        reorder_by: ReorderMode::None,
        arithmetic: false,
        kmer_reference: false,
        debruijn: false,
        kmer_size: 20,
        quality_compressor: QualityCompressor::Bsc,
        sequence_compressor: SequenceCompressor::Bsc,
        header_compressor: HeaderCompressor::Bsc,
        bsc_static: false,
        chunked: false,
    };
    fqz::compression::compress(&compress_args2).unwrap();

    // Decompress both
    let output1 = temp_path.join("output1.fastq");
    let output2 = temp_path.join("output2.fastq");

    let decompress_args1 = DecompressArgs {
        input: archive1,
        output: vec![output1.clone()],
        working_dir: temp_path.clone(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args1).unwrap();

    let decompress_args2 = DecompressArgs {
        input: archive2,
        output: vec![output2.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        range: None,
    };
    fqz::compression::decompress(&decompress_args2).unwrap();

    // Both should produce identical output
    let content1 = fs::read_to_string(&output1).unwrap();
    let content2 = fs::read_to_string(&output2).unwrap();
    assert_eq!(content1, content2, "Multi-threaded compression should produce same output as single-threaded");
}
