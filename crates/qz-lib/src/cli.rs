use std::path::PathBuf;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReorderMode {
    /// Sort within each 5M-record chunk (fast, bounded memory, local ordering only)
    Local,
    /// Two-pass bucket sort across entire file (better ordering, bounded memory)
    Global,
}

#[derive(Clone)]
pub struct CompressConfig {
    /// Input FASTQ file(s) (one for single-end, two for paired-end)
    pub input: Vec<PathBuf>,
    /// Output QZ archive file
    pub output: PathBuf,
    /// Working directory for temporary files
    pub working_dir: PathBuf,
    /// Number of threads (0 = auto-detect)
    pub threads: usize,
    /// Do not preserve quality values
    pub no_quality: bool,
    /// Input files are gzipped FASTQ
    pub gzipped: bool,
    /// Input is FASTA format (no quality scores)
    pub fasta: bool,
    /// Quality compression mode
    pub quality_mode: QualityMode,
    /// Enable delta encoding for sequences (experimental)
    pub delta_encoding: bool,
    /// Enable run-length encoding for homopolymers (experimental)
    pub rle_encoding: bool,
    /// Enable positional quality modeling (experimental)
    pub quality_modeling: bool,
    /// Enable quality delta encoding between adjacent reads (experimental)
    pub quality_delta: bool,
    /// Enable zstd dictionary training
    pub dict_training: bool,
    /// Dictionary size in KB
    pub dict_size: usize,
    /// Compression level for legacy zstd mode (1-22)
    pub compression_level: i32,
    /// Enable arithmetic coding for sequences and qualities (experimental)
    pub arithmetic: bool,
    /// Enable de Bruijn graph compression for sequences
    pub debruijn: bool,
    /// K-mer size for de Bruijn / reference compression (0 = auto)
    pub kmer_size: usize,
    /// Quality score compressor
    pub quality_compressor: QualityCompressor,
    /// Sequence compressor
    pub sequence_compressor: SequenceCompressor,
    /// Header compressor
    pub header_compressor: HeaderCompressor,
    /// Use BSC static coder instead of adaptive
    pub bsc_static: bool,
    /// Use 2-bit sequence encoding (4 bases per byte + N-mask bitmap)
    pub twobit: bool,
    /// Use template-based header encoding
    pub header_template: bool,
    /// Reorder reads by content similarity
    pub reorder: Option<ReorderMode>,
    /// Reverse-complement canonicalization
    pub rc_canon: bool,
    /// Prepend a syncmer-derived hint byte before each read's sequence
    pub sequence_hints: bool,
    /// Inline delta encoding against cached similar reads
    pub sequence_delta: bool,
    /// Two-pass pattern routing for BWT clustering
    pub factorize: bool,
    /// Local reordering with delta encoding
    pub local_reorder: bool,
    /// Ultra compression with optional level (1-5, 0=auto)
    pub ultra: Option<u8>,
    /// Deprecated: use ultra level 2 instead
    pub fast_ultra: bool,
    /// Number of reads per quality_ctx sub-block (default 500K)
    pub quality_ctx_block_size: usize,
}

impl Default for CompressConfig {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            output: PathBuf::new(),
            working_dir: PathBuf::from("."),
            threads: 0,
            no_quality: false,
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
            arithmetic: false,
            debruijn: false,
            kmer_size: 0,
            quality_compressor: QualityCompressor::Bsc,
            sequence_compressor: SequenceCompressor::Bsc,
            header_compressor: HeaderCompressor::Bsc,
            bsc_static: false,
            twobit: false,
            header_template: false,
            rc_canon: false,
            reorder: None,
            sequence_hints: false,
            sequence_delta: false,
            factorize: false,
            local_reorder: false,
            ultra: None,
            fast_ultra: false,
            quality_ctx_block_size: 500_000,
        }
    }
}

#[derive(Clone)]
pub struct DecompressConfig {
    /// Input QZ archive
    pub input: PathBuf,
    /// Output FASTQ file(s)
    pub output: Vec<PathBuf>,
    /// Working directory for temporary files
    pub working_dir: PathBuf,
    /// Number of threads
    pub num_threads: usize,
    /// Output gzipped FASTQ
    pub gzipped: bool,
    /// Gzip compression level (0-9)
    pub gzip_level: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QualityMode {
    /// Lossless quality preservation
    Lossless,
    /// Illumina 8-level binning
    IlluminaBin,
    /// Illumina 4-level binning (more aggressive)
    Illumina4,
    /// Binary thresholding
    Binary,
    /// QVZ lossy compression
    Qvz,
    /// Discard quality scores entirely
    Discard,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QualityCompressor {
    /// Legacy zstd (faster, ~20% larger)
    Zstd,
    /// BSC/BWT (best compression, default)
    Bsc,
    /// OpenZL format-aware compression
    OpenZl,
    /// fqzcomp context-modeled compression
    Fqzcomp,
    /// Context-adaptive range coding
    QualityCtx,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequenceCompressor {
    /// Legacy zstd with 2-bit encoding + N-mask
    Zstd,
    /// BSC on raw ASCII sequences (best compression, default)
    Bsc,
    /// OpenZL format-aware compression
    OpenZl,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeaderCompressor {
    /// Legacy template + zstd
    Zstd,
    /// BSC on raw headers (best compression, default)
    Bsc,
    /// OpenZL format-aware compression
    OpenZl,
}

pub fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}
