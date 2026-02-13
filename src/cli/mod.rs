use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "qz")]
#[command(author = "QZ Contributors")]
#[command(version = "0.1.0")]
#[command(about = "Columnar FASTQ compression using BSC/BWT encoding", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Compress FASTQ files
    Compress(CompressArgs),
    /// Decompress QZ archives
    Decompress(DecompressArgs),
}

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum ReorderMode {
    /// Sort within each 5M-record chunk (fast, bounded memory, local ordering only)
    Local,
    /// Two-pass bucket sort across entire file (better ordering, bounded memory)
    Global,
}

#[derive(Parser)]
pub struct CompressArgs {
    /// Input FASTQ file(s) (one for single-end, two for paired-end)
    #[arg(short, long, value_name = "FILE", required = true)]
    pub input: Vec<PathBuf>,

    /// Output QZ archive file
    #[arg(short, long, value_name = "FILE", required = true)]
    pub output: PathBuf,

    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    pub working_dir: PathBuf,

    /// Number of threads (0 = auto-detect, default: auto)
    #[arg(short = 't', long, default_value = "0")]
    pub threads: usize,

    /// Do not preserve quality values
    #[arg(long)]
    pub no_quality: bool,

    /// Input files are gzipped FASTQ
    #[arg(short, long)]
    pub gzipped: bool,

    /// Input is FASTA format (no quality scores)
    #[arg(long)]
    pub fasta: bool,

    /// Quality compression mode
    #[arg(short, long, value_enum, default_value = "lossless")]
    pub quality_mode: QualityMode,

    /// Enable delta encoding for sequences (experimental)
    #[arg(long)]
    pub delta_encoding: bool,

    /// Enable run-length encoding for homopolymers (experimental)
    #[arg(long)]
    pub rle_encoding: bool,

    /// Enable positional quality modeling (experimental)
    #[arg(long)]
    pub quality_modeling: bool,

    /// Enable quality delta encoding between adjacent reads (experimental)
    #[arg(long)]
    pub quality_delta: bool,

    /// Enable zstd dictionary training (only used with --quality-compressor zstd).
    /// WARNING: Only beneficial for very large datasets (>1M reads) or multi-file scenarios.
    /// The dictionary overhead (64KB+) outweighs gains for typical single-file compression.
    #[arg(long, hide = true)]
    pub dict_training: bool,

    /// Dictionary size in KB (only with --dict-training and --quality-compressor zstd)
    #[arg(long, default_value = "64", hide = true)]
    pub dict_size: usize,

    /// Compression level for legacy zstd mode (1-22). Ignored when using BSC (default).
    #[arg(long, short = 'l', default_value = "3", hide = true)]
    pub compression_level: i32,

    /// Enable arithmetic coding for sequences and qualities (experimental, 6-8x compression)
    #[arg(long)]
    pub arithmetic: bool,

    /// Enable de Bruijn graph compression for sequences (patent-safe, RC-aware)
    #[arg(long)]
    pub debruijn: bool,

    /// K-mer size for de Bruijn / reference compression (0 = auto-select based on depth, range: 9-31)
    #[arg(long, default_value = "0")]
    pub kmer_size: usize,

    /// Quality score compressor
    #[arg(long, value_enum, default_value = "bsc")]
    pub quality_compressor: QualityCompressor,

    /// Sequence compressor
    #[arg(long, value_enum, default_value = "bsc")]
    pub sequence_compressor: SequenceCompressor,

    /// Header compressor
    #[arg(long, value_enum, default_value = "bsc")]
    pub header_compressor: HeaderCompressor,

    /// Use BSC static coder instead of adaptive (for benchmarking: slightly worse compression, slightly faster)
    #[arg(long)]
    pub bsc_static: bool,

    /// Use chunked streaming mode for lower memory usage (reads in 5M-record chunks with pipelined I/O)
    #[arg(long)]
    pub chunked: bool,

    /// Use 2-bit sequence encoding (4 bases per byte + N-mask bitmap, may improve compression ratio)
    #[arg(long)]
    pub twobit: bool,

    /// Use template-based header encoding (extracts common prefix, delta-encodes coordinates)
    #[arg(long)]
    pub header_template: bool,

    /// Reorder reads by content similarity for better compression (destroys original read order).
    /// Modes: 'local' = sort within 5M-record chunks (fast, bounded memory);
    /// 'global' = two-pass bucket sort across entire file (better ordering, bounded memory).
    #[arg(long, value_enum)]
    pub reorder: Option<ReorderMode>,

    /// Reverse-complement canonicalization: store each read as the lexicographically smaller
    /// of {read, revcomp(read)} plus a 1-bit strand flag. Doubles effective k-mer overlap
    /// for BWT compression without reordering.
    #[arg(long)]
    pub rc_canon: bool,

    /// Prepend a syncmer-derived hint byte before each read's sequence to aid BWT clustering.
    /// Preserves read order while giving BSC's BWT a grouping signal for similar reads.
    #[arg(long)]
    pub sequence_hints: bool,

    /// Inline delta encoding: encode reads as deltas against cached similar reads.
    /// Keeps everything in one BSC block; matching bases become 0x00 for better BWT compression.
    #[arg(long)]
    pub sequence_delta: bool,

    /// Factorized sequence compression: two-pass inverted-index delta encoding
    /// with separate metadata stream. Uses syncmer-based matching within each chunk
    /// to find similar reads and encode them as deltas, without reordering.
    #[arg(long)]
    pub factorize: bool,
}

impl Default for CompressArgs {
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
            chunked: false,
            twobit: false,
            header_template: false,
            rc_canon: false,
            reorder: None,
            sequence_hints: false,
            sequence_delta: false,
            factorize: false,
        }
    }
}

#[derive(Parser)]
pub struct DecompressArgs {
    /// Input QZ archive
    #[arg(short, long, value_name = "FILE", required = true)]
    pub input: PathBuf,

    /// Output FASTQ file(s)
    #[arg(short, long, value_name = "FILE", required = true)]
    pub output: Vec<PathBuf>,

    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    pub working_dir: PathBuf,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = num_cpus())]
    pub num_threads: usize,

    /// Output gzipped FASTQ
    #[arg(short, long)]
    pub gzipped: bool,

    /// Gzip compression level (0-9)
    #[arg(long, default_value = "6")]
    pub gzip_level: u32,

}

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
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

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum QualityCompressor {
    /// Legacy zstd (faster, ~20% larger)
    Zstd,
    /// BSC/BWT (best compression, default)
    Bsc,
    /// OpenZL format-aware compression (Meta)
    #[value(name = "openzl")]
    OpenZl,
    /// fqzcomp context-modeled compression with mean-quality reordering (~8% smaller than BSC)
    Fqzcomp,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum SequenceCompressor {
    /// Legacy zstd with 2-bit encoding + N-mask (faster, ~9% larger)
    Zstd,
    /// BSC on raw ASCII sequences (best compression, default)
    Bsc,
    /// OpenZL format-aware compression (Meta)
    #[value(name = "openzl")]
    OpenZl,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum HeaderCompressor {
    /// Legacy template + zstd (~2x larger)
    Zstd,
    /// BSC on raw headers (best compression, default)
    Bsc,
    /// OpenZL format-aware compression (Meta)
    #[value(name = "openzl")]
    OpenZl,
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

