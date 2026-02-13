use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "fqz")]
#[command(author = "FQZ Contributors")]
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
    /// Decompress FQZ archives
    Decompress(DecompressArgs),
    /// Benchmark all compression strategies on each FASTQ stream (headers, sequences, qualities)
    Benchmark(BenchmarkArgs),
}

#[derive(Parser)]
pub struct CompressArgs {
    /// Input FASTQ file(s) (one for single-end, two for paired-end)
    #[arg(short, long, value_name = "FILE", required = true)]
    pub input: Vec<PathBuf>,

    /// Output FQZ archive file
    #[arg(short, long, value_name = "FILE", required = true)]
    pub output: PathBuf,

    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    pub working_dir: PathBuf,

    /// Number of threads (default: number of CPU cores)
    #[arg(short = 't', long, default_value_t = num_cpus())]
    pub num_threads: usize,

    /// Enable read reordering for better compression (experimental, may have patent implications)
    #[arg(long)]
    pub allow_reordering: bool,

    /// Patent-safe read reordering strategy for better compression.
    /// Flowcell: Sort by tile/X/Y coordinates (best for quality correlation).
    /// GC: Sort by GC content bins. Length: Sort by read length.
    /// Lexicographic: Alphabetical sort. Smart: Multi-level metadata sort.
    #[arg(long, value_enum, default_value = "none")]
    pub reorder_by: ReorderMode,

    /// Do not preserve quality values
    #[arg(long)]
    pub no_quality: bool,

    /// Do not preserve read IDs
    #[arg(long)]
    pub no_ids: bool,

    /// Enable long read mode (for reads > 511 bases)
    #[arg(long)]
    pub long_mode: bool,

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

    /// DEPRECATED: BSC is now the default compressor (kept for backward compatibility)
    #[arg(long, hide = true)]
    pub use_zstd: bool,

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

    /// Number of threads to use (0 = auto-detect). Parallelizes stream compression.
    #[arg(long, short = 'j', default_value = "0")]
    pub threads: usize,

    /// Enable arithmetic coding for sequences and qualities (experimental, 6-8x compression)
    #[arg(long)]
    pub arithmetic: bool,

    /// Enable k-mer reference compression for sequences (6-8x compression, patent-safe)
    #[arg(long)]
    pub kmer_reference: bool,

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
}

#[derive(Parser)]
pub struct DecompressArgs {
    /// Input FQZ archive
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

    /// Decompress only reads in range [start, end]
    #[arg(long, value_names = ["START", "END"], num_args = 2)]
    pub range: Option<Vec<u64>>,
}

#[derive(Parser)]
pub struct BenchmarkArgs {
    /// Input FASTQ file to benchmark against
    #[arg(short, long, value_name = "FILE", required = true)]
    pub input: PathBuf,

    /// Input file is gzipped FASTQ
    #[arg(short, long)]
    pub gzipped: bool,

    /// Maximum number of reads to load (default: all)
    #[arg(long)]
    pub max_reads: Option<usize>,

    /// Only benchmark a specific stream (headers, sequences, qualities)
    #[arg(long)]
    pub stream: Option<String>,
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
pub enum ReorderMode {
    /// No reordering (preserve input order)
    None,
    /// Sort by flowcell tile/X/Y coordinates (patent-safe, best for quality correlation)
    Flowcell,
    /// Sort by GC content percentage (patent-safe)
    Gc,
    /// Sort by read length (patent-safe)
    Length,
    /// Sort alphabetically by sequence (patent-safe)
    Lexicographic,
    /// Multi-level metadata sort: flowcell + GC + length (patent-safe)
    Smart,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum QualityCompressor {
    /// Legacy zstd (faster, ~20% larger)
    Zstd,
    /// BSC/BWT (best compression, default)
    Bsc,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum SequenceCompressor {
    /// Legacy zstd with 2-bit encoding + N-mask (faster, ~9% larger)
    Zstd,
    /// BSC on raw ASCII sequences (best compression, default)
    Bsc,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum HeaderCompressor {
    /// Legacy template + zstd (~2x larger)
    Zstd,
    /// BSC on raw headers (best compression, default)
    Bsc,
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

