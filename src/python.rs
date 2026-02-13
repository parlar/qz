use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::path::PathBuf;

use crate::cli::{
    CompressArgs, DecompressArgs, HeaderCompressor, QualityCompressor, QualityMode,
    SequenceCompressor,
};
use crate::compression;

/// Compress a FASTQ file to QZ format
///
/// Args:
///     input_file: Path to input FASTQ file (str or list of str for paired-end)
///     output_file: Path to output QZ archive
///     working_dir: Working directory for temporary files (default: ".")
///     num_threads: Number of threads to use (0 = auto, default: auto)
///     no_quality: Do not preserve quality values (default: False)
///     quality_mode: Quality mode - "lossless", "illumina-bin", "illumina4", "binary", "qvz", "discard" (default: "lossless")
///     gzipped: Input files are gzipped (default: False)
///     fasta: Input is FASTA format (default: False)
///     quality_compressor: Quality compressor - "bsc", "zstd", or "openzl" (default: "bsc")
///     sequence_compressor: Sequence compressor - "bsc", "zstd", or "openzl" (default: "bsc")
///     header_compressor: Header compressor - "bsc", "zstd", or "openzl" (default: "bsc")
///     chunked: Use chunked streaming mode for lower memory (default: False)
///     bsc_static: Use BSC static coder instead of adaptive (default: False)
///     arithmetic: Enable arithmetic coding (experimental, default: False)
///     debruijn: Enable de Bruijn graph compression (experimental, default: False)
///
/// Returns:
///     None
///
/// Example:
///     >>> import qz
///     >>> qz.compress("reads.fastq", "compressed.qz")
///     >>> qz.compress(["R1.fastq", "R2.fastq"], "paired.qz")
#[pyfunction]
#[pyo3(signature = (
    input_file,
    output_file,
    working_dir = None,
    num_threads = None,
    no_quality = false,
    quality_mode = "lossless",
    gzipped = false,
    fasta = false,
    quality_compressor = "bsc",
    sequence_compressor = "bsc",
    header_compressor = "bsc",
    chunked = false,
    bsc_static = false,
    arithmetic = false,
    debruijn = false
))]
fn compress(
    input_file: PathOrList,
    output_file: String,
    working_dir: Option<String>,
    num_threads: Option<usize>,
    no_quality: bool,
    quality_mode: &str,
    gzipped: bool,
    fasta: bool,
    quality_compressor: &str,
    sequence_compressor: &str,
    header_compressor: &str,
    chunked: bool,
    bsc_static: bool,
    arithmetic: bool,
    debruijn: bool,
) -> PyResult<()> {
    let input_paths = match input_file {
        PathOrList::Single(path) => vec![PathBuf::from(path)],
        PathOrList::List(paths) => paths.into_iter().map(PathBuf::from).collect(),
    };

    let qmode = parse_quality_mode(quality_mode)?;
    let qcomp = parse_compressor::<QualityCompressor>(quality_compressor, "quality_compressor")?;
    let scomp = parse_compressor::<SequenceCompressor>(sequence_compressor, "sequence_compressor")?;
    let hcomp = parse_compressor::<HeaderCompressor>(header_compressor, "header_compressor")?;

    let args = CompressArgs {
        input: input_paths,
        output: PathBuf::from(output_file),
        working_dir: PathBuf::from(working_dir.unwrap_or_else(|| ".".to_string())),
        threads: num_threads.unwrap_or(0),
        no_quality,
        gzipped,
        fasta,
        quality_mode: qmode,
        quality_compressor: qcomp,
        sequence_compressor: scomp,
        header_compressor: hcomp,
        chunked,
        bsc_static,
        arithmetic,
        debruijn,
        ..CompressArgs::default()
    };

    compression::compress(&args).map_err(|e| PyRuntimeError::new_err(format!("{:?}", e)))
}

/// Decompress a QZ archive to FASTQ format
///
/// Args:
///     input_file: Path to input QZ archive
///     output_file: Path to output FASTQ file (str or list of str for paired-end)
///     working_dir: Working directory for temporary files (default: ".")
///     num_threads: Number of threads to use (default: auto)
///     gzipped: Output gzipped FASTQ (default: False)
///     gzip_level: Gzip compression level 0-9 (default: 6)
///
/// Returns:
///     None
///
/// Example:
///     >>> import qz
///     >>> qz.decompress("compressed.qz", "reads.fastq")
///     >>> qz.decompress("compressed.qz", "reads.fastq.gz", gzipped=True)
#[pyfunction]
#[pyo3(signature = (
    input_file,
    output_file,
    working_dir = None,
    num_threads = None,
    gzipped = false,
    gzip_level = 6
))]
fn decompress(
    input_file: String,
    output_file: PathOrList,
    working_dir: Option<String>,
    num_threads: Option<usize>,
    gzipped: bool,
    gzip_level: u32,
) -> PyResult<()> {
    let output_paths = match output_file {
        PathOrList::Single(path) => vec![PathBuf::from(path)],
        PathOrList::List(paths) => paths.into_iter().map(PathBuf::from).collect(),
    };

    let args = DecompressArgs {
        input: PathBuf::from(input_file),
        output: output_paths,
        working_dir: PathBuf::from(working_dir.unwrap_or_else(|| ".".to_string())),
        num_threads: num_threads.unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(8)
        }),
        gzipped,
        gzip_level,
    };

    compression::decompress(&args).map_err(|e| PyRuntimeError::new_err(format!("{:?}", e)))
}

/// Get version information
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Helper enum to accept either a string or list of strings
#[derive(FromPyObject)]
enum PathOrList {
    Single(String),
    List(Vec<String>),
}

fn parse_quality_mode(mode: &str) -> PyResult<QualityMode> {
    match mode {
        "lossless" => Ok(QualityMode::Lossless),
        "illumina-bin" => Ok(QualityMode::IlluminaBin),
        "illumina4" => Ok(QualityMode::Illumina4),
        "binary" => Ok(QualityMode::Binary),
        "qvz" => Ok(QualityMode::Qvz),
        "discard" => Ok(QualityMode::Discard),
        _ => Err(PyRuntimeError::new_err(format!(
            "Invalid quality mode: '{}'. Must be one of: lossless, illumina-bin, illumina4, binary, qvz, discard",
            mode
        ))),
    }
}

fn parse_compressor<T: FromStrCompressor>(s: &str, param_name: &str) -> PyResult<T> {
    T::from_str_py(s, param_name)
}

trait FromStrCompressor: Sized {
    fn from_str_py(s: &str, param_name: &str) -> PyResult<Self>;
}

impl FromStrCompressor for QualityCompressor {
    fn from_str_py(s: &str, param_name: &str) -> PyResult<Self> {
        match s {
            "bsc" => Ok(QualityCompressor::Bsc),
            "zstd" => Ok(QualityCompressor::Zstd),
            "openzl" => Ok(QualityCompressor::OpenZl),
            _ => Err(PyRuntimeError::new_err(format!(
                "Invalid {}: '{}'. Must be 'bsc', 'zstd', or 'openzl'", param_name, s
            ))),
        }
    }
}

impl FromStrCompressor for SequenceCompressor {
    fn from_str_py(s: &str, param_name: &str) -> PyResult<Self> {
        match s {
            "bsc" => Ok(SequenceCompressor::Bsc),
            "zstd" => Ok(SequenceCompressor::Zstd),
            "openzl" => Ok(SequenceCompressor::OpenZl),
            _ => Err(PyRuntimeError::new_err(format!(
                "Invalid {}: '{}'. Must be 'bsc', 'zstd', or 'openzl'", param_name, s
            ))),
        }
    }
}

impl FromStrCompressor for HeaderCompressor {
    fn from_str_py(s: &str, param_name: &str) -> PyResult<Self> {
        match s {
            "bsc" => Ok(HeaderCompressor::Bsc),
            "zstd" => Ok(HeaderCompressor::Zstd),
            "openzl" => Ok(HeaderCompressor::OpenZl),
            _ => Err(PyRuntimeError::new_err(format!(
                "Invalid {}: '{}'. Must be 'bsc', 'zstd', or 'openzl'", param_name, s
            ))),
        }
    }
}

/// Python module definition
#[pymodule]
fn qz(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compress, m)?)?;
    m.add_function(wrap_pyfunction!(decompress, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;

    // Add version as module attribute
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
