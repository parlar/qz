use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::path::PathBuf;

use crate::cli::{CompressArgs, DecompressArgs, QualityMode};
use crate::compression;

/// Compress a FASTQ file to COIL format
///
/// Args:
///     input_file: Path to input FASTQ file (or list of 2 files for paired-end)
///     output_file: Path to output COIL archive
///     working_dir: Working directory for temporary files (default: ".")
///     num_threads: Number of threads to use (default: auto)
///     no_quality: Do not preserve quality values (default: False)
///     no_ids: Do not preserve read IDs (default: False)
///     quality_mode: Quality compression mode - "lossless", "illumina-bin", or "binary" (default: "lossless")
///     gzipped: Input files are gzipped (default: False)
///
/// Returns:
///     None
///
/// Example:
///     >>> import coil
///     >>> coil.compress("reads.fastq", "compressed.coil")
#[pyfunction]
#[pyo3(signature = (
    input_file,
    output_file,
    working_dir = None,
    num_threads = None,
    no_quality = false,
    no_ids = false,
    quality_mode = "lossless",
    gzipped = false
))]
fn compress(
    input_file: PathOrList,
    output_file: String,
    working_dir: Option<String>,
    num_threads: Option<usize>,
    no_quality: bool,
    no_ids: bool,
    quality_mode: &str,
    gzipped: bool,
) -> PyResult<()> {
    // Parse input files
    let input_paths = match input_file {
        PathOrList::Single(path) => vec![PathBuf::from(path)],
        PathOrList::List(paths) => paths.into_iter().map(PathBuf::from).collect(),
    };

    // Parse quality mode
    let qmode = match quality_mode {
        "lossless" => QualityMode::Lossless,
        "illumina-bin" => QualityMode::IlluminaBin,
        "binary" => QualityMode::Binary,
        "qvz" => QualityMode::Qvz,
        _ => {
            return Err(PyRuntimeError::new_err(format!(
                "Invalid quality mode: {}. Must be one of: lossless, illumina-bin, binary, qvz",
                quality_mode
            )))
        }
    };

    let args = CompressArgs {
        input: input_paths,
        output: PathBuf::from(output_file),
        working_dir: PathBuf::from(working_dir.unwrap_or_else(|| ".".to_string())),
        num_threads: num_threads.unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(8)
        }),
        allow_reordering: false,  // Not exposed in API for patent safety
        no_quality,
        no_ids,
        long_mode: false,
        gzipped,
        fasta: false,
        quality_mode: qmode,
        delta_encoding: false,
        rle_encoding: false,
    };

    compression::compress(&args).map_err(|e| PyRuntimeError::new_err(format!("{}", e)))
}

/// Decompress a COIL archive to FASTQ format
///
/// Args:
///     input_file: Path to input COIL archive
///     output_file: Path to output FASTQ file (or list of 2 files for paired-end)
///     working_dir: Working directory for temporary files (default: ".")
///     num_threads: Number of threads to use (default: auto)
///     gzipped: Output gzipped FASTQ (default: False)
///     gzip_level: Gzip compression level 0-9 (default: 6)
///
/// Returns:
///     None
///
/// Example:
///     >>> import coil
///     >>> coil.decompress("compressed.coil", "reads.fastq")
///     >>> # Decompress to gzipped output
///     >>> coil.decompress("compressed.coil", "reads.fastq.gz", gzipped=True)
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
    // Parse output files
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
        range: None,
    };

    compression::decompress(&args).map_err(|e| PyRuntimeError::new_err(format!("{}", e)))
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

/// Python module definition
#[pymodule]
fn coil(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compress, m)?)?;
    m.add_function(wrap_pyfunction!(decompress, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;

    // Add version as module attribute
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
