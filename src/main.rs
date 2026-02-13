use anyhow::Result;
use clap::Parser;
use tracing::info;

mod cli;
mod compression;
mod encoding;
mod io;

use cli::{BenchmarkArgs, Cli, Commands};

fn main() -> Result<()> {
    // Initialize logging
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();

    // Show version banner (hide with FQZ_NO_BANNER=1)
    if std::env::var("FQZ_NO_BANNER").is_err() {
        eprintln!("FQZ v{} - Columnar FASTQ compression", env!("CARGO_PKG_VERSION"));
        eprintln!("Compression: BSC/BWT-based columnar encoding");
        eprintln!();
    }

    match cli.command {
        Commands::Compress(args) => {
            info!("Starting compression...");
            compression::compress(&args)?;
            info!("Compression complete!");
        }
        Commands::Decompress(args) => {
            info!("Starting decompression...");
            compression::decompress(&args)?;
            info!("Decompression complete!");
        }
        Commands::Benchmark(args) => {
            info!("Starting benchmark...");
            run_benchmark(&args)?;
            info!("Benchmark complete!");
        }
    }

    Ok(())
}

fn run_benchmark(args: &BenchmarkArgs) -> Result<()> {
    use io::FastqReader;

    info!("Loading FASTQ data from: {:?}", args.input);
    let mut reader = FastqReader::from_path(&args.input, false)?;
    let mut records = Vec::new();

    while let Some(record) = reader.next()? {
        records.push(record);
        if let Some(max) = args.max_reads {
            if records.len() >= max {
                break;
            }
        }
    }

    info!("Loaded {} reads", records.len());

    let report = compression::benchmark::run_benchmark(&records, args.stream.as_deref())?;
    report.print();

    Ok(())
}
