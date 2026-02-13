use anyhow::Result;
use clap::Parser;
use tracing::info;

mod cli;
mod compression;
mod io;

use cli::{Cli, Commands};

fn main() -> Result<()> {
    // Initialize logging
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();

    // Show version banner (hide with QZ_NO_BANNER=1)
    if std::env::var("QZ_NO_BANNER").is_err() {
        eprintln!("QZ v{} - Columnar FASTQ compression", env!("CARGO_PKG_VERSION"));
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
    }

    Ok(())
}
