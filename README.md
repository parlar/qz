# QZ

FASTQ compression tool using columnar encoding and BSC/BWT.

QZ splits FASTQ records into three streams (headers, sequences, qualities) and compresses each independently using [libbsc](https://github.com/IlyaGrebnov/libbsc)'s block-sorting transform and adaptive coding. Read order is preserved.

## Features

- Columnar encoding of headers, sequences, and quality scores
- BSC/BWT compression with adaptive QLFC (default) or static QLFC (`--bsc-static`)
- Lossless roundtrip (bit-exact reconstruction)
- Lossy quality modes: Illumina 8-level binning, 4-level binning, binary thresholding, discard
- Chunked streaming mode (`--chunked`) for constant-memory compression of large files
- Multithreaded compression and decompression
- Gzipped FASTQ input/output


## Requirements

- Rust nightly (edition 2024)
- C++ compiler (for libbsc)
- Git (to clone the libbsc submodule)

## Building

```bash
git clone <repo-url> qz
cd qz
git clone https://github.com/IlyaGrebnov/libbsc.git third_party/libbsc
rustup install nightly
rustup run nightly cargo build --release
# Binary is at target/release/qz
```

## Usage

### Compress

```bash
qz compress -i reads.fastq -o reads.qz
qz compress -i reads.fastq.gz -o reads.qz --gzipped
qz compress -i reads.fastq -o reads.qz --quality-mode illumina-bin
qz compress -i reads.fastq -o reads.qz --quality-mode discard
qz compress -i large.fastq -o large.qz --chunked
qz compress -i reads.fastq -o reads.qz --bsc-static
```

### Decompress

```bash
qz decompress -i reads.qz -o reads.fastq
qz decompress -i reads.qz -o reads.fastq.gz --gzipped
```

## Quality modes

| Mode | Flag | Description |
|------|------|-------------|
| Lossless | `--quality-mode lossless` | Bit-exact preservation (default) |
| Illumina 8-level | `--quality-mode illumina-bin` | Bins to 8 representative values |
| Illumina 4-level | `--quality-mode illumina4` | 4-level binning |
| Binary | `--quality-mode binary` | Threshold at Q20 |
| Discard | `--quality-mode discard` | Drop all quality scores |

## Compression modes

| Mode | Flag | Notes |
|------|------|-------|
| BSC adaptive | *(default)* | Two-pass QLFC coder: scans data to build frequency models, then encodes. Best ratio. |
| BSC static | `--bsc-static` | Single-pass QLFC with fixed frequency model. Slightly faster, ~5% larger. |
| Zstd | `--sequence-compressor zstd` etc. | Legacy mode, per-stream override. |

## Environment variables

| Variable | Effect |
|----------|--------|
| `QZ_NO_BANNER=1` | Suppress the version banner on stderr |
| `RUST_LOG=debug` | Enable debug logging |

## Benchmark

100M single-end WGS reads (150bp, ERR3239334, 33.8 GB FASTQ, 72 threads):

| Tool | Mode | Ratio | Wall time | Peak memory | Preserves order |
|------|------|-------|-----------|-------------|-----------------|
| QZ | BSC adaptive (default) | 15.40x | 2m 28s | 6.1 GB | Yes |
| QZ | BSC static | 14.65x | 2m 31s | 6.1 GB | Yes |
| SPRING | default (reorders reads) | 16.52x | 6m 18s | 23.5 GB | No |

SPRING achieves ~7% better compression by reordering reads, at the cost of not preserving input order. QZ preserves order, compresses in less time, and uses less memory.

## Architecture

```
FASTQ input
    ├─ Headers   → BSC (BWT + LZP + QLFC)
    ├─ Sequences → BSC (BWT + LZP + QLFC)
    └─ Qualities → BSC (BWT + LZP + QLFC)
         ↓
    QZ archive (single file)
```

Each stream is split into 25 MB blocks, compressed independently. The archive stores block offsets, record count, read lengths, and quality binning metadata.

## License

See LICENSE file.
