# QZ

Order-preserving FASTQ compression using columnar stream separation and block-sorting transforms.

## Overview

QZ decomposes each FASTQ record into three streams — headers, DNA sequences, and quality scores — and compresses each independently. The primary compression engine is [libbsc](https://github.com/IlyaGrebnov/libbsc) (Grebnov), which applies the Burrows-Wheeler Transform followed by Quantized Local Frequency Coding (QLFC). Separating the streams allows the BWT to operate on large, statistically homogeneous blocks: DNA sequences from the same organism share extensive k-mer content, and the BWT clusters these shared substrings, producing long runs amenable to entropy coding. On raw DNA sequences, this achieves ~1.85 bits per base without read reordering.

For quality scores, QZ implements a context-adaptive range coder (`quality_ctx`) that conditions on read position, previous quality value, a stability flag, and the local DNA sequence context. This exploits the strong positional and sequential correlations in Illumina quality profiles that block-sorting alone cannot capture, reducing quality stream size by ~7% relative to BSC.

All streams are split into blocks (25 MB for BSC, 500K reads for quality_ctx) and compressed in parallel via [rayon](https://github.com/rayon-rs/rayon). Input is read in 2.5M-record chunks with pipelined I/O. Read order is preserved.

## Results

10M reads, 150 bp, whole-genome sequencing (ERR3239334), 3,492 MB uncompressed. 72-core machine, sequential runs, byte-identical roundtrip verified via MD5.

| Tool | Size (MB) | Ratio | Compress | Comp RAM | Decompress | Dec RAM |
|------|-----------|-------|----------|----------|------------|---------|
| **QZ default** | 435 | **8.03x** | 17.6 s | 5.9 GB | 14.7 s | 6.7 GB |
| **QZ ultra 3** | 416 | **8.39x** | 34.2 s | 13.8 GB | 21.3 s | 8.5 GB |
| SPRING | 431 | 8.10x | 65.2 s | 9.9 GB | 16.9 s | 10.0 GB |
| bzip2 -9 | 542 | 6.44x | 168.5 s | 7.3 MB | 88.2 s | 4.5 MB |
| pigz -9 | 695 | 5.02x | 10.0 s | 22.5 MB | 8.0 s | 1.7 MB |

SPRING was run without `-r` (read order preserved). Raw timing data in [`benchmarks/results_10m_optimized.txt`](benchmarks/results_10m_optimized.txt).

## Method

### Stream separation

Each FASTQ record contributes to three independent byte streams:

```
FASTQ record
    ├─ Header   (@SRR... instrument/run/flowcell/...)
    ├─ Sequence  (ACGTN...)
    └─ Quality   (Phred+33 ASCII)
```

Headers, sequences, and qualities are concatenated separately, then each stream is split into 25 MB blocks and compressed in parallel using libbsc's adaptive BWT + QLFC pipeline. The three streams are compressed concurrently via `rayon::join`.

### Quality score compression

Quality scores exhibit strong positional correlation (systematic 3' degradation), sequential correlation (adjacent scores are similar), and sequence-dependent effects (certain motifs produce characteristic quality profiles). Generic block-sorting captures some of this structure, but a context model that conditions directly on these features achieves better compression.

**Compressor selection.** QZ automatically selects between two quality compressors based on input characteristics:

| Condition | Compressor | Rationale |
|-----------|-----------|-----------|
| Lossless, >= 100K reads | `quality_ctx` | Adaptive context models require sufficient data to converge; above this threshold they outperform BSC by ~7% |
| Lossless, < 100K reads | BSC | Insufficient data for context model convergence |
| Lossy modes | BSC | Reduced alphabet (3–7 bits) compresses efficiently under BWT without context modeling |

**Context model.** `quality_ctx` is an LZMA-style forward range coder with 160,000 adaptive contexts. For each quality symbol, the context is the cross-product of:

| Feature | States | Description |
|---------|--------|-------------|
| Position bin | 64 | `floor(pos / 5)`, captures positional quality curve |
| Previous quality | 50 | Phred score of preceding base (0–49) |
| Stability | 2 | `|q[i-1] - q[i-2]| <= 2` (stable) vs. changing |
| Base pair | 25 | `prev_base * 5 + cur_base` (A/C/G/T/N x A/C/G/T/N) |

Each context maintains a Laplace-smoothed frequency table over the observed quality alphabet (~20 symbols for Illumina), updated on every symbol and rescaled (halving with minimum frequency 1) when the total exceeds 2^20. The range coder uses 64-bit carry propagation and renormalizes below 2^24.

On 150 bp Illumina WGS, `quality_ctx` achieves ~0.15 bits per quality score. Quality data accounts for ~6% of total archive size in lossless mode.

**Parallelism.** For large inputs, quality scores are partitioned into 500K-read sub-blocks, each encoded independently and in parallel. The sub-block format is self-describing: `[num_blocks: u32][block_len: u32, blob]...`.

**Lossy quantization.** Before compression, quality scores may be quantized: Illumina 8-level binning maps ~40 Phred values to 8 representatives (3 bits/symbol); discard mode replaces all scores with a constant. Quantized scores are bit-packed to the minimum width before compression.

**Relation to prior work.** The context model draws on two key contributions:

- *fqzcomp* (Bonfield, 2013; [paper](https://doi.org/10.1093/bioinformatics/btac010)) demonstrated that conditioning on previous quality values and a running delta counter substantially improves quality compression. fqzcomp uses an adaptive arithmetic coder with up to 16 configurable context bits and was adopted as the quality codec in CRAM 3.1. QZ integrates fqzcomp via FFI through [htscodecs](https://github.com/samtools/htscodecs) as an alternative backend.

- *ENANO* (Dufort y Álvarez et al., 2020; [paper](https://doi.org/10.1093/bioinformatics/btaa551)) introduced DNA sequence context (6-base sliding window) and multi-level error confidence tracking into quality modeling, achieving strong results on nanopore data with ~32K contexts.

QZ adapts both ideas for Illumina short reads: a 2-base sequence window (vs. ENANO's 6) and a binary stability flag (vs. ENANO's 4-tier confidence), yielding a larger context space (160K vs. 32K) that converges faster on the narrower Illumina quality distribution.

| | fqzcomp | ENANO | QZ |
|---|---------|-------|----|
| Target data | General (CRAM 3.1) | Nanopore | Illumina short reads |
| Entropy coder | Adaptive arithmetic | Range coder (Shelwien) | LZMA-style range coder |
| Contexts | Up to 2^16 | ~32K | 160K |
| Position | Raw or LUT-quantized | Implicit via sequence window | Binned by 5 (64 bins) |
| Previous quality | Last 2, shifted+OR | Last 2, quantized difference | Last 1, direct (0–49) |
| Sequence context | None | 6-base sliding window | Previous + current base (5x5) |
| Stability/delta | Running delta counter | 4-tier error confidence | Binary stable/changing |
| Self-describing | Yes | No | No |

### Ultra mode

Ultra mode (`--ultra [1-5]`) increases compression by processing 5M-record chunks, enabling BSC with LZP (Lempel-Ziv Prediction before BWT), and compressing 2 chunks in parallel. Levels trade speed and memory for ratio:

| Level | RAM | Description |
|-------|-----|-------------|
| 1 | ~8 GB | Fast |
| 3 | ~17 GB | High compression |
| 5 | ~17 GB | Maximum compression |
| auto | varies | Selects highest level fitting available RAM |

## Installation

### Requirements

- Rust nightly (edition 2024)
- C++ compiler with OpenMP support (for libbsc)

### Build

```bash
git clone <repo-url> qz && cd qz
git clone https://github.com/IlyaGrebnov/libbsc.git third_party/libbsc
rustup install nightly
rustup run nightly cargo build --release
# Binary: target/release/qz
```

### Python bindings

```bash
cd crates/qz-python
pip install maturin
maturin develop --release
```

## Usage

### Compress

```bash
qz compress -i reads.fastq -o reads.qz                          # lossless (default)
qz compress -i reads.fastq.gz -o reads.qz                       # gzipped input (auto-detected)
qz compress -i reads.fastq -o reads.qz --quality-mode illumina-bin  # lossy 8-level binning
qz compress -i reads.fastq -o reads.qz --quality-mode discard   # discard quality scores
qz compress -i reads.fastq -o reads.qz --ultra 3                # ultra compression, level 3
qz compress -i seqs.fasta -o seqs.qz --fasta                    # FASTA input
```

### Decompress

Decompression auto-detects all settings from the archive header.

```bash
qz decompress -i reads.qz -o reads.fastq
qz decompress -i reads.qz -o reads.fastq.gz --gzipped           # gzipped output
```

### CLI reference

**Compress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input FASTQ file (gzipped auto-detected) | required |
| `-o, --output FILE` | Output QZ archive | required |
| `-w, --working-dir PATH` | Working directory for temp files | `.` |
| `-t, --threads N` | Thread count (0 = auto) | auto |
| `--fasta` | Input is FASTA format | off |
| `-q, --quality-mode MODE` | `lossless`, `illumina-bin`, or `discard` | `lossless` |
| `--no-quality` | Equivalent to `--quality-mode discard` | off |
| `--ultra [LEVEL]` | Ultra compression (1–5, or omit for auto) | off |

**Decompress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input QZ archive | required |
| `-o, --output FILE` | Output FASTQ file | required |
| `-w, --working-dir PATH` | Working directory for temp files | `.` |
| `-t, --threads N` | Thread count | auto |
| `--gzipped` | Output gzipped FASTQ | off |
| `--gzip-level N` | Gzip level (0–9) | `6` |

### Python

```python
import qz

qz.compress("reads.fastq", "reads.qz")
qz.compress("reads.fastq", "reads.qz", ultra=3)
qz.compress("reads.fastq", "reads.qz", quality_mode="illumina-bin")

qz.decompress("reads.qz", "reads.fastq")
qz.decompress("reads.qz", "reads.fastq.gz", gzipped=True)

print(qz.version())
```

### Library (Rust)

```rust
use qz_lib::cli::{CompressConfig, QualityMode};
use qz_lib::compression;

let config = CompressConfig {
    input: vec!["reads.fastq".into()],
    output: "reads.qz".into(),
    quality_mode: QualityMode::Lossless,
    ..CompressConfig::default()
};
compression::compress(&config)?;
```

## Project Structure

```
qz/
├── Cargo.toml                     workspace root
├── third_party/
│   ├── libbsc/                    libbsc (BWT + QLFC engine)
│   └── htscodecs/                 htscodecs (fqzcomp quality codec)
├── crates/
│   ├── qz-lib/                    core library (all algorithms, no CLI deps)
│   │   ├── src/
│   │   │   ├── compression/
│   │   │   │   ├── mod.rs             module hub, constants, archive I/O helpers
│   │   │   │   ├── compress_impl.rs   compression orchestration + streaming
│   │   │   │   ├── decompress_impl.rs decompression + archive header parsing
│   │   │   │   ├── codecs.rs          per-stream compress/decompress functions
│   │   │   │   ├── bsc.rs             libbsc FFI bindings, block-parallel BSC
│   │   │   │   ├── quality_ctx.rs     context-adaptive range coder for qualities
│   │   │   │   ├── ultra.rs           ultra mode (reorder + HARC delta encoding)
│   │   │   │   ├── template.rs        de Bruijn graph template sequence compression
│   │   │   │   └── ...                20+ codec/algorithm modules
│   │   │   └── io/fastq.rs        FASTQ/FASTA reader/writer
│   │   ├── build.rs               compiles libbsc + htscodecs via cc
│   │   └── tests/                 43 roundtrip integration tests
│   ├── qz-cli/                    CLI binary (Clap) → produces `qz` executable
│   ├── qz-python/                 Python bindings (PyO3/maturin)
│   └── qz-bench/                  development benchmark binaries
├── benchmarks/                    benchmark scripts and results
└── real_data/                     test data (not tracked)
```

## Testing

```bash
rustup run nightly cargo test --release
```

131 tests total: 88 unit tests (per-module codec roundtrips) and 43 integration tests covering lossless, lossy quality modes, ultra mode, FASTA, arithmetic/de Bruijn/delta/RLE encodings, error handling on corrupt archives, and edge cases.

## Environment Variables

| Variable | Effect |
|----------|--------|
| `QZ_NO_BANNER=1` | Suppress version banner on stderr |
| `RUST_LOG=debug` | Enable debug logging (tracing) |

## Acknowledgments

- [libbsc](https://github.com/IlyaGrebnov/libbsc) by Ilya Grebnov — block-sorting compression (BWT + LZP + QLFC). Apache 2.0.
- [libsais](https://github.com/IlyaGrebnov/libsais) by Ilya Grebnov — suffix array construction used by libbsc's BWT. Apache 2.0.
- [fqzcomp](https://github.com/jkbonfield/fqzcomp) by James Bonfield — context-modeled quality compression; quality codec in CRAM 3.1. Integrated via [htscodecs](https://github.com/samtools/htscodecs). BSD.
  - Bonfield, J.K. and Mahoney, M.V. (2013). Compression of FASTQ and SAM format sequencing data. *PLoS ONE*, 8(3):e59190.
- [ENANO](https://github.com/guilledufort/EnanoFASTQ) by Guillermo Dufort y Álvarez et al. — sequence-aware quality modeling with stability tracking.
  - Dufort y Álvarez, G. et al. (2020). ENANO: encoder for nanopore FASTQ files. *Bioinformatics*, 36(16):4506–4507.

## License

See LICENSE file.
