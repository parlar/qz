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
| **QZ default** | 435 | **8.03x** | 17.4 s | 5.8 GB | 13.8 s | 6.9 GB |
| **QZ ultra 1** | 426 | **8.21x** | 28.4 s | 7.6 GB | 14.2 s | 8.7 GB |
| **QZ ultra 3** | 416 | **8.39x** | 37.1 s | 13.8 GB | 21.6 s | 8.5 GB |
| **QZ ultra 5** | 416 | **8.39x** | 31.4 s | 14.0 GB | 1:02.5 | 8.5 GB |
| SPRING | 431 | 8.10x | 1:01.4 | 11.9 GB | 15.4 s | 10.0 GB |
| bzip2 -9 | 542 | 6.44x | 2:47.8 | 7.3 MB | 1:26.6 | 4.5 MB |
| pigz -9 | 695 | 5.02x | 9.7 s | 20.8 MB | 7.9 s | 1.7 MB |

SPRING was run without `-r` (read order preserved). Raw timing data in [`benchmarks/results_10m_all_ultra.txt`](benchmarks/results_10m_all_ultra.txt).

## Architecture

### System overview

```
                            ┌──────────────────────────────────────────────┐
                            │              QZ Compression                  │
                            │                                              │
  FASTQ input ──►  FastqReader ──►  Chunk loop (2.5M records) ─────────────┤
  (file, .gz,      (buffered,       │                                      │
   or stdin)        auto-detect     │  ┌────────── rayon::join ─────────┐  │
                    gzip)           │  │                                │  │
                                    │  │  Headers ──► BSC blocks ──┐    │  │
                                    │  │                           │    │  │
                                    │  │  ┌── rayon::join ──────┐  │    │  │
                                    │  │  │                     │  │    │  │
                                    │  │  │ Sequences ► BSC ──┐ │  │    │  │
                                    │  │  │                   │ │  │    │  │
                                    │  │  │ Qualities ► BSC ──┤ │  │    │  │
                                    │  │  │  or quality_ctx   │ │  │    │  │
                                    │  │  └───────────────────┘ │  │    │  │
                                    │  └────────────────────────┘  │    │  │
                                    │                              │    │  │
                                    │         Compressed blocks ◄──┘    │  │
                                    │              │                    │  │
                                    │              ▼                    │  │
                                    │     Temp files or memory          │  │
                                    └──────────────┬────────────────────┘  │
                                                   │                       │
                                                   ▼                       │
                                          Archive assembly                 │
                                          (header + stream blocks)         │
                                                   │                       │
                                                   ▼                       │
                                              .qz archive ───────────────► │
                                          (file or stdout)                 │
                            └──────────────────────────────────────────────┘


                            ┌──────────────────────────────────────────────┐
                            │            QZ Decompression                  │
                            │                                              │
  .qz archive ──►  Parse archive header ──► Spawn 3 decompressor threads   │
  (file or           (magic, version,        │                             │
   stdin→tmpfile)     stream offsets)        │  Thread 1: Headers ──┐      │
                                             │    BSC decompress    │      │
                                             │    (batch=8 blocks)  │      │
                                             │                      │      │
                                             │  Thread 2: Seqs ─────┤      │
                                             │    BSC decompress    │      │
                                             │                      │      │
                                             │  Thread 3: Quals ────┤      │
                                             │    BSC decompress    │      │
                                             │                      │      │
                                             │     bounded channels │      │
                                             │     (capacity = 2)   │      │
                                             │          │           │      │
                                             │          ▼           │      │
                                             │  Main thread:        │      │
                                             │    Reconstruct       │      │
                                             │    FASTQ records     │      │
                                             │    from 3 streams    │      │
                                             │          │           │      │
                                             │          ▼           │      │
                                             │  FASTQ output ──────►│      │
                                             │  (file, .gz,         │      │
                                             │   or stdout)         │      │
                            └──────────────────────────────────────────────┘
```

### Compression pipeline

Compression proceeds in five stages: reading, stream building, parallel compression, block accumulation, and archive assembly.

**Stage 1: Chunked reading.** The FASTQ reader (`FastqReader`) reads records in chunks of 2.5M (default mode) or 5M (ultra mode). Input may be a file, gzipped file (auto-detected via magic bytes), or stdin. Each chunk yields a `Vec<FastqRecord>` where records store raw bytes (`Vec<u8>`) to avoid UTF-8 validation overhead. Reading is pipelined: the main thread reads the next chunk while a background `std::thread::scope` compresses the current one.

```
Main thread:     ┃ Read chunk 0 ┃ Read chunk 1 ┃ Read chunk 2 ┃ ...
                 ┃              ┃              ┃              ┃
Compress thread: ┃              ┃ Compress 0   ┃ Compress 1   ┃ ...
                                ╰──overlapped──╯
```

**Stage 2: Stream building.** Each chunk's records are split into three byte streams by `records_to_streams()`:

| Stream | Format | Content |
|--------|--------|---------|
| Headers | `[varint(len), raw bytes]...` | `@SRR...` identifier lines |
| Sequences | `[varint(len), bases]...` | `ACGTN...` DNA bases |
| Qualities | `[varint(len), packed bytes]...` | Phred+33 scores, optionally bit-packed |

When all reads have the same length (common for Illumina), per-read varint framing is omitted and the constant length is stored once in the archive header, saving ~1 byte per read.

**Stage 3: Parallel compression.** The three streams are compressed concurrently using nested `rayon::join`:

```
rayon::join(
    ║ Headers ──► BSC (25 MB blocks, adaptive QLFC, no LZP)
    ║
    ║ rayon::join(
    ║     ║ Sequences ──► BSC (25 MB blocks)
    ║     ║ Qualities ──► BSC or quality_ctx (500K-read sub-blocks)
    ║ )
)
```

Each stream is split into 25 MB blocks that are compressed independently via `rayon::par_iter`. BSC uses the Burrows-Wheeler Transform followed by adaptive Quantized Local Frequency Coding (QLFC). LZP (Lempel-Ziv Prediction) is disabled in the default path because the BWT already captures repeating patterns; LZP adds 5–10% runtime with negligible compression gain on genomic data.

Compressed blocks are shrunk via `shrink_to_fit()` immediately after compression. Without this, BSC's output buffer allocation (input size + header per block) would waste ~35 GB across 1400 blocks.

**Stage 4: Block accumulation.** Compressed blocks are accumulated either in memory (`Vec<Vec<u8>>` for small inputs) or streamed to temp files (for large inputs or reorder mode). Temp files use a RAII cleanup guard (`TmpCleanup`) that deletes files on drop, including on panic or error.

**Stage 5: Archive assembly.** The archive is written sequentially: v2 header, then header blocks, sequence blocks, quality blocks. No seeking is required, so output can go to stdout. The archive format is:

```
┌──────────────────────────────────────────────────────────────────┐
│  V2 prefix (8 bytes)                                             │
│  ┌──────────┬─────────┬──────────┬────────────────┐              │
│  │ "QZ"     │ ver=02  │ rsvd=00  │ header_size u32│              │
│  └──────────┴─────────┴──────────┴────────────────┘              │
│                                                                  │
│  Header body (variable length)                                   │
│  ┌────────────────────────────────────────────────────────┐      │
│  │ encoding_type, flags, quality_binning                  │      │
│  │ quality_compressor, sequence_compressor, header_comp   │      │
│  │ num_reads (u64), stream lengths (u64 x 4)              │      │
│  │ const_seq_len, const_qual_len (if flag set)            │      │
│  └────────────────────────────────────────────────────────┘      │
│                                                                  │
│  Stream data                                                     │
│  ┌────────────────────────────────────────────────────────┐      │
│  │ Headers:   num_blocks u32, [len u32, data]...          │      │
│  │ Sequences: num_blocks u32, [len u32, data]...          │      │
│  │ Qualities: num_blocks u32, [len u32, data]...          │      │
│  └────────────────────────────────────────────────────────┘      │
└──────────────────────────────────────────────────────────────────┘
```

The header is self-describing: all compressor types, encoding modes, and stream lengths are recorded, so decompression requires no external metadata.

### Decompression pipeline

Decompression uses a streaming architecture with three parallel decompressor threads feeding a single reconstruction thread through bounded channels.

**Header parsing.** The archive header is read and validated (magic bytes `QZ`, version `0x02`). Stream offsets are computed from the recorded lengths: headers start at `data_offset`, sequences at `data_offset + headers_len`, qualities at `data_offset + headers_len + sequences_len`.

**Parallel decompression.** Three background threads are spawned via `std::thread::scope`, each responsible for one stream:

1. Seek to stream offset in the archive file
2. Read `num_blocks` from the stream header
3. In batches of 8: read block headers and data sequentially, decompress the batch in parallel via `rayon::par_iter`, send decompressed blocks through a bounded `SyncSender` channel

The channel capacity is 2 blocks per stream. This backpressures the decompressor threads when the reconstruction thread falls behind, bounding peak memory to ~300 MB (3 streams x ~4 blocks x 25 MB + output buffer).

**Record reconstruction.** The main thread reads from all three channels through `ChannelStreamBuffer` wrappers that provide varint/byte-slice reading over the channel. For each of the `num_reads` records, it reads the header length and bytes, sequence length and bytes, and quality length and packed bytes, unpacks qualities to ASCII, and writes the four FASTQ lines (`@header\nseq\n+\nqual\n`).

**Stdin input.** Since the decompressor needs to seek to three different offsets in the archive, stdin input is spooled to a temp file first, then decompressed normally.

**Gzip output.** When `--gzipped` is requested, output is piped through `gzp::ParCompress` for parallel gzip compression (multi-threaded block compression into a multi-member gzip stream).

### BSC block compression

[libbsc](https://github.com/IlyaGrebnov/libbsc) applies the Burrows-Wheeler Transform (via [libsais](https://github.com/IlyaGrebnov/libsais)) to sort all rotations of the input block, clustering repeated substrings. The BWT output is then compressed by adaptive QLFC, which models local symbol frequencies. The block size is 25 MB.

QZ compiles libbsc and libsais from source via `build.rs` with `-O3 -march=native` and OpenMP support. The complete libbsc pipeline is available: BWT, LZP (Lempel-Ziv Prediction), QLFC (static and adaptive), sort transform, and preprocessing filters. QZ drives it through Rust FFI bindings (`bsc.rs`).

```
Stream (75 MB)
    ├─ Block 0 (25 MB) ──► BWT ──► QLFC ──► compressed (rayon worker 1)
    ├─ Block 1 (25 MB) ──► BWT ──► QLFC ──► compressed (rayon worker 2)
    └─ Block 2 (25 MB) ──► BWT ──► QLFC ──► compressed (rayon worker 3)
```

The default configuration uses adaptive QLFC without LZP (`compress_parallel_adaptive_no_lzp`). Adaptive QLFC learns symbol frequencies during encoding, improving compression by 0.5–0.7% over the static variant. LZP is disabled because the BWT already captures the repeating k-mer structure in genomic data; LZP adds 5–10% runtime with negligible compression benefit.

**Two-level threading model.** libbsc has its own OpenMP-based multithreading for parallelizing the BWT within a single block. QZ adds a second level of parallelism above it using rayon for inter-block parallelism. These two levels interact carefully:

```
                    ┌───────────────────────────────────────┐
                    │       rayon thread pool               │
                    │       (inter-block parallelism)       │
                    │                                       │
                    │  Worker 1: Block 0                    │
                    │    └─ bsc_compress(FASTMODE)          │
                    │       └─ BWT (single-threaded)        │
                    │       └─ QLFC adaptive                │
                    │                                       │
                    │  Worker 2: Block 1                    │
                    │    └─ bsc_compress(FASTMODE)          │
                    │       └─ BWT (single-threaded)        │
                    │       └─ QLFC adaptive                │
                    │                                       │
                    │  Worker N: Block N                    │
                    │    └─ ...                             │
                    └───────────────────────────────────────┘

                    vs.

                    ┌───────────────────────────────────────┐
                    │       Single block, MT mode           │
                    │       (intra-block parallelism)       │
                    │                                       │
                    │  bsc_compress(FASTMODE|MULTITHREADING)│
                    │    └─ BWT via libsais_bwt_omp()       │
                    │       └─ up to 16 OpenMP threads      │
                    │    └─ QLFC adaptive                   │
                    └───────────────────────────────────────┘
```

libbsc supports two feature flags passed per call:

| Flag | Effect |
|------|--------|
| `FASTMODE` | Enables fast decompression; skips expensive data-type detection |
| `MULTITHREADING` | Enables OpenMP parallelism within BWT (libsais) |

The default compression path uses **FASTMODE only**: each 25 MB block gets a single-threaded BWT, but rayon runs many blocks concurrently. This avoids thread explosion (N rayon workers x M OpenMP threads) and gives better throughput than intra-block parallelism for typical genomic workloads with many blocks.

The MT path (`FASTMODE | MULTITHREADING`) is available for single large blocks where intra-block parallelism matters more. When used, OpenMP threads are capped at 12 via `omp_set_num_threads()`.

**libbsc modification.** QZ patches one file in the libbsc source: `bwt.cpp`. Upstream libbsc caps the thread count passed to `libsais_bwt_omp()` at 8. QZ raises this cap to 16 to better utilize high-core machines:

```c
// upstream:  numThreads > 8 ? 8 : numThreads
// QZ:        numThreads > 16 ? 16 : numThreads
```

This affects both `libsais_bwt_aux_omp()` and `libsais_bwt_omp()` calls. On a 72-core system, the effective thread count per BWT call is `min(omp_get_max_threads() / omp_get_num_threads(), 16)`, allowing better scaling when the MT path is used for large individual blocks.

### Quality score compression

QZ automatically selects the quality compressor based on input:

| Condition | Compressor | Rationale |
|-----------|-----------|-----------|
| Lossless, >= 100K reads | `quality_ctx` | Context model outperforms BSC by ~7% once converged |
| Lossless, < 100K reads | BSC | Insufficient data for context model convergence |
| Lossy modes | BSC | Reduced alphabet (3–7 bits) compresses well under BWT |

**Context model.** `quality_ctx` is an LZMA-style forward range coder with 160,000 adaptive contexts. Each quality symbol is coded in a context formed by the cross-product of:

| Feature | States | Description |
|---------|--------|-------------|
| Position bin | 64 | `floor(pos / 5)`, captures positional quality curve |
| Previous quality | 50 | Phred score of preceding base (0–49) |
| Stability | 2 | `\|q[i-1] - q[i-2]\| <= 2` (stable) vs. changing |
| Base pair | 25 | `prev_base * 5 + cur_base` (A/C/G/T/N x A/C/G/T/N) |

Each context maintains a Laplace-smoothed frequency table (~20 symbols for Illumina), updated per symbol and rescaled when the total exceeds 2^20. For large inputs, quality scores are partitioned into 500K-read sub-blocks compressed independently in parallel.

On 150 bp Illumina WGS, `quality_ctx` achieves ~0.15 bits per quality score. Quality data accounts for ~6% of total archive size in lossless mode.

**Lossy quantization.** Quality scores may be quantized before compression: Illumina 8-level binning maps ~40 Phred values to 8 representatives (3 bits/symbol); discard mode replaces all scores with a constant. Quantized scores are bit-packed to the minimum width.

**Relation to prior work.** The context model draws on fqzcomp (Bonfield, 2013; [paper](https://doi.org/10.1093/bioinformatics/btac010)), which demonstrated quality-conditioned arithmetic coding and was adopted as the CRAM 3.1 quality codec, and ENANO (Dufort y Alvarez et al., 2020; [paper](https://doi.org/10.1093/bioinformatics/btaa551)), which introduced DNA sequence context and stability tracking for nanopore data. QZ adapts both for Illumina short reads with a 2-base window and binary stability flag, yielding 160K contexts that converge fast on the narrow Illumina quality distribution.

### Ultra mode

Ultra mode (`--ultra [1-5]`) increases compression by reordering reads to group similar sequences together before compression. Levels control chunk size and parallelism:

| Level | Chunk size | Parallel chunks | Quality sub-block | RAM |
|-------|-----------|-----------------|-------------------|-----|
| 1 | 1M | 4 | 250K | ~8 GB |
| 3 | 5M | 2 | 500K | ~14 GB |
| 5 | 10M | 1 | 1M | ~14 GB |
| auto | varies | varies | varies | fits available RAM |

**Reordering strategy.** Reads are grouped by sequence similarity using center-hash grouping: two 32-bit hashes are computed from the central 32 bases of each read, and reads sharing a hash are placed adjacent in the output. Singletons (unique hashes) remain in input order to preserve BWT-friendly locality. The permutation is stored in the archive so the decompressor can restore original order.

### Memory management

Memory is the primary constraint for high-throughput genomic compression. QZ uses several strategies to keep peak usage bounded:

- **Sequential stream compression.** Within each chunk, headers, sequences, and qualities are compressed sequentially (not all three simultaneously). This prevents 70+ concurrent BWT allocations that would consume ~14 GB.
- **Pipelined I/O.** Reading the next chunk overlaps with compressing the current one, hiding I/O latency without doubling memory.
- **Block shrinking.** `shrink_to_fit()` on each BSC output block releases the unused allocation headroom (BSC allocates output = input + header).
- **Temp file accumulation.** For large inputs, compressed blocks are streamed to disk rather than held in memory. A RAII drop guard ensures cleanup on success, error, or panic.
- **Bounded decompression channels.** Decompressor threads send blocks through channels with capacity 2, preventing unbounded memory growth when the writer is slower than decompression.

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

### Piping (stdin/stdout)

Use `-` for `-i` or `-o` to read from stdin or write to stdout:

```bash
cat reads.fastq | qz compress -i - -o reads.qz                  # compress from stdin
qz compress -i reads.fastq -o - > reads.qz                      # compress to stdout
qz decompress -i reads.qz -o - > reads.fastq                    # decompress to stdout
cat reads.qz | qz decompress -i - -o reads.fastq                # decompress from stdin
cat reads.fastq | qz compress -i - -o - | qz decompress -i - -o - > out.fastq  # full pipe
```

All log output goes to stderr, so stdout remains clean for piped data. Decompression from stdin spools the archive to a temp file first (the decompressor needs to seek to stream offsets).

### CLI reference

**Compress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input FASTQ file (gzipped auto-detected, `-` for stdin) | required |
| `-o, --output FILE` | Output QZ archive (`-` for stdout) | required |
| `-w, --working-dir PATH` | Working directory for temp files | `.` |
| `-t, --threads N` | Thread count (0 = auto) | auto |
| `--fasta` | Input is FASTA format | off |
| `-q, --quality-mode MODE` | `lossless`, `illumina-bin`, or `discard` | `lossless` |
| `--no-quality` | Equivalent to `--quality-mode discard` | off |
| `--ultra [LEVEL]` | Ultra compression (1–5, or omit for auto) | off |

**Decompress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input QZ archive (`-` for stdin) | required |
| `-o, --output FILE` | Output FASTQ file (`-` for stdout) | required |
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
│   ├── libbsc/                    libbsc (BWT + QLFC), compiled via build.rs
│   │   └── libbsc/bwt/bwt.cpp     ← patched: BWT thread cap 8→16
│   └── htscodecs/                 htscodecs (fqzcomp quality codec), unmodified
├── crates/
│   ├── qz-lib/                    core library (all algorithms, no CLI deps)
│   │   ├── src/
│   │   │   ├── compression/
│   │   │   │   ├── mod.rs             archive format, I/O helpers, public API
│   │   │   │   ├── compress_impl.rs   chunked compression orchestrator
│   │   │   │   ├── decompress_impl.rs streaming decompression + header parsing
│   │   │   │   ├── codecs.rs          per-stream compress/decompress dispatch
│   │   │   │   ├── bsc.rs             libbsc FFI, block-parallel BSC, threading
│   │   │   │   ├── quality_ctx.rs     context-adaptive range coder (160K contexts)
│   │   │   │   ├── ultra.rs           ultra mode (reorder + quality_ctx)
│   │   │   │   ├── columnar.rs        quality binning + bit-packing
│   │   │   │   ├── fqzcomp.rs         htscodecs FFI for fqzcomp quality codec
│   │   │   │   ├── header_col.rs      columnar header compression
│   │   │   │   ├── n_mask.rs          2-bit DNA encoding + N-bitmap
│   │   │   │   ├── dna_utils.rs       k-mer hashing, reverse complement
│   │   │   │   └── ...                additional codec modules
│   │   │   ├── io/fastq.rs        FASTQ/FASTA reader (buffered, gzip, stdin)
│   │   │   └── cli.rs             CompressConfig, DecompressConfig (no Clap)
│   │   ├── build.rs               compiles libbsc + htscodecs as static C/C++ libs
│   │   └── tests/                 38 roundtrip integration tests
│   ├── qz-cli/                    CLI binary (Clap) → produces `qz` executable
│   ├── qz-python/                 Python bindings (PyO3/maturin)
│   └── qz-bench/                  development benchmark binaries
├── benchmarks/                    benchmark scripts and results
└── real_data/                     test data (not tracked)
```

### C/C++ dependencies

`build.rs` compiles two C/C++ libraries as static archives linked into the final binary:

**libbsc** — Block-sorting compressor. All source files are compiled with `-O3 -march=native -std=c++11 -fopenmp`. The `LIBBSC_OPENMP_SUPPORT` and `LIBSAIS_OPENMP` defines enable OpenMP-parallel BWT/suffix-array construction. One source file is patched (see [BSC block compression](#bsc-block-compression) above).

**htscodecs** — Only `fqzcomp_qual.c` and `utils.c` are compiled (not the full library). These provide the fqzcomp quality compression algorithm, available as an alternative quality backend. No modifications to upstream source.

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
