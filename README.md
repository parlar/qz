# QZ

High-performance FASTQ compression using columnar encoding and BSC/BWT.


** do not use in production **

QZ is under active development and the documentation below in not fully accurate. The code is currently messy but will be cleaned up shortly. Bear with me. It is released as open source anyway as it is.


QZ splits FASTQ records into three independent streams — headers, sequences, and quality scores — and compresses each using [libbsc](https://github.com/IlyaGrebnov/libbsc)'s block-sorting transform with adaptive coding. All three streams compress in parallel. Read order is preserved.

## Why QZ?

Most FASTQ compressors treat reads as opaque byte streams and hand them to a general-purpose compressor like gzip or zstd. This ignores the structure of the data: headers follow predictable templates, DNA sequences are drawn from a 4–5 letter alphabet with strong local correlation, and quality scores have their own distinct statistical profile. Mixing all three into a single stream forces the compressor to model very different distributions at once, leaving compression ratio on the table.

QZ takes a columnar approach instead. Each FASTQ record is split into its three constituent streams, and each stream is compressed independently with the algorithm best suited to it. The key insight is that the Burrows-Wheeler Transform — the same algorithm behind bzip2 — is remarkably effective on genomic data when applied to large, homogeneous blocks. DNA sequences from the same organism share vast amounts of k-mer content, and the BWT naturally clusters these shared substrings together, creating long runs that compress extremely well.

The engine behind this is [libbsc](https://github.com/IlyaGrebnov/libbsc), a block-sorting compression library by **Ilya Grebnov**. libbsc pairs a fast BWT implementation (built on Grebnov's [libsais](https://github.com/IlyaGrebnov/libsais) suffix array library) with Lempel-Ziv Prediction (LZP) and Quantized Local Frequency Coding (QLFC), a combination that consistently outperforms zstd on structured genomic data. On raw DNA sequences, BSC achieves ~1.85 bits per base — close to the theoretical limit for order-preserving compression without read reordering.

QZ wraps libbsc in a parallel pipeline: each data stream is split into 25 MB blocks that are compressed independently via rayon, giving near-linear scaling across cores. The result is a compressor that approaches the ratio of research tools like SPRING while running an order of magnitude faster and using a fraction of the memory — without sacrificing original read order.


## Features

- **Columnar encoding** — headers, sequences, and qualities compressed independently
- **BSC/BWT compression** with adaptive QLFC (default) or static QLFC
- **Lossless roundtrip** — bit-exact reconstruction of the original FASTQ
- **Lossy quality modes** — Illumina 8-level binning, 4-level, binary, QVZ, or full discard
- **Ultra mode** — maximum compression with read reordering and de Bruijn graph encoding
- **Paired-end support** — compress and decompress read pairs
- **Chunked streaming** (`--chunked`) — constant-memory compression for arbitrarily large files
- **Multithreaded** — parallel compression and decompression via rayon
- **Gzipped I/O** — read gzipped FASTQ input, write gzipped FASTQ output
- **FASTA support** — compress FASTA files (no quality scores)

## Requirements

- **Rust nightly** (edition 2024)
- **C++ compiler** with OpenMP support (for libbsc)
- **Git** (to clone the libbsc dependency)

## Building

```bash
git clone <repo-url> qz
cd qz
git clone https://github.com/IlyaGrebnov/libbsc.git third_party/libbsc
rustup install nightly
rustup run nightly cargo build --release
# Binary: target/release/qz
```

## Usage

### Default Mode

The default mode gives the best balance of speed and compression. No extra flags needed.

```bash
# Lossless compression
qz compress -i reads.fastq -o reads.qz

# Gzipped input
qz compress -i reads.fastq.gz -o reads.qz --gzipped

# Paired-end
qz compress -i R1.fastq -i R2.fastq -o reads.qz

# FASTA input (no quality scores)
qz compress -i seqs.fasta -o seqs.qz --fasta

# Lossy quality binning
qz compress -i reads.fastq -o reads.qz --quality-mode illumina-bin

# Drop quality scores entirely
qz compress -i reads.fastq -o reads.qz --quality-mode discard
# or equivalently: --no-quality

# Chunked streaming for large files (constant memory)
qz compress -i large.fastq -o large.qz --chunked

# Limit threads
qz compress -i reads.fastq -o reads.qz -t 8
```

# Auto-select ultra level
qz compress -i reads.fastq -o reads.qz --ultra

# Specific ultra level
qz compress -i reads.fastq -o reads.qz --ultra 3
```

### Decompress

Decompression auto-detects all settings from the archive — no flags needed to match what was used during compression.

```bash
qz decompress -i reads.qz -o reads.fastq

# Gzipped output
qz decompress -i reads.qz -o reads.fastq.gz --gzipped --gzip-level 6

# Paired-end
qz decompress -i reads.qz -o R1.fastq -o R2.fastq
```

## Compress Options Reference

### Input / Output

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input FASTQ file(s). One for single-end, two for paired-end. | required |
| `-o, --output FILE` | Output QZ archive. | required |
| `-w, --working-dir PATH` | Working directory for temp files. | `.` |
| `-t, --threads N` | Number of threads (0 = auto-detect). | auto |
| `--gzipped` | Input is gzipped FASTQ. | off |
| `--fasta` | Input is FASTA format (no quality scores). | off |

### Quality

| Flag | Description | Default |
|------|-------------|---------|
| `--quality-mode MODE` | Quality compression mode (see table below). | `lossless` |
| `--no-quality` | Discard quality scores. | off |
| `--quality-compressor COMP` | Quality compressor: `bsc`, `zstd`, `openzl`, `fqzcomp`, `quality-ctx`. | `bsc` |
| `--quality-delta` | Delta encoding between adjacent reads (experimental). | off |
| `--quality-modeling` | Positional quality modeling (experimental). | off |

### Sequences

| Flag | Description | Default |
|------|-------------|---------|
| `--sequence-compressor COMP` | Sequence compressor: `bsc`, `zstd`, `openzl`. | `bsc` |
| `--twobit` | 2-bit encoding (4 bases/byte) + N-mask bitmap. | off |
| `--rc-canon` | Reverse-complement canonicalization with strand flag. | off |
| `--sequence-hints` | Syncmer-derived hint byte before each read for BWT clustering. | off |
| `--sequence-delta` | Inline delta encoding against cached similar reads. | off |
| `--delta-encoding` | Delta encoding for sequences (experimental). | off |
| `--rle-encoding` | Run-length encoding for homopolymers (experimental). | off |
| `--debruijn` | De Bruijn graph compression (patent-safe, RC-aware). | off |
| `--kmer-size N` | K-mer size for de Bruijn compression (0 = auto, range 9–31). | `0` |

### Headers

| Flag | Description | Default |
|------|-------------|---------|
| `--header-compressor COMP` | Header compressor: `bsc`, `zstd`, `openzl`. | `bsc` |
| `--header-template` | Template-based header encoding (extract prefix, delta-encode coordinates). | off |

### Compression Strategy

| Flag | Description | Default |
|------|-------------|---------|
| `--bsc-static` | Use static QLFC instead of adaptive (~5 % larger, slightly faster). | off |
| `--chunked` | Streaming mode: 5 M-record chunks, constant memory. | off |
| `--ultra [LEVEL]` | Ultra compression (levels 1–5, or auto). | off |
| `--reorder {local\|global}` | Reorder reads by similarity. `local` = within chunks; `global` = two-pass sort. Destroys original order. | off |
| `--local-reorder` | Group similar reads by center-hash with delta encoding. Restores order on decompress. | off |
| `--factorize` | Two-pass pattern routing for BWT clustering. | off |
| `--arithmetic` | Arithmetic coding for sequences and qualities (experimental). | off |

## Decompress Options Reference

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input QZ archive. | required |
| `-o, --output FILE` | Output FASTQ file(s). | required |
| `-w, --working-dir PATH` | Working directory for temp files. | `.` |
| `-t, --num-threads N` | Number of threads. | auto |
| `--gzipped` | Output gzipped FASTQ. | off |
| `--gzip-level N` | Gzip compression level (0–9). | `6` |

## Quality Modes

| Mode | Flag | Description |
|------|------|-------------|
| Lossless | `--quality-mode lossless` | Bit-exact preservation (default) |
| Illumina 8-level | `--quality-mode illumina-bin` | Bins to 8 representative values |
| Illumina 4-level | `--quality-mode illumina4` | 4-level binning |
| Binary | `--quality-mode binary` | Threshold at Q20 |
| QVZ | `--quality-mode qvz` | QVZ-style quantization |
| Discard | `--quality-mode discard` | Drop all quality scores |

## Quality Score Compression

Quality scores are the most compressible stream in a FASTQ file, but also the most nuanced. A raw Phred score is a small integer (typically 0–41 for Illumina), and adjacent scores within a read are highly correlated — quality tends to degrade toward the 3' end, with occasional dips at specific positions. Across reads, the same positions tend to produce similar quality profiles. A generic compressor like BSC captures some of this structure through the BWT, but it has no awareness that position 50 of one read is related to position 50 of another, or that a quality of 37 following a 38 is far more likely than a quality of 5 following a 38.

The idea of exploiting this structure with explicit context modeling was pioneered by **James Bonfield**'s [fqzcomp](https://github.com/jkbonfield/fqzcomp) and later refined by **Guillermo Dufort y Álvarez** et al. in [ENANO](https://github.com/guilledufort/EnanoFASTQ). QZ's quality compressor builds on the insights from both, while making different trade-offs suited to short-read Illumina data.

### fqzcomp

[fqzcomp](https://github.com/jkbonfield/fqzcomp) (Bonfield, 2013) was one of the first compressors to show that quality scores could be compressed far more effectively by conditioning on context rather than treating them as raw bytes. It uses an adaptive arithmetic coder with up to 16 bits of context assembled from four configurable sources:

- **Previous quality values** — the last two quality scores, shifted and bitwise-ORed together, optionally remapped through a lookup table
- **Position along the read** — the base offset from the 5' end, optionally quantized via a position lookup table
- **Running delta** — a counter that increments each time a quality value differs from the previous one, capturing whether the signal is smooth or volatile. As Bonfield observed, this running variance has better discriminative power than raw position for predicting quality behavior
- **Selector bits** — arbitrary bits copied from the data stream into the context for further discretization

The context construction is written into the compressed stream itself, allowing the decoder to reconstruct it without external parameters. fqzcomp offers multiple compression strategies (fast to high) and was incorporated into the [CRAM 3.1](https://doi.org/10.1093/bioinformatics/btac010) file format as its quality codec. QZ integrates fqzcomp via FFI through the [htscodecs](https://github.com/samtools/htscodecs) library as an alternative quality backend (`--quality-compressor fqzcomp`).

### ENANO

[ENANO](https://doi.org/10.1093/bioinformatics/btaa551) (Dufort y Álvarez et al., 2020) was designed for nanopore FASTQ files, where the quality alphabet is much larger (94 Phred values vs. ~40 for Illumina) and the scores are noisier. ENANO's key contribution was a richer context model for quality scores that combines:

- **DNA sequence neighborhood** — a sliding window of recent base calls (default length 6), encoded as a rolling hash. This captures the fact that certain sequence motifs produce systematically different quality profiles on nanopore
- **Previous quality values** — the last two Phred scores, quantized and combined into a difference context that encodes both magnitude and direction of change
- **Running prediction error** — per-context running averages of quality scores, with the error between predicted and actual quality fed back to classify contexts by confidence level. This creates a form of adaptive stability tracking: contexts where the model predicts well use different coding tables from contexts where quality is volatile
- **Error magnitude averaging** — accumulated error magnitudes are binned into confidence levels (4 tiers), giving the model a second-order measure of local predictability

These components combine into ~32,768 distinct contexts, each driving its own adaptive arithmetic coder (based on Eugene Shelwien's range coder). ENANO achieves >24% improvement over pigz and ~6% over SPRING on nanopore data.

### QZ's approach (`quality-ctx`)

QZ's quality compressor is a custom LZMA-style forward range coder that draws from both fqzcomp and ENANO but is tuned for fixed-length Illumina short reads. For each quality score, the context is a combination of:

- **Read position** (binned by 5) — captures the characteristic quality curve along the read, similar to fqzcomp's position context but with coarser binning to limit context space
- **Previous quality score** (0–49 Phred) — direct conditioning on the last quality value, as in both fqzcomp and ENANO
- **Stability flag** — a single bit indicating whether quality is stable (|q[i-1] − q[i-2]| ≤ 2) or changing. This is a simplified version of ENANO's multi-level error confidence, reduced to a binary signal that works well for the narrower Illumina quality range
- **DNA base context** — the cross-product of the previous and current sequence base (5 × 5 = 25 states, including N). Like ENANO's sequence neighborhood but limited to a 2-base window, reflecting that Illumina quality artifacts are more localized than nanopore's

This gives ~160,000 distinct contexts (64 position bins × 50 quality levels × 2 stability states × 25 base pairs), each with its own Laplace-smoothed adaptive frequency table that updates on every symbol. On typical Illumina WGS data, this achieves ~0.15 bits per base — about 7% smaller than BSC on the same data. The compressor is enabled automatically in default mode for datasets with fixed-length reads and at least 100,000 records.

### How the three approaches differ

| | fqzcomp | ENANO | QZ (`quality-ctx`) |
|---|---------|-------|--------------------|
| **Target data** | General (CRAM 3.1) | Nanopore | Illumina short reads |
| **Entropy coder** | Adaptive arithmetic | Range coder (Shelwien) | LZMA-style range coder |
| **Context bits** | Up to 16 (configurable) | ~15 (32K contexts) | ~18 (160K contexts) |
| **Position** | Raw or LUT-quantized | Implicit via sequence window | Binned by 5 (64 bins) |
| **Previous quality** | Last 2, shifted+OR'd | Last 2, quantized difference | Last 1, direct (0–49) |
| **Sequence context** | None | Sliding window (6 bases) | Previous + current base (5×5) |
| **Stability/delta** | Running delta counter | Multi-level error confidence | Binary stable/changing flag |
| **Prediction** | None | Per-context running average | None |
| **Read length** | Variable | Variable | Fixed only |
| **Self-describing** | Yes (context spec in stream) | No | No |

The main insight QZ borrows from fqzcomp is the idea that conditioning on quality context (not just position) is essential. From ENANO, QZ takes the stability/delta concept and the use of DNA sequence as a quality predictor. Where QZ simplifies relative to both is in using fewer context levels (binary stability instead of 4-tier, 2-base window instead of 6) — this works because Illumina quality profiles are more regular than nanopore, and the narrower context space converges faster on smaller datasets.

### Lossy quality modes

When lossless quality preservation is not required, QZ can bin quality scores before compression, dramatically reducing the alphabet size and improving ratios across all backends. Illumina 8-level binning maps the full Phred range to 8 representative values that preserve variant calling accuracy, while more aggressive modes (4-level, binary, discard) trade quality fidelity for size.

## Architecture

```
FASTQ input
    ├─ Headers   → BSC (BWT + LZP + adaptive QLFC)
    ├─ Sequences → BSC (BWT + LZP + adaptive QLFC)
    └─ Qualities → BSC (BWT + LZP + adaptive QLFC)
         ↓
    QZ archive (single file with metadata + block offsets)
```

Each stream is split into **25 MB blocks** that are compressed independently and in parallel. The archive stores block offsets, record count, read lengths, and quality metadata for self-contained decompression.

### Chunked Streaming

For files that exceed available memory, `--chunked` mode reads in 5 M-record chunks and pipelines I/O with compression. Temporary files keep memory constant regardless of input size — peak usage is ~6 GB for 100 M reads on 72 threads.

## Environment Variables

| Variable | Effect |
|----------|--------|
| `QZ_NO_BANNER=1` | Suppress the version banner on stderr |
| `RUST_LOG=debug` | Enable debug logging via tracing |

## Testing

```bash
rustup run nightly cargo test --release
```

Integration tests perform full compress-decompress roundtrips across quality modes, compressor combinations, chunked mode, paired-end, and more.

## Acknowledgments

- [libbsc](https://github.com/IlyaGrebnov/libbsc) by Ilya Grebnov — block-sorting compression library (Apache 2.0)
- [libsais](https://github.com/IlyaGrebnov/libsais) by Ilya Grebnov — suffix array construction used by libbsc's BWT (Apache 2.0)
- [fqzcomp](https://github.com/jkbonfield/fqzcomp) by James Bonfield — quality score context modeling, integrated via [htscodecs](https://github.com/samtools/htscodecs) (BSD)
- [ENANO](https://github.com/guilledufort/EnanoFASTQ) by Guillermo Dufort y Álvarez et al. — delta-stability context and sequence-aware quality modeling ([paper](https://doi.org/10.1093/bioinformatics/btaa551))

## License

See LICENSE file.
