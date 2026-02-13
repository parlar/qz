# QZ - High-Performance FASTQ Compression

## Build & Environment
- Edition 2024 requires **nightly Rust**: `rustup install nightly`
- C++ dependency: `third_party/libbsc/` (clone from https://github.com/IlyaGrebnov/libbsc.git)
- Build: `cargo build --release` (build.rs compiles libbsc as static C++ lib, linked via FFI)
- Release profile: LTO "fat", 1 codegen unit, panic="abort", stripped

## Architecture Overview

```
FASTQ input
    ├─ Headers   → Raw + BSC (BWT + LZP + adaptive QLFC)
    ├─ Sequences → Raw + BSC (block-parallel, 25 MB blocks)
    └─ Qualities → BSC (with optional delta/modeling/dictionary)
         ↓
    QZ archive (single file with metadata + block offsets)
```

All three streams compress **in parallel** via `rayon::join` (headers parallel with nested seq+qual).

### Key Files
- `src/cli/mod.rs` - Clap-based CLI: compress, decompress subcommands
- `src/compression/mod.rs` - Main orchestrator (~2400 lines), compress() and decompress()
- `src/compression/bsc.rs` - FFI bindings to libbsc; `compress_parallel` / `decompress_parallel` split into 25 MB blocks
- `src/io/fastq.rs` - `FastqRecord { id, sequence, quality: Option<String> }`
- `tests/integration_test.rs` - 11 roundtrip tests (~760 lines)

### Module Visibility
- **Private** (internal, access via `super::`): columnar, delta, n_mask, quality, quality_delta, quality_model, read_id, rle, zstd_dict
- **Public**: arithmetic_quality, arithmetic_sequence, bsc, debruijn, dna_utils, paired_end

## Key Design Decisions

### BSC is the default compressor for all streams
- BSC adaptive (BWT + LZP + adaptive QLFC) beats Zstd for structured genomic data
- Headers: Raw+BSC is ~2x better than template+Zstd
- Sequences: Raw+BSC gives ~4.6x (order-preserving); ~1.85 bits/base is near theoretical limit without reordering
- BSC threading: ONLY use `LIBBSC_FEATURE_FASTMODE` (never `MULTITHREADING`) -- rayon handles parallelism

### Block-parallel BSC
- Data split into 25 MB blocks, each compressed/decompressed independently via rayon
- `compress_parallel` / `decompress_parallel` in bsc.rs handle block splitting + reassembly
- `compress_parallel_adaptive()` uses adaptive QLFC (slightly better, ~0.5-0.7%)

### Chunked streaming mode (`--chunked`)
- Reads in 5M-record chunks, BSC-compresses each chunk's streams into 25 MB blocks
- Pipelined: reads next chunk on main thread while compressing previous in `std::thread::scope`
- Streams compressed SEQUENTIALLY (headers -> seq -> qual) to limit BSC working memory
- Temp files (`.qz_chunked_{h,s,q}.tmp`) keep memory constant regardless of input size
- TmpCleanup drop guard ensures cleanup on success, error, or panic
- Peak memory: ~6 GB for 100M reads with 72 threads

### Memory lessons
- BSC allocates output buffer = input_size + header. Call `shrink_to_fit()` on compressed blocks
- Avoid `rayon::join` across all 3 streams simultaneously (70+ concurrent BWT = 14 GB). Process streams sequentially in chunked mode
- Don't accumulate `Vec<Vec<u8>>` blocks in memory for large inputs -- use temp files

### Read order is preserved (no reordering by default)
- SPRING reorders reads for ~5% better compression but destroys original order
- QZ keeps reads in input order for downstream tool compatibility

## Archive Format
```
[encoding_type: 1B] [arithmetic_mode: 1B]
[read_lengths (if arithmetic)]
[quality_binning: 1B] [quality_compressor: 1B] [sequence_compressor: 1B] [header_compressor: 1B]
[quality_model (if enabled)] [quality_delta: 1B] [quality_dict (if enabled)]
[template_prefix + has_comment]
[num_reads: 8B] [headers_len: 8B] [sequences_len: 8B] [nmasks_len: 8B] [qualities_len: 8B]
[headers data] [sequences data] [nmasks data] [qualities data]
```

Compressor codes: zstd=0, bsc=1. Header compressor: zstd=0 (template+zstd), bsc=1 (raw+bsc).

## CompressArgs Gotcha
When adding new fields to `CompressArgs`, you must update ALL constructors in `tests/integration_test.rs` (12+ instances, some nested in for loops with different indentation). Search for `CompressArgs {` to find them all.

## Compression Performance (1M reads, 150bp WGS)
| Tool | Ratio | Compress | Decompress |
|------|-------|----------|------------|
| QZ (BSC default) | 7.67x | 3.4s | 2.4s |
| SPRING (reorders) | 8.06x | 104s | 34s |
| gzip | 4.76x | 32s | 1.8s |

## Test Data
If available in `real_data/`:
- `NA12878_exome_50k.fastq` - 50K reads, 101bp exome
- `ERR3239334_1.1m.fastq` - 1M reads, 150bp WGS
- `ERR3239334_1.5m.fastq` - ~5M reads, 150bp WGS
