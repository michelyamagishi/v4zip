# v4zip

Lossless FASTA/FASTQ compressor based on Klein four-group ($V_4$) and dihedral ($D_4$) DNA symmetry.

`v4zip` preserves every input byte bit-for-bit (verified via CRC32 and SHA-256), including headers, letter case, ambiguity symbols, line wrapping, LF/CRLF endings, and FASTQ quality scores. It is completely species-agnostic and learns all sequence statistics online.

## Highlights (v3.1)

- **CRISPR & Palindrome Profiling + Hotspot Indexing**: Detects palindromic dyad self-symmetry and CRISPR periodic repeats (`20..120` bp), adding collision-free direct-indexed palindrome/reverse-complement experts (`rc8_`, `rc9_`, `rc10_`), a 4-way `MatchModelV3` (`k = 13` CRISPR/palindrome scale), and a 4-way `RepeatDnaCodec` with a 64-entry MDL-gated Hotspot Index Table.
- **Two Orthogonal Partitioning Strategies**: Supports `--split SIZE` (exact-size volumes `.001`, `.002`, ... with 0 bytes compression loss) and `--chunk-size SIZE` (parallel multi-core FASTA/FASTQ boundary-aligned `.partXXXX.v4z` chunks + `.v4zm` manifest).
- **100% Backward Compatibility & Hardened Decoder**: Writes **Format 4** while decoding Formats **1, 2 (`v2.1`), 3 (`v3.0`), and 4 (`v3.1`)** bit-for-bit with anti-OOM and anti-Zstd-bomb caps.

## Documentation

To keep this README concise, detailed topics are organized in [`docs/`](docs/):

- **[User Guide (Use Cases, Genomes, Complexity & Transfer)](docs/user-guide.md)** — Decision matrix and practical examples by genome type, size, complexity, and data transfer channel.
- **[Architecture & Compression Methods](docs/architecture.md)** — `GenomeProfiler`, 37-expert `PredictorV3`, `RepeatDnaCodec` ($D_4$ + CRISPR/Palindrome Hotspot Index Table), and `QualityModel` v4.
- **[Partitioning Strategies (`--split` & `--chunk-size`)](docs/partitioning.md)** — Exact volume splitting vs. parallel biological chunked manifests (`.v4zm`).
- **[Archive Formats & Compatibility](docs/formats-and-compatibility.md)** — Format 1/2/3/4 compatibility matrix and decoder hardening.
- **[Benchmarks](docs/benchmarks.md)** — Full `v2.1` vs. `v3.0` vs. `v3.1` tables, synthetic CRISPR/palindrome locus results, and historical comparisons.

## Quick Benchmark Summary (`v2.1` vs. `v3.0` vs. `v3.1`)

| Dataset | Input | `v4zip 2.1` | `v4zip 3.0` | **`v4zip 3.1`** | Gain vs. `v2.1` | Gain vs. `v3.0` |
|---|---:|---:|---:|---:|---:|---:|
| `m_genitalium.fna` | 587,389 B | 130,851 B | 130,000 B | **129,973 B** | -0.67% | -0.02% |
| `e_coli_k12.fna` | 4,699,745 B | 1,098,231 B | 1,094,427 B | **1,093,948 B** | -0.39% | -0.04% |
| `err000589_10k.fastq` | 1,457,148 B | 437,403 B | 405,047 B | **404,839 B** | **-7.44%** | -0.05% |
| `mg_5x_snp1.fna` (5 strains) | 2,936,980 B | ~654,000 B | 191,532 B | **189,808 B** | **~3.45x smaller** | **-0.90%** |
| **Total (4 datasets)** | **9,681,262 B** | **2,320,485 B** | **1,821,006 B** | **1,818,568 B** | **-21.63%** | **-0.13%** |

See [`docs/benchmarks.md`](docs/benchmarks.md) for full timings, synthetic CRISPR/palindrome benchmarks (`-11.04%`), and partitioning measurements.

## Building

Requires C++17, CMake 3.16+, `libzstd`, and a 64-bit POSIX system (Linux or macOS).

```bash
# Linux (Debian/Ubuntu): sudo apt-get install build-essential cmake libzstd-dev
# macOS (Homebrew):      brew install cmake zstd

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

## Usage

```bash
# Single-file FASTA / FASTQ
./build/v4zip compress input.fasta output.v4z [--verbose]
./build/v4zip decompress output.v4z restored.fasta [--max-output BYTES] [--verbose]

./build/v4zip fastq compress input.fastq output.v4zq [--verbose]
./build/v4zip fastq decompress output.v4zq restored.fastq [--verbose]

# Exact-size volume split (0% compression loss)
./build/v4zip compress genome.fna genome.v4z --split 20M --verbose
./build/v4zip decompress genome.v4z.001 restored.fna

# Parallel biological chunked manifest (.v4zm)
./build/v4zip fastq compress reads.fastq reads.v4zm --chunk-size 32M --verbose
./build/v4zip fastq decompress reads.v4zm restored.fastq
```

Existing outputs are never overwritten. Default decompression limit is `1 GiB` (`--max-output` overrides it).

## Citation

```bibtex
@book{yamagishi2017mathematical,
  title     = {Mathematical Grammar of Biology},
  author    = {Yamagishi, Michel Eduardo Beleza},
  year      = {2017},
  publisher = {Springer International Publishing},
  series    = {SpringerBriefs in Mathematics},
  doi       = {10.1007/978-3-319-62689-5},
  isbn      = {978-3-319-62689-5}
}
```

## License

MIT
