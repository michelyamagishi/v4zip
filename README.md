# v4zip

Lossless FASTA/FASTQ compressor based on Klein four-group (V4) symmetry.
Version 3.0 introduces biological genome profiling with specialized compression
pipelines, multi-core stream/Zstd parallelism, 4-ary APM probability calibration,
canonical V4 long-range repeat indexing, and an 8-expert adaptive FASTQ quality
range coder. It preserves every input byte losslessly (verified bit-exact via
SHA-256), including headers, case, ambiguity symbols, line wrapping, LF/CRLF
endings, and quality scores.

## v3.0 Design Premises & Architecture

Premises & architecture enhancements by **Julio Fabio** (`contato@juliofabio.com`),
building on the Klein V4 mathematical foundation:

1. **Compression Ratio First (Terabyte-Scale Readiness):** storage efficiency is
   the primary objective; models prioritize minimal bit/base output over raw
   single-core throughput, with dedicated handling for multi-strain / mapped
   repeat-heavy collections.
2. **Biological Profiling & Strategy Dispatch (`GenomeProfiler`):** instead of
   forcing one generic pipeline on all inputs, a fast pre-pass measures GC ratio,
   3-base codon periodicity, homopolymer run density, record fragmentation, and
   canonical V4 16-mer repeat density to route each DNA stream to a specialized
   strategy (`coding_small`, `coding_large`, `short_reads`, `repeat_heavy`,
   `low_complex`) recorded in the v3 stream header.
3. **Specialized Pipelines:**
   - **`repeat_heavy` (`RepeatDnaCodec`):** canonical Klein V4 LZ77 self-index
     (forward, complement, reverse, reverse-complement) + Zstd-19 match tokens +
     `PredictorV3` residual coding for non-matching bases/SNPs.
   - **`PredictorV3`:** 34-expert log-linear/linear mixture (dense orders `0..9`,
     5 Klein V4 `OrbitTableV3` orders `{10, 12, 15, 18, 22}`, codon phase `2/5/8`
     + wobble, skip/gapped, run-length, and `MatchModelV3` for `k = {9/11, 16, 24}`)
     followed by 4-ary Adaptive Probability Map (APM) post-calibration.
   - **`QualityModel` v4:** 8-expert adaptive range coder for FASTQ quality strings
     (`o0`, `o1`, collision-free `o2` `64x64`, 2-way `o3`, exact `pos x p1`,
     `pos x p1/2 x p2/2`, `run`, `mean`).
4. **Multi-Core Execution:** independent archive streams (`dna`, `quality`,
   `identifiers`/`separators`, layout) are encoded/decoded concurrently via
   `std::async`, paired with multi-threaded `zstd` level 19 (`ZSTD_c_nbWorkers`).
5. **Strict Lossless & Backward Compatibility:** the v3.0 encoder writes format 3;
   the decoder restores formats 1, 2, and 3 bit-exactly.

## Building

Requires a C++17 compiler, CMake 3.16+, libzstd development headers/library,
and a 64-bit POSIX system. Run the commands from the directory containing
`CMakeLists.txt`. This source package builds only the command-line program;
Python and separate test or development directories are not required.

Compression ratio is the priority. Model tables use approximately **3.94 GiB**;
file buffers require additional RAM. The maximum total process memory measured
on the benchmark inputs below was **4.54 GiB**. Larger inputs may require more.

### Linux (GCC)

On Debian or Ubuntu, install the dependencies:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake libzstd-dev
```

Build and install into a local directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
cmake --install build --prefix "$PWD/install"
```

### macOS (Clang)

Install Apple's Command Line Tools if they are not already installed:

```bash
xcode-select --install
```

With [Homebrew](https://brew.sh/) installed, install
[CMake](https://formulae.brew.sh/formula/cmake) and
[zstd](https://formulae.brew.sh/formula/zstd), then build:

```bash
brew install cmake zstd
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$(brew --prefix zstd)"
cmake --build build --parallel 2
cmake --install build --prefix "$PWD/install"
```

The executable is `build/v4zip`; the installed copy is `install/bin/v4zip`.
Linux and macOS (Apple Silicon / AppleClang 21) were verified locally.

## Usage

```bash
# FASTA
./build/v4zip compress input.fasta output.v4z
./build/v4zip decompress output.v4z restored.fasta

# FASTQ
./build/v4zip fastq compress input.fastq output.v4zq
./build/v4zip fastq decompress output.v4zq restored.fastq
```

Options: `--verbose`, `--help`, `--version`, and `--max-output BYTES` for
decompression. The default restored-file limit is 1 GiB. Existing output files
are never overwritten. The previous `-o` and `--stats` options are not supported.

FASTA and four-line FASTQ are supported, including empty records and missing
final newlines. FASTQ sequence and quality lines must have equal byte lengths;
wrapped FASTQ is not supported. Quality scores are preserved losslessly.

The encoder writes rewrite format 3; the decoder reads rewrite formats 1, 2, and 3.
Archives from the original upstream v4zip 1.x legacy tool require the original executable.

## Benchmarks

### v3.0 vs. v2.1 (Apple M5, macOS, AppleClang 21, SHA-256 verified)

| Dataset | Input (bytes) | v4zip 2.1 | v4zip 3.0 | Gain | Enc / Dec |
|---|---:|---:|---:|---:|---:|
| `m_genitalium.fna` | 587,389 | 130,851 B | **130,000 B** | -851 B (-0.65%) | 1.12s / 1.10s |
| `e_coli_k12.fna` | 4,699,745 | 1,098,231 B | **1,094,427 B** | -3,804 B (-0.35%) | 9.89s / 10.27s |
| `err000589_10k.fastq` | 1,457,148 | 437,403 B | **405,047 B** | **-32,356 B (-7.40%)** | **0.85s / 0.75s** |
| `mg_5x_snp1.fna` (5 strains, 1% SNPs) | 2,936,980 | ~654,000 B | **191,532 B** | **~3.4x smaller** | **1.51s / 1.34s** |

### Historical v2.1 Benchmarks (2026-11-10)

The v4zip 2.1 measurements are from 2026-11-10: three serial repetitions on each
complete input, with byte-exact restoration in every repetition. Hardware:
Intel Core i7-10700, Linux x86-64, GCC 15.2.0, libzstd 1.5.7, one fixed CPU.
The previously measured competitor results are from 2026-11-05 on the same inputs.

### FASTA (bits per ACGT base)

| Genome | v4zip 2.1 | GeCo3-l9 (DNA only) | gzip-9 | xz-9 |
|---|---:|---:|---:|---:|
| M. genitalium | 1.8046 | 1.7942 | 2.2858 | 2.1136 |
| E. coli K-12 | 1.8928 | 1.8831 | 2.3769 | 2.1745 |
| S. cerevisiae | 1.9411 | 1.8166 | 2.5278 | 2.2000 |
| H. sapiens chr22 | 1.6669 | 1.6008 | 2.2792 | 1.9666 |
| C. elegans | 1.8960 | 1.7017 | 2.5284 | 2.1404 |
| H. sapiens chr1 | 1.6244 | 1.5684 | 2.2695 | 1.9018 |

The v4zip, gzip, and xz columns count complete archives, including headers and
formatting. The GeCo3 column restores only the uppercase ACGT projection, so it
represents a different restoration task. Lower bits/base means a smaller archive
within the same restoration task. GeCo3 used level 9; gzip used 1.14 and xz 5.8.3.

### FASTQ (Illumina 159bp, 621,630 reads, 231.5 MiB)

| Compressor | Archive size (MiB) | Ratio |
|---|---:|---:|
| v4zip 2.1 | 17.50 | 13.23x |
| gzip-9 | 39.89 | 5.80x |
| bzip2-9 | 33.57 | 6.90x |
| xz-9 | 26.16 | 8.85x |
| zstd-19 | 30.70 | 7.54x |

All FASTQ entries restore the complete original file. MiB means 2^20 bytes;
ratio is original file bytes divided by archive bytes. Competitor versions:
gzip 1.14, bzip2 1.0.8, xz 5.8.3, and zstd 1.5.7.

Across all seven complete inputs, aggregate archives decreased from 117,686,700
bytes with the 2.0 rewrite to 100,126,639 bytes with 2.1 (**14.92% smaller**).
No individual input grew. These are measured corpus results, not a claim of
optimal compression for every dataset.

## Testing

The reduced package has no automated test suite. After building, run these
manual round-trip checks from the project directory, using filenames that do
not already exist:

```bash
printf '>example\nAcGTNN\n' > example.fasta
./build/v4zip compress example.fasta example.v4z
./build/v4zip decompress example.v4z restored.fasta
cmp example.fasta restored.fasta

printf '@example\nAcGT\n+details\n!I~#\n' > example.fastq
./build/v4zip fastq compress example.fastq example.v4zq
./build/v4zip fastq decompress example.v4zq restored.fastq
cmp example.fastq restored.fastq
```

Each `cmp` command finishes silently when the restored file is identical.

## Citation

If you use v4zip in your research, please cite:

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
