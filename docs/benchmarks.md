# v4zip Benchmarks

All `v3.1` and `v3.0` measurements were run on **Apple M5, macOS, AppleClang 21** (`-DCMAKE_BUILD_TYPE=Release`) with **100% bit-exact lossless restoration verified via SHA-256**.

## 1. v3.1 (CRISPR & Palindromes, Format 4) vs. v3.0 vs. v2.1

| Dataset | Input (bytes) | `v4zip 2.1` | `v4zip 3.0` | **`v4zip 3.1`** | Gain vs. `v2.1` | Gain vs. `v3.0` | Enc / Dec |
|---|---:|---:|---:|---:|---:|---:|---:|
| `m_genitalium.fna` | 587,389 | 130,851 B | 130,000 B | **129,973 B** | -878 B (-0.67%) | -27 B (-0.02%) | 1.55s / 1.60s |
| `e_coli_k12.fna` | 4,699,745 | 1,098,231 B | 1,094,427 B | **1,093,948 B** | -4,283 B (-0.39%) | -479 B (-0.04%) | 12.80s / 12.27s |
| `err000589_10k.fastq` | 1,457,148 | 437,403 B | 405,047 B | **404,839 B** | **-32,564 B (-7.44%)** | -208 B (-0.05%) | **0.86s / 0.82s** |
| `mg_5x_snp1.fna` (5 strains, 1% SNPs) | 2,936,980 | ~654,000 B | 191,532 B | **189,808 B** | **-464,192 B (~3.45x smaller)** | **-1,724 B (-0.90%)** | **1.76s / 1.64s** |
| **Total (4 datasets)** | **9,681,262** | **2,320,485 B** | **1,821,006 B** | **1,818,568 B** | **-501,917 B (-21.63%)** | **-2,438 B (-0.13%)** | — |

## 2. Synthetic CRISPR & Palindromic Dyad Locus (120,001 bp)

A synthetic genomic locus containing 64 canonical 16-bp CRISPR direct repeats and palindromic dyads under all 4 dihedral $D_4$ transformations (`I`, `R`, `C`, `RC`) interspersed with unique 24-bp spacers:

| Codec | DNA Stream Size | Bits / Base | Gain | Lossless |
|---|---:|---:|---:|:---:|
| Baseline (`PredictorV3` / LZ `MIN_MATCH=18`) | 24,723 B | 1.6482 bpb | — | 100% |
| **v3.1 CRISPR & Palindrome Hotspot Index Table** | **21,994 B** | **1.4663 bpb** | **-11.04% (-2,729 B)** | **100%** |

## 3. Partitioning E2E Verification (`--split` & `--chunk-size`)

| Mode & Dataset | Output Parts | Total Size | Enc Time | SHA-256 |
|---|---|---:|---:|:---:|
| `--split 50K` on `m_genitalium.fna` | 3 volumes (`51,200`, `51,200`, `27,573` B) | **129,973 B** (0 B loss) | 1.35s | EXACT |
| `--chunk-size 800K` on `mg_5x_snp1.fna` | 3 chunks + `153 B` `.v4zm` manifest | 435,245 B | 1.73s | EXACT |
| `--chunk-size 400K` on `err000589_10k.fastq` | 4 chunks (4-line aligned) + `.v4zm` | 411,275 B | **0.48s** (~1.8× faster) | EXACT |

## 4. Historical v2.1 Benchmarks (2026-11-10)

Measured on Intel Core i7-10700, Linux x86-64, GCC 15.2.0, libzstd 1.5.7, one fixed CPU.

### FASTA (bits per ACGT base)

| Genome | v4zip 2.1 | GeCo3-l9 (DNA only) | gzip-9 | xz-9 |
|---|---:|---:|---:|---:|
| M. genitalium | 1.8046 | 1.7942 | 2.2858 | 2.1136 |
| E. coli K-12 | 1.8928 | 1.8831 | 2.3769 | 2.1745 |
| S. cerevisiae | 1.9411 | 1.8166 | 2.5278 | 2.2000 |
| H. sapiens chr22 | 1.6669 | 1.6008 | 2.2792 | 1.9666 |
| C. elegans | 1.8960 | 1.7017 | 2.5284 | 2.1404 |
| H. sapiens chr1 | 1.6244 | 1.5684 | 2.2695 | 1.9018 |

### FASTQ (Illumina 159bp, 621,630 reads, 231.5 MiB)

| Compressor | Archive size (MiB) | Ratio |
|---|---:|---:|
| v4zip 2.1 | 17.50 | 13.23x |
| gzip-9 | 39.89 | 5.80x |
| bzip2-9 | 33.57 | 6.90x |
| xz-9 | 26.16 | 8.85x |
| zstd-19 | 30.70 | 7.54x |
