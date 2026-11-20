# v4zip User Guide — Use Cases, Genome Types, Complexity & Data Transfer

This guide explains how to choose the right `v4zip` command, flags, and partitioning strategy depending on your **genome type**, **file size**, **sequence complexity**, and **transfer/storage constraints**.

---

## 1. Decision Matrix (Quick Selection)

| Use Case / Genome Type | Typical Input Size | Sequence Complexity | Recommended Command | Why |
|---|---|---|---|---|
| **Viral, Plasmid, Mitochondrial, Small Bacteria** *(e.g. M. genitalium, SARS-CoV-2)* | `< 1.5 Mbp` (< 2 MiB) | High coding density, low repeat rate | `v4zip compress in.fna out.v4z` | Auto-dispatches `coding_small` (`PredictorV3` compact APM); single file is already small (< 500 KB). |
| **Bacteria & Fungi** *(e.g. E. coli, S. cerevisiae)* | `1.5 Mbp – 50 Mbp` | Coding-rich, moderate repeats & palindromes | `v4zip compress in.fna out.v4z` | Auto-dispatches `coding_large` (`PredictorV3` 37 experts + `rc8/9/10` + 4-way `k=13` match). |
| **Multi-Strain Pan-Genomes & CRISPR/Transposon Loci** *(e.g. 5–100 bacterial strains)* | `2 Mbp – 500 Mbp` | High canonical 16-mer, CRISPR (`20..120` bp) or palindromic dyad repeats | `v4zip compress pangenome.fna out.v4z` *(or add `--split 20M` for transfer)* | Auto-dispatches `repeat_heavy` (`RepeatDnaCodec` 4-way $D_4$ LZ + 64-entry Hotspot Index Table); achieves 3×–20× compression over single-strain baselines. |
| **Low-Complexity / Tandem-Repeat Sequences** *(centromeres, telomeres, homopolymers)* | Any | `run_rate >= 0.25` | `v4zip compress seq.fna out.v4z` | Auto-dispatches `low_complex` (boosted run-length smoothing and fast local adaptation). |
| **Illumina / Short-Read FASTQ** *(medium runs)* | `1 MiB – 500 MiB` | Fragmented short reads + quality scores | `v4zip fastq compress in.fastq out.v4zq` | Auto-dispatches `short_reads` (`k0 = 9`) + 8-expert `QualityModel` v4. |
| **Large FASTQ Runs / Multi-Core Speedup** | `> 200 MiB` (up to tens of GiB) | High read count (4-line records) | `v4zip fastq compress in.fastq out.v4zm --chunk-size 32M` | Compresses up to 4 biological 4-line-aligned chunks concurrently (~1.8×–3× faster) and emits standalone `.partXXXX.v4z` files. |
| **Human / Mammalian Chromosomes or Assemblies** *(e.g. H. sapiens chr1 / whole genome)* | `50 Mbp – 3.2 Gbp` | Soft-masked repeats (`a/c/g/t`), large size | Per-chromosome: `v4zip compress chr1.fna chr1.v4z --split 25M`<br>Whole multi-record genome: `v4zip compress hg38.fna hg38.v4zm --chunk-size 64M` | Respects `1 GiB` default safety cap (`--max-output` on decompress) and avoids single-stream RAM/position limits (`2^30` bases per stream). |
| **Strict Transfer Limits** *(Email <= 20 MB, flaky FTP, multipart upload)* | Any | Any — **zero compression loss required** | `v4zip compress in.fna out.v4z --split 20M` | Global compression ratio is **100% preserved (0 B loss)** while producing exact-size `.001`, `.002`, ... volumes. |

---

## 2. Detailed Examples by Genome Type

### Case A: Small Bacterial or Viral Genome (`< 1.5 Mbp`)
- **Profile:** Compact coding sequence, low repeat density.
- **Complexity:** $O(N)$ streaming; ~1.5s for 587 KB (*M. genitalium*), ~500 MiB RAM.

```bash
# Compress
./build/v4zip compress m_genitalium.fna m_genitalium.v4z --verbose

# Decompress and verify bit-exact restoration
./build/v4zip decompress m_genitalium.v4z restored.fna --verbose
cmp m_genitalium.fna restored.fna
```

### Case B: Standard Bacterial / Eukaryotic Genome (`1.5 Mbp – 250 Mbp`)
- **Profile:** Large coding genome with local reverse-complement and palindromic motifs (*E. coli K-12*, *S. cerevisiae*, *C. elegans*, human chromosomes `chr1`–`chr22`).
- **Complexity:** $O(N)$ streaming; ~12s for 4.7 MB (*E. coli K-12*).

```bash
./build/v4zip compress e_coli_k12.fna e_coli_k12.v4z --verbose
./build/v4zip decompress e_coli_k12.v4z restored_ecoli.fna --verbose
```

### Case C: Pan-Genomes, Multi-Strain Collections & CRISPR Arrays
- **Profile:** Multiple strains of the same species (with 0.1%–5% SNPs) or genomes rich in CRISPR direct-repeat arrays (`20..120` bp period) and palindromic dyads.
- **What happens internally:** `GenomeProfiler` detects `repeat_16_rate >= 0.30`, `crispr_array_rate >= 0.08`, or `palindrome_dyad_rate >= 0.18` and routes to `RepeatDnaCodec` (4-way $D_4$ self-index + 64-entry Hotspot Index Table + `PredictorV3` for SNP residuals).
- **Important tip:** If you want **maximum cross-strain compression**, keep all strains in one global stream (use `--split` if you need smaller files for transfer, NOT `--chunk-size`, because `--chunk-size` isolates each chunk and prevents later chunks from referencing earlier strains!).

```bash
# Maximum pan-genome compression + 10 MiB transfer volumes (0 bytes lost)
./build/v4zip compress pan_genome_50strains.fna pan.v4z --split 10M --verbose

# Restore from the first volume (.001)
./build/v4zip decompress pan.v4z.001 restored_pan.fna --verbose
```

### Case D: Short-Read FASTQ Sequencing Runs (Single-File vs. Parallel Chunks)
- **Profile:** Four-line Illumina/short-read FASTQ files where quality scores dominate raw size.
- **Note:** FASTQ sequence and quality lines must have equal byte lengths (wrapped multi-line FASTQ is not supported).

```bash
# Option 1: Single file (best ratio for < 100 MiB FASTQ)
./build/v4zip fastq compress reads.fastq reads.v4zq --verbose
./build/v4zip fastq decompress reads.v4zq restored.fastq --verbose

# Option 2: Parallel 32 MiB chunks (faster multi-core encoding + resumable transfer)
./build/v4zip fastq compress big_reads.fastq big_reads.v4zm --chunk-size 32M --verbose
./build/v4zip fastq decompress big_reads.v4zm restored_big.fastq --verbose
```

### Case E: Human / Multi-GiB Genomes (Exceeding the Default `1 GiB` Safety Cap)
- **Why `--max-output` matters:** By default, `v4zip decompress` refuses to produce more than **1 GiB (`1073741824` bytes)** to protect against zip-bomb attacks.
- **Position index limit:** A single `PredictorV3` stream indexes matches up to $2^{30}$ bases (~1.07 Gbp). For a complete 3.1 Gbp human FASTA (`hg38.fna`) containing multiple chromosomes (`>chr1`, `>chr2`, ...), use **`--chunk-size 64M`** so each chromosome/chunk stays well within the index window and compresses in parallel.

```bash
# Compress a 3.1 GiB multi-chromosome human assembly in parallel 64 MiB chunks
./build/v4zip compress hg38.fna hg38.v4zm --chunk-size 64M --verbose

# Decompress passing --max-output (e.g. 4 GiB = 4294967296 bytes)
./build/v4zip decompress hg38.v4zm restored_hg38.fna --max-output 4294967296 --verbose
```

---

## 3. Complexity, Memory & Performance Guide

| Pipeline | Time Complexity | Peak RAM (Approx.) | Throughput Behavior |
|---|---|---|---|
| **`PredictorV3`** (`coding_small`, `coding_large`, `short_reads`, `low_complex`) | $O(N)$ | ~520 MiB + input/output buffers | CPU-bound on 37-expert arithmetic mixing (~0.35–0.55 MB/s per DNA stream). |
| **`RepeatDnaCodec`** (`repeat_heavy`) | $O(N)$ | ~520 MiB peak *(160 MiB 4-way hash table is freed before `PredictorV3` runs on residuals)* | Much faster on high-repeat pan-genomes (~1.7 MB/s on `mg_5x_snp1.fna`). |
| **`--split SIZE`** | Same as single file | Same as single file | Zero compression overhead; writes exact `SIZE` volumes sequentially. |
| **`--chunk-size SIZE`** | $O(N / \min(4, \text{chunks}))$ wall-clock | Up to `4 ×` active chunk models in parallel batches | ~1.8×–3.5× faster wall-clock time on multi-core CPUs; slight ratio trade-off (~0.2–1.5%) only if chunks are too small (`< 16M`). |

---

## 4. Data Transfer & Integrity Best Practices

### Choosing `--split` vs. `--chunk-size` for Transfer

- **Use `--split SIZE`** when:
  - You are sending via **email** (`--split 15M` or `--split 20M`) or **FTP / HTTP** with strict file-size caps.
  - You are compressing a **pan-genome** where cross-record redundancy matters most (0% compression loss).
  - You have a **single gigantic FASTA record** (e.g., a single 250 Mbp chromosome `chr1.fna`), since `--chunk-size` only splits at `>` record headers.
- **Use `--chunk-size SIZE`** when:
  - You are transferring a **multi-GiB FASTQ or multi-chromosome FASTA** over unreliable networks and want each `.partXXXX.v4z` file to be **independently valid, verifiable, and re-transmittable**.
  - You want **multi-core speedup** during compression.

### Transfer Verification Checklist

1. **Email / FTP with `--split`**:
   - Transfer all generated volumes: `archive.v4z.001`, `archive.v4z.002`, ..., `archive.v4z.NNN` into the same directory.
   - Run `v4zip decompress archive.v4z.001 restored.fna`. If any volume is missing or corrupted by 1 bit, the archive header/stream CRC32 fails closed immediately.
2. **Object Storage / Resumable FTP with `--chunk-size`**:
   - Transfer the manifest `dataset.v4zm` and all `dataset.part0001.v4z` .. `dataset.partNNNN.v4z` into the same directory.
   - Run `v4zip [fastq] decompress dataset.v4zm restored.fastq`. `v4zip` verifies every part's raw size and CRC32 as well as the global CRC32 before publishing the final output atomically.
3. **Never-Overwrite Guarantee**:
   - `v4zip` refuses to overwrite any existing destination file (including `.001`, `.partXXXX.v4z`, or `.v4zm`). Remove or rename old outputs explicitly before re-running.
