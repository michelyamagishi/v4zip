# v4zip Architecture & Compression Methods

`v4zip` is a lossless FASTA/FASTQ compressor grounded in **Klein four-group ($V_4$) / dihedral group ($D_4$)** DNA symmetries. Rather than applying a single monolithic model to every input, it profiles biological structure in a fast pre-pass and dispatches each DNA stream to a specialized lossless pipeline.

## 1. Biological Profiling (`GenomeProfiler`)

Implemented in [`src/detail/genome_profiler.hpp`](../src/detail/genome_profiler.hpp), the pre-pass inspects universal sequence statistics in $O(N)$ time (with stride-4 canonical 16-mer sampling):

- **GC content & homopolymer run density** (`run_rate`)
- **3-base codon reading-frame periodicity** (`codon_bias`)
- **Record fragmentation** (`avg_record_len`)
- **Canonical 16-mer repeat density** (`repeat_16_rate`)
- **Palindromic dyad symmetry rate** (`palindrome_dyad_rate`): detects mirror (`rev(w) == w`) and reverse-complement hairpin/stem-loop (`(rev(w) ^ mask) == w`) self-symmetry on 16-mers.
- **CRISPR array periodicity rate** (`crispr_array_rate`): detects quasi-periodic canonical 16-mer repeats within the biological CRISPR direct-repeat + spacer window (`20..120` bp).

### Strategy Dispatch Table

| Strategy | Trigger Condition | Primary Pipeline |
|---|---|---|
| `low_complex` | `run_rate >= 0.25` | `PredictorV3` (high run-length smoothing) |
| `repeat_heavy` | `repeat_16_rate >= 0.30` or `crispr_array_rate >= 0.08` or `palindrome_dyad_rate >= 0.18` | `RepeatDnaCodec` (4-way $D_4$ LZ + Hotspot Index Table + `PredictorV3` residuals) |
| `short_reads` | Multi-record and `avg_record_len <= 250` | `PredictorV3` (`k0 = 9`, fast local weight decay) |
| `coding_small` | `< 1.5 Mbp` | `PredictorV3` (compact APM) |
| `coding_large` | `>= 1.5 Mbp` | `PredictorV3` (full 3072-context APM) |

## 2. `PredictorV3` (37-Expert Log-Linear/Linear Mixture + 4-ary APM)

Implemented in [`src/detail/predictor_v3.hpp`](../src/detail/predictor_v3.hpp):

- **Dense orders `0..9`** (`10` experts)
- **Klein $V_4$ `OrbitTableV3` orders `{10, 12, 15, 18, 22}` × 3 views** (`raw`, `rc`, `pooled` = `15` experts)
- **Codon phase (`pos % 3`) + wobble experts** (`4` experts)
- **Sparse skip1, skip2, gapped & run-length experts** (`4` experts)
- **`MatchModelV3` 4-way multi-scale repeat matcher** (`1` expert): tracks 4 recent positions per canonical $D_4$ bucket across 4 scales `{k0, 13, 16, 24}` (where `k = 13` captures CRISPR/palindrome core motifs).
- **Direct-indexed collision-free Palindrome / Reverse-Complement tables (`rc8_`, `rc9_`, `rc10_`)** (`3` experts): share reverse-complement evidence on orders 7, 8, and 9 with zero hash collisions and no double-counting on self-complementary palindromic dyads (`w == rc`).
- **4-ary Adaptive Probability Map (APM)**: post-calibrates the mixture distribution before arithmetic range coding.

## 3. `RepeatDnaCodec` (4-Way $D_4$ LZ + CRISPR/Palindrome Hotspot Index Table)

Implemented in [`src/detail/repeat_codec.hpp`](../src/detail/repeat_codec.hpp):

- **4-Way Canonical $D_4$ Hash Index**: indexes up to 4 recent positions per canonical 16-mer (`forward`, `reverse`, `complement`, `reverse-complement`), preventing bucket eviction across multi-strain pan-genomes and inverted repeats.
- **CRISPR & Palindrome Hotspot Index Table (Format 4)**:
  - Pre-scans palindromic dyad 16-mers (`count >= 3` MDL gate, capped at `64` entries).
  - Emits compact 2-byte dictionary tokens (`kind = 5 + transform`, `varint(dict_idx)`) for short (`16..17` bp) repeats below the standard LZ `MIN_MATCH = 18` threshold.
  - Frees the 160 MiB hash table (`table.clear(); table.shrink_to_fit()`) before allocating `PredictorV3` for residual literals.
  - Enforces fail-closed decoder caps (`dict_size <= 64`, `dict[d] <= 0xFFFFFFFF`, `tok_raw_size <= expected_dna * 2 + 64`).

## 4. `QualityModel` v4 & Multi-Core Stream Execution

- **`QualityModel`** ([`src/detail/quality_coder.hpp`](../src/detail/quality_coder.hpp)): 8-expert adaptive range coder for FASTQ quality scores (`o0`, `o1`, collision-free `o2` `64x64`, 2-way `o3`, `pos x p1`, `pos x p1/2 x p2/2`, `run`, `mean`).
- **Stream Parallelism** ([`src/archive2.cpp`](../src/archive2.cpp)): `dna`, `quality`, `identifiers`/`separators`, and layout streams are compressed/decompressed concurrently via `std::async`, paired with multi-threaded `zstd` level 19.
