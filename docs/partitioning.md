# Partitioning Strategies (`--split` & `--chunk-size`)

For very large genomes (multi-GiB human assemblies, metagenomes, or deep FASTQ runs) or constrained transfer channels (email, FTP, object-storage multipart upload), `v4zip` provides two orthogonal, lossless partitioning strategies implemented in [`src/detail/partition.hpp`](../src/detail/partition.hpp) and [`src/main.cpp`](../src/main.cpp).

The two flags are **mutually exclusive** on `compress` (`v4zip` fails closed if both are supplied).

## Strategy 1: Exact Volume Split (`--split SIZE`)

Compresses the entire input as a single global archive (preserving **100% of global compression ratio, 0 bytes lost**) and slices the resulting `.v4z` bitstream into exact-size volumes: `OUTPUT.001`, `OUTPUT.002`, ..., `OUTPUT.NNN`.

- **Best for:** Maximum compression ratio + strict file-size limits (email attachments, FTP, FAT32/USB).
- **Size suffixes:** Decimal bytes or `K` (1024), `M` (1 MiB), `G` (1 GiB).

```bash
# Compress into exact 20 MiB volumes
./build/v4zip compress genome.fna genome.v4z --split 20M --verbose

# Decompress automatically detects .001, .002, ...
./build/v4zip decompress genome.v4z.001 restored.fna --verbose
```

## Strategy 2: Biological Chunked Manifest (`--chunk-size SIZE`)

Splits the input along **lossless biological record boundaries** — `>` headers for FASTA and exact **4-line read blocks** for FASTQ — and compresses chunks in parallel (batches of up to 4 concurrent `std::async` workers) into standalone `.part0001.v4z`, `.part0002.v4z`, ... archives indexed by a `V4ZM1` manifest (`.v4zm`).

- **Best for:** Multi-GiB genomes, multi-core compression speedup, independent per-chunk verification/retransmission.
- **Note:** The output path must end with `.v4zm`.

```bash
# Compress FASTQ in parallel 32 MiB 4-line-aligned chunks
./build/v4zip fastq compress reads.fastq reads.v4zm --chunk-size 32M --verbose

# Decompress automatically reads manifest, verifies every part CRC32 + global CRC32
./build/v4zip fastq decompress reads.v4zm restored.fastq --verbose
```

### `V4ZM1` Manifest Structure

```text
V4ZM1 <fastq_flag:0|1> <total_raw_bytes> <total_raw_crc32> <part_count>
<part_filename_1> <raw_bytes_1> <raw_crc32_1>
<part_filename_2> <raw_bytes_2> <raw_crc32_2>
...
```

## Comparison

| Property | `--split SIZE` | `--chunk-size SIZE` |
|---|---|---|
| Compression ratio vs. single file | **Identical (0 B loss)** | Slight boundary overhead (~0.2–1.5% on small chunks; negligible at `>= 32M`) |
| Output file sizes | Exact `SIZE` (last volume `<= SIZE`) | Determined by compressed chunk size (~4× smaller than `SIZE`) |
| Multi-core chunk parallelism | Stream-level only | **Up to 4 chunks concurrently** |
| Standalone parts | No (requires all volumes) | **Yes** (each `.partXXXX.v4z` is a valid standalone `.v4z` archive) |
| Integrity checks | Global archive header + stream CRC32 | Per-part CRC32 + global raw CRC32 + path-traversal guard |
