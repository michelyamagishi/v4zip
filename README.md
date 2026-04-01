# v4zip

Lossless DNA/FASTQ compressor based on Klein four-group (V4) symmetry. Exploits Chargaff's second parity rule to share statistics across reverse-complement contexts, achieving better compression than general-purpose tools on genomic data.

## Building

Requires a C++17 compiler and CMake 3.14+.

### Linux (GCC)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

### macOS (Clang)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(sysctl -n hw.ncpu)
```

### Windows (MSVC)

```powershell
cmake -B build
cmake --build build --config Release
```

The binary is at `build/v4zip` (Linux/macOS) or `build\Release\v4zip.exe` (Windows).

## Usage

```bash
# FASTA
v4zip compress input.fasta output.v4z
v4zip decompress output.v4z restored.fasta

# FASTQ
v4zip fastq compress input.fastq output.v4zq
v4zip fastq decompress output.v4zq restored.fastq
```

Options: `-o N` (context order, default 16), `--verbose`, `--stats` (FASTQ stream breakdown).

## Benchmarks

### FASTA (bits per ACGT base)

| Genome | v4zip | GeCo3-l9 | gzip-9 | xz-9 |
|--------|-------|----------|--------|------|
| M. genitalium (580 KB) | 1.83 | 1.79 | 2.29 | 2.11 |
| E. coli K-12 (4.6 MB) | 1.94 | 1.88 | 2.38 | 2.17 |
| S. cerevisiae (12 MB) | 1.91 | 1.82 | 2.40 | 2.13 |
| H. sapiens chr22 (34 MB) | 1.78 | 1.60 | 2.28 | 1.97 |
| C. elegans (100 MB) | 1.85 | 1.18 | 2.53 | 2.14 |
| H. sapiens chr1 (249 MB) | 1.78 | 1.57 | 2.30 | 1.93 |

### FASTQ (Illumina 150bp, 621K reads, 231 MB)

| Compressor | Size | Ratio |
|------------|------|-------|
| v4zip | 34.6 MB | 6.7x |
| gzip-9 | 39.8 MB | 5.8x |
| bzip2-9 | 33.5 MB | 6.9x |
| xz-9 | 26.1 MB | 8.8x |
| zstd-19 | 31.1 MB | 7.4x |

## Testing

```bash
cd build && ctest --output-on-failure
```

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
