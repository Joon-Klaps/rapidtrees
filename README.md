# rust-python-tree-distances CLI

Compute pairwise tree distances (Robinson–Foulds, weighted RF, Kuhner–Felsenstein) from BEAST/NEXUS `.trees` files and write a labeled distance matrix.

## Performance ZIKA dataset (283 taxa, 4000 trees, ~8M comparisons)
- **Robinson-Foulds (RF)**: ~3.5s total → **~2.3M comparisons/sec**
- **Weighted RF**: ~3.5s total → **~2.3M comparisons/sec**
- **Kuhner-Felsenstein (KF)**: ~3.5s total → **~2.3M comparisons/sec**

## Install / Build

Requirements: Rust toolchain (stable). Then build the binary:

> [!NOTE]
> Install the Rust toolchain from https://rustup.rs/ if you don't have it yet.
>
> ```bash
> curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
> ```
>

Clone the repository:
```bash
git clone https://github.com/Joon-Klaps/rust-python-tree-distances.git
```

Build the project:
```bash
cd rust-python-tree-distances
cargo build --release
```

The executable will be at `target/release/rust-python-tree-distances`.

## Usage

```bash
# if rust-python-tree-distances is not in your PATH, use the full path, e.g. ./target/release/rust-python-tree-distances
rust-python-tree-distances \
  --input <path/to/file.trees> \
  --output <path/to/output.tsv[.gz]> \
  [--burnin-trees <N>] \
  [--burnin-states <STATE>] \
  [--use-real-taxa] \
  [--metric rf|weighted|kf] \
  [-q|--quiet]
```

Flags and options:

- `-i, --input <INPUT>`: Path to BEAST `.trees` (NEXUS) file.
- `-o, --output <OUTPUT>`: Output path for the TSV distance matrix. If the path ends with `.gz` it will be gzip-compressed. Use `-` to write to stdout (uncompressed).
- `-t, --burnin-trees <N>`: Drop the first N trees (default: 0).
- `-s, --burnin-states <STATE>`: Keep only trees with `STATE_ > STATE` (default: 0).
- `--use-real-taxa`: Map numeric taxon IDs to labels using the TRANSLATE block if present.
- `--metric <rf|weighted|kf>`: Choose the distance metric (default: `rf`), weighted referring to weighted RF.
- `-q, --quiet`: Suppress progress messages on stdout. Errors still go to stderr.

The output is a square TSV matrix where both the header row and the first column are tree names, constructed as `<file_basename>_tree_STATE<state>`. When writing to stdout (`-o -`), the matrix is printed to stdout, allowing easy piping.

## Examples

- Compute RF matrix and write to gzipped file:

```bash
rust-python-tree-distances \
  -i tests/data/hiv1.trees \
  -o out/hiv1_rf.tsv.gz \
  --metric rf

# Reading in beast 0.003s
# Read in 162 taxons for 21 trees
# Creating tree bit snapshots 0.002s
# Determining distances using RF for 210 combinations
# Determining distances using RF 0.000s
# Writing to output 0.000s
```

- Apply burn-in by tree count and state:

```bash
rust-python-tree-distances \                                                                                                      2 ↵
  -i tests/data/hiv1.trees \
  -o out/hiv1_rf.tsv \
  -t 2 \

# Reading in beast 0.003s
# Read in 162 taxons for 19 trees
# Creating tree bit snapshots 0.001s
# Determining distances using RF for 171 combinations
# Determining distances using RF 0.000s
# Writing to output 0.000s
```

## Performance notes

- Trees are parsed once. Bitset snapshots are built once and reused for pairwise comparisons. Parallelism is provided by `rayon`.
- Weighted RF and KF produce floating-point matrices; RF produces integer matrices.

## Benchmarks

We provide a benchmark suite to evaluate both memory usage and runtime performance for varying dataset sizes.

To run the benchmark:

```bash
cargo bench --bench memory_time_benchmark
```

This will output a table showing estimated memory usage, wall-clock time, and CPU time for pairwise distance calculations.

Various output (on a MacBook Pro M1):

| Taxa (N) | Trees (T) | Combinations | Est. Memory | Actual Memory | Wall Time | CPU Time |
|----------|-----------|--------------|-------------|---------------|-----------|----------|
| 10       | 100       | 10.0K        | 51.56 KB    | 64.00 KB      | 391.71 µs | 1.77 ms  |
| 10       | 1000      | 1.0M         | 515.62 KB   | 448.00 KB     | 1.16 ms   | 9.52 ms  |
| 10       | 10000     | 100.0M       | 5.04 MB     | 5.02 MB       | 81.86 ms  | 757.66 ms |
| 10       | 100000    | 10.0B        | 50.35 MB    | 52.42 MB      | 8.17 s (est) | 1.27 min (est) |
| 100      | 100       | 10.0K        | 481.25 KB   | 464.00 KB     | 230.54 µs | 1.45 ms  |
| 100      | 1000      | 1.0M         | 4.70 MB     | 1.17 MB       | 42.00 ms  | 147.22 ms |
| 100      | 10000     | 100.0M       | 47.00 MB    | 38.42 MB      | 1.78 s    | 14.75 s  |
| 100      | 100000    | 10.0B        | 469.97 MB   | 451.47 MB     | 2.77 min (est) | 25.29 min (est) |
| 500      | 100       | 10.0K        | 4.59 MB     | 80.00 KB      | 1.61 ms   | 14.18 ms |
| 500      | 1000      | 1.0M         | 45.90 MB    | 25.95 MB      | 165.58 ms | 1.46 s   |
| 500      | 10000     | 100.0M       | 458.98 MB   | 399.44 MB     | 19.75 s   | 2.87 min |
| 500      | 100000    | 10.0B        | 4.48 GB     | 4.51 GB       | 31.75 min (est) | 294.98 min (est) |
| 1000     | 100       | 10.0K        | 15.27 MB    | 5.30 MB       | 4.69 ms   | 44.25 ms |
| 1000     | 1000      | 1.0M         | 152.71 MB   | 126.77 MB     | 660.03 ms | 5.30 s   |
| 1000     | 10000     | 100.0M       | 1.49 GB     | 1.44 GB       | 1.13 min  | 9.72 min |
| 1000     | 100000    | 10.0B        | 14.91 GB    | 9.41 GB       | 111.42 min (est) | 983.70 min (est) |
| 2000     | 100       | 10.0K        | 54.94 MB    | 34.56 MB      | 23.11 ms  | 182.28 ms |
| 2000     | 1000      | 1.0M         | 549.44 MB   | 480.02 MB     | 2.89 s    | 19.48 s  |
| 2000     | 10000     | 100.0M       | 5.37 GB     | 5.10 GB       | 3.92 min  | 34.76 min |
| 2000     | 100000    | 10.0B        | 53.66 GB    | Skipped (>30GB) | -         | -        |
| 5000     | 100       | 10.0K        | 316.63 MB   | 299.31 MB     | 270.13 ms | 2.09 s   |
| 5000     | 1000      | 1.0M         | 3.09 GB     | 3.06 GB       | 27.61 s   | 4.05 min |
| 5000     | 10000     | 100.0M       | 30.92 GB    | Skipped (>30GB) | -         | -        |
| 5000     | 100000    | 10.0B        | 309.21 GB   | Skipped (>30GB) | -         | -        |

*Note: Times for large T are extrapolated from a smaller sample of comparisons.*

## Troubleshooting

- If no trees are parsed, verify the input is a valid NEXUS `.trees` file and adjust `--burnin-*` settings.
- Use `-q` when writing to stdout and piping to other tools to suppress timing messages.
- For gzipped output, ensure the output filename ends with `.gz`.

---

## Python API

The package also provides Python bindings for easy integration into Python workflows.

### Installation

From source (requires Rust):

```bash
# Using pip (recommended)
pip install -e .

# Or using maturin directly
pip install maturin
maturin develop --release
```

### Quick Start

```python
import rust_python_tree_distances as rtd

# Compute Robinson-Foulds distances
tree_names, rf_matrix = rtd.pairwise_rf(
    paths=["file1.trees", "file2.trees"],
    burnin_trees=10,  # Skip first 10 trees from each file
    burnin_states=0,   # Skip trees with STATE < 0
    use_real_taxa=True # Use TRANSLATE block if available, set to true if multiple files are provided
)

# Compute Weighted RF distances (considers branch lengths)
tree_names, wrf_matrix = rtd.pairwise_weighted_rf(
    paths=["file1.trees"],
    burnin_trees=10
)

# Compute Kuhner-Felsenstein distances
tree_names, kf_matrix = rtd.pairwise_kf(
    paths=["file1.trees"],
    burnin_trees=10
)

# Output is a list of tree names and a 2D distance matrix
print(f"Computed distances for {len(tree_names)} trees")
print(f"RF distance between tree 0 and 1: {rf_matrix[0][1]}")
```
