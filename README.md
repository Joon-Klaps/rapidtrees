# 🌲⚡ rapidtrees

<p align="center">
  <em>Blazingly fast pairwise phylogenetic tree distance calculations — Robinson–Foulds, Weighted RF, and Kuhner–Felsenstein — powered by Rust</em>
</p>

<p align="center">
  <a href="https://crates.io/crates/rapidtrees"><img src="https://img.shields.io/crates/v/rapidtrees.svg?style=flat-square&logo=rust&logoColor=white" alt="Crates.io" /></a>
  <a href="https://pypi.org/project/rapidtrees"><img src="https://img.shields.io/pypi/v/rapidtrees.svg?style=flat-square&logo=python&logoColor=white" alt="PyPI" /></a>
  <a href="https://github.com/Joon-Klaps/rapidtrees/actions/workflows/ci.yml"><img src="https://img.shields.io/github/actions/workflow/status/Joon-Klaps/rapidtrees/ci.yml?style=flat-square&logo=github&label=CI" alt="CI" /></a>
  <a href="https://codecov.io/gh/Joon-Klaps/rapidtrees"><img src="https://img.shields.io/codecov/c/github/Joon-Klaps/rapidtrees?style=flat-square&logo=codecov" alt="Coverage" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square" alt="License: MIT" /></a>
</p>

<p align="center">
  <a href="#%EF%B8%8F-overview">Overview</a> •
  <a href="#-installing">Installing</a> •
  <a href="#-usage">Usage</a> •
  <a href="#-python-api">Python API</a> •
  <a href="#-snap-format">Snap Format</a> •
  <a href="#%EF%B8%8F-benchmarks">Benchmarks</a>
</p>

---

## 🗺️ Overview

`rapidtrees` computes pairwise tree distances from [BEAST](https://beast.community/)/NEXUS `.trees` files or from precomputed `.snap` files and writes a labeled distance matrix. Three metrics are supported:

| Metric             | Flag                | Output  | Description                                   |
| ------------------ | ------------------- | ------- | --------------------------------------------- |
| Robinson–Foulds    | `--metric rf`       | integer | Symmetric difference of bipartitions          |
| Weighted RF        | `--metric weighted` | float   | Branch-length-weighted bipartition difference |
| Kuhner–Felsenstein | `--metric kf`       | float   | Euclidean distance on branch lengths          |

### ✨ Why `rapidtrees`?

- 🦀 **Rust core** — zero-overhead bitset operations with a cache-friendly memory layout
- 🔀 **Parallel by default** — powered by [`rayon`](https://github.com/rayon-rs/rayon), automatically scales across all cores
- 🐍 **Python bindings** — drop into any Python/NumPy workflow via [`PyO3`](https://pyo3.rs/)
- 📦 **No Rust toolchain required** — pre-built wheels on PyPI for Linux, macOS, and Windows
- 🗜️ **Gzip output** — stream directly to `.tsv.gz` without a separate compression step

### 🚀 Performance

Benchmarked on a ZIKA dataset (283 taxa · 4 000 trees · ~8 M comparisons):

| Metric             | Total time | Throughput                 |
| ------------------ | ---------- | -------------------------- |
| Robinson-Foulds    | ~3.5 s     | **~2.3 M comparisons/sec** |
| Weighted RF        | ~3.5 s     | **~2.3 M comparisons/sec** |
| Kuhner-Felsenstein | ~3.5 s     | **~2.3 M comparisons/sec** |

---

## 🔧 Installing

### 🐍 Python (PyPI) — _recommended_

Pre-built wheels for Linux, macOS, and Windows. No Rust toolchain needed.

```bash
pip install rapidtrees
```

### 🦀 CLI (crates.io)

Install the standalone command-line binary. Requires the [Rust toolchain](https://rustup.rs/).

```bash
cargo install rapidtrees
```

### 🛠️ From source

```bash
git clone https://github.com/Joon-Klaps/rapidtrees.git
cd rapidtrees
cargo build --release   # CLI binary → target/release/rapidtrees
pip install -e .        # Python bindings (requires Rust toolchain + maturin)
```

---

## 💻 Usage

```bash
rapidtrees \
  (--input <path/to/file.trees> | --snap-input <path/to/file.snap>) \
  --output <path/to/output.tsv[.gz]> \
  [--burnin-trees <N>] \
  [--burnin-states <STATE>] \
  [--use-real-taxa] \
  [--metric rf|weighted|kf] \
  [-q|--quiet]
```

| Flag                          | Description                                                        |
| ----------------------------- | ------------------------------------------------------------------ |
| `-i, --input <INPUT>`         | Path to BEAST `.trees` (NEXUS) file                                |
| `--snap-input <SNAP_INPUT>`   | Path to `.snap` file (currently supports only `--metric rf`)       |
| `-o, --output <OUTPUT>`       | Output path. Use `.gz` suffix for gzip compression; `-` for stdout |
| `-t, --burnin-trees <N>`      | Drop the first N trees (default: `0`)                              |
| `-s, --burnin-states <STATE>` | Keep only trees with `STATE > STATE` (default: `0`)                |
| `--use-real-taxa`             | Map numeric taxon IDs via the TRANSLATE block                      |
| `--metric <rf\|weighted\|kf>` | Distance metric (default: `rf`)                                    |
| `-q, --quiet`                 | Suppress progress messages (errors still go to stderr)             |

The output is a **square TSV matrix** where both the header row and first column contain tree names formatted as `<file_basename>_tree_STATE<state>`. Use `-o -` to write to stdout for easy piping.

### 💡 Examples

<details open>
<summary><strong>Compute RF matrix → gzipped file</strong></summary>

```bash
rapidtrees \
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

</details>

<details>
<summary><strong>Apply burn-in by tree count</strong></summary>

```bash
rapidtrees \
  -i tests/data/hiv1.trees \
  -o out/hiv1_rf.tsv \
  -t 2

# Reading in beast 0.003s
# Read in 162 taxons for 19 trees
# Creating tree bit snapshots 0.001s
# Determining distances using RF for 171 combinations
# Determining distances using RF 0.000s
# Writing to output 0.000s
```

</details>

<details>
<summary><strong>Compute RF matrix directly from a snapshot file</strong></summary>

```bash
rapidtrees \
  --snap-input out/hiv1.snap \
  -o out/hiv1_rf_from_snap.tsv.gz \
  --metric rf

# Read snap with 21 trees and 2193 bipartitions in 0.001s
# Determined distances using RF in 0.000s
# Writing to out/hiv1_rf_from_snap.tsv.gz in 0.000s
```

</details>

---

## 🐍 Python API

`rapidtrees` ships Python bindings for seamless integration into Python workflows.

```python
import rapidtrees as rtd

# Robinson-Foulds distances
tree_names, rf_matrix = rtd.pairwise_rf(
    paths=["file1.trees", "file2.trees"],
    burnin_trees=10,    # skip first 10 trees per file
    burnin_states=0,    # skip trees with STATE < 0
    use_real_taxa=True  # use TRANSLATE block when merging multiple files
)

# Weighted RF (considers branch lengths)
tree_names, wrf_matrix = rtd.pairwise_weighted_rf(
    paths=["file1.trees"],
    burnin_trees=10
)

# Kuhner-Felsenstein distances
tree_names, kf_matrix = rtd.pairwise_kf(
    paths=["file1.trees"],
    burnin_trees=10
)

print(f"Computed distances for {len(tree_names)} trees")
print(f"RF distance between tree 0 and 1: {rf_matrix[0][1]}")
```

---

## 📦 Snap Format

`rapidtrees` can export tree snapshots to a compressed binary `.snap` file for downstream analyses (ESS computation, ASDSF, convergence diagnostics) on HPC clusters without re-parsing the original `.trees` files. The CLI can also compute RF distance matrices directly from `.snap` files via `--snap-input`.

### What is a snapshot?

A **tree snapshot** is a compact bitset representation of a phylogenetic tree. Each bipartition (split) is encoded as a bitset over leaf indices. The full set of snapshots for a tree collection captures everything needed for RF-family distance computations and convergence diagnostics — without storing the original Newick strings.

> **Note:** Snapshots are not human-readable and are not intended for general interchange. They are an internal format optimized for fast distance calculations and cannot be convert to it's original newick-style format.

### File layout

A `.snap` file is a **gzip-compressed** binary stream with the following sections in order:

```
┌─────────────────────────────────────────────────────┐
│  HEADER                                              │
│  4 bytes  magic        "SNAP" (0x534E4150)          │
│  1 byte   version      format version (currently 1) │
│  8 bytes  n_trees      u64 LE — number of trees     │
│  8 bytes  n_taxa       u64 LE — number of leaf taxa │
│  8 bytes  n_bip        u64 LE — number of unique    │
│                         bipartitions across all trees│
├─────────────────────────────────────────────────────┤
│  TAXA NAMES                                          │
│  For each of n_taxa:                                 │
│    4 bytes  length     u32 LE — byte length of name │
│    N bytes  name       UTF-8 string                  │
├─────────────────────────────────────────────────────┤
│  TREE NAMES                                          │
│  For each of n_trees:                                │
│    4 bytes  length     u32 LE — byte length of name │
│    N bytes  name       UTF-8 string                  │
├─────────────────────────────────────────────────────┤
│  PRESENCE MATRIX                                     │
│  n_trees × n_bip bytes, row-major uint8             │
│  presence[i][j] = 1 if bipartition j is in tree i  │
│                   0 otherwise                        │
└─────────────────────────────────────────────────────┘
```

Bipartition column order is **deterministic**: columns are sorted in ascending `Bitset` order (lexicographic over `u64` words, i.e. by leaf-index bit pattern), so the same tree set always produces the same column indices regardless of parse order.

### What the presence matrix gives you

The presence matrix is sufficient for all major convergence diagnostics:

| Diagnostic              | What you need                    |
| ----------------------- | -------------------------------- |
| **Pseudo ESS**          | `presence.mean(axis=0)` per chain |
| **ASDSF**               | Per-chain split frequencies       |
| **Fréchet ESS**         | RF distances via XOR of rows      |
| **WRF / KF distances**  | Branch lengths *(future section)* |

Note: `sum(presence[i] XOR presence[j]) == RF(tree_i, tree_j)` exactly.

### Loading in Python

```python
import rapidtrees
import numpy as np

tree_names, taxa_names, n_bip, pres_bytes = rapidtrees.load_snap("run1.snap")
n = len(tree_names)
presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n, n_bip).copy()

# Global split frequencies (for Pseudo ESS / ASDSF)
global_freq = presence.mean(axis=0)

# RF distance between any two trees — no recomputation needed
rf_01 = int((presence[0] ^ presence[1]).sum())
```

---



## ⏱️ Benchmarks

Benchmarks were run on a MacBook Pro M1. Trees are parsed **once** and bitset snapshots are reused across all pairwise comparisons. Parallelism is provided by [`rayon`](https://github.com/rayon-rs/rayon) — no manual thread management needed.

<details>
<summary><strong>Show full benchmark table</strong></summary>

| Taxa (N) | Trees (T) | Combinations | Est. Memory | Actual Memory   | Wall Time        | CPU Time         |
| -------- | --------- | ------------ | ----------- | --------------- | ---------------- | ---------------- |
| 10       | 100       | 10.0K        | 51.56 KB    | 64.00 KB        | 391.71 µs        | 1.77 ms          |
| 10       | 1000      | 1.0M         | 515.62 KB   | 448.00 KB       | 1.16 ms          | 9.52 ms          |
| 10       | 10000     | 100.0M       | 5.04 MB     | 5.02 MB         | 81.86 ms         | 757.66 ms        |
| 10       | 100000    | 10.0B        | 50.35 MB    | 52.42 MB        | 8.17 s (est)     | 1.27 min (est)   |
| 100      | 100       | 10.0K        | 481.25 KB   | 464.00 KB       | 230.54 µs        | 1.45 ms          |
| 100      | 1000      | 1.0M         | 4.70 MB     | 1.17 MB         | 42.00 ms         | 147.22 ms        |
| 100      | 10000     | 100.0M       | 47.00 MB    | 38.42 MB        | 1.78 s           | 14.75 s          |
| 100      | 100000    | 10.0B        | 469.97 MB   | 451.47 MB       | 2.77 min (est)   | 25.29 min (est)  |
| 500      | 100       | 10.0K        | 4.59 MB     | 80.00 KB        | 1.61 ms          | 14.18 ms         |
| 500      | 1000      | 1.0M         | 45.90 MB    | 25.95 MB        | 165.58 ms        | 1.46 s           |
| 500      | 10000     | 100.0M       | 458.98 MB   | 399.44 MB       | 19.75 s          | 2.87 min         |
| 500      | 100000    | 10.0B        | 4.48 GB     | 4.51 GB         | 31.75 min (est)  | 294.98 min (est) |
| 1000     | 100       | 10.0K        | 15.27 MB    | 5.30 MB         | 4.69 ms          | 44.25 ms         |
| 1000     | 1000      | 1.0M         | 152.71 MB   | 126.77 MB       | 660.03 ms        | 5.30 s           |
| 1000     | 10000     | 100.0M       | 1.49 GB     | 1.44 GB         | 1.13 min         | 9.72 min         |
| 1000     | 100000    | 10.0B        | 14.91 GB    | 9.41 GB         | 111.42 min (est) | 983.70 min (est) |
| 2000     | 100       | 10.0K        | 54.94 MB    | 34.56 MB        | 23.11 ms         | 182.28 ms        |
| 2000     | 1000      | 1.0M         | 549.44 MB   | 480.02 MB       | 2.89 s           | 19.48 s          |
| 2000     | 10000     | 100.0M       | 5.37 GB     | 5.10 GB         | 3.92 min         | 34.76 min        |
| 2000     | 100000    | 10.0B        | 53.66 GB    | Skipped (>30GB) | -                | -                |
| 5000     | 100       | 10.0K        | 316.63 MB   | 299.31 MB       | 270.13 ms        | 2.09 s           |
| 5000     | 1000      | 1.0M         | 3.09 GB     | 3.06 GB         | 27.61 s          | 4.05 min         |
| 5000     | 10000     | 100.0M       | 30.92 GB    | Skipped (>30GB) | -                | -                |
| 5000     | 100000    | 10.0B        | 309.21 GB   | Skipped (>30GB) | -                | -                |

</details>

> **Note:** Weighted RF and KF produce floating-point matrices; RF produces integer matrices.

---

## 🔍 Troubleshooting

- **No trees parsed?** Verify the input is a valid NEXUS `.trees` file and adjust `--burnin-*` settings.
- **Piping to other tools?** Use `-q` to suppress timing messages on stdout.
- **Gzipped output not working?** Ensure the output filename ends with `.gz`.

---

## ⚖️ License

`rapidtrees` is provided under the [MIT License](LICENSE).
