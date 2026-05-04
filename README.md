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

Benchmarked on a ZIKA dataset (283 taxa · 4 000 trees · ~8 M comparisons) (non-gzipped output):

| Metric             | Total time | Throughput                 |
| ------------------ | ---------- | -------------------------- |
| Robinson-Foulds    | ~2.7 s     | **~3.0 M comparisons/sec** |
| Weighted RF        | ~3.4 s     | **~2.3 M comparisons/sec** |
| Kuhner-Felsenstein | ~3.3 s     | **~2.3 M comparisons/sec** |

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

#### Prerequisites
- [Rust toolchain](https://rustup.rs/) — for building the Rust core
- [pixi](https://pixi.sh/) — for managing Python and R dependencies

#### Setup

```bash
git clone https://github.com/Joon-Klaps/rapidtrees.git
cd rapidtrees
# Set up environment — installs Python, R (with phangorn), and builds the package
pixi install
```

#### Development tasks

```bash
# Run Python API tests (includes R/phangorn cross-validation)
pixi run test-python

# Run Rust unit tests
pixi run test-rust
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

`rapidtrees` exposes four functions from its Rust core. All accept a Python **iterator** of newick strings, keeping memory constant regardless of tree count.

| Function | Returns |
| --- | --- |
| `pairwise_rf_from_newick_iter` | `(names, bytes)` — RF matrix as flat `uint32` bytes, row-major |
| `pairwise_rf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes)` — RF + bipartition presence matrix |
| `pairwise_wrf_from_newick_iter` | `(names, list[float])` — Weighted RF, flat row-major |
| `pairwise_kf_from_newick_iter` | `(names, list[float])` — Kuhner-Felsenstein, flat row-major |

```python
import rapidtrees as rtd
import numpy as np

trees = [
    "(A:0.1,(B:0.1,C:0.1):0.1);",
    "(A:0.1,(C:0.1,B:0.1):0.1);",
    "((A:0.1,B:0.1):0.1,C:0.1);",
]
names = ["t1", "t2", "t3"]

tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
    names, iter(trees), [{}], [0, 0, 0]
)
rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(len(tree_names), -1)
```

For BEAST `.trees` files, translate maps, the snapshot API, and multi-file usage see **[docs/python-api.md](docs/python-api.md)**.

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

### Presence matrix from Python

Snap files are written and read by the CLI only. From Python, use
`pairwise_rf_with_snapshots_from_newick_iter` to obtain the same presence
matrix in memory without writing a file — see the [Python API section](#-python-api) above.

```python
import rapidtrees as rtd
import numpy as np

tree_names, rf_bytes, leaf_names, n_bip, pres_bytes = (
    rtd.pairwise_rf_with_snapshots_from_newick_iter(
        names, iter(newicks), translate_maps, map_indices
    )
)
n = len(tree_names)
presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n, n_bip).copy()

# Global split frequencies (for Pseudo ESS / ASDSF)
global_freq = presence.mean(axis=0)

# RF distance between any two trees — no recomputation needed
rf_01 = int((presence[0].astype(int) ^ presence[1].astype(int)).sum())
```

---



## ⏱️ Benchmarks

Benchmarks were run on a MacBook Pro M1. Trees are parsed **once** and bitset snapshots are reused across all pairwise comparisons. Parallelism is provided by [`rayon`](https://github.com/rayon-rs/rayon) — no manual thread management needed.

<details>
<summary><strong>Show full benchmark table</strong></summary>

> Output of `cargo bench --bench memory_time_benchmark`

| Taxa (N) | Trees (T) | Combinations | Est. Memory | Actual Memory | Wall Time | CPU Time |
|----------|-----------|--------------|-------------|---------------|-----------|----------|
| 10       | 100       | 10.0K        | 13.96 KB    | 17.87 KB      | 172.29 µs | 733.00 µs |
| 10       | 1000      | 1.0M         | 139.65 KB   | 168.95 KB     | 5.94 ms   | 12.99 ms |
| 10       | 10000     | 100.0M       | 1.36 MB     | 1.64 MB       | 717.43 ms | 1.36 s   |
| 10       | 100000    | 10.0B        | 13.64 MB    | 16.40 MB      | 3.53 min  | 4.04 min |
| 100      | 100       | 10.0K        | 133.59 KB   | 136.09 KB     | 244.42 µs | 1.75 ms  |
| 100      | 1000      | 1.0M         | 1.30 MB     | 1.21 MB       | 17.49 ms  | 50.99 ms |
| 100      | 10000     | 100.0M       | 13.05 MB    | 11.95 MB      | 1.08 s    | 4.81 s   |
| 100      | 100000    | 10.0B        | 130.46 MB   | 119.41 MB     | 4.29 min  | 9.86 min |
| 500      | 100       | 10.0K        | 706.84 KB   | 713.93 KB     | 2.02 ms   | 3.08 ms  |
| 500      | 1000      | 1.0M         | 6.90 MB     | 5.89 MB       | 27.77 ms  | 216.87 ms |
| 500      | 10000     | 100.0M       | 69.03 MB    | 57.84 MB      | 2.82 s    | 21.30 s  |
| 500      | 100000    | 10.0B        | 690.27 MB   | 577.28 MB     | 7.44 min  | 37.13 min |
| 1000     | 100       | 10.0K        | 1.50 MB     | 1.51 MB       | 873.75 µs | 5.63 ms  |
| 1000     | 1000      | 1.0M         | 15.03 MB    | 11.86 MB      | 51.22 ms  | 420.94 ms |
| 1000     | 10000     | 100.0M       | 150.32 MB   | 115.30 MB     | 5.12 s    | 42.92 s  |
| 1000     | 100000    | 10.0B        | 1.47 GB     | 1.12 GB       | 11.26 min | 74.74 min |
| 2000     | 100       | 10.0K        | 3.50 MB     | 3.51 MB       | 1.14 ms   | 9.92 ms  |
| 2000     | 1000      | 1.0M         | 35.01 MB    | 24.15 MB      | 93.73 ms  | 838.01 ms |
| 2000     | 10000     | 100.0M       | 350.06 MB   | 230.59 MB     | 9.42 s    | 1.38 min |
| 2000     | 100000    | 10.0B        | 3.42 GB     | 2.24 GB       | 18.57 min | 141.11 min |
| 5000     | 100       | 10.0K        | 14.30 MB    | 12.43 MB      | 61.36 ms  | 83.85 ms |
| 5000     | 1000      | 1.0M         | 143.02 MB   | 63.97 MB      | 224.40 ms | 2.09 s   |
| 5000     | 10000     | 100.0M       | 1.40 GB     | 579.40 MB     | 22.45 s   | 3.44 min |             |

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
