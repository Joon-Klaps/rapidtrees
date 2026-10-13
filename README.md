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
  <a href="#-presence-matrix">Presence Matrix</a> •
  <a href="#%EF%B8%8F-benchmarks">Benchmarks</a>
</p>

---

## 🗺️ Overview

`rapidtrees` computes pairwise tree distances from [BEAST](https://beast.community/)/NEXUS `.trees` files or plain multi-tree Newick files, and writes a labeled distance matrix. Three metrics are supported:

| Metric             | Flag                | Output  | Description                                   |
| ------------------ | ------------------- | ------- | --------------------------------------------- |
| Robinson–Foulds    | `--metric rf`       | integer | Symmetric difference of bipartitions          |
| Weighted RF        | `--metric weighted` | float   | Branch-length-weighted bipartition difference |
| Kuhner–Felsenstein | `--metric kf`       | float   | Euclidean distance on branch lengths          |

### ✨ Why `rapidtrees`?

- 🦀 **Rust core** — splits as 128-bit fingerprints, compared as `u32` IDs in a cache-friendly layout
- 🔀 **Parallel by default** — powered by [`rayon`](https://github.com/rayon-rs/rayon), automatically scales across all cores
- 🐍 **Python bindings** — drop into any Python/NumPy workflow via [`PyO3`](https://pyo3.rs/)
- 📦 **No Rust toolchain required** — pre-built wheels on PyPI for Linux, macOS, and Windows
- 🗜️ **Gzip output** — stream directly to `.tsv.gz` without a separate compression step

### 🚀 Performance

Measured on an Apple M1 Max laptop (10 cores, 64 GB RAM) with rapidtrees 0.11.0 from crates.io, on four published BEAST posteriors from the BEAST workshop tutorials ([YFV](https://beast.community/workshop_continuous_diffusion_yfv), [WNV](https://beast.community/workshop_continuous_diffusion_wnv), [batRABV](https://beast.community/workshop_discrete_diffusion_rabv), [EBOV](https://beast.community/workshop_discrete_diffusion_ebov)). Wall-clock time for the whole run, reading the NEXUS file included, with all 10 cores and with one thread in brackets; median of three runs after a warm-up, matrix written to `/dev/null`:

| Dataset | Taxa · trees | Pairs   | Robinson-Foulds | Weighted RF     | Kuhner-Felsenstein | Peak memory |
| ------- | ------------ | ------- | --------------- | --------------- | ------------------ | ----------- |
| YFV     | 65 · 5 001   | 12.5 M  | 0.12 s (0.56 s) | 0.41 s (2.38 s) | 0.37 s (2.21 s)    | 393 MiB     |
| WNV     | 104 · 1 001  | 0.50 M  | 0.02 s (0.07 s) | 0.04 s (0.15 s) | 0.04 s (0.14 s)    | 85 MiB      |
| batRABV | 372 · 1 001  | 0.50 M  | 0.05 s (0.15 s) | 0.07 s (0.27 s) | 0.06 s (0.27 s)    | 157 MiB     |
| EBOV    | 836 · 1 167  | 0.68 M  | 0.07 s (0.30 s) | 0.11 s (0.70 s) | 0.10 s (0.54 s)    | 159 MiB     |

On YFV that is about **100 M Robinson-Foulds comparisons per second**. Command: `RAYON_NUM_THREADS=<n> rapidtrees --input <file>.trees --metric <rf|weighted|kf> --output /dev/null --use-real-taxa --quiet`, no burn-in, unrooted.

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
  --input <path/to/file.trees> \
  --output <path/to/output.tsv[.gz]> \
  [--burnin-trees <N>] \
  [--burnin-states <STATE>] \
  [--use-real-taxa] \
  [--metric rf|weighted|kf] \
  [-q|--quiet]
```

| Flag                          | Description                                                        |
| ----------------------------- | ------------------------------------------------------------------ |
| `-i, --input <INPUT>`         | Path to a NEXUS `.trees` or plain Newick file (auto-detected)      |
| `-o, --output <OUTPUT>`       | Output path. Use `.gz` suffix for gzip compression; `-` for stdout |
| `-t, --burnin-trees <N>`      | Drop the first N trees (default: `0`)                              |
| `-s, --burnin-states <STATE>` | Keep only trees with `STATE > STATE` (default: `0`)                |
| `--use-real-taxa`             | Map numeric taxon IDs via the TRANSLATE block                      |
| `--metric <rf\|weighted\|kf>` | Distance metric (default: `rf`)                                    |
| `-q, --quiet`                 | Suppress progress messages (errors still go to stderr)             |

**Input formats.** `--input` takes either a NEXUS/BEAST trees file or a plain Newick file with one tree per line; the format is sniffed from the contents, not the extension. Newick files name nothing, carry no `TRANSLATE` block and no `STATE_` labels, so their trees are named after the line they start on, and `--burnin-states` / `--use-real-taxa` have nothing to act on (`--burnin-trees` still works).

The output is a **square TSV matrix** where both the header row and first column contain tree names: `<file_basename>_<tree_name>` for NEXUS (e.g. `hiv1_STATE_10000`) and `<file_basename>_line<n>` for Newick (e.g. `hiv1_line3`). Use `-o -` to write to stdout for easy piping.

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

---


### 🎛️ Tuning the kernels

Every metric keeps a dense column only for splits held by a large enough share of the trees, and gives the rarer ones posting lists. The defaults were tuned on simulated posteriors. `RAPIDTREES_RF_DENSE_SHARE` (default `0.03`) and `RAPIDTREES_WEIGHTED_DENSE_SHARE` (default `0.25`) override them, for re-tuning on other data or hardware. They change the speed only, never a distance.


## 🐍 Python API

`rapidtrees` exposes six functions and the `ProgressCounter` class from its Rust core. All six functions accept a Python **iterator** of newick strings, keeping memory constant regardless of tree count, and share the signature `(names, newick_iter, translate_maps, map_indices, rooted=False, progress=None)`.

| Function | Returns |
| --- | --- |
| `pairwise_rf_from_newick_iter` | `(names, bytes)` — RF matrix as flat `uint32` bytes, row-major |
| `pairwise_rf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — RF matrix + presence matrix + clade bitmasks |
| `pairwise_wrf_from_newick_iter` | `(names, list[float])` — Weighted RF, flat row-major |
| `pairwise_wrf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — wRF matrix + branch-length matrix + clade bitmasks |
| `pairwise_kf_from_newick_iter` | `(names, list[float])` — Kuhner-Felsenstein, flat row-major |
| `pairwise_kf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — KF matrix + branch-length matrix + clade bitmasks |

```python
import rapidtrees as rtd
import numpy as np

trees = [
    "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
    "((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);",
    "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.2);",
]
names = ["t1", "t2", "t3"]

# Newick strings already use real taxon names, so pass one empty translate map for every tree
tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
    names, iter(trees), [{}], [0] * len(trees)
)
rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(len(tree_names), -1)
# [[0 2 0]
#  [2 0 2]
#  [0 2 0]]
```

For BEAST `.trees` files, translate maps, the snapshot API, and multi-file usage see **[docs/python-api.md](docs/python-api.md)**.

---

## 📦 Presence Matrix

Alongside a distance matrix, `rapidtrees` can hand back the **presence matrix**: which bipartition appears in which tree. That is everything the convergence diagnostics need (ESS, ASDSF, split frequencies), and it comes out of the same single parse as the distances.

> **Splits are identified by a 128-bit fingerprint,** not by comparing leaf sets — that is what keeps building a tree linear in its taxon count rather than quadratic. Two distinct splits are merged if their fingerprints collide, with probability about `e²/2¹²⁹` for `e` distinct splits in the run: `1.5 × 10⁻²³` at a hundred million splits, some nineteen orders of magnitude below the rate at which the machine's own RAM flips a bit unnoticed. Every run prints its own `e` and bound.

Bipartition column order is **deterministic**: columns are sorted in ascending `Bitset` order (lexicographic over `u64` words, i.e. by leaf-index bit pattern), so the same tree set always produces the same column indices regardless of parse order.

### What the presence matrix gives you

The presence matrix is sufficient for all major convergence diagnostics:

| Diagnostic              | What you need                    |
| ----------------------- | -------------------------------- |
| **Pseudo ESS**          | `presence.mean(axis=0)` per chain |
| **ASDSF**               | Per-chain split frequencies       |
| **Fréchet ESS**         | RF distances via XOR of rows      |
| **WRF / KF distances**  | Branch lengths — use `pairwise_wrf_with_snapshots_from_newick_iter` or `pairwise_kf_with_snapshots_from_newick_iter` |

Note: `sum(presence[i] XOR presence[j]) == RF(tree_i, tree_j)` exactly.

### Presence matrix from Python

Use `pairwise_rf_with_snapshots_from_newick_iter` to get the presence matrix
alongside the RF distances, from one parse — see the
[Python API section](#-python-api) above.

```python
import math

import numpy as np
import pandas as pd
import rapidtrees as rtd

# `trees` and `names` are the ones from the Python API example above
tree_names, rf_bytes, leaf_names, n_bip, pres_bytes, bip_clade_bytes = (
    rtd.pairwise_rf_with_snapshots_from_newick_iter(
        names, iter(trees), [{}], [0] * len(trees)
    )
)
n = len(tree_names)
presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n, n_bip).copy()

# Global split frequencies (for Pseudo ESS / ASDSF)
global_freq = presence.mean(axis=0)

# RF distance between any two trees — no recomputation needed
rf_01 = int((presence[0].astype(int) ^ presence[1].astype(int)).sum())

# Named presence matrix — decode clade bitmasks for column labels
bytes_per_bip = math.ceil(len(leaf_names) / 8)
bip_arr  = np.frombuffer(bip_clade_bytes, dtype=np.uint8).reshape(n_bip, bytes_per_bip)
bip_bool = np.unpackbits(bip_arr, axis=1, bitorder="little")[:, :len(leaf_names)]
col_labels = ["|".join(leaf_names[i] for i in np.where(bip_bool[j])[0]) for j in range(n_bip)]
df = pd.DataFrame(presence, index=tree_names, columns=col_labels)
```

---



## ⏱️ Benchmarks

Benchmarks were run on a MacBook Pro M1. Trees are parsed **once** and their interned split IDs are reused across all pairwise comparisons. Parallelism is provided by [`rayon`](https://github.com/rayon-rs/rayon) — no manual thread management needed.

<details>
<summary><strong>Show full benchmark table</strong></summary>

> Output of `cargo bench --bench memory_time_benchmark`

| Taxa (N) | Trees (T) | Combinations | Actual Memory | Wall Time | CPU Time |
|----------|-----------|--------------|---------------|-----------|----------|
| 10       | 100       | 5.0K         | ████░░░░░░ 22.19 KB | ░░░░░░░░░░ 96.31 µs | 393.50 µs |
| 10       | 1000      | 500.0K       | █████░░░░░ 215.55 KB | ░░░░░░░░░░ 911.52 µs | 4.63 ms  |
| 10       | 10000     | 50.0M        | ██████░░░░ 2.10 MB | ░░░░░░░░░░ 41.71 ms | 132.37 ms |
| 100      | 100       | 5.0K         | █████░░░░░ 206.12 KB | ░░░░░░░░░░ 167.12 µs | 834.50 µs |
| 100      | 1000      | 500.0K       | ██████░░░░ 1.94 MB | ░░░░░░░░░░ 1.27 ms | 5.02 ms  |
| 100      | 10000     | 50.0M        | ███████░░░ 19.27 MB | ░░░░░░░░░░ 42.22 ms | 131.66 ms |
| 500      | 100       | 5.0K         | ██████░░░░ 1.00 MB | ░░░░░░░░░░ 151.79 µs | 818.00 µs |
| 500      | 1000      | 500.0K       | ███████░░░ 9.60 MB | ░░░░░░░░░░ 1.08 ms | 6.33 ms  |
| 500      | 10000     | 50.0M        | ████████░░ 95.61 MB | ░░░░░░░░░░ 43.19 ms | 151.38 ms |
| 1000     | 100       | 5.0K         | ██████░░░░ 2.01 MB | ░░░░░░░░░░ 118.50 µs | 683.00 µs |
| 1000     | 1000      | 500.0K       | ███████░░░ 19.20 MB | ░░░░░░░░░░ 1.44 ms | 8.53 ms  |
| 1000     | 10000     | 50.0M        | ████████░░ 191.03 MB | ░░░░░░░░░░ 51.35 ms | 137.48 ms |
| 2000     | 100       | 5.0K         | ██████░░░░ 4.04 MB | ░░░░░░░░░░ 216.77 µs | 1.09 ms  |
| 2000     | 1000      | 500.0K       | ███████░░░ 38.39 MB | ░░░░░░░░░░ 1.45 ms | 9.02 ms  |
| 2000     | 10000     | 50.0M        | ████████░░ 381.89 MB | ░░░░░░░░░░ 36.07 ms | 127.28 ms |
| 5000     | 100       | 5.0K         | ███████░░░ 10.24 MB | ░░░░░░░░░░ 195.81 µs | 981.50 µs |
| 5000     | 1000      | 500.0K       | ████████░░ 96.09 MB | ░░░░░░░░░░ 1.58 ms | 9.22 ms  |
| 5000     | 10000     | 50.0M        | █████████░ 954.57 MB | ░░░░░░░░░░ 39.01 ms | 152.59 ms |
| 10000    | 100       | 5.0K         | ███████░░░ 20.56 MB | ░░░░░░░░░░ 267.69 µs | 1.65 ms  |
| 10000    | 1000      | 500.0K       | ████████░░ 192.24 MB | ░░░░░░░░░░ 2.10 ms | 13.65 ms |
| 10000    | 10000     | 50.0M        | █████████░ 1.86 GB | ░░░░░░░░░░ 48.65 ms | 224.52 ms |
| 20000    | 100       | 5.0K         | ███████░░░ 41.28 MB | ░░░░░░░░░░ 1.37 ms | 3.38 ms  |
| 20000    | 1000      | 500.0K       | ████████░░ 384.62 MB | ░░░░░░░░░░ 3.40 ms | 21.16 ms |
| 20000    | 10000     | 50.0M        | █████████░ 3.73 GB | ░░░░░░░░░░ 113.57 ms | 352.03 ms |
</details>

---

## 🔍 Troubleshooting

- **No trees parsed?** Verify the input is a valid NEXUS `.trees` or Newick file and adjust `--burnin-*` settings.
- **Piping to other tools?** Use `-q` to suppress timing messages on stdout.
- **Gzipped output not working?** Ensure the output filename ends with `.gz`.

---

## ⚖️ License

`rapidtrees` is provided under the [MIT License](LICENSE).
