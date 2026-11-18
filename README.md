# 🌲⚡ rapidtrees

<p align="center">
  <em>Blazingly fast pairwise phylogenetic tree distance calculations — Robinson–Foulds, Weighted RF, and Kuhner–Felsenstein — powered by Rust</em>
</p>

<p align="center">
  <a href="https://crates.io/crates/rapidtrees"><img src="https://img.shields.io/crates/v/rapidtrees.svg?style=flat-square&logo=rust&logoColor=white" alt="Crates.io" /></a>
  <a href="https://pypi.org/project/rapidtrees"><img src="https://img.shields.io/pypi/v/rapidtrees.svg?style=flat-square&logo=python&logoColor=white" alt="PyPI" /></a>
  <a href="https://www.npmjs.com/package/rapidtrees"><img src="https://img.shields.io/npm/v/rapidtrees.svg?style=flat-square&logo=npm&logoColor=white" alt="npm" /></a>
  <a href="https://github.com/Joon-Klaps/rapidtrees/actions/workflows/ci.yml"><img src="https://img.shields.io/github/actions/workflow/status/Joon-Klaps/rapidtrees/ci.yml?style=flat-square&logo=github&label=CI" alt="CI" /></a>
  <a href="https://codecov.io/gh/Joon-Klaps/rapidtrees"><img src="https://img.shields.io/codecov/c/github/Joon-Klaps/rapidtrees?style=flat-square&logo=codecov" alt="Coverage" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square" alt="License: MIT" /></a>
  <a href="https://doi.org/10.5281/zenodo.23811107"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23811107.svg" alt="DOI" /></a>
</p>

<p align="center">
  <a href="#%EF%B8%8F-overview">Overview</a> •
  <a href="#-installing">Installing</a> •
  <a href="#-command-line">Command line</a> •
  <a href="#-from-code">From code</a> •
  <a href="#-documentation">Documentation</a>
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
- 🐍 **Python, CLI and WebAssembly** — one kernel behind a NumPy-friendly Python API, a command-line tool and an npm package
- 📦 **No Rust toolchain required** — pre-built wheels on PyPI for Linux, macOS, and Windows
- 🗜️ **Gzip output** — stream directly to `.tsv.gz` without a separate compression step

Two distinct splits are merged only if their 128-bit fingerprints collide, with probability about `e²/2¹²⁹` for `e` distinct splits in the run: `1.5 × 10⁻²³` at a hundred million splits. Every run prints its own `e` and bound.

### 🚀 Performance

Measured on an Apple M1 Max laptop (10 cores, 64 GB RAM) with rapidtrees 0.11.0 from crates.io, on four published BEAST posteriors from the BEAST workshop tutorials ([YFV](https://beast.community/workshop_continuous_diffusion_yfv), [WNV](https://beast.community/workshop_continuous_diffusion_wnv), [batRABV](https://beast.community/workshop_discrete_diffusion_rabv), [EBOV](https://beast.community/workshop_discrete_diffusion_ebov)). Wall-clock time for the whole run, reading the NEXUS file included, with all 10 cores and with one thread in brackets; median of three runs after a warm-up, matrix written to `/dev/null`:

| Dataset | Taxa · trees | Pairs   | Robinson-Foulds | Weighted RF     | Kuhner-Felsenstein | Peak memory |
| ------- | ------------ | ------- | --------------- | --------------- | ------------------ | ----------- |
| YFV     | 65 · 5 001   | 12.5M  | 0.12s (0.56s) | 0.41s (2.38s) | 0.37s (2.21s)    | 393 MiB     |
| WNV     | 104 · 1 001  | 0.50M  | 0.02s (0.07s) | 0.04s (0.15s) | 0.04s (0.14s)    | 85 MiB      |
| batRABV | 372 · 1 001  | 0.50M  | 0.05s (0.15s) | 0.07s (0.27s) | 0.06s (0.27s)    | 157 MiB     |
| EBOV    | 836 · 1 167  | 0.68M  | 0.07s (0.30s) | 0.11s (0.70s) | 0.10s (0.54s)    | 159 MiB     |

Command: `rapidtrees --input <file>.trees --metric <rf|weighted|kf> --output /dev/null --use-real-taxa --quiet`, no burn-in, unrooted. Synthetic scaling results are in [benches/README.md](https://github.com/Joon-Klaps/rapidtrees/blob/master/benches/README.md).

---

## 🔧 Installing

```bash
pip install rapidtrees     # Python, pre-built wheels, no Rust toolchain needed
cargo install rapidtrees   # command-line tool, needs the Rust toolchain
npm install rapidtrees     # WebAssembly build for browsers and Node
```

To build from source and run the tests, see [CONTRIBUTING.md](https://github.com/Joon-Klaps/rapidtrees/blob/master/CONTRIBUTING.md).

---

## 💻 Command line

```bash
rapidtrees -i tests/data/hiv1.trees -o out/hiv1_rf.tsv.gz --metric rf
```

| Flag                          | Description                                                        |
| ----------------------------- | ------------------------------------------------------------------ |
| `-i, --input <INPUT>`         | Path to a NEXUS `.trees` or plain Newick file (auto-detected)      |
| `-o, --output <OUTPUT>`       | Output path. Use `.gz` suffix for gzip compression; `-` for stdout |
| `-t, --burnin-trees <N>`      | Drop the first N trees (default: `0`)                              |
| `-s, --burnin-states <STATE>` | Drop trees sampled at MCMC state `STATE` or before (default: `0`)  |
| `--use-real-taxa`             | Map numeric taxon IDs via the TRANSLATE block                      |
| `--metric <rf\|weighted\|kf>` | Distance metric (default: `rf`)                                    |
| `-q, --quiet`                 | Suppress progress messages (errors still go to stderr)             |

`--input` takes either a NEXUS/BEAST trees file or a plain Newick file with one tree per line; the format is sniffed from the contents, not the extension. Newick files carry no `TRANSLATE` block and no `STATE_` labels, so `--burnin-states` and `--use-real-taxa` have nothing to act on there (`--burnin-trees` still works).

The output is a **square TSV matrix** whose header row and first column hold the tree names: `<file_basename>_<tree_name>` for NEXUS (e.g. `hiv1_STATE_10000`) and `<file_basename>_line<n>` for Newick (e.g. `hiv1_line3`). Progress messages go to stdout, so pass `-q` when writing the matrix to stdout with `-o -`.

---

## 🧩 From code

**Python.** Every function takes an iterator of Newick strings, so memory stays flat however many trees you feed it:

```python
import numpy as np
import rapidtrees as rtd

trees = [
    "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
    "((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);",
    "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.2);",
]
names = ["t1", "t2", "t3"]

# Newick strings already use real taxon names, so pass one empty translate map for every tree
tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(names, iter(trees), [{}], [0] * len(trees))
rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(len(tree_names), -1)
# [[0 2 0]
#  [2 0 2]
#  [0 2 0]]
```

**Rust.** `io::read_tree_file` reads a NEXUS or Newick file, applies burn-in and returns the tree names with a `Snapshots` ready for the pairwise kernels:

```rust
use rapidtrees::io::{ReadOptions, read_tree_file};

let mut opts = ReadOptions::default();
opts.burnin_trees = 100;
let (names, snaps) = read_tree_file("posterior.trees", &opts)?;
let rf = snaps.into_pairwise_rf(None); // row-major, names.len() × names.len()
```

**JavaScript.** `pairwise_rf`, `pairwise_wrf` and `pairwise_kf` take the text of a tree file and return the names and the matrix. The build is single-threaded and does not stream.

---

## 📚 Documentation

| Document | Covers |
| --- | --- |
| [Python API](https://github.com/Joon-Klaps/rapidtrees/blob/master/docs/python-api.md) | All six functions, BEAST files and translate maps, progress reporting, the presence and branch-length matrices for convergence diagnostics |
| [Rust API (docs.rs)](https://docs.rs/rapidtrees) | `io::read_trees`, `io::read_tree_file`, `Snapshots` and the distance kernels |
| [JavaScript API](https://github.com/Joon-Klaps/rapidtrees/blob/master/docs/javascript-api.md) | The npm package: all functions and return types, browser and Node usage, memory limits, building it locally |
| [rapidtrees for Dummies](https://github.com/Joon-Klaps/rapidtrees/blob/master/docs/rapidtrees-for-dummies.md) | How it works, from Newick string to distance matrix, and the kernel tuning knobs |
| [Benchmarks](https://github.com/Joon-Klaps/rapidtrees/blob/master/benches/README.md) | Synthetic memory and time scaling |
| [Contributing](https://github.com/Joon-Klaps/rapidtrees/blob/master/CONTRIBUTING.md) | Building from source, tests, releases |

---

## ⚖️ License

`rapidtrees` is provided under the [MIT License](LICENSE).
