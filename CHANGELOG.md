# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
This project uses release names based on random words from [codenamegenerator.com](https://www.codenamegenerator.com/):
    - PREFIX: Microsoft Corperation
    - DICTIONARY: Snakes

## [Unreleased]
- **Speed (RF/WRF/KF):** All three pairwise metrics now run on a dense, contiguous-row backend instead of the per-pair sorted-merge, sharing one parallel upper-triangle driver.
  - **RF** is reformulated as `RF(i,j) = aᵢ + aⱼ − 2·G[i][j]`, where the shared-split count `G` is `popcount(Pᵢ & Pⱼ)` over packed `u64` bit-rows. Splits present in every tree (the "universal" columns, which always include the pendant edges) cancel exactly in RF and are dropped before packing, shrinking the rows. Results are byte-identical to the old merge. On posterior-like samples (many near-identical trees) this is ~15× faster.
  - **KF** uses the same shape with branch lengths: `KF(i,j) = sqrt(‖wᵢ‖² + ‖wⱼ‖² − 2·⟨wᵢ,wⱼ⟩)`, where the shared term is the dot product of two length rows — the weighted echo of RF's shared-split popcount. A near-zero clamp guards the square root; results match the old merge within floating-point tolerance (verified to <1e-12 on real 162-taxon trees).
  - **WRF** (`Σ |wᵢ − wⱼ|`) has no inner-product shortcut — an absolute difference doesn't factor into a matrix product — so it sweeps the two dense length rows directly. Same column order as the merge, so results are effectively identical, with better cache behaviour and vectorisation.
  - The earlier dense-vs-merge gate is removed: rapidtrees targets similar-tree MCMC sets where the dense backend always wins, so the metrics commit to it unconditionally (the unused merge functions are deleted). The non-export WRF/KF Python paths also release the bipartition table before the O(n²) loop, matching the RF path's memory saving. ([#14](https://github.com/Joon-Klaps/rapidtrees/pull/14))
- **API (Rust):** Unified each metric's two `Snapshots` methods into one. `pairwise_rf`/`pairwise_wrf`/`pairwise_kf` now take a `progress: Option<&AtomicUsize>` argument (`None` to skip progress tracking, `Some(counter)` to report it), and the separate `pairwise_*_counted` methods are removed. This is a breaking change for direct Rust callers only — **the Python API is unchanged** (the `pairwise_*_from_newick_iter` bindings keep identical signatures). ([#14](https://github.com/Joon-Klaps/rapidtrees/pull/14))
- **Memory:** Reduced the RAM footprint of `Snapshots`, targeting the large-tip regime (≤5000 trees, 10,000+ tips). Two changes: (1) removed the never-read `bipartition_index` field, which stored a second full copy of every unique bipartition bitset (~`U × ceil(tips/64) × 8` bytes; GB-scale at high tip counts) — interning now dedups via a `hashbrown::HashTable<u32>` that keeps exactly one copy of each bitset; (2) RF-only paths no longer store per-tree branch lengths, and the non-export RF path also releases the bipartition table before the pairwise loop. No Python API change. ([#13](https://github.com/Joon-Klaps/rapidtrees/pull/13))

## [0.6.0] - Powershell Cobra (2026-05-28)

- Added `rapidtrees.ProgressCounter` — a lock-free, Python-facing handle for observing live progress of a pairwise distance call. Instantiate one, pass it as `progress=` to any of the six `pairwise_*_from_newick_iter` functions, and read `.value()` / `.total()` / `.fraction()` from another Python thread while the call blocks. The GIL is released for the duration of the rayon loop only when a counter is supplied; passing `None` (the default) preserves the original zero-overhead, GIL-held behaviour. Reading is a single atomic load (~1 ns) — no callback, no GIL acquisition from Rust, no monitor thread. ([#12](https://github.com/Joon-Klaps/rapidtrees/pull/12))
- CLI: when stderr is a TTY and `--quiet` is not set, the `rapidtrees` binary now renders a live progress bar to stderr during pairwise distance computation (e.g. `[█████████░░░░░] 64.2% (3211/5000 pairs, ETA 4.1s)`). Piped/redirected runs are unaffected — the bar auto-suppresses when stderr is not a terminal.

## [0.5.1] - 0.5.1 - Azure Mamba (2026-05-12)

- wRF and KF distances now match phangorn ([#9](https://github.com/Joon-Klaps/rapidtrees/pull/9)).
- include tests for load_beast_raw in io.rs ([#10](https://github.com/Joon-Klaps/rapidtrees/pull/10)).
- Add `pairwise_wrf_with_snapshots_from_newick_iter` and `pairwise_kf_with_snapshots_from_newick_iter`. Instead of absence-presence matrix for every bipartition, `branch_length_bytes` field is a flat `float64` matrix of shape `(n_trees, n_bip)` where entry `[i, j]` is the branch length of edge `j` in tree `i` (0.0 if absent). ([#11](https://github.com/Joon-Klaps/rapidtrees/pull/11))

## [0.5.0] - Obsidian Sidewinder (2026-05-11)

- Add `bipartition_clade_bytes: bytes` in the return value of `pairwise_rf_with_snapshots_from_newick_iter` ([#8](https://github.com/Joon-Klaps/rapidtrees/pull/8)). The new field is a packed bitmask buffer of shape `(n_bip, ceil(n_leaves/8))` — decode with `np.unpackbits(..., bitorder='little')`. This representation is more compact (up to 16× smaller for large trees) and has zero Python object overhead.


## [0.4.0] - Xbox Treeboa (2026-05-01)

- Interned snapshots for faster RF by @hongsamL in https://github.com/Joon-Klaps/rapidtrees/pull/4
- Export tree snapshots by @Joon-Klaps in https://github.com/Joon-Klaps/rapidtrees/pull/3

## [0.3.0] - Gale Sidewinder (2026-03-03)

Release with API changes for compatibility with `treetracer`

## [0.2.1] - Slate King Cobra (2025-11-16)

### Changed

#### Performance Improvements 🚀

- **Major algorithmic optimization**: Replaced HashMap/HashSet-based intersection with sorted merge algorithm for distance calculations
  - Direct array indexing instead of HashMap lookups (10-20× faster for branch length access)
  - Performance: ~2.3M tree comparisons/second on ZIKA dataset (283 taxa, 4000 trees, ~8M comparisons)
  - All three distance metrics (RF, Weighted RF, KF) now complete in ~3.5 seconds, and are deterministic

- **Optimized data structures**:
  - `TreeSnapshot` now uses `Vec<Bitset>` instead of `HashSet<Bitset>` for partitions (sorted for merge)
  - Parallel vectors for partitions and branch lengths instead of `HashMap<Bitset, f64>`
  - Switched to `FxHashMap` from `rustc-hash` crate for remaining hash operations (faster than std HashMap)

- **Added `#[inline]` annotations** to critical distance calculation functions for better compiler optimization

#### Bug Fixes 🐛

- Fixed BEAST tree name parsing in `io.rs`:
  - Previously only extracted state numbers, now properly extracts full tree names
  - Handles various BEAST tree name formats (e.g., `classic2_STATE_968940000`, `STATE_8500000`)

### Added

- Updated pre-commit configuration

## [0.1.0] - Aero Boa (Initial Release)

Initial release with basic functionality for computing phylogenetic tree distances:

- Robinson-Foulds distance
- Weighted Robinson-Foulds distance
- Kuhner-Felsenstein distance
- BEAST/NEXUS file format support
- Python API via maturin

