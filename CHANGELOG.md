# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
This project uses release names based on random words from [codenamegenerator.com](https://www.codenamegenerator.com/):
    - PREFIX: Microsoft Corperation
    - DICTIONARY: Snakes


## [0.5.1] - (YYYY-MM-DD)

- wRF and KF distances now match phangorn ([#9](https://github.com/Joon-Klaps/rapidtrees/pull/9)).

## [0.5.0] - (2026-05-11)

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

