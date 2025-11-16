# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
This project uses release names based on random words from [codenamegenerator.com](https://www.codenamegenerator.com/):
    - PREFIX: Microsoft Corperation
    - DICTIONARY: Snakes

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

