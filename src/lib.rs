//! Crate root: lightweight module orchestration and public re-exports.
//!
//! Modules:
//! - `distances`: single-pair and pairwise distance functions.
//! - `io`: reading and parsing BEAST/NEXUS tree files.
//! - `bitset`: compact bitset representation for tree partitions.
//! - `snapshot`: tree snapshot and interned snapshot types (crate-internal).
//! - `api`: Python bindings via `pyo3` (gated behind "python" feature).

pub mod bitset;
pub mod distances;
pub mod io;
pub(crate) mod parallel;
pub(crate) mod snapshot;

#[cfg(feature = "python")]
pub mod api;
#[cfg(feature = "python")]
pub(crate) mod progress;

// In-crate Newick parser used in place of `phylotree` for the wasm build. Only
// compiled when it is actually the active backend (wasm without phylotree), plus
// under `test` so a native test can cross-check it against phylotree.
#[cfg(any(all(feature = "wasm", not(feature = "phylotree")), test))]
#[path = "phylotree-patch.rs"]
pub(crate) mod phylotree_patch;

// WebAssembly bindings. The `compute_distances_core` helper is also compiled
// under `test` (native) so it can be unit-tested without a browser; only the
// `#[wasm_bindgen]` wrapper inside requires the `wasm` feature.
#[cfg(any(feature = "wasm", test))]
mod wasm;

pub use bitset::Bitset;
#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::{Snapshot, Snapshots};
