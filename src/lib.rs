//! Crate root: lightweight module orchestration and public re-exports.
//!
//! This crate does one thing: pairwise phylogenetic tree distances, fast. It is
//! consumed by a CLI (`main.rs`), by Python through PyO3 (`api.rs`), and — as a
//! plain Rust dependency — by the `treetracer-core` workspace member, which
//! holds the browser bindings and the diagnostics built on top of a distance
//! matrix. Nothing wasm-specific lives here; the crate is merely kept
//! wasm-*compatible*, which is what the phylotree fork, the `parallel` feature
//! gate and the pure-Rust dependency choices are for.
//!
//! Modules:
//! - `distances`: single-pair and pairwise distance functions.
//! - `io`: reading and parsing BEAST/NEXUS tree files.
//! - `bitset`: compact bitset representation for tree partitions.
//! - `snapshot`: tree snapshot and interned snapshot types (crate-internal).
//! - `par`: rayon-or-sequential shim; public so dependents share one policy.
//! - `api`: Python bindings via `pyo3` (gated behind "python" feature).

pub mod bitset;
pub mod distances;
pub mod io;
/// Parallelism shim. Public so that dependent crates get the same
/// rayon-or-sequential behaviour this crate's `parallel` feature selects,
/// rather than each picking its own and diverging on wasm32.
pub mod par;
pub(crate) mod snapshot;

#[cfg(feature = "python")]
pub mod api;
#[cfg(feature = "python")]
pub(crate) mod progress;

pub use bitset::Bitset;
#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::{Snapshot, Snapshots};
