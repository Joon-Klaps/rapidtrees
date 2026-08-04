//! Crate root: lightweight module orchestration and public re-exports.
//!
//! Modules:
//! - `distances`: single-pair and pairwise distance functions.
//! - `io`: reading and parsing BEAST/NEXUS tree files.
//! - `bitset`: compact bitset representation for tree partitions.
//! - `snapshot`: tree snapshot and interned snapshot types (crate-internal).
//! - `api`: Python bindings via `pyo3` (gated behind "python" feature).

pub mod bitset;
pub mod clades;
pub mod distances;
pub mod ess;
pub mod hipstr;
pub mod io;
pub mod layout;
pub mod mds;
pub(crate) mod par;
pub(crate) mod snapshot;
pub mod stats;

#[cfg(feature = "python")]
pub mod api;
#[cfg(feature = "python")]
pub(crate) mod progress;

/// Browser bindings via `wasm-bindgen` (gated behind the "wasm" feature).
#[cfg(feature = "wasm")]
pub mod wasm;

pub use bitset::Bitset;
#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::{Snapshot, Snapshots};
