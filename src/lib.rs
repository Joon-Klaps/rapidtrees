//! Crate root: lightweight module orchestration and public re-exports.
//!
//! Modules:
//! - `distances`: single-pair and pairwise distance functions.
//! - `io`: reading and parsing BEAST/NEXUS tree files.
//! - `bitset`: compact bitset representation for tree partitions.
//! - `gpu_layout`: GPU buffer layouts and WGSL sources (no GPU dependency).
//! - `snapshot`: tree snapshot and interned snapshot types (crate-internal).
//! - `par`: rayon-or-sequential shim; public so dependents share one policy.
//! - `api`: Python bindings via `pyo3` (gated behind "python" feature).

pub mod bitset;
pub mod distances;
/// GPU buffer layouts and WGSL kernel sources. Like the `pub` items in `io`,
/// this is treetracer-web's API: it owns the `wgpu` device, this owns the data
/// layout the shaders read.
pub mod gpu_layout;
pub mod io;
/// Parallelism shim. Rayon vs sequential iterators.
pub mod par;
pub(crate) mod snapshot;

#[cfg(feature = "python")]
pub mod api;
#[cfg(feature = "python")]
pub(crate) mod progress;

pub use bitset::Bitset;
#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::Snapshots;
