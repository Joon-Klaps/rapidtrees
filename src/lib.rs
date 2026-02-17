//! Crate root: lightweight module orchestration and public re-exports.
//!
//! Modules:
//! - `distances`: generic distance trait + RF / cluster affinity implementations.
//! - `io`: reading and parsing BEAST/NEXUS tree files.
//! - `bitset`: compact bitset representation for tree partitions.
//! - `snapshot`: tree snapshot for efficient distance calculations.
//! - `api`: Python bindings via `pyo3` (gated behind "python" feature).
//!
//! Public API kept stable by re-exporting key items from the new modules.

pub mod bitset;
pub mod distances;
pub mod io;
pub mod snapshot;

#[cfg(feature = "python")]
pub mod api;

// Re-export frequently used types & functions
pub use bitset::Bitset;
pub use io::{read_beast_trees, rename_leaf_nodes, strip_beast_annotations};
#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::TreeSnapshot;
