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
pub(crate) mod snapshot;

#[cfg(feature = "gpu")]
pub(crate) mod gpu;

#[cfg(feature = "python")]
pub mod api;
#[cfg(feature = "python")]
pub(crate) mod progress;

pub use bitset::Bitset;
#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::{Snapshot, Snapshots};

/// Human-readable label for the GPU adapter rapidtrees would use (e.g.
/// `Tesla P100-SXM2-16GB (Vulkan, DiscreteGpu)`), or `None` if no compatible GPU
/// is available and computation would fall back to the CPU.
///
/// Forces the lazy GPU context to initialise on first call. Lets callers confirm
/// up front whether the GPU will actually be used instead of inferring it from
/// timings. Only present when the crate is built with the `gpu` feature.
#[cfg(feature = "gpu")]
pub fn gpu_adapter_label() -> Option<String> {
    gpu::adapter_label()
}
