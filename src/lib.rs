//! Crate root: lightweight module orchestration and public re-exports.
//!
//! Modules:
//! - `distances`: single-pair and pairwise distance functions.
//! - `io`: reading and parsing BEAST/NEXUS tree files.
//! - `snapshot`: tree snapshot and interned snapshot types (crate-internal).
//! - `par`: rayon-or-sequential shim; public so dependents share one policy.
//! - `cpu`: runtime choice of instruction set for the pairwise row loops.
//! - `api`: Python bindings via `pyo3` (gated behind "python" feature).

pub(crate) mod cpu;
pub mod distances;
pub mod io;
/// Parallelism shim. Rayon vs sequential iterators.
pub mod par;
pub(crate) mod snapshot;

#[cfg(feature = "python")]
pub mod api;
#[cfg(feature = "python")]
pub(crate) mod progress;

#[cfg(feature = "cli")]
pub use io::write_matrix_tsv;
pub use snapshot::{Retain, Snapshots};

/// The instruction-set level the pairwise row loops run at on this machine:
/// `"baseline"`, `"x86-64-v2"`, `"x86-64-v3"` or `"x86-64-v4"`.
///
/// The best level the CPU supports, chosen at runtime, so a prebuilt wheel or
/// a plain `cargo install` still uses POPCNT, AVX2 or AVX-512 where the CPU
/// has them. Setting `RAPIDTREES_CPU` to one of these names caps it. Off
/// x86-64 this is always `"baseline"`, which on aarch64 already includes NEON.
pub fn cpu_level() -> &'static str {
    cpu::Level::current().name()
}
