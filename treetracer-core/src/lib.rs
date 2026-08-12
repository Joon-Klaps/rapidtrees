//! Crate root: the browser-facing half of treetracer.
//!
//! `rapidtrees` answers one question — how far apart are two trees — and this
//! crate is everything built on top of that answer for the web UI: projecting a
//! distance matrix, summarising clade support, estimating ESS, laying trees out
//! for drawing, and the wasm-bindgen surface that hands it all to JavaScript.
//!
//! The split is by *audience*, not by subject. `rapidtrees` is consumed by
//! Python (PyO3) and by a CLI; nothing here is reachable from either. Keeping
//! the two apart means a change to a projection or a layout cannot break a
//! published Python wheel, and `rapidtrees` no longer compiles wasm-bindgen it
//! never uses.
//!
//! Build with:
//!
//! ```text
//! wasm-pack build treetracer-core --target web
//! ```
//!
//! Modules:
//! - `clades`: clade frequencies and maximum clade credibility trees.
//! - `ess`: pseudo-ESS over a distance matrix.
//! - `hipstr`: HIPSTR summary-tree construction.
//! - `layout`: tree coordinates for drawing.
//! - `mds`: principal coordinates analysis of the distance matrix.
//! - `stats`: small shared statistical helpers.
//! - `wasm`: the wasm-bindgen boundary — keep it thin.

pub mod clades;
pub mod ess;
pub mod hipstr;
pub mod layout;
pub mod mds;
pub mod stats;

pub mod wasm;
