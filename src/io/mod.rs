//! Reading BEAST/NEXUS and plain-Newick tree files, and writing distance matrices.
//!
//! - `annotations`: stripping BEAST `[&...]` annotations from a Newick string.
//! - `nexus`: NEXUS primitives — tree headers, tree lines, the `TRANSLATE` block.
//! - `format`: telling NEXUS from plain Newick.
//! - `reader`: the streaming reader that hands trees out one at a time.
//! - `load`: [`load_beast_trees`], a whole file into [`crate::Snapshots`].
//! - `write`: [`write_matrix_tsv`] (behind the `cli` feature).
//!
//! Everything public here is re-exported flat, so the paths stay
//! `rapidtrees::io::*`. `strip_beast_annotations`, `extract_name_state` and
//! `parse_taxon_block` are treetracer-web's API.

mod annotations;
mod format;
mod load;
mod nexus;
mod reader;
#[cfg(feature = "cli")]
mod write;

pub use annotations::strip_beast_annotations;
pub use format::{TreeFormat, detect_format};
pub use load::load_beast_trees;
pub use nexus::{extract_name_state, parse_taxon_block};
#[cfg(feature = "cli")]
pub use write::write_matrix_tsv;

/// The 21-tree BEAST posterior every loader test reads.
#[cfg(test)]
fn hiv2_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hiv2.trees")
}

/// The same 21 trees as `hiv2.trees`, written one per line as plain Newick.
#[cfg(test)]
fn hiv2_newick_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hiv2.newick")
}
