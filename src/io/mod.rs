//! Reading BEAST/NEXUS and plain-Newick tree files, and writing distance matrices.
//!
//! - `annotations`: stripping BEAST `[&...]` annotations from a Newick string.
//! - `nexus`: NEXUS primitives — tree headers, tree lines, the `TRANSLATE` block.
//! - `format`: telling NEXUS from plain Newick.
//! - `reader`: the streaming reader that hands trees out one at a time.
//! - `load`: [`read_trees`] and [`read_tree_file`], a whole input into
//!   [`crate::Snapshots`], with their [`ReadOptions`] and [`ReadTreesError`].
//! - `write`: [`write_matrix_tsv`] (behind the `cli` feature).
//!
//! Everything public here is re-exported flat, so the paths stay
//! `rapidtrees::io::*`. `strip_beast_annotations`, `extract_name_state`,
//! `nexus_tree_lines` and `parse_taxon_block` are what readers outside this
//! crate build on, so their signatures are published surface.

mod annotations;
mod format;
mod load;
mod nexus;
mod reader;
#[cfg(feature = "cli")]
mod write;

pub use annotations::strip_beast_annotations;
pub use format::{TreeFormat, detect_format};
// Re-exporting a deprecated item counts as using it.
#[allow(deprecated)]
pub use load::load_beast_trees;
pub use load::{MIN_TREES, ReadOptions, ReadTreesError, read_tree_file, read_trees};
pub use nexus::{extract_name_state, nexus_tree_lines, parse_taxon_block};
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
