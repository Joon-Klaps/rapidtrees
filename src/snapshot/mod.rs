//! Tree snapshot types and bulk-collection constructors.
//!
//! # What is a bipartition?
//! Each internal branch in a tree divides the leaves into two groups.
//! ```text
//!      root
//!     /    \
//!   {A,B}  {C,D}  ← This branch creates split {A,B}|{C,D}
//! ```
//! We store one canonical side per split (the side not containing leaf 0).
//!
//! # The pipeline
//! One tree at a time, and never more than a chunk of them alive at once:
//!
//! ```text
//!   newick ──parse──▶ PhyloTree ──[build]──▶ Snapshot ──[intern]──▶ InternSnap
//!                                            (per tree,             (u32 split
//!                                             dropped after)         IDs, kept)
//! ```
//!
//! - [`fingerprint`] — the run-wide tables that make one tree's splits
//!   comparable to another's.
//! - [`build`] — one `PhyloTree` folded into one [`Snapshot`] of fingerprinted
//!   edges.
//! - [`intern`] — those edges deduplicated across every tree into `u32` IDs.
//! - [`export`] — flat byte buffers of the result for the Python side.
//!
//! Pairwise distances then run on the `u32` IDs alone; see [`crate::distances`].

mod build;
mod clades;
mod export;
mod fingerprint;
mod intern;

use build::{Part, Snapshot};
use fingerprint::{build_leaf_index, taxon_labels};
use intern::Interner;

pub(crate) use intern::InternSnap;

use crate::distances::{Backend, Distances};
use crate::par::*;
use clades::CladeTable;
use phylotree::tree::Tree as PhyloTree;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};

/// A bulk collection of tree snapshots in an interned split-ID representation.
///
/// All trees in the set share a single bipartition table: each unique split is
/// assigned a `u32` ID once. Pairwise distance computation then works on
/// `Vec<u32>` instead of leaf sets, which is faster and more cache-friendly.
///
/// # Construction
/// - [`Snapshots::from_newick_iter`] — main constructor for BEAST/NEXUS inputs
/// - [`Snapshots::from_newicks`] — convenience constructor for plain Newick slices
///
/// # Pairwise distances
/// - [`Snapshots::pairwise_rf`], [`Snapshots::pairwise_wrf`], [`Snapshots::pairwise_kf`]
#[derive(Debug)]
pub struct Snapshots {
    pub(crate) snapshots: Vec<InternSnap>,
    /// The canonical leaf set of split ID `i`. Empty when the caller asked not
    /// to build it — see [`Retain`].
    pub(crate) clades: CladeTable,
    pub words_per_bitset: usize,
    /// Alphabetically sorted taxon names shared by all trees in this set.
    pub leaf_names: Vec<String>,
}

/// What to keep besides the split IDs themselves.
///
/// Neither is read by the distance kernels, and both cost real work: branch
/// lengths are `~n_bip × 8` bytes per tree, and materialising a bipartition is
/// an `O(subtree)` walk per distinct split. A path that computes a matrix and
/// exports nothing wants both off.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Retain {
    /// Per-edge branch lengths. Required by WRF and KF; dead weight for RF.
    pub(crate) lengths: bool,
    /// The canonical leaf set per distinct split. Required only to export a
    /// bipartition table to Python or to order export columns.
    pub(crate) bipartitions: bool,
}

impl Retain {
    /// Keep everything — what the public constructors use.
    pub(crate) fn everything() -> Self {
        Self {
            lengths: true,
            bipartitions: true,
        }
    }
}

impl Snapshots {
    /// Build a `Snapshots` collection from a lazy iterator of `(newick, translate_map)` pairs.
    ///
    /// Each `newick` may contain BEAST-format `[&...]` annotations — they are stripped
    /// automatically. The `translate_map` is applied to rename leaf labels (pass an empty
    /// map for plain Newick files with no BEAST translate block).
    ///
    /// All trees must share the same leaf set; an error is returned if any tree differs.
    ///
    /// # Errors
    /// Returns `Err(String)` if any newick fails to parse or leaf sets are inconsistent.
    pub fn from_newick_iter<'a>(
        entries: impl IntoIterator<Item = (&'a str, &'a HashMap<String, String>)>,
        rooted: bool,
    ) -> Result<Self, String> {
        Self::from_newick_iter_opts(entries, rooted, Retain::everything())
    }

    /// Like [`Snapshots::from_newick_iter`], but lets internal callers skip work
    /// they will not read. See [`Retain`].
    ///
    /// The public constructor retains everything, so the Rust/Python API is
    /// unaffected.
    pub(crate) fn from_newick_iter_opts<'a>(
        entries: impl IntoIterator<Item = (&'a str, &'a HashMap<String, String>)>,
        rooted: bool,
        retain: Retain,
    ) -> Result<Self, String> {
        let entries: Vec<_> = entries.into_iter().collect();
        if entries.is_empty() {
            return Ok(Self::empty());
        }

        // Parse the first tree to establish the reference leaf set.
        let (first_newick, first_translate) = entries[0];
        let first_tree = parse_and_rename(first_newick, first_translate, 0)?;

        //sanity check: ensure all leaf names are unique within the first tree
        let first_leaves = first_tree.get_leaves();
        if first_leaves.len() != first_leaves.iter().collect::<HashSet<_>>().len() {
            return Err(
                "Trees have duplicate leaf names. All leaf names must be unique.".to_string(),
            );
        }

        let mut sorted_leaf_names: Vec<String> = first_leaves
            .iter()
            .filter_map(|&id| first_tree.get(&id).ok()?.name.clone())
            .collect();
        sorted_leaf_names.sort_unstable();
        sorted_leaf_names.dedup();

        // One label set and one name → bit table for the whole run: fingerprints
        // are only comparable across trees if every tree draws from the same
        // tables, and deriving them per tree would repeat the same sort on the
        // same names once per tree.
        let labels = taxon_labels(sorted_leaf_names.len());
        let leaf_index = build_leaf_index(&sorted_leaf_names);

        let first_snap = Snapshot::from_tree(&first_tree, rooted, &labels, &leaf_index)
            .map_err(|e| format!("Failed to snapshot tree at index 0: {e}"))?;
        drop(first_tree);

        // Bound how many raw snapshots are alive at once. Holding every tree's
        // un-interned parts simultaneously is the dominant memory cost at
        // construction — it can dwarf the deduplicated result and OOM the process.
        // Estimate one snapshot's raw bytes from the first tree (all trees share
        // the leaf set, so this is representative), parse the rest in chunks sized
        // to ~`CHUNK_TARGET_BYTES`, and fold each chunk into the interner — freeing
        // its parts — before parsing the next.
        const CHUNK_TARGET_BYTES: usize = 256 * 1024 * 1024;
        let per_snap_bytes = first_snap.parts.len() * size_of::<Part>()
            + first_snap.leaf_order.len() * size_of::<u32>();
        let chunk = (CHUNK_TARGET_BYTES / per_snap_bytes.max(1)).clamp(1, 4096);

        let mut interner = Interner::new(
            entries.len(),
            first_snap.words,
            sorted_leaf_names.len(),
            rooted,
            retain,
        );
        interner.push(first_snap);

        let mut base = 1usize; // tree 0 is already interned
        for chunk_entries in entries[1..].chunks(chunk) {
            // Parse this chunk in parallel, validating leaf sets.
            let raw: Vec<Snapshot> = chunk_entries
                .par_iter()
                .enumerate()
                .map(|(k, &(newick, translate))| {
                    let i = base + k;
                    let tree = parse_and_rename(newick, translate, i)?;
                    check_leaf_set(&tree, &leaf_index, i)?;
                    Snapshot::from_tree(&tree, rooted, &labels, &leaf_index)
                        .map_err(|e| format!("Failed to snapshot tree at index {i}: {e}"))
                })
                .collect::<Result<_, _>>()?;

            // Sequential fold: each raw snapshot is dropped right after it is
            // interned, so peak stays near the deduplicated footprint.
            for snap in raw {
                interner.push(snap);
            }
            base += chunk_entries.len();
        }

        Ok(interner.finish(sorted_leaf_names))
    }

    /// Build a `Snapshots` collection from a slice of plain Newick strings.
    ///
    /// No BEAST annotation stripping or taxon renaming is performed — the strings
    /// must already be in standard Newick format.
    ///
    /// # Errors
    /// Returns `Err(String)` if any newick fails to parse or leaf sets differ.
    pub fn from_newicks(newicks: &[&str], rooted: bool) -> Result<Self, String> {
        let empty: HashMap<String, String> = HashMap::new();
        Self::from_newick_iter(newicks.iter().map(|&n| (n, &empty)), rooted)
    }

    /// Number of trees in this collection.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Returns `true` if the collection contains no trees.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// How many distinct splits this collection contains — the `e` in the
    /// `e²/2¹²⁹` fingerprint-collision bound.
    ///
    /// Derived from the split IDs rather than the bipartition table, so it is
    /// correct even on the paths that never materialise one.
    pub fn n_distinct_splits(&self) -> usize {
        self.snapshots
            .iter()
            .filter_map(|s| s.split_ids.last().copied())
            .max()
            .map_or(0, |max_id| max_id as usize + 1)
    }

    /// Compute all pairwise Robinson–Foulds distances as a symmetric n×n matrix.
    ///
    /// Pass `Some(counter)` to track progress: after each row `i` finishes, the
    /// counter is bumped by `n - i - 1`, reaching `n*(n-1)/2` when done. Pass
    /// `None` to skip the (negligible) counter work entirely.
    pub fn pairwise_rf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<u32> {
        self.pairwise_rf_with(progress, Backend::Auto).matrix
    }

    /// Compute all pairwise Weighted Robinson–Foulds distances as a symmetric n×n matrix.
    ///
    /// See [`Self::pairwise_rf`] for the `progress` argument.
    pub fn pairwise_wrf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<f64> {
        self.pairwise_wrf_with(progress, Backend::Auto).matrix
    }

    /// Compute all pairwise Kuhner–Felsenstein distances as a symmetric n×n matrix.
    ///
    /// See [`Self::pairwise_rf`] for the `progress` argument.
    pub fn pairwise_kf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<f64> {
        self.pairwise_kf_with(progress, Backend::Auto).matrix
    }

    /// [`Self::pairwise_rf`] with the backend forced, reporting the kernel that
    /// ran.
    ///
    /// [`Backend::Auto`] is what the plain entry points use and what any caller
    /// should want; the explicit variants exist so a benchmark can pin one, and
    /// the returned [`Distances::kernel`] is how a caller logs `Auto`'s choice.
    /// Both backends return identical matrices.
    pub fn pairwise_rf_with(
        &self,
        progress: Option<&std::sync::atomic::AtomicUsize>,
        backend: Backend,
    ) -> Distances<u32> {
        crate::distances::distance_rf(self, progress, backend)
    }

    /// [`Self::pairwise_wrf`] with the backend forced. See [`Self::pairwise_rf_with`].
    pub fn pairwise_wrf_with(
        &self,
        progress: Option<&std::sync::atomic::AtomicUsize>,
        backend: Backend,
    ) -> Distances<f64> {
        crate::distances::distance_wrf(self, progress, backend)
    }

    /// [`Self::pairwise_kf`] with the backend forced. See [`Self::pairwise_rf_with`].
    pub fn pairwise_kf_with(
        &self,
        progress: Option<&std::sync::atomic::AtomicUsize>,
        backend: Backend,
    ) -> Distances<f64> {
        crate::distances::distance_kf(self, progress, backend)
    }

    fn empty() -> Self {
        Self {
            snapshots: Vec::new(),
            clades: CladeTable::new(),
            words_per_bitset: 0,
            leaf_names: Vec::new(),
        }
    }
}

/// Confirm one tree carries exactly the run's taxa — each one once, none
/// missing, none unknown.
///
/// Checked against `leaf_index` rather than by building this tree's own
/// `HashSet<String>`: the set comparison is the same, but the set costs one
/// `String` allocation per leaf *per tree*, which is the run's whole taxon list
/// cloned once for every tree in the file. The `seen` vector is `n` bytes and is
/// what catches a name repeated within one tree — a plain count would let
/// `(A,A,B)` pass against `{A,B,C}`.
fn check_leaf_set(
    tree: &PhyloTree,
    leaf_index: &FxHashMap<&str, usize>,
    index: usize,
) -> Result<(), String> {
    let mismatch = || {
        format!(
            "Tree {index} has a different leaf set than tree 0. All trees must share the same taxa."
        )
    };

    let mut seen = vec![false; leaf_index.len()];
    for leaf_id in tree.get_leaves() {
        let node = tree
            .get(&leaf_id)
            .map_err(|e| format!("Tree {index}: {e}"))?;
        let name = node.name.as_deref().ok_or_else(|| {
            format!("Tree {index} has an unnamed leaf. All leaves must be named.")
        })?;
        let &bit = leaf_index.get(name).ok_or_else(mismatch)?;
        if std::mem::replace(&mut seen[bit], true) {
            return Err(format!(
                "Tree {index} repeats the leaf name {name:?}. All leaf names must be unique."
            ));
        }
    }

    if seen.iter().all(|&s| s) {
        Ok(())
    } else {
        Err(mismatch())
    }
}

/// Strip BEAST annotations, parse the Newick, and apply taxon renaming.
///
/// `index` only feeds the parse-error message so it points at the offending tree.
fn parse_and_rename(
    newick: &str,
    translate: &HashMap<String, String>,
    index: usize,
) -> Result<PhyloTree, String> {
    let clean = crate::io::strip_beast_annotations(newick);
    let mut tree = PhyloTree::from_newick(&clean)
        .map_err(|e| format!("Failed to parse newick at index {index}: {e}"))?;
    crate::io::rename_leaf_nodes(&mut tree, translate);
    Ok(tree)
}

#[cfg(test)]
mod tests;
