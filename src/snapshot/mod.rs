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
//!   newick ──[newick]──▶ Snapshot ──[intern]──▶ InternSnap
//!                        (per tree,             (u32 split
//!                         dropped after)         IDs, kept)
//! ```
//!
//! - [`fingerprint`] — the run-wide tables that make one tree's splits
//!   comparable to another's.
//! - [`newick`] — one tree's text read straight into a [`Snapshot`] of
//!   fingerprinted edges, with no tree built in between.
//! - [`build`] — the [`Snapshot`] and [`Part`] types themselves.
//! - [`intern`] — those edges deduplicated across every tree into `u32` IDs.
//! - [`export`] — flat byte buffers of the result for the Python side.
//!
//! Pairwise distances then run on the `u32` IDs alone; see [`crate::distances`].

mod build;
mod clades;
mod export;
mod fingerprint;
mod intern;
mod newick;
mod rooted_facts;

use build::{Part, Snapshot};
use fingerprint::{build_leaf_index, taxon_labels};
use intern::Interner;
use rooted_facts::{RawRootedFacts, RootedFactsStore};

pub(crate) use intern::InternSnap;

use crate::par::*;
use clades::CladeTable;
use std::collections::HashMap;

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
    /// How many trees hold each split, indexed by split ID.
    pub(crate) split_counts: Vec<u32>,
    pub words_per_bitset: usize,
    /// Alphabetically sorted taxon names shared by all trees in this set.
    pub leaf_names: Vec<String>,
    /// Optional tree-aligned node heights and directly observed rooted splits.
    /// Existing distance and snapshot-export paths leave this absent.
    rooted_facts: Option<RootedFactsStore>,
}

/// What to keep besides the split IDs themselves.
///
/// None is read by the distance kernels, and each costs real work: branch
/// lengths are `~n_bip × 8` bytes per tree, and materialising a bipartition is
/// an `O(subtree)` walk per distinct split. Rooted facts add a second linear
/// tree traversal plus compact per-node storage. A path that computes a matrix
/// and exports nothing wants all three off.
#[derive(Debug, Clone, Copy)]
pub struct Retain {
    /// Per-edge branch lengths. Required by WRF and KF; dead weight for RF.
    pub lengths: bool,
    /// The canonical leaf set per distinct split. Required only to export a
    /// bipartition table to Python or to order export columns.
    pub bipartitions: bool,
    /// Exact node heights and observed child splits for rooted summary trees.
    /// Existing public constructors and distance APIs deliberately leave this
    /// off; the dedicated rooted-facts endpoint will opt in.
    pub rooted_facts: bool,
}

impl Retain {
    /// Keep the pre-existing branch-length and bipartition exports.
    ///
    /// Rooted facts are a separate optional feature, so public constructors
    /// retain their exact historical validation, memory use, and behavior.
    pub fn everything() -> Self {
        Self {
            lengths: true,
            bipartitions: true,
            rooted_facts: false,
        }
    }

    /// Keep only what a distance matrix needs
    pub fn for_distances(weighted: bool) -> Self {
        Self {
            lengths: weighted,
            bipartitions: false,
            rooted_facts: false,
        }
    }
}

impl Snapshots {
    /// Build a `Snapshots` collection from a lazy iterator of `(newick, translate_map)` pairs.
    ///
    /// Each `newick` may contain BEAST-format `[&...]` annotations, or any other
    /// `[...]` comment; the reader skips them. The `translate_map` is applied to
    /// rename leaf labels (pass an empty map for plain Newick files with no BEAST
    /// translate block).
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

    /// Like [`Snapshots::from_newick_iter`], but lets a caller skip work it will
    /// not read. See [`Retain`].
    ///
    /// [`Snapshots::from_newick_iter`] retains everything, so callers that do
    /// not opt in are unaffected.
    pub fn from_newick_iter_opts<'a>(
        entries: impl IntoIterator<Item = (&'a str, &'a HashMap<String, String>)>,
        rooted: bool,
        retain: Retain,
    ) -> Result<Self, String> {
        if retain.rooted_facts && !rooted {
            return Err("Rooted facts require rooted snapshot mode.".to_string());
        }

        let entries: Vec<_> = entries.into_iter().collect();
        if entries.is_empty() {
            return Ok(Self::empty(retain.rooted_facts));
        }

        // Tree 0 defines the run's taxa, so its names are read first and on
        // their own: every tree's leaf check, tree 0's included, needs the very
        // table they are about to build. An unnamed leaf fails inside
        // `leaf_names`; a repeated one is caught here.
        let (first_newick, first_translate) = entries[0];
        let mut sorted_leaf_names = newick::leaf_names(first_newick, first_translate)?;
        let n_leaves = sorted_leaf_names.len();
        sorted_leaf_names.sort_unstable();
        sorted_leaf_names.dedup();
        if sorted_leaf_names.len() != n_leaves {
            return Err(
                "Tree 0 has duplicate leaf names. All leaf names must be unique.".to_string(),
            );
        }

        // One label set and one name → bit table for the whole run: fingerprints
        // are only comparable across trees if every tree draws from the same
        // tables, and deriving them per tree would repeat the same sort on the
        // same names once per tree.
        let labels = taxon_labels(sorted_leaf_names.len());
        let leaf_index = build_leaf_index(&sorted_leaf_names);
        let run = newick::RunTables {
            leaf_index: &leaf_index,
            labels: &labels,
            total: labels.iter().fold(0, |acc, &label| acc ^ label),
            rooted,
            require_explicit_lengths: retain.rooted_facts,
        };

        let first_snap = newick::snapshot(first_newick, first_translate, 0, &run)?;
        let first_facts = collect_raw_rooted_facts(&first_snap, 0, retain.rooted_facts)?;

        // Bound how many raw snapshots are alive at once. Holding every tree's
        // un-interned parts simultaneously is the dominant memory cost at
        // construction — it can dwarf the deduplicated result and OOM the process.
        // Estimate one snapshot's raw bytes from the first tree (all trees share
        // the leaf set, so this is representative), parse the rest in chunks sized
        // to ~`CHUNK_TARGET_BYTES`, and fold each chunk into the interner — freeing
        // its parts — before parsing the next.
        const CHUNK_TARGET_BYTES: usize = 256 * 1024 * 1024;
        let per_snap_bytes = estimated_raw_snapshot_bytes(&first_snap, first_facts.as_ref());
        let chunk = (CHUNK_TARGET_BYTES / per_snap_bytes.max(1)).clamp(1, 4096);

        let mut interner = Interner::new(
            entries.len(),
            first_snap.words,
            sorted_leaf_names.len(),
            rooted,
            retain,
        );
        interner.push(first_snap);
        if let Some(facts) = first_facts {
            interner
                .push_rooted_facts(facts)
                .map_err(|e| format!("Failed to intern rooted facts for tree at index 0: {e}"))?;
        }

        let mut base = 1usize; // tree 0 is already interned
        for chunk_entries in entries[1..].chunks(chunk) {
            // Parse this chunk in parallel, validating leaf sets.
            let raw: Vec<(Snapshot, Option<RawRootedFacts>)> = chunk_entries
                .par_iter()
                .enumerate()
                .map(|(k, &(newick, translate))| {
                    let i = base + k;
                    let snap = newick::snapshot(newick, translate, i, &run)?;
                    let facts = collect_raw_rooted_facts(&snap, i, retain.rooted_facts)?;
                    Ok::<_, String>((snap, facts))
                })
                .collect::<Result<_, _>>()?;

            // Sequential fold: each raw snapshot is dropped right after it is
            // interned, so peak stays near the deduplicated footprint.
            for (offset, (snap, facts)) in raw.into_iter().enumerate() {
                interner.push(snap);
                if let Some(facts) = facts {
                    let i = base + offset;
                    interner.push_rooted_facts(facts).map_err(|e| {
                        format!("Failed to intern rooted facts for tree at index {i}: {e}")
                    })?;
                }
            }
            base += chunk_entries.len();
        }

        Ok(interner.finish(sorted_leaf_names))
    }

    /// Build a `Snapshots` collection from a slice of plain Newick strings.
    ///
    /// No taxon renaming is performed, so leaves must carry their taxon names.
    /// `[...]` comments, BEAST annotations included, are skipped as in
    /// [`Self::from_newick_iter`].
    ///
    /// # Errors
    /// Returns `Err(String)` if any newick fails to parse or leaf sets differ.
    pub fn from_newicks(newicks: &[&str], rooted: bool) -> Result<Self, String> {
        let empty: HashMap<String, String> = HashMap::new();
        Self::from_newick_iter(newicks.iter().map(|&n| (n, &empty)), rooted)
    }

    /// Number of trees in this collection.
    pub fn len(&self) -> usize {
        debug_assert!(
            self.rooted_facts
                .as_ref()
                .is_none_or(|facts| facts.len() == self.snapshots.len()),
            "rooted facts must stay aligned with snapshot rows"
        );
        self.snapshots.len()
    }

    /// Returns `true` if the collection contains no trees.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// How many distinct splits this collection contains — the `e` in the
    /// `e²/2¹²⁹` fingerprint-collision bound.
    ///
    /// Read off the per-split tree counts rather than the bipartition table, so
    /// it is correct even on the paths that never materialise one.
    pub fn n_distinct_splits(&self) -> usize {
        self.split_counts.len()
    }

    /// Compute all pairwise Robinson–Foulds distances as a symmetric n×n matrix.
    ///
    /// Pass `Some(counter)` to track progress: after each row `i` finishes, the
    /// counter is bumped by `n - i - 1`, reaching `n*(n-1)/2` when done. Pass
    /// `None` to skip the (negligible) counter work entirely.
    pub fn pairwise_rf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<u32> {
        crate::distances::distance_rf(self, progress)
    }

    /// Compute all pairwise Weighted Robinson–Foulds distances as a symmetric n×n matrix.
    ///
    /// See [`Self::pairwise_rf`] for the `progress` argument.
    pub fn pairwise_wrf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<f64> {
        crate::distances::distance_wrf(self, progress)
    }

    /// Compute all pairwise Kuhner–Felsenstein distances as a symmetric n×n matrix.
    ///
    /// See [`Self::pairwise_rf`] for the `progress` argument.
    pub fn pairwise_kf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<f64> {
        crate::distances::distance_kf(self, progress)
    }

    fn empty(retain_rooted_facts: bool) -> Self {
        Self {
            snapshots: Vec::new(),
            clades: CladeTable::new(),
            split_counts: Vec::new(),
            words_per_bitset: 0,
            leaf_names: Vec::new(),
            rooted_facts: retain_rooted_facts.then(RootedFactsStore::default),
        }
    }
}

/// Approximate heap bytes held by one parsed snapshot until it is interned.
///
/// The original estimate covers the ordinary snapshot vectors. Optional raw
/// rooted facts are added only on the facts-enabled path, keeping the existing
/// chunk sizing unchanged otherwise.
fn estimated_raw_snapshot_bytes(
    snapshot: &Snapshot,
    rooted_facts: Option<&RawRootedFacts>,
) -> usize {
    snapshot.parts.len() * size_of::<Part>()
        + snapshot.leaf_order.len() * size_of::<u32>()
        + rooted_facts.map_or(0, RawRootedFacts::estimated_heap_bytes)
}

/// Collect the optional sidecar from the direct parser's rooted snapshot.
fn collect_raw_rooted_facts(
    snapshot: &Snapshot,
    index: usize,
    retain: bool,
) -> Result<Option<RawRootedFacts>, String> {
    if !retain {
        return Ok(None);
    }
    RawRootedFacts::from_snapshot(snapshot)
        .map(Some)
        .map_err(|e| format!("Failed to collect rooted facts for tree at index {index}: {e}"))
}
#[cfg(test)]
mod rooted_facts_tests;
#[cfg(test)]
mod tests;
