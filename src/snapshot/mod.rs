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

use build::{Part, Snapshot};
use fingerprint::{build_leaf_index, taxon_labels};
use intern::Interner;

pub(crate) use intern::InternSnap;

use crate::par::*;
use clades::CladeTable;
use std::collections::HashMap;

/// Bytes of tree text and raw parts parsed per chunk; see
/// [`Snapshots::from_newick_iter_opts`].
///
/// The chunk is held whole until it is interned, so below about 50 000 taxa
/// it, not the split table, sets peak construction memory. 16 MiB is where
/// shrinking it stops paying: a smaller chunk saves little more, and from
/// about 10^5 taxa a chunk is already at its floor of one tree per thread.
const CHUNK_TARGET_BYTES: usize = 16 * 1024 * 1024;

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
}

/// What to keep besides the split IDs themselves.
///
/// Neither is read by the distance kernels, and both cost real work: branch
/// lengths are `~n_bip × 8` bytes per tree, and materialising a bipartition is
/// an `O(subtree)` walk per distinct split. A path that computes a matrix and
/// exports nothing wants both off.
#[derive(Debug, Clone, Copy)]
pub struct Retain {
    /// Per-edge branch lengths. Required by WRF and KF; dead weight for RF.
    pub lengths: bool,
    /// The canonical leaf set per distinct split. Required only to export a
    /// bipartition table to Python or to order export columns.
    pub bipartitions: bool,
}

impl Retain {
    /// Keep everything — what the public constructors use.
    pub fn everything() -> Self {
        Self {
            lengths: true,
            bipartitions: true,
        }
    }

    /// Keep only what a distance matrix needs
    pub fn for_distances(weighted: bool) -> Self {
        Self {
            lengths: weighted,
            bipartitions: false,
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
    /// not read, and takes the Newick text as anything that derefs to `str`.
    ///
    /// Entries are pulled lazily, one chunk of trees at a time, and a chunk's
    /// text is dropped as soon as it is parsed. A caller that hands over owned
    /// strings from a lazy source (a file being read, a Python iterator)
    /// therefore never holds more than one chunk of text. See [`Retain`] for
    /// what can be skipped; [`Snapshots::from_newick_iter`] retains everything,
    /// so callers that do not opt in are unaffected.
    pub fn from_newick_iter_opts<'a, T>(
        entries: impl IntoIterator<Item = (T, &'a HashMap<String, String>)>,
        rooted: bool,
        retain: Retain,
    ) -> Result<Self, String>
    where
        T: AsRef<str> + Send + Sync,
    {
        Self::from_chunks(entries.into_iter(), rooted, retain, CHUNK_TARGET_BYTES)
    }

    /// [`Self::from_newick_iter_opts`] with the chunk budget given.
    ///
    /// A chunk is `budget` bytes of text and raw parts, estimated from tree 0
    /// (every tree shares its leaf set, so it is representative), rounded to a
    /// whole number of rounds of the thread pool: never fewer trees than
    /// threads, or wide trees would parse on a few cores while the rest idle.
    fn from_chunks<'a, T>(
        mut entries: impl Iterator<Item = (T, &'a HashMap<String, String>)>,
        rooted: bool,
        retain: Retain,
        budget: usize,
    ) -> Result<Self, String>
    where
        T: AsRef<str> + Send + Sync,
    {
        let Some((first, first_translate)) = entries.next() else {
            return Ok(Self::empty());
        };
        let first_newick = first.as_ref();

        // Tree 0 defines the run's taxa, so its names are read first and on
        // their own: every tree's leaf check, tree 0's included, needs the very
        // table they are about to build. An unnamed leaf fails inside
        // `leaf_names`; a repeated one is caught here.
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
        };

        let first_snap = newick::snapshot(first_newick, first_translate, 0, &run)?;
        let per_tree = first_snap.parts.len() * size_of::<Part>()
            + first_snap.leaf_order.len() * size_of::<u32>()
            + first_newick.len();
        let threads = current_num_threads().max(1);
        let rounds = (budget / (per_tree.max(1) * threads)).clamp(1, (4096 / threads).max(1));
        let chunk = rounds * threads;
        drop(first);

        let mut interner = Interner::new(
            entries.size_hint().0 + 1,
            first_snap.words,
            sorted_leaf_names.len(),
            rooted,
            retain,
        );
        interner.push_chunk(vec![first_snap])?;

        let mut base = 1usize; // tree 0 is already interned
        loop {
            let batch: Vec<_> = entries.by_ref().take(chunk).collect();
            if batch.is_empty() {
                break;
            }
            // Parse this chunk in parallel, validating leaf sets. Its text is
            // not needed past the parse.
            let raw: Vec<Snapshot> = batch
                .par_iter()
                .enumerate()
                .map(|(k, (newick, translate))| {
                    newick::snapshot(newick.as_ref(), translate, base + k, &run)
                })
                .collect::<Result<_, _>>()?;
            base += batch.len();
            drop(batch);

            // The chunk's raw parts are dropped as soon as it is interned, so
            // peak stays near the deduplicated footprint.
            interner.push_chunk(raw)?;
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

    /// A collection of no trees.
    pub(crate) fn empty() -> Self {
        Self {
            snapshots: Vec::new(),
            clades: CladeTable::new(),
            split_counts: Vec::new(),
            words_per_bitset: 0,
            leaf_names: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests;
