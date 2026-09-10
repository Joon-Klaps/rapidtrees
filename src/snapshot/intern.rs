//! Bipartition deduplication: many trees' splits collapse into one table of
//! unique splits, and each tree becomes a sorted list of `u32` IDs into it.
//!
//! This is what lets pairwise distance work on integers: a comparison is a
//! single integer compare, and a large run's working set fits in L2 rather
//! than DRAM.

use super::Retain;
use super::Snapshots;
use super::build::Snapshot;
use super::clades::CladeTable;
use super::fingerprint::Fingerprint;
use hashbrown::HashTable;
use rustc_hash::FxBuildHasher;
use std::hash::BuildHasher;

/// One tree's bipartitions in interned form (private implementation detail of [`Snapshots`]).
///
/// `split_ids` is sorted ascending so the RF sorted-merge can run on integers.
/// `lengths[i]` is the branch length of `split_ids[i]` (parallel arrays).
#[derive(Debug, Clone)]
pub(crate) struct InternSnap {
    pub(crate) split_ids: Vec<u32>,
    pub(crate) lengths: Vec<f64>,
}

/// Incremental bipartition interner.
///
/// Deduplicates bipartitions into split IDs using a `HashTable<u32>` that stores
/// only the integer IDs — each unique leaf set is kept exactly once, in
/// `clades`. Snapshots are folded in one at a time via [`Interner::push`],
/// which consumes each raw `Snapshot` so its parts are freed immediately.
/// This keeps construction peak memory near the deduplicated footprint instead
/// of holding every tree's raw parts at once.
///
/// IDs are assigned in first-seen order, so feeding snapshots in tree-index
/// order yields exactly the same clade ordering and `split_ids` as
/// interning the whole `Vec<Snapshot>` at once.
pub(super) struct Interner {
    hasher: FxBuildHasher,
    table: HashTable<u32>,
    /// `(canonical fingerprint, smaller side's cardinality)` per split ID —
    /// what a candidate is matched against.
    keys: Vec<(Fingerprint, u32)>,
    clades: CladeTable,
    snapshots: Vec<InternSnap>,
    words: usize,
    num_leaves: usize,
    rooted: bool,
    retain: Retain,
}

impl Interner {
    pub(super) fn new(
        n_trees: usize,
        words: usize,
        num_leaves: usize,
        rooted: bool,
        retain: Retain,
    ) -> Self {
        Self {
            hasher: FxBuildHasher,
            table: HashTable::new(),
            keys: Vec::new(),
            clades: CladeTable::new(),
            snapshots: Vec::with_capacity(n_trees),
            words,
            num_leaves,
            rooted,
            retain,
        }
    }

    /// Intern one raw snapshot. When `store_lengths` is `false`, branch lengths
    /// are dropped (RF-only paths never read them).
    ///
    /// A split is matched on its 128-bit fingerprint *and* the cardinality of
    /// its smaller side, which is equal for both sides of a bipartition and
    /// rules out a slice of the collision space for one integer compare. Only a
    /// fingerprint the table has never seen records a leaf set, so that
    /// `O(subtree)` walk runs once per unique split rather than once per
    /// occurrence — and not at all when [`Retain::bipartitions`] is off.
    pub(super) fn push(&mut self, snap: Snapshot) {
        // Bind disjoint fields to locals so the closures below can borrow the
        // table mutably while reading `keys`/`hasher`.
        let hasher = &self.hasher;
        let table = &mut self.table;
        let keys = &mut self.keys;
        let clades = &mut self.clades;
        let (num_leaves, rooted, retain) = (self.num_leaves, self.rooted, self.retain);

        let mut paired: Vec<(u32, f64)> = snap
            .parts
            .iter()
            .map(|part| {
                let candidate = (part.key, part.size.min(num_leaves as u32 - part.size));
                let hash = hasher.hash_one(part.key);
                let id = match table.find(hash, |&id| keys[id as usize] == candidate) {
                    Some(&id) => id,
                    None => {
                        let new_id = keys.len() as u32;
                        keys.push(candidate);
                        if retain.bipartitions {
                            snap.push_canonical(part, rooted, clades);
                        }
                        // register the new ID in the hash table
                        table.insert_unique(hash, new_id, |&id| {
                            hasher.hash_one(keys[id as usize].0)
                        });
                        new_id
                    }
                };
                (id, part.length)
            })
            .collect();

        paired.sort_unstable_by_key(|&(id, _)| id);
        // On RF-only paths skip materialising the lengths column entirely
        // instead of unzipping then discarding it.
        let (split_ids, lengths): (Vec<u32>, Vec<f64>) = if self.retain.lengths {
            paired.into_iter().unzip()
        } else {
            (paired.into_iter().map(|(id, _)| id).collect(), Vec::new())
        };
        self.snapshots.push(InternSnap { split_ids, lengths });
    }

    pub(super) fn finish(self, leaf_names: Vec<String>) -> Snapshots {
        Snapshots {
            snapshots: self.snapshots,
            clades: self.clades,
            words_per_bitset: self.words,
            leaf_names,
        }
    }
}
