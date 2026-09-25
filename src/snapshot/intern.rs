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
use super::rooted_facts::{
    InternedRootedFacts, ROOT_ID, RawCladeRef, RawRootedFacts, RootedFactsBuilder,
};
use hashbrown::HashTable;
use rustc_hash::{FxBuildHasher, FxHashMap};
use std::hash::BuildHasher;

/// One tree's bipartitions in interned form (private implementation detail of [`Snapshots`]).
///
/// `split_ids` follows the order of the tree's parts, not ID order: the kernels
/// and the export index by ID, so nothing downstream needs it sorted.
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
    /// How many trees hold each split ID.
    counts: Vec<u32>,
    clades: CladeTable,
    snapshots: Vec<InternSnap>,
    rooted_facts: Option<RootedFactsBuilder>,
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
            counts: Vec::new(),
            clades: CladeTable::new(),
            snapshots: Vec::with_capacity(n_trees),
            rooted_facts: retain
                .rooted_facts
                .then(|| RootedFactsBuilder::with_capacity(n_trees)),
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
        let counts = &mut self.counts;
        let clades = &mut self.clades;
        let (num_leaves, rooted, retain) = (self.num_leaves, self.rooted, self.retain);

        let split_ids: Vec<u32> = snap
            .parts
            .iter()
            .map(|part| {
                let candidate = (part.key, part.size.min(num_leaves as u32 - part.size));
                let hash = hasher.hash_one(part.key);
                match table.find(hash, |&id| keys[id as usize] == candidate) {
                    Some(&id) => {
                        counts[id as usize] += 1;
                        id
                    }
                    None => {
                        let new_id = keys.len() as u32;
                        keys.push(candidate);
                        counts.push(1);
                        if retain.bipartitions {
                            snap.push_canonical(part, rooted, clades);
                        }
                        // register the new ID in the hash table
                        table.insert_unique(hash, new_id, |&id| {
                            hasher.hash_one(keys[id as usize].0)
                        });
                        new_id
                    }
                }
            })
            .collect();

        // On RF-only paths skip materialising the lengths column entirely.
        let lengths = if retain.lengths {
            snap.parts.iter().map(|part| part.length).collect()
        } else {
            Vec::new()
        };
        self.snapshots.push(InternSnap { split_ids, lengths });
    }

    /// Resolve one tree's optional rooted facts after its ordinary snapshot.
    ///
    /// [`Interner::push`] has already assigned every clade ID for this tree.
    /// This method therefore performs lookup only: it cannot create IDs or
    /// alter first-seen ordering. The completed facts are appended to a
    /// separate tree-aligned sidecar, leaving [`InternSnap`] unchanged.
    pub(super) fn push_rooted_facts(&mut self, raw: RawRootedFacts) -> Result<(), String> {
        if !self.rooted {
            return Err("cannot intern rooted facts in unrooted mode".to_string());
        }
        let sidecar_rows = self
            .rooted_facts
            .as_ref()
            .ok_or_else(|| "rooted-facts retention is not enabled".to_string())?
            .len();
        if self.snapshots.len() != sidecar_rows + 1 {
            return Err(format!(
                "rooted-facts row {sidecar_rows} has no matching newly interned snapshot"
            ));
        }

        let RawRootedFacts {
            nodes,
            root_height,
            splits,
        } = raw;
        let tree_split_ids = &self.snapshots[sidecar_rows].split_ids;
        if nodes.len() != tree_split_ids.len() {
            return Err(format!(
                "rooted facts contain {} non-root nodes but snapshot row {sidecar_rows} contains {} clades",
                nodes.len(),
                tree_split_ids.len()
            ));
        }
        let expected_splits = self.num_leaves.saturating_sub(1);
        if splits.len() != expected_splits {
            return Err(format!(
                "rooted facts contain {} splits but a binary {}-taxon tree requires {expected_splits}",
                splits.len(),
                self.num_leaves
            ));
        }
        let expected_root_splits = usize::from(expected_splits > 0);
        let root_splits = splits.iter().filter(|split| split.parent.is_none()).count();
        if root_splits != expected_root_splits {
            return Err(format!(
                "rooted facts contain {root_splits} root splits; expected {expected_root_splits}"
            ));
        }

        // The direct parser emits both `Snapshot::parts` and raw node facts in
        // the same postorder. Master intentionally keeps interned IDs in that
        // order (they are no longer sorted), so validate the alignment without
        // sorting either side and retain heights in matching order.
        let mut id_by_clade = FxHashMap::default();
        let mut node_heights = Vec::with_capacity(nodes.len());
        for (node, &expected_id) in nodes.into_iter().zip(tree_split_ids) {
            let id = self.resolve_rooted_clade(node.clade)?;
            if id != expected_id {
                return Err(format!(
                    "rooted facts are not aligned with snapshot clades for row {sidecar_rows}"
                ));
            }
            if id_by_clade.insert(node.clade, id).is_some() {
                return Err(format!(
                    "rooted facts contain a duplicate clade in row {sidecar_rows}"
                ));
            }
            node_heights.push(node.height);
        }

        let resolve_in_tree = |clade: RawCladeRef| -> Result<u32, String> {
            id_by_clade.get(&clade).copied().ok_or_else(|| {
                format!("rooted split references a clade absent from snapshot row {sidecar_rows}")
            })
        };

        let mut resolved_splits = Vec::with_capacity(splits.len());
        for split in splits {
            let parent = match split.parent {
                Some(parent) => resolve_in_tree(parent)?,
                None => ROOT_ID,
            };
            let mut children = [
                resolve_in_tree(split.children[0])?,
                resolve_in_tree(split.children[1])?,
            ];
            children.sort_unstable();
            resolved_splits.push([parent, children[0], children[1]]);
        }

        let builder = self
            .rooted_facts
            .as_mut()
            .ok_or_else(|| "rooted-facts retention was disabled during interning".to_string())?;
        let mut split_ids = Vec::with_capacity(resolved_splits.len());
        for split in resolved_splits {
            split_ids.push(builder.intern_split(split)?);
        }
        split_ids.sort_unstable();
        if split_ids.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(format!(
                "rooted facts contain a duplicate observed split in row {sidecar_rows}"
            ));
        }
        builder.push(InternedRootedFacts {
            node_heights,
            root_height,
            split_ids,
        });
        Ok(())
    }

    /// Find the existing global ID for one rooted raw-clade reference.
    fn resolve_rooted_clade(&self, clade: RawCladeRef) -> Result<u32, String> {
        let num_leaves = u32::try_from(self.num_leaves)
            .map_err(|_| "rooted facts cannot represent more than u32::MAX taxa".to_string())?;
        if clade.size == 0 || clade.size >= num_leaves {
            return Err(format!(
                "invalid non-root clade size {} for a {num_leaves}-taxon tree",
                clade.size
            ));
        }

        let candidate = (clade.key, clade.size.min(num_leaves - clade.size));
        let hash = self.hasher.hash_one(clade.key);
        self.table
            .find(hash, |&id| self.keys[id as usize] == candidate)
            .copied()
            .ok_or_else(|| {
                format!(
                    "rooted clade 0x{:032x}/{} was not interned by its snapshot",
                    clade.key, clade.size
                )
            })
    }

    pub(super) fn finish(self, leaf_names: Vec<String>) -> Snapshots {
        debug_assert!(
            self.rooted_facts
                .as_ref()
                .is_none_or(|facts| facts.len() == self.snapshots.len()),
            "rooted facts must stay aligned with snapshot rows"
        );
        let rooted_facts = self.rooted_facts.map(RootedFactsBuilder::finish);
        Snapshots {
            snapshots: self.snapshots,
            clades: self.clades,
            split_counts: self.counts,
            words_per_bitset: self.words,
            leaf_names,
            rooted_facts,
        }
    }
}
