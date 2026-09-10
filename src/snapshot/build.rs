//! One tree in, one [`Snapshot`] out.
//!
//! Two iterative passes over a flat arena fold a 128-bit XOR fingerprint and a
//! leaf run per node, so a node costs `O(1)` and a tree `Θ(n)`.
//!
//! XOR is what makes canonicalisation free: a set and its complement satisfy
//! `label(A′) = total ^ label(A)`, so a bipartition's canonical form is
//! `min(h, h ^ total)`, computed without touching either leaf set. Leaves are
//! numbered in traversal order, which keeps every subtree a contiguous run of
//! `leaf_order` — that run is how [`Snapshot::push_canonical`] recovers a
//! split's leaf set later, once per *distinct* split.

use super::clades::CladeTable;
use super::fingerprint::Fingerprint;
use phylotree::tree::{Node, Tree as PhyloTree, TreeError};
use rustc_hash::FxHashMap;

/// One edge of a tree, identified by a fingerprint instead of a leaf set.
///
/// `first`/`size` locate the subtree's leaves as a contiguous run of
/// [`Snapshot::leaf_order`], which is what lets the canonical leaf set be
/// recovered later — once per *unique* split — with the tree already dropped.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Part {
    /// Canonical fingerprint: `min(h, h ^ total)` for an internal bipartition,
    /// so both of its sides give the same value without either leaf set being
    /// touched. Pendant edges and rooted clades keep the raw subtree
    /// fingerprint.
    pub(super) key: Fingerprint,
    /// Where this subtree's leaves start in [`Snapshot::leaf_order`].
    pub(super) first: u32,
    /// Leaves in the subtree, before canonicalisation. Load-bearing: the
    /// pendant and near-full filters test it, so they stay integer comparisons
    /// rather than counting a leaf set.
    pub(super) size: u32,
    /// Length of the edge above the node.
    pub(super) length: f64,
}

/// One tree's bipartitions, on the way into a [`Snapshots`] collection.
///
/// Short-lived: built per tree, handed straight to the interner, dropped. That
/// is what keeps construction peak memory near the deduplicated footprint.
/// Taxon names live once on [`Snapshots`], not per tree.
///
/// # Canonicalization
/// Each internal bipartition has two complementary sides. The fingerprint
/// canonicalises them for free by taking `min(h, h ^ total)`; the leaf set the
/// interner records is still the side without leaf 0, so the exported
/// bipartition table is unchanged. Pendant edges are stored as-is.
#[derive(Debug, Clone)]
pub(crate) struct Snapshot {
    /// All edges (internal bipartitions + pendant edges), filtered and sorted
    /// by fingerprint.
    pub(crate) parts: Vec<Part>,

    /// Leaf bit indices in traversal order. Every subtree occupies a
    /// contiguous run, so `leaf_order[first..first + size]` is a part's leaves.
    pub(crate) leaf_order: Vec<u32>,

    /// `u64` words needed to bit-pack one leaf set, for the export format.
    pub(crate) words: usize,
}

impl Snapshot {
    /// Extract a snapshot from an already-parsed tree.
    ///
    /// # Algorithm
    /// 1. One DFS from the root, folding a 128-bit fingerprint and a leaf count
    ///    per node
    /// 2. Drop the trivial splits and canonicalize by fingerprint
    /// 3. Sort by fingerprint, merging any duplicate
    ///
    /// `leaf_index` gives each taxon its bit position and `labels` holds one
    /// entry per taxon at that same position. Both are properties of the *run*,
    /// not of this tree: fingerprints are only comparable across trees if every
    /// tree draws from the same tables — see [`build_leaf_index`] and
    /// [`taxon_labels`].
    ///
    /// # Errors
    /// Returns `TreeError` if the tree is empty, malformed, has unnamed leaves,
    /// or carries a taxon absent from `leaf_index`.
    pub(crate) fn from_tree(
        tree: &PhyloTree,
        rooted_mode: bool,
        labels: &[Fingerprint],
        leaf_index: &FxHashMap<&str, usize>,
    ) -> Result<Self, TreeError> {
        let num_leaves = leaf_index.len();

        // Step 1: fingerprint every node
        let root_id = tree.get_root()?;
        let (leaf_order, parts) =
            Self::collect_partitions(tree, root_id, leaf_index, labels, rooted_mode)?;

        // Steps 2-3
        // rooted_mode=false: bipartitions (canonicalized, deduped, trivial-filtered)
        // rooted_mode=true:  clades (raw subtree fingerprints, just sorted)
        let parts = Self::canonicalize_partitions(parts, num_leaves, rooted_mode);

        Ok(Snapshot {
            parts,
            leaf_order,
            words: num_leaves.div_ceil(64),
        })
    }

    /// Every node's subtree fingerprint and leaf run, paired with its edge.
    ///
    /// Two iterative passes over one flat arena indexed by node id. The DFS
    /// resolves each node once and keeps the reference, so the passes and the
    /// branch-length read that follow them cost no further lookups.
    ///
    /// The root is skipped — it creates no bipartition. Missing branch lengths
    /// are treated as 0.0.
    fn collect_partitions(
        tree: &PhyloTree,
        root_id: usize,
        leaf_index: &FxHashMap<&str, usize>,
        labels: &[Fingerprint],
        rooted_mode: bool,
    ) -> Result<(Vec<u32>, Vec<Part>), TreeError> {
        let mut order: Vec<(usize, &Node)> = Vec::with_capacity(tree.size());
        let mut stack = vec![root_id];
        while let Some(id) = stack.pop() {
            let node = tree.get(&id)?;
            order.push((id, node));
            stack.extend(node.children.iter().copied());
        }

        let mut acc = vec![Acc::default(); tree.size()];
        let mut leaf_order = Vec::with_capacity(leaf_index.len());

        // Pre-order: each leaf takes the next slot, which is what keeps a
        // subtree's leaves adjacent.
        for &(id, node) in order.iter().filter(|(_, n)| n.children.is_empty()) {
            let name = node.name.as_deref().ok_or(TreeError::UnnamedLeaves)?;
            let bit = *leaf_index.get(name).ok_or(TreeError::DifferentTipIndices)?;
            acc[id] = Acc {
                fp: labels[bit],
                first: leaf_order.len() as u32,
                size: 1,
            };
            leaf_order.push(bit as u32);
        }

        // Post-order: an internal node is the XOR of its children.
        for &(id, node) in order.iter().rev().filter(|(_, n)| !n.children.is_empty()) {
            let folded = node.children.iter().fold(Acc::empty(), |a, &c| Acc {
                fp: a.fp ^ acc[c].fp,
                first: a.first.min(acc[c].first),
                size: a.size + acc[c].size,
            });
            acc[id] = folded;
        }

        let total = acc[root_id].fp;
        let parts = order
            .into_iter()
            .filter(|&(id, _)| id != root_id)
            .map(|(id, node)| {
                let Acc { fp, first, size } = acc[id];
                Part {
                    // A pendant keeps its raw fingerprint: its complement is
                    // filtered out below, so nothing else can claim the key,
                    // and at two taxa `min(h, h ^ total)` would collapse the
                    // two pendants onto each other.
                    key: if rooted_mode || size == 1 {
                        fp
                    } else {
                        fp.min(fp ^ total)
                    },
                    first,
                    size,
                    length: node.parent_edge.unwrap_or(0.0),
                }
            })
            .collect();

        Ok((leaf_order, parts))
    }

    /// Drop the trivial parts, then sort by fingerprint and merge duplicates.
    ///
    /// Three cases, decided by the raw subtree size, which the DFS already
    /// accumulated:
    ///   size == 1              → pendant edge, kept as-is.
    ///   size >= num_leaves - 1 → the complement of a pendant; dropped, since
    ///                            the pendant itself is kept above.
    ///   otherwise              → internal bipartition. `key` already carries
    ///                            `min(h, h ^ total)`, so both sides of a split
    ///                            arrive equal and there is nothing to flip.
    ///
    /// Sorting compares a 16-byte key, independent of taxon count.
    fn canonicalize_partitions(
        mut parts: Vec<Part>,
        num_leaves: usize,
        rooted_mode: bool,
    ) -> Vec<Part> {
        if rooted_mode {
            // Rooted clade mode: raw subtree fingerprints as clades.
            // No canonicalization, no trivial filter, no dedup.
            // Both root children are distinct clades; L-1 leaf clades are valid.
            parts.sort_unstable_by_key(|p| p.key);
            return parts;
        }

        parts.retain(|p| p.size == 1 || (p.size as usize) < num_leaves - 1);
        parts.sort_unstable_by_key(|p| p.key);

        // Deduplicate identical fingerprints (e.g. the root bipartition in a
        // rooted binary tree, where both root children canonicalize to the same
        // split). Merge by summing branch lengths.
        parts.dedup_by(|later, kept| {
            later.key == kept.key && {
                kept.length += later.length;
                true
            }
        });

        parts
    }

    /// Append one part's canonical leaf set to `table`.
    ///
    /// The canonical side is the one without leaf 0, and taking the complement
    /// costs no set arithmetic: leaves are numbered in traversal order, so a
    /// subtree owns a contiguous run of `leaf_order` and everything *else* is
    /// simply the rest of the array.
    ///
    /// A pendant edge is kept as-is, leaf 0 or not — its complement was filtered
    /// out, so no other occurrence can disagree about which side is canonical.
    pub(super) fn push_canonical(&self, part: &Part, rooted_mode: bool, table: &mut CladeTable) {
        let (lo, hi) = (
            part.first as usize,
            part.first as usize + part.size as usize,
        );
        let run = &self.leaf_order[lo..hi];

        if !rooted_mode && part.size > 1 && run.contains(&0) {
            table.push(
                self.leaf_order[..lo]
                    .iter()
                    .chain(&self.leaf_order[hi..])
                    .copied(),
            );
        } else {
            table.push(run.iter().copied());
        }
    }
}

/// Per-node DFS accumulator: a fingerprint and the subtree's leaf run.
#[derive(Debug, Clone, Copy, Default)]
struct Acc {
    fp: Fingerprint,
    pub(super) first: u32,
    pub(super) size: u32,
}

impl Acc {
    /// Identity for the children fold: XOR's zero, and a `first` that any real
    /// child wins the `min` against.
    fn empty() -> Self {
        Acc {
            fp: 0,
            first: u32::MAX,
            size: 0,
        }
    }
}
