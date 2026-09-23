//! What one tree becomes on its way into a [`Snapshots`](super::Snapshots) collection.
//!
//! [`super::newick`] builds these straight from the text.
//!
//! XOR is what makes canonicalisation free: a set and its complement satisfy
//! `label(A′) = total ^ label(A)`, so a bipartition's canonical form is
//! `min(h, h ^ total)`, computed without touching either leaf set. Leaves are
//! numbered in traversal order, which keeps every subtree a contiguous run of
//! `leaf_order` — that run is how [`Snapshot::push_canonical`] recovers a
//! split's leaf set later, once per *distinct* split.

use super::clades::CladeTable;
use super::fingerprint::Fingerprint;

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

/// One tree's bipartitions, on the way into a [`Snapshots`](super::Snapshots) collection.
///
/// Short-lived: built per tree, handed straight to the interner, dropped. That
/// is what keeps construction peak memory near the deduplicated footprint.
/// Taxon names live once on [`Snapshots`](super::Snapshots), not per tree.
///
/// # Canonicalization
/// Each internal bipartition has two complementary sides. The fingerprint
/// canonicalises them for free by taking `min(h, h ^ total)`; the leaf set the
/// interner records is still the side without leaf 0, so the exported
/// bipartition table is unchanged. Pendant edges are stored as-is.
#[derive(Debug, Clone)]
pub(crate) struct Snapshot {
    /// All edges (internal bipartitions + pendant edges), filtered,
    /// canonicalised and free of duplicates. The reader emits them in
    /// post-order; nothing downstream depends on the order.
    pub(crate) parts: Vec<Part>,

    /// Leaf bit indices in traversal order. Every subtree occupies a
    /// contiguous run, so `leaf_order[first..first + size]` is a part's leaves.
    pub(crate) leaf_order: Vec<u32>,

    /// `u64` words needed to bit-pack one leaf set, for the export format.
    pub(crate) words: usize,
}

impl Snapshot {
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
