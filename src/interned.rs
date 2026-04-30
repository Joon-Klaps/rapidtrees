//! Interned tree snapshots: every distinct bipartition stored once globally,
//! each tree represented as a sorted `Vec<u32>` of split IDs.
//!
//! # Why this exists
//!
//! `TreeSnapshot` (in [`crate::snapshot`]) stores each tree's bipartitions as a
//! sorted `Vec<Bitset>`. At small taxa counts that's fine, but at the
//! 1000–2000-taxa scale a single bitset is 16–32 `u64` words and `Bitset::cmp`
//! becomes a multi-word `memcmp`. Worse, the snapshot vec dwarfs L3 cache:
//! 5000 trees × ~600 KB ≈ 3 GB of pointer-chased data, so the RF inner loop
//! spends most of its time waiting on DRAM.
//!
//! Interning replaces every `Bitset` in every snapshot with a `u32` ID into a
//! single global table. The RF inner loop then runs on `Vec<u32>` — single
//! integer compare per cmp, full working set fits in L1. Empirically this is
//! the load-bearing optimization for large-taxa workloads.
//!
//! # Algorithm
//!
//! 1. Walk every snapshot's `parts` list, assigning each unique bipartition a
//!    fresh `u32` ID via `FxHashMap<Bitset, u32>`.
//! 2. For each tree, collect its `(split_id, branch_length)` pairs and **sort
//!    by split_id** (paired so lengths stay aligned). The sorted-merge RF
//!    algorithm requires a consistent total order on splits across trees;
//!    numeric u32 ordering is the cheapest possible.
//! 3. RF becomes the same merge as before, but on `&[u32]` rather than
//!    `&[Bitset]` — equality and ordering are now O(1) instead of O(words).
//!
//! Output is bit-identical to [`crate::distances::pairwise_rf_matrix`] —
//! interning is a representation change, not an algorithmic one. The included
//! `rf_matches_legacy` test gates this invariant.

use crate::bitset::Bitset;
use crate::distances::rf_from_snapshots;
use crate::snapshot::TreeSnapshot;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

/// One snapshot in interned form.
///
/// `split_ids` is sorted ascending so the RF sorted-merge can run on integers.
/// `lengths[i]` is the branch length of `split_ids[i]` (parallel arrays).
#[derive(Debug, Clone)]
pub struct InternedSnapshot {
    pub split_ids: Vec<u32>,
    pub lengths: Vec<f64>,
    pub root_children_ids: Vec<u32>,
    pub rooted: bool,
}

/// A collection of interned snapshots backed by a single global bipartition
/// table.
///
/// `bipartitions[i]` is the `Bitset` for split ID `i`. `bipartition_index`
/// maps each `Bitset` back to its ID — useful when a caller wants to look up
/// split membership from outside.
#[derive(Debug)]
pub struct InternedSnapshots {
    pub snapshots: Vec<InternedSnapshot>,
    pub bipartitions: Vec<Bitset>,
    pub bipartition_index: FxHashMap<Bitset, u32>,
    pub words_per_bitset: usize,
    pub leaf_names: Vec<String>,
}

impl InternedSnapshots {
    /// Build the interning table from a previously-computed `Vec<TreeSnapshot>`.
    ///
    /// Walks every snapshot's `parts` list, assigning each new bipartition a
    /// fresh u32 ID. Snapshots come out with `split_ids` sorted ascending and
    /// `lengths` re-aligned to match.
    ///
    /// `leaf_names` is just stored verbatim; it indexes the bits inside each
    /// `Bitset` (alphabetically-sorted taxon order, established by
    /// `TreeSnapshot::from_tree`).
    pub fn from_snapshots(snaps: Vec<TreeSnapshot>, leaf_names: Vec<String>) -> Self {
        let words = snaps.first().map(|s| s.words).unwrap_or(0);
        let mut index: FxHashMap<Bitset, u32> = FxHashMap::default();
        let mut bipartitions: Vec<Bitset> = Vec::new();

        let interned: Vec<InternedSnapshot> = snaps
            .into_iter()
            .map(|snap| {
                let split_ids: Vec<u32> = snap
                    .parts
                    .iter()
                    .map(|b| match index.get(b) {
                        Some(&id) => id,
                        None => {
                            let id = bipartitions.len() as u32;
                            index.insert(b.clone(), id);
                            bipartitions.push(b.clone());
                            id
                        }
                    })
                    .collect();

                // Re-pair so lengths follow split_ids through the sort.
                let mut paired: Vec<(u32, f64)> =
                    split_ids.into_iter().zip(snap.lengths).collect();
                paired.sort_unstable_by_key(|&(id, _)| id);
                let (split_ids, lengths): (Vec<u32>, Vec<f64>) = paired.into_iter().unzip();

                // Root children's bitsets are guaranteed to be in `index` by
                // construction (they were inserted alongside `parts` when
                // present, and trivial leaf clades may not be there in unrooted
                // mode — those silently drop, mirroring the existing
                // root-children semantics in TreeSnapshot).
                let root_children_ids: Vec<u32> = snap
                    .root_children
                    .iter()
                    .filter_map(|b| index.get(b).copied())
                    .collect();

                InternedSnapshot {
                    split_ids,
                    lengths,
                    root_children_ids,
                    rooted: snap.rooted,
                }
            })
            .collect();

        Self {
            snapshots: interned,
            bipartitions,
            bipartition_index: index,
            words_per_bitset: words,
            leaf_names,
        }
    }
}

/// Robinson–Foulds distance between two interned snapshots.
///
/// Identical structure to [`crate::distances::rf_from_snapshots`] but the
/// merge runs on `Vec<u32>` instead of `Vec<Bitset>` — every cmp is a single
/// u32 compare instead of a multi-word memcmp.
#[inline]
pub fn rf_from_interned(a: &InternedSnapshot, b: &InternedSnapshot) -> usize {
    let mut inter = 0;
    let mut i = 0;
    let mut j = 0;
    while i < a.split_ids.len() && j < b.split_ids.len() {
        match a.split_ids[i].cmp(&b.split_ids[j]) {
            std::cmp::Ordering::Equal => {
                inter += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
    a.split_ids.len() + b.split_ids.len() - 2 * inter
}

/// Compute all pairwise RF distances over an [`InternedSnapshots`].
///
/// Mirrors the structure of [`crate::distances::pairwise_rf_matrix`] —
/// parallel upper-triangle walk via rayon, mirrored into the lower triangle.
pub fn pairwise_rf_matrix_interned(snaps: &InternedSnapshots) -> Vec<Vec<usize>> {
    let n = snaps.snapshots.len();
    let pairs: Vec<(usize, usize, usize)> = (0..n)
        .into_par_iter()
        .flat_map_iter(move |i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| {
            let dist = rf_from_interned(&snaps.snapshots[i], &snaps.snapshots[j]);
            (i, j, dist)
        })
        .collect();

    let mut matrix = vec![vec![0usize; n]; n];
    for (i, j, dist) in pairs {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }
    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distances::pairwise_rf_matrix;
    use phylotree::tree::Tree as PhyloTree;

    // First three trees from the treedist fixture (same as
    // distances::pairwise_tests). Known RF: [0,1]=4, [0,2]=2, [1,2]=2.
    const T0: &str = "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";
    const T1: &str = "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);";
    const T2: &str = "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";

    fn three_snaps() -> Vec<TreeSnapshot> {
        [T0, T1, T2]
            .iter()
            .map(|nwk| {
                let tree = PhyloTree::from_newick(nwk).unwrap();
                TreeSnapshot::from_tree(&tree, false).unwrap()
            })
            .collect()
    }

    #[test]
    fn rf_matches_legacy() {
        let snaps = three_snaps();
        let leaf_names: Vec<String> =
            ('A'..='J').map(|c| c.to_string()).collect();

        let legacy = pairwise_rf_matrix(&snaps);
        let interned =
            InternedSnapshots::from_snapshots(snaps, leaf_names);
        let new = pairwise_rf_matrix_interned(&interned);

        assert_eq!(legacy, new, "interned RF must match legacy element-wise");

        // Sanity: known values from the treedist fixture.
        assert_eq!(new[0][1], 4);
        assert_eq!(new[0][2], 2);
        assert_eq!(new[1][2], 2);
    }

    #[test]
    fn split_ids_are_sorted_per_snapshot() {
        let snaps = three_snaps();
        let leaf_names: Vec<String> =
            ('A'..='J').map(|c| c.to_string()).collect();
        let interned = InternedSnapshots::from_snapshots(snaps, leaf_names);
        for snap in &interned.snapshots {
            for w in snap.split_ids.windows(2) {
                assert!(w[0] < w[1], "split_ids must be strictly ascending (no dup IDs in a tree)");
            }
        }
    }

    #[test]
    fn lengths_remain_aligned() {
        let snaps = three_snaps();
        let leaf_names: Vec<String> =
            ('A'..='J').map(|c| c.to_string()).collect();
        let interned = InternedSnapshots::from_snapshots(snaps, leaf_names);
        for snap in &interned.snapshots {
            assert_eq!(snap.split_ids.len(), snap.lengths.len());
        }
    }

    /// Suppress the "unused import" warning since rf_from_snapshots is brought
    /// in for parity context but only directly invoked by `pairwise_rf_matrix`
    /// inside `crate::distances`. The compile-time check is enough.
    #[allow(dead_code)]
    fn _unused_proof(_: fn(&TreeSnapshot, &TreeSnapshot) -> usize) {}
    #[allow(dead_code)]
    const _PROOF: fn() = || _unused_proof(rf_from_snapshots);
}
