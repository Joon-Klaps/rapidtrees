//! Tree distance functions.
//!
//! # Public API
//! - [`rf_distance`], [`wrf_distance`], [`kf_distance`] — compare two standalone
//!   [`crate::snapshot::Snapshot`]s (use bitset sorted-merge; panics if leaf sets differ).
//!
//! # Bulk pairwise computation
//! Call [`crate::snapshot::Snapshots::pairwise_rf`] etc. — these use the fast
//! interned-integer path internally via the `pub(crate)` helpers in this module.

use crate::snapshot::{Snapshot, Snapshots};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(test)]
use itertools::Itertools;

// ─── Public: Snapshot-based distance functions ────────────────────────────────

/// Robinson–Foulds distance between two snapshots.
///
/// Uses a sorted-merge on `Vec<Bitset>`. RF = |A| + |B| − 2|A ∩ B|.
///
/// # Panics
/// Panics if `a` and `b` have different leaf sets — that is a caller error.
pub fn rf_distance(a: &Snapshot, b: &Snapshot) -> usize {
    assert_eq!(
        a.leaf_names, b.leaf_names,
        "rf_distance: leaf sets differ — cannot compare snapshots from different taxa"
    );
    let mut inter = 0;
    let mut i = 0;
    let mut j = 0;
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Equal => {
                inter += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
    a.parts.len() + b.parts.len() - 2 * inter
}

/// Weighted Robinson–Foulds distance between two snapshots.
///
/// For each bipartition:
/// - Shared: contributes `|len_a − len_b|`.
/// - Only in A: contributes `len_a`.
/// - Only in B: contributes `len_b`.
///
/// # Panics
/// Panics if `a` and `b` have different leaf sets.
pub fn wrf_distance(a: &Snapshot, b: &Snapshot) -> f64 {
    assert_eq!(
        a.leaf_names, b.leaf_names,
        "wrf_distance: leaf sets differ — cannot compare snapshots from different taxa"
    );
    let mut dist = 0.0;
    let mut i = 0;
    let mut j = 0;
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Equal => {
                dist += (a.lengths[i] - b.lengths[j]).abs();
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                dist += a.lengths[i];
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                dist += b.lengths[j];
                j += 1;
            }
        }
    }
    dist += a.lengths[i..].iter().sum::<f64>();
    dist += b.lengths[j..].iter().sum::<f64>();
    dist
}

/// Kuhner–Felsenstein (Branch Score) distance between two snapshots.
///
/// Like WRF but uses squared differences: `sqrt(Σ (len_a − len_b)²)`.
/// More sensitive to large branch-length differences; a Euclidean metric in
/// branch-length space.
///
/// # Panics
/// Panics if `a` and `b` have different leaf sets.
pub fn kf_distance(a: &Snapshot, b: &Snapshot) -> f64 {
    assert_eq!(
        a.leaf_names, b.leaf_names,
        "kf_distance: leaf sets differ — cannot compare snapshots from different taxa"
    );
    let mut sum_sq = 0.0;
    let mut i = 0;
    let mut j = 0;
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Equal => {
                let d = a.lengths[i] - b.lengths[j];
                sum_sq += d * d;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                sum_sq += a.lengths[i] * a.lengths[i];
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                sum_sq += b.lengths[j] * b.lengths[j];
                j += 1;
            }
        }
    }
    sum_sq += a.lengths[i..].iter().map(|&l| l * l).sum::<f64>();
    sum_sq += b.lengths[j..].iter().map(|&l| l * l).sum::<f64>();
    sum_sq.sqrt()
}

// ─── pub(crate): bulk pairwise paths used by Snapshots::pairwise_rf/wrf/kf ──
//
// All three metrics lay every tree's splits out as one contiguous row — packed
// bits for RF, branch lengths for WRF/KF — and then sweep the upper triangle
// once. RF and KF share the same shape: "what each tree has on its own, minus
// twice what the pair has in common", where "in common" is a popcount for RF
// and a dot product for KF. WRF is the odd one out (an absolute difference has
// no shared-term shortcut), so it just sweeps the two length rows directly.

/// Fill a symmetric `n × n` distance matrix in parallel.
///
/// Every distance shows up twice — the matrix is a mirror image across its
/// diagonal — so we only ever compute the upper triangle (one rayon task per
/// row) and then copy each value across to its mirror slot. The diagonal stays
/// at zero: a tree's distance to itself. `cell(i, j)` returns the distance for
/// the tree pair `(i, j)`.
///
/// If `progress` is given, it's bumped by that row's pair count the moment the
/// row finishes, so another thread can watch the total climb to `n*(n-1)/2`.
fn fill_symmetric<T, F>(n: usize, progress: Option<&AtomicUsize>, cell: F) -> Vec<T>
where
    T: Copy + Default + Send,
    F: Fn(usize, usize) -> T + Sync,
{
    let mut matrix = vec![T::default(); n * n];

    matrix.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        for (j, slot) in row.iter_mut().enumerate().skip(i + 1) {
            *slot = cell(i, j);
        }
        if let Some(counter) = progress {
            counter.fetch_add(n.saturating_sub(i + 1), Ordering::Relaxed);
        }
    });

    for i in 0..n {
        for j in (i + 1)..n {
            matrix[j * n + i] = matrix[i * n + j];
        }
    }

    matrix
}

/// How many distinct splits exist across all the trees.
///
/// Split IDs are handed out densely from 0 and each tree's list is sorted, so
/// the biggest ID anywhere is the largest last element; +1 turns that into a
/// count. We read it from `split_ids` (never the bitset table), so it still
/// works on the paths that free the bitsets to save memory.
pub(crate) fn n_distinct_splits(snaps: &Snapshots) -> usize {
    snaps
        .snapshots
        .iter()
        .filter_map(|s| s.split_ids.last().copied())
        .max()
        .map(|max_id| max_id as usize + 1)
        .unwrap_or(0)
}

/// Borrow row `i` of a flat, row-major matrix whose rows are `stride` wide.
#[inline]
fn row_slice<T>(flat: &[T], i: usize, stride: usize) -> &[T] {
    &flat[i * stride..][..stride]
}

// ─── pub(crate): dense bit-packed popcount RF (the hot path) ────────────────

/// Work out the packed bit slot for every split, shared by the CPU (u64) and GPU
/// (u32) RF packers so both agree on which splits are dropped and where each lands.
///
/// We count how many trees each split appears in; a split present in all `n` trees
/// is an "everywhere" split that cancels out of RF, so it gets no slot. Every kept
/// split is handed a fresh, packed-together bit slot. Returns `(bit_slot, kept,
/// everywhere)` where `bit_slot[id]` is the slot index (or `u32::MAX` if dropped),
/// `kept` is the number of slots, and `everywhere` is the number dropped.
pub(crate) fn split_bit_layout(snaps: &Snapshots) -> (Vec<u32>, usize, usize) {
    let n = snaps.snapshots.len();
    let n_splits = n_distinct_splits(snaps);

    let mut tree_count = vec![0u32; n_splits];
    for snap in &snaps.snapshots {
        for &id in &snap.split_ids {
            tree_count[id as usize] += 1;
        }
    }

    let mut bit_slot = vec![u32::MAX; n_splits];
    let mut kept = 0u32;
    for (id, &count) in tree_count.iter().enumerate() {
        if count < n as u32 {
            bit_slot[id] = kept;
            kept += 1;
        }
    }
    let kept = kept as usize;
    let everywhere = n_splits - kept;
    (bit_slot, kept, everywhere)
}

/// Compute every pairwise Robinson–Foulds distance with bit-packed popcounts.
///
/// The trick: the RF distance between two trees is just "how many splits each
/// one has, minus twice the splits they share." We write that as
///
/// ```text
///   RF(i, j) = aᵢ + aⱼ − 2 · (splits that i and j share)
/// ```
///
/// and get the shared count by AND-ing the two trees' presence bitmasks and
/// counting the set bits (`popcount`). That's much faster than walking two
/// sorted lists.
///
/// One extra simplification: a split that shows up in *every* tree (always the
/// pendant/leaf edges, plus anything the whole set happens to agree on) adds the
/// same amount to both `a` values and to the shared count, so it cancels out of
/// RF entirely. We drop those "everywhere" splits before packing, which makes
/// the bitmasks shorter and the popcounts cheaper without changing any distance.
pub(crate) fn pairwise_rf_packed(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<usize> {
    let n = snaps.snapshots.len();
    if n == 0 {
        return Vec::new();
    }

    let (bit_slot, kept, everywhere) = split_bit_layout(snaps);
    let words = kept.div_ceil(64);

    // Build one bitmask row per tree: a set bit means "this tree has that split."
    let mut packed = vec![0u64; n * words];
    if words > 0 {
        packed
            .par_chunks_mut(words)
            .zip(&snaps.snapshots)
            .for_each(|(row, snap)| {
                for &id in &snap.split_ids {
                    let slot = bit_slot[id as usize];
                    if slot != u32::MAX {
                        let slot = slot as usize;
                        row[slot >> 6] |= 1u64 << (slot & 63);
                    }
                }
            });
    }

    // Splits per tree once the everywhere-splits are ignored. Every tree holds
    // all of them, so it's just each tree's total minus that fixed count.
    let kept_per_tree: Vec<usize> = snaps
        .snapshots
        .iter()
        .map(|snap| snap.split_ids.len() - everywhere)
        .collect();

    fill_symmetric(n, progress, |i, j| {
        let row_i = row_slice(&packed, i, words);
        let row_j = row_slice(&packed, j, words);
        let shared: usize = row_i
            .iter()
            .zip(row_j)
            .map(|(&x, &y)| (x & y).count_ones() as usize)
            .sum();
        kept_per_tree[i] + kept_per_tree[j] - 2 * shared
    })
}

// ─── pub(crate): dense branch-length WRF / KF (the weighted hot paths) ──────

/// Lay every tree's branch lengths out as one long row of `f64`s, where column
/// `k` holds the length of split `k` (0 if that tree doesn't have the split).
///
/// This is the weighted-metric twin of the RF bit-packing: instead of one bit
/// per split we store the actual branch length, so the inner loops below can
/// read two trees' lengths straight from contiguous memory. Returns the flat
/// `n × n_splits` matrix and `n_splits` (the row stride).
///
/// Unlike RF, we keep *every* split column: an "everywhere" split still has a
/// different length in each tree, so it doesn't cancel out of a weighted score.
pub(crate) fn dense_length_rows(snaps: &Snapshots) -> (Vec<f64>, usize) {
    let n = snaps.snapshots.len();
    let n_splits = n_distinct_splits(snaps);

    let mut rows = vec![0.0f64; n * n_splits];
    if n_splits > 0 {
        rows.par_chunks_mut(n_splits)
            .zip(&snaps.snapshots)
            .for_each(|(row, snap)| {
                for (&id, &length) in snap.split_ids.iter().zip(&snap.lengths) {
                    row[id as usize] = length;
                }
            });
    }
    (rows, n_splits)
}

/// Compute every pairwise Weighted Robinson–Foulds distance on dense length rows.
///
/// WRF just adds up, split by split, how far apart the two trees' branch lengths
/// are: `Σ |lenᵢ − lenⱼ|`. A split missing from one tree counts as length 0
/// there, which the dense layout hands us for free. Walking the two rows side by
/// side — no sorted-merge, no per-split branching — is friendlier to the CPU and
/// vectorises well.
///
/// WRF is the one metric that can't reuse RF/KF's "shared term" shortcut: an
/// absolute difference doesn't factor into a single matrix product, so this
/// stays an honest per-pair sweep — just a faster-laid-out one.
pub(crate) fn pairwise_wrf_dense(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    let n = snaps.snapshots.len();
    if n == 0 {
        return Vec::new();
    }

    let (rows, n_splits) = dense_length_rows(snaps);

    fill_symmetric(n, progress, |i, j| {
        let row_i = row_slice(&rows, i, n_splits);
        let row_j = row_slice(&rows, j, n_splits);
        row_i.iter().zip(row_j).map(|(&a, &b)| (a - b).abs()).sum()
    })
}

/// Compute every pairwise Kuhner–Felsenstein distance with the same
/// "selfᵢ + selfⱼ − 2·shared" trick RF uses, adapted for branch lengths.
///
/// KF is the straight-line (Euclidean) distance between two trees' branch
/// lengths: `sqrt(Σ (lenᵢ − lenⱼ)²)`. Multiplying out the bracket rewrites the
/// sum as
///
/// ```text
///   Σ lenᵢ²  +  Σ lenⱼ²  −  2 · Σ (lenᵢ · lenⱼ)
/// ```
///
/// The first two pieces are each tree's own "self length" (`norm_sq`), which we
/// work out once up front. The last piece is the dot product of the two length
/// rows — the weighted echo of RF's shared-split popcount. Tiny rounding error
/// can nudge the bracket a hair below zero for near-identical trees, so we clamp
/// at zero before taking the square root.
pub(crate) fn pairwise_kf_dense(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    let n = snaps.snapshots.len();
    if n == 0 {
        return Vec::new();
    }

    let (rows, n_splits) = dense_length_rows(snaps);

    // Each tree's own squared length, summed once so the pair loop only needs
    // the dot product between the two rows.
    let norm_sq: Vec<f64> = (0..n)
        .map(|i| row_slice(&rows, i, n_splits).iter().map(|&l| l * l).sum())
        .collect();

    fill_symmetric(n, progress, |i, j| {
        let row_i = row_slice(&rows, i, n_splits);
        let row_j = row_slice(&rows, j, n_splits);
        let dot: f64 = row_i.iter().zip(row_j).map(|(&a, &b)| a * b).sum();
        (norm_sq[i] + norm_sq[j] - 2.0 * dot).max(0.0).sqrt()
    })
}

// ─── pub: GPU-aware dispatch ──────────────────────────────────────────────────
//
// These are the GPU-or-CPU entry points used by api.rs (Python) and main.rs
// (CLI). When the `gpu` feature is not compiled in, or no adapter is available,
// or n < GPU_THRESHOLD, they fall through to the CPU path silently.

/// Generate a `dispatch_*` entry point that tries the GPU when `use_gpu` is true
/// and the `gpu` feature is compiled in, falling back to the CPU path otherwise.
macro_rules! gpu_dispatch {
    ($(#[$meta:meta])* $name:ident -> $ret:ty, gpu = $gpu:path, cpu = $cpu:path) => {
        $(#[$meta])*
        pub fn $name(
            snaps: &Snapshots,
            progress: Option<&std::sync::atomic::AtomicUsize>,
            use_gpu: bool,
        ) -> $ret {
            #[cfg(feature = "gpu")]
            if use_gpu && let Some(result) = $gpu(snaps) {
                return result;
            }
            let _ = use_gpu;
            $cpu(snaps, progress)
        }
    };
}

gpu_dispatch! {
    dispatch_rf -> Vec<usize>, gpu = crate::gpu::try_pairwise_rf, cpu = pairwise_rf_packed
}

gpu_dispatch! {
    dispatch_wrf -> Vec<f64>, gpu = crate::gpu::try_pairwise_wrf, cpu = pairwise_wrf_dense
}

gpu_dispatch! {
    dispatch_kf -> Vec<f64>, gpu = crate::gpu::try_pairwise_kf, cpu = pairwise_kf_dense
}

/// Twelve 10-taxon trees from the PHYLIP treedist reference suite.
///
/// Ground-truth RF, WRF, and KF distances verified against
/// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
#[cfg(test)]
const TREEDIST_TREES: [&str; 12] = [
    "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
];

#[test]
// Robinson–Foulds distances according to
// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
fn robinson_foulds_treedist() {
    let trees = TREEDIST_TREES;
    let rfs = [
        vec![0, 4, 2, 10, 10, 10, 10, 10, 10, 10, 2, 10],
        vec![4, 0, 2, 10, 8, 10, 8, 10, 8, 10, 2, 10],
        vec![2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
        vec![10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
        vec![10, 8, 10, 2, 0, 4, 2, 4, 2, 2, 10, 4],
        vec![10, 10, 10, 2, 4, 0, 2, 2, 4, 2, 10, 2],
        vec![10, 8, 10, 4, 2, 2, 0, 4, 2, 4, 10, 4],
        vec![10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
        vec![10, 8, 10, 4, 2, 4, 2, 2, 0, 4, 10, 2],
        vec![10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
        vec![2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
        vec![10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
    ];

    trees.iter().combinations(2).for_each(|pair| {
        let (t0, t1) = (pair[0], pair[1]);
        let (i0, i1) = (
            trees.iter().position(|&t| t == *t0).unwrap(),
            trees.iter().position(|&t| t == *t1).unwrap(),
        );
        let s0 = Snapshot::from_newick(t0, false).unwrap();
        let s1 = Snapshot::from_newick(t1, false).unwrap();
        assert_eq!(
            rf_distance(&s0, &s1),
            rfs[i0][i1],
            "RF mismatch at [{i0}, {i1}]"
        );
    });
}

#[test]
// Weighted Robinson–Foulds distances according to
// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
fn weighted_robinson_foulds_treedist() {
    let trees = TREEDIST_TREES;
    let expected = [
        [
            0.,
            0.4,
            0.2,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.9999999999999999,
        ],
        [
            0.4,
            0.,
            0.2,
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.2,
            0.9999999999999999,
        ],
        [
            0.2,
            0.2,
            0.,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.9999999999999999,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.2,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.9999999999999999,
            0.2,
        ],
        [
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.2,
            0.,
            0.4,
            0.2,
            0.4,
            0.2,
            0.2,
            0.9999999999999999,
            0.4,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.4,
            0.,
            0.2,
            0.2,
            0.4,
            0.2,
            0.9999999999999999,
            0.2,
        ],
        [
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.4,
            0.2,
            0.2,
            0.,
            0.4,
            0.2,
            0.4,
            0.9999999999999999,
            0.4,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.2,
            0.2,
            0.9999999999999999,
            0.,
        ],
        [
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.4,
            0.2,
            0.4,
            0.2,
            0.2,
            0.,
            0.4,
            0.9999999999999999,
            0.2,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.2,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.9999999999999999,
            0.2,
        ],
        [
            0.2,
            0.2,
            0.,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.9999999999999999,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.2,
            0.2,
            0.9999999999999999,
            0.,
        ],
    ];

    trees.iter().combinations(2).for_each(|pair| {
        let (t0, t1) = (pair[0], pair[1]);
        let (i0, i1) = (
            trees.iter().position(|&t| t == *t0).unwrap(),
            trees.iter().position(|&t| t == *t1).unwrap(),
        );
        let s0 = Snapshot::from_newick(t0, false).unwrap();
        let s1 = Snapshot::from_newick(t1, false).unwrap();
        assert!(
            (wrf_distance(&s0, &s1) - expected[i0][i1]).abs() <= f64::EPSILON,
            "WRF mismatch at [{i0}, {i1}]"
        );
    });
}

#[test]
// Branch score distances according to
// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
fn kuhner_felsenstein_treedist() {
    let trees = TREEDIST_TREES;
    let expected = [
        [
            0.,
            0.2,
            0.14142135623730953,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.316227766016838,
        ],
        [
            0.2,
            0.,
            0.14142135623730953,
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.14142135623730953,
            0.316227766016838,
        ],
        [
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.316227766016838,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.14142135623730953,
            0.,
            0.2,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.14142135623730953,
            0.316227766016838,
            0.2,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.2,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.2,
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.2,
            0.14142135623730953,
            0.2,
            0.316227766016838,
            0.2,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.316227766016838,
            0.,
        ],
        [
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.2,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.2,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.316227766016838,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.316227766016838,
            0.,
        ],
    ];

    trees.iter().combinations(2).for_each(|pair| {
        let (t0, t1) = (pair[0], pair[1]);
        let (i0, i1) = (
            trees.iter().position(|&t| t == *t0).unwrap(),
            trees.iter().position(|&t| t == *t1).unwrap(),
        );
        let s0 = Snapshot::from_newick(t0, false).unwrap();
        let s1 = Snapshot::from_newick(t1, false).unwrap();
        assert!(
            (kf_distance(&s0, &s1) - expected[i0][i1]).abs() <= f64::EPSILON,
            "KF mismatch at [{i0}, {i1}]"
        );
    });
}

#[cfg(test)]
mod tests {
    use super::{TREEDIST_TREES, kf_distance, rf_distance, wrf_distance};
    use crate::snapshot::{Snapshot, Snapshots};

    const T0: &str = "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";
    const T1: &str = "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);";
    const T2: &str = "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";

    fn three_snapshots() -> Snapshots {
        Snapshots::from_newicks(&[T0, T1, T2], false).unwrap()
    }

    #[test]
    #[allow(clippy::erasing_op)] // Keep O*n for clarity of the flattened matrix indexing
    fn rf_known_values() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len(); // n = 3
        let mat = snaps.pairwise_rf(None);
        assert_eq!(mat[0 * n + 1], 4, "RF(T0,T1)");
        assert_eq!(mat[0 * n + 2], 2, "RF(T0,T2)");
        assert_eq!(mat[n + 2], 2, "RF(T1,T2)");
        for i in 0..n {
            assert_eq!(mat[i * n + i], 0, "diagonal [{i}]");
            for j in 0..n {
                assert_eq!(mat[i * n + j], mat[j * n + i], "RF symmetry [{i}][{j}]");
            }
        }
    }

    #[test]
    fn wrf_symmetric_zero_diagonal() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_wrf(None);

        for (i, row) in mat.chunks(n).enumerate() {
            assert_eq!(row[i], 0.0, "WRF diagonal [{i}]");
            for (j, v) in row.iter().enumerate() {
                assert!(
                    (v - mat[j * n + i]).abs() < f64::EPSILON,
                    "WRF symmetry [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn kf_symmetric_zero_diagonal() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_kf(None);
        for (i, row) in mat.chunks(n).enumerate() {
            assert_eq!(row[i], 0.0, "KF diagonal [{i}]");
            for (j, v) in row.iter().enumerate() {
                assert!(
                    (v - mat[j * n + i]).abs() < f64::EPSILON,
                    "KF symmetry [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn kf_differs_from_wrf() {
        let snaps = three_snapshots();
        let wrf = snaps.pairwise_wrf(None);
        let kf = snaps.pairwise_kf(None);
        let n = snaps.snapshots.len();
        let any_different = (0..n)
            .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
            .any(|(i, j)| (kf[i * n + j] - wrf[i * n + j]).abs() > f64::EPSILON);
        assert!(any_different, "KF and WRF should produce different values");
    }

    #[test]
    fn metrics_satisfy_triangle_inequality_rf() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_rf(None);

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let d_ij = mat[i * n + j] as f64;
                    let d_jk = mat[j * n + k] as f64;
                    let d_ik = mat[i * n + k] as f64;
                    // Add a small epsilon for floats (WRF/KF)
                    assert!(
                        d_ik <= d_ij + d_jk + 1e-10,
                        "Triangle inequality failed for indices {}, {}, {}",
                        i,
                        j,
                        k
                    );
                }
            }
        }
    }

    #[test]
    fn metrics_satisfy_triangle_inequality_wrf() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_wrf(None);

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let d_ij = mat[i * n + j];
                    let d_jk = mat[j * n + k];
                    let d_ik = mat[i * n + k];
                    // Add a small epsilon for floats (WRF/KF)
                    assert!(
                        d_ik <= d_ij + d_jk + 1e-10,
                        "Triangle inequality failed for indices {}, {}, {}",
                        i,
                        j,
                        k
                    );
                }
            }
        }
    }

    #[test]
    fn metrics_satisfy_triangle_inequality_kf() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_kf(None);

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let d_ij = mat[i * n + j];
                    let d_jk = mat[j * n + k];
                    let d_ik = mat[i * n + k];
                    // Add a small epsilon for floats (WRF/KF)
                    assert!(
                        d_ik <= d_ij + d_jk + 1e-10,
                        "Triangle inequality failed for indices {}, {}, {}",
                        i,
                        j,
                        k
                    );
                }
            }
        }
    }

    #[test]
    fn split_ids_sorted_per_snapshot() {
        let snaps = three_snapshots();
        for snap in &snaps.snapshots {
            for w in snap.split_ids.windows(2) {
                assert!(w[0] < w[1], "split_ids must be strictly ascending");
            }
        }
    }

    #[test]
    fn lengths_aligned_with_split_ids() {
        let snaps = three_snapshots();
        for snap in &snaps.snapshots {
            assert_eq!(snap.split_ids.len(), snap.lengths.len());
        }
    }

    fn treedist_snapshots() -> Snapshots {
        Snapshots::from_newicks(&TREEDIST_TREES, false).unwrap()
    }

    fn treedist_singles() -> Vec<Snapshot> {
        TREEDIST_TREES
            .iter()
            .map(|t| Snapshot::from_newick(t, false).unwrap())
            .collect()
    }

    // The bulk `pairwise_*` paths pack every tree into one matrix and sweep it,
    // while the public per-pair `*_distance` functions walk two trees directly.
    // The per-pair functions are pinned to PHYLIP ground truth by the
    // `*_treedist` tests above, so matching them pins the bulk paths too.

    #[test]
    fn pairwise_rf_matches_per_pair_ground_truth() {
        let snaps = treedist_snapshots();
        let singles = treedist_singles();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_rf(None);
        for i in 0..n {
            for j in 0..n {
                assert_eq!(
                    mat[i * n + j],
                    rf_distance(&singles[i], &singles[j]),
                    "bulk RF disagrees with per-pair RF at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn pairwise_wrf_matches_per_pair_ground_truth() {
        let snaps = treedist_snapshots();
        let singles = treedist_singles();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_wrf(None);
        for i in 0..n {
            for j in 0..n {
                let expected = wrf_distance(&singles[i], &singles[j]);
                assert!(
                    (mat[i * n + j] - expected).abs() < 1e-12,
                    "bulk WRF {} disagrees with per-pair WRF {expected} at [{i}][{j}]",
                    mat[i * n + j],
                );
            }
        }
    }

    #[test]
    fn pairwise_kf_matches_per_pair_ground_truth() {
        let snaps = treedist_snapshots();
        let singles = treedist_singles();
        let n = snaps.snapshots.len();
        let mat = snaps.pairwise_kf(None);
        for i in 0..n {
            for j in 0..n {
                let expected = kf_distance(&singles[i], &singles[j]);
                // KF's Gram reformulation (selfᵢ + selfⱼ − 2·dot) reorders the
                // summation, so it agrees within floating-point tolerance, not
                // bit-for-bit.
                assert!(
                    (mat[i * n + j] - expected).abs() < 1e-9,
                    "bulk KF {} disagrees with per-pair KF {expected} at [{i}][{j}]",
                    mat[i * n + j],
                );
            }
        }
    }

    #[test]
    fn pairwise_rf_identical_trees_all_zero() {
        // Every split is present in all 5 trees → every column is universal and
        // gets dropped, leaving zero-width bit rows. RF must still be 0 here.
        let t = TREEDIST_TREES[0];
        let snaps = Snapshots::from_newicks(&[t, t, t, t, t], false).unwrap();
        let mat = snaps.pairwise_rf(None);
        assert_eq!(mat.len(), 25);
        assert!(
            mat.iter().all(|&d| d == 0),
            "identical trees must have RF 0 everywhere"
        );
    }

    #[test]
    fn pairwise_kf_identical_trees_all_zero() {
        // Identical trees make `selfᵢ + selfⱼ − 2·dot` algebraically zero, but
        // rounding can nudge it slightly negative — the clamp must keep the
        // square root from producing NaN.
        let t = TREEDIST_TREES[0];
        let snaps = Snapshots::from_newicks(&[t, t, t], false).unwrap();
        let mat = snaps.pairwise_kf(None);
        assert!(
            mat.iter().all(|&d| d == 0.0),
            "identical trees must have KF 0 everywhere (clamp must avoid NaN)"
        );
    }
}
