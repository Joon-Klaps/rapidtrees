//! Tree distance functions.
//!
//! # Public API
//! - [`rf_distance`], [`wrf_distance`], [`kf_distance`] — compare two standalone
//!   [`crate::snapshot::Snapshot`]s (use bitset sorted-merge; panics if leaf sets differ).
//!
//! # Bulk pairwise computation
//! Call [`crate::snapshot::Snapshots::pairwise_rf`] etc. — these use the fast
//! interned-integer path internally via the `pub(crate)` helpers in this module.

use crate::snapshot::{InternSnap, Snapshot, Snapshots};
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

// ─── pub(crate): fast interned-integer paths used by Snapshots::pairwise_* ──

/// RF distance on two interned snapshots (integer sorted-merge, no leaf check).
#[inline]
pub(crate) fn rf_distance_fast(a: &InternSnap, b: &InternSnap) -> usize {
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

/// WRF distance on two interned snapshots.
#[inline]
pub(crate) fn wrf_distance_fast(a: &InternSnap, b: &InternSnap) -> f64 {
    let mut dist = 0.0;
    let mut i = 0;
    let mut j = 0;
    while i < a.split_ids.len() && j < b.split_ids.len() {
        match a.split_ids[i].cmp(&b.split_ids[j]) {
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

/// KF distance on two interned snapshots.
#[inline]
pub(crate) fn kf_distance_fast(a: &InternSnap, b: &InternSnap) -> f64 {
    let mut sum_sq = 0.0;
    let mut i = 0;
    let mut j = 0;
    while i < a.split_ids.len() && j < b.split_ids.len() {
        match a.split_ids[i].cmp(&b.split_ids[j]) {
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

/// Generic symmetric pairwise matrix computation over all trees in a [`Snapshots`].
///
/// Computes `metric(a, b)` for every upper-triangle pair in parallel via rayon,
/// then mirrors into the lower triangle. The diagonal stays at `T::default()`.
pub(crate) fn pairwise_symmetric<T, F>(snaps: &Snapshots, metric: F) -> Vec<T>
where
    T: Copy + Default + Send,
    F: Fn(&InternSnap, &InternSnap) -> T + Sync,
{
    pairwise_symmetric_counted(snaps, metric, None)
}

/// Same as [`pairwise_symmetric`], but bumps `progress` by `n - i - 1` after every
/// row `i` completes. `progress = None` is equivalent to [`pairwise_symmetric`].
///
/// The increment is the number of upper-triangle pairs produced by that row, so
/// when all rows have finished `progress` equals `n*(n-1)/2`.
pub(crate) fn pairwise_symmetric_counted<T, F>(
    snaps: &Snapshots,
    metric: F,
    progress: Option<&AtomicUsize>,
) -> Vec<T>
where
    T: Copy + Default + Send,
    F: Fn(&InternSnap, &InternSnap) -> T + Sync,
{
    let n = snaps.snapshots.len();

    // 1. Single contiguous allocation
    let mut matrix = vec![T::default(); n * n];

    // 2. Safely divide the flat array into mutable chunks (rows)
    matrix.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        for (j, dist) in row.iter_mut().enumerate().take(n).skip(i + 1) {
            *dist = metric(&snaps.snapshots[i], &snaps.snapshots[j]);
        }
        if let Some(counter) = progress {
            // Pairs produced by this row: indices (i, i+1), ..., (i, n-1).
            counter.fetch_add(n.saturating_sub(i + 1), Ordering::Relaxed);
        }
    });

    // 3. Sequential mirroring
    for i in 0..n {
        for j in (i + 1)..n {
            matrix[j * n + i] = matrix[i * n + j];
        }
    }

    matrix
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
    use crate::snapshot::Snapshots;

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
        let mat = snaps.pairwise_rf();
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
        let mat = snaps.pairwise_wrf();

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
        let mat = snaps.pairwise_kf();
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
        let wrf = snaps.pairwise_wrf();
        let kf = snaps.pairwise_kf();
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
        let mat = snaps.pairwise_rf();

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
        let mat = snaps.pairwise_wrf();

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
        let mat = snaps.pairwise_kf();

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
}
