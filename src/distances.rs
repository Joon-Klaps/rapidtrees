//! Pairwise tree distance backends for [`crate::snapshot::Snapshots`].
//!
//! All three metrics have the form `selfᵢ + selfⱼ − 2·shared`, where `shared` is
//! a popcount (RF), a running minimum (WRF), or a dot product (KF). That shape
//! lets each drop the splits that cannot affect `shared` before the O(n²) sweep:
//! RF drops splits held by every tree, WRF/KF drop splits held by only one. Both
//! filters are pure optimisations — disabling either changes no distance.
//!
//! Two kernels implement that shape, chosen per call by [`Backend::Auto`]:
//! dense rows swept word by word, or a two-pointer merge over the sorted
//! split-ID lists. Dense costs `⌈U/64⌉` words per pair regardless of how much
//! the trees share, where `U` is the collection's distinct-split count — cheap
//! on a posterior.
//!
//! There is no per-pair entry point. To compare two trees, build a two-tree
//! `Snapshots` and read the off-diagonal cell.

use crate::par::*;
use crate::snapshot::Snapshots;
use std::cmp::Ordering::{Equal, Greater, Less};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Which per-pair kernel to run.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum Backend {
    /// Choose from row density and the dense matrix's memory footprint.
    #[default]
    Auto,
    /// Dense presence/length rows, swept word by word.
    Dense,
    /// Two-pointer merge over the sorted split-ID lists.
    Sparse,
}

/// Refuse a dense matrix larger than this and merge instead.
const DENSE_BUDGET_BYTES: u64 = 5 << 30;

/// Bytes a dense `trees × width` matrix of 8-byte cells would take.
///
/// In `u64` because on wasm32 the `usize` product overflows well inside the
/// sizes this guard exists to refuse.
#[inline]
fn matrix_bytes(trees: usize, width: usize) -> u64 {
    trees as u64 * width as u64 * 8
}

/// The kernel a call actually ran. [`Backend::Auto`] resolves to one of these.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Kernel {
    /// Dense rows swept word by word.
    Dense,
    /// Two-pointer merge over the sorted split-ID lists.
    Sparse,
}

impl std::fmt::Display for Kernel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Kernel::Dense => "dense",
            Kernel::Sparse => "sparse",
        })
    }
}

/// A pairwise distance matrix and the kernel that produced it.
///
/// Row-major `n × n`, symmetric, zero diagonal. Both kernels return identical
/// matrices, so `kernel` is reporting only — nothing should branch on it.
#[derive(Clone, Debug)]
pub struct Distances<T> {
    /// The distances themselves.
    pub matrix: Vec<T>,
    /// Which kernel [`Backend::Auto`] resolved to for this call.
    pub kernel: Kernel,
}

/// Resolve `backend` to the kernel that will run.
///
/// `sweep` is the dense per-pair cost in merge-step equivalents and
/// `occupancy` is the kept tree × split incidences, so `occupancy / trees` is
/// the mean row length the merge would walk. One decision per call, never per
/// pair.
fn choose_kernel(
    backend: Backend,
    sweep: usize,
    bytes: u64,
    occupancy: u64,
    trees: usize,
) -> Kernel {
    let dense = match backend {
        Backend::Dense => true,
        Backend::Sparse => false,
        Backend::Auto => bytes <= DENSE_BUDGET_BYTES && sweep as u64 * trees as u64 <= occupancy,
    };
    if dense { Kernel::Dense } else { Kernel::Sparse }
}

/// Fold the values at every split ID both ascending lists hold.
///
/// This two-pointer walk is what makes the sparse backends independent of the
/// collection's distinct-split count.
#[inline]
fn merge_shared<T: Default + std::ops::AddAssign>(
    a: &[u32],
    b: &[u32],
    mut at: impl FnMut(usize, usize) -> T,
) -> T {
    let (mut i, mut j, mut acc) = (0, 0, T::default());
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Less => i += 1,
            Greater => j += 1,
            Equal => {
                acc += at(i, j);
                i += 1;
                j += 1;
            }
        }
    }
    acc
}

/// Fill a symmetric `n × n` matrix from `cell(i, j)`, one rayon task per row.
///
/// Only the upper triangle is computed; the diagonal stays at `T::default()`.
/// `progress` is bumped by each row's pair count as that row finishes.
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

    // Mirror into the lower triangle in tiles: a plain row-read/column-write
    // sweep puts every write on its own cache line at large `n`.
    const TILE: usize = 64;
    for i0 in (0..n).step_by(TILE) {
        for j0 in (i0..n).step_by(TILE) {
            for i in i0..(i0 + TILE).min(n) {
                for j in j0.max(i + 1)..(j0 + TILE).min(n) {
                    matrix[j * n + i] = matrix[i * n + j];
                }
            }
        }
    }

    matrix
}

/// Number of distinct splits across all trees.
///
/// Read from `split_ids` (dense from 0, sorted per tree) rather than the bitset
/// table, so it still works on paths that free the bitsets.
fn n_distinct_splits(snaps: &Snapshots) -> usize {
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

/// How many trees hold each split, indexed by split ID.
fn split_tree_counts(snaps: &Snapshots) -> Vec<u32> {
    let mut counts = vec![0u32; n_distinct_splits(snaps)];
    for snap in &snaps.snapshots {
        for &id in &snap.split_ids {
            counts[id as usize] += 1;
        }
    }
    counts
}

/// Give every split `keep` accepts a packed column index.
///
/// Returns `(column_of, n_columns, occupancy)`; `column_of[id]` is `u32::MAX`
/// if dropped and `occupancy` is the surviving tree × split incidences, which
/// is what [`choose_kernel`] weighs the dense sweep against. Packing the survivors
/// together is what shrinks the per-pair sweep.
fn assign_columns(counts: &[u32], keep: impl Fn(u32) -> bool) -> (Vec<u32>, usize, u64) {
    let mut column_of = vec![u32::MAX; counts.len()];
    let (mut n_columns, mut occupancy) = (0u32, 0u64);
    for (id, &count) in counts.iter().enumerate() {
        if keep(count) {
            column_of[id] = n_columns;
            n_columns += 1;
            occupancy += u64::from(count);
        }
    }
    (column_of, n_columns as usize, occupancy)
}

// ─── Robinson–Foulds ────────────────────────────────────────────────────────

/// `RF(i, j) = aᵢ + aⱼ − 2·popcount(rowᵢ & rowⱼ)` over presence bit-rows, or the
/// same count from a merge of the two sorted split-ID lists.
///
/// Splits held by *every* tree add equally to both `a` values and to the shared
/// count, so they cancel exactly and are dropped before packing. They cancel in
/// the merge too, which is why it can walk the raw ID lists with no filter.
pub(crate) fn distance_rf(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    backend: Backend,
) -> Distances<u32> {
    let n = snaps.snapshots.len();
    if n == 0 {
        return Distances {
            matrix: Vec::new(),
            kernel: Kernel::Dense,
        };
    }

    let counts = split_tree_counts(snaps);
    let n_splits = counts.len();
    let (bit_slot, kept, occupancy) = assign_columns(&counts, |count| count < n as u32);
    let everywhere = n_splits - kept;
    let words = kept.div_ceil(64);

    // A word carries 64 columns and the merge branches per split, so dense
    // stays ahead to ~16 words per kept split; past that the budget usually
    // decides anyway. Both constants are measured, not derived.
    let kernel = choose_kernel(
        backend,
        words.div_ceil(16),
        matrix_bytes(n, words),
        occupancy,
        n,
    );
    if kernel == Kernel::Sparse {
        let snapshots = &snaps.snapshots;
        let matrix = fill_symmetric(n, progress, |i, j| {
            let (a, b) = (&snapshots[i].split_ids, &snapshots[j].split_ids);
            (a.len() + b.len()) as u32 - 2 * merge_shared(a, b, |_, _| 1u32)
        });
        return Distances { matrix, kernel };
    }

    // One bitmask row per tree: a set bit means "this tree has that split".
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

    // Every tree holds all the everywhere-splits, so `a` is just its total minus
    // that fixed count.
    let kept_per_tree: Vec<u32> = snaps
        .snapshots
        .iter()
        .map(|snap| (snap.split_ids.len() - everywhere) as u32)
        .collect();

    let matrix = fill_symmetric(n, progress, |i, j| {
        let row_i = row_slice(&packed, i, words);
        let row_j = row_slice(&packed, j, words);
        let shared: u32 = row_i
            .iter()
            .zip(row_j)
            .map(|(&x, &y)| (x & y).count_ones())
            .sum();
        kept_per_tree[i] + kept_per_tree[j] - 2 * shared
    });
    Distances { matrix, kernel }
}

// ─── weighted metrics (WRF, KF) ─────────────────────────────────────────────

/// Lay out branch lengths over only the splits at least two trees hold.
///
/// A split held by one tree alone can never be shared, so it gets no column —
/// its `term(length)` folds into that tree's `unique_self` instead. On diverse
/// sets this cuts ~48 000 splits to ~2 600 (~19×); on similar sets it is nearly
/// a no-op. "Everywhere" splits are kept: unlike in RF, they do not cancel out
/// of a weighted score.
///
/// Returns the flat `n × stride` length matrix and, per tree, `Σ term(length)`
/// over the splits no other tree holds.
fn shared_length_rows(
    snaps: &Snapshots,
    column_of: &[u32],
    stride: usize,
    term: impl Fn(f64) -> f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = snaps.snapshots.len();
    let mut rows = vec![0.0f64; n * stride];
    if stride > 0 {
        rows.par_chunks_mut(stride)
            .zip(&snaps.snapshots)
            .for_each(|(row, snap)| {
                for (&id, &length) in snap.split_ids.iter().zip(&snap.lengths) {
                    let col = column_of[id as usize];
                    if col != u32::MAX {
                        row[col as usize] = length;
                    }
                }
            });
    }

    let unique_self = snaps
        .snapshots
        .iter()
        .map(|snap| {
            snap.split_ids
                .iter()
                .zip(&snap.lengths)
                .filter(|&(&id, _)| column_of[id as usize] == u32::MAX)
                .map(|(_, &l)| term(l))
                .sum()
        })
        .collect();

    (rows, unique_self)
}

/// `finish(selfᵢ + selfⱼ − 2·Σ overlap)` — the shape WRF and KF share.
///
/// `term` maps a length to its `self` contribution, `overlap` is the per-split
/// shared term, `finish` is applied last.
///
/// `self` is summed in the same order `overlap` walks — column order for the
/// dense backend, split-ID order for the merge — so identical trees cancel to
/// exactly 0.0 either way. The clamp stops rounding from handing `finish` a
/// negative.
fn weighted_distances(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    backend: Backend,
    term: impl Fn(f64) -> f64 + Sync,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Distances<f64> {
    let n = snaps.snapshots.len();
    if n == 0 {
        return Distances {
            matrix: Vec::new(),
            kernel: Kernel::Dense,
        };
    }

    let counts = split_tree_counts(snaps);
    let (column_of, stride, occupancy) = assign_columns(&counts, |count| count >= 2);

    // A dense f64 column costs roughly a quarter of a merge step.
    let kernel = choose_kernel(
        backend,
        stride.div_ceil(4),
        matrix_bytes(n, stride),
        occupancy,
        n,
    );
    if kernel == Kernel::Sparse {
        let snapshots = &snaps.snapshots;
        let self_total: Vec<f64> = snapshots
            .iter()
            .map(|snap| snap.lengths.iter().map(|&l| term(l)).sum())
            .collect();

        let matrix = fill_symmetric(n, progress, |i, j| {
            let (a, b) = (&snapshots[i], &snapshots[j]);
            let shared: f64 = merge_shared(&a.split_ids, &b.split_ids, |x, y| {
                overlap(a.lengths[x], b.lengths[y])
            });
            finish((self_total[i] + self_total[j] - 2.0 * shared).max(0.0))
        });
        return Distances { matrix, kernel };
    }

    let (rows, unique_self) = shared_length_rows(snaps, &column_of, stride, &term);
    let self_total: Vec<f64> = (0..n)
        .map(|i| {
            row_slice(&rows, i, stride)
                .iter()
                .map(|&l| term(l))
                .sum::<f64>()
                + unique_self[i]
        })
        .collect();

    let matrix = fill_symmetric(n, progress, |i, j| {
        let row_i = row_slice(&rows, i, stride);
        let row_j = row_slice(&rows, j, stride);
        let shared_term: f64 = row_i.iter().zip(row_j).map(|(&a, &b)| overlap(a, b)).sum();
        finish((self_total[i] + self_total[j] - 2.0 * shared_term).max(0.0))
    });
    Distances { matrix, kernel }
}

/// `WRF(i, j) = Σ lenᵢ + Σ lenⱼ − 2·Σ min(lenᵢ, lenⱼ)`.
///
/// The `min` form follows from `|a − b| = a + b − 2·min(a, b)`. Assumes
/// non-negative branch lengths; missing lengths parse as 0.0.
pub(crate) fn distance_wrf(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    backend: Backend,
) -> Distances<f64> {
    weighted_distances(snaps, progress, backend, |l| l, f64::min, |d| d)
}

/// `KF(i, j) = sqrt(Σ lenᵢ² + Σ lenⱼ² − 2·Σ lenᵢ·lenⱼ)` — Euclidean distance in
/// branch-length space.
pub(crate) fn distance_kf(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    backend: Backend,
) -> Distances<f64> {
    weighted_distances(snaps, progress, backend, |l| l * l, |a, b| a * b, f64::sqrt)
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

    let snaps = Snapshots::from_newicks(&trees, false).unwrap();
    let n = snaps.len();
    let mat = snaps.pairwise_rf(None);

    for (i, want_row) in rfs.iter().enumerate() {
        for (j, &want) in want_row.iter().enumerate() {
            assert_eq!(mat[i * n + j], want, "RF mismatch at [{i}, {j}]");
        }
    }
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

    let snaps = Snapshots::from_newicks(&trees, false).unwrap();
    let n = snaps.len();
    let mat = snaps.pairwise_wrf(None);

    // The shared-term rewrite reassociates the sum, so this lands a few ulp off
    // the published value. Distances are ~1.0, so 1e-12 is still far tighter
    // than any real disagreement.
    for (i, want_row) in expected.iter().enumerate() {
        for (j, &want) in want_row.iter().enumerate() {
            assert!(
                (mat[i * n + j] - want).abs() <= 1e-12,
                "WRF mismatch at [{i}, {j}]: got {}, want {want}",
                mat[i * n + j]
            );
        }
    }
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

    let snaps = Snapshots::from_newicks(&trees, false).unwrap();
    let n = snaps.len();
    let mat = snaps.pairwise_kf(None);

    // Tolerance as in `weighted_robinson_foulds_treedist`.
    for (i, want_row) in expected.iter().enumerate() {
        for (j, &want) in want_row.iter().enumerate() {
            assert!(
                (mat[i * n + j] - want).abs() <= 1e-12,
                "KF mismatch at [{i}, {j}]: got {}, want {want}",
                mat[i * n + j]
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Backend, Kernel, TREEDIST_TREES};
    use crate::snapshot::{InternSnap, Snapshots};
    use std::collections::{BTreeMap, BTreeSet};

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

    /// Test oracle: all three metrics for one pair, straight from the
    /// definitions over the union of both trees' splits (absent = length 0).
    ///
    /// Deliberately naive — no filtering, no `selfᵢ + selfⱼ − 2·shared` rewrite —
    /// so it cannot reproduce a backend bug by construction.
    fn reference_distances(a: &InternSnap, b: &InternSnap) -> (usize, f64, f64) {
        let lengths_by_split = |s: &InternSnap| -> BTreeMap<u32, f64> {
            s.split_ids
                .iter()
                .copied()
                .zip(s.lengths.iter().copied())
                .collect()
        };
        let (ma, mb) = (lengths_by_split(a), lengths_by_split(b));

        let (mut rf, mut wrf, mut sum_sq) = (0usize, 0.0, 0.0);
        for id in ma
            .keys()
            .chain(mb.keys())
            .copied()
            .collect::<BTreeSet<u32>>()
        {
            let x = ma.get(&id).copied().unwrap_or(0.0);
            let y = mb.get(&id).copied().unwrap_or(0.0);
            if ma.contains_key(&id) != mb.contains_key(&id) {
                rf += 1;
            }
            wrf += (x - y).abs();
            sum_sq += (x - y) * (x - y);
        }
        (rf, wrf, sum_sq.sqrt())
    }

    /// Assert all three metrics agree with [`reference_distances`] on `snaps`,
    /// under both the dense and the sparse kernel.
    fn assert_matches_reference(snaps: &Snapshots, ctx: &str) {
        for backend in [Backend::Dense, Backend::Sparse] {
            let ctx = format!("{ctx}, {backend:?}");
            assert_cells(
                snaps,
                &snaps.pairwise_rf_with(None, backend).matrix,
                &snaps.pairwise_wrf_with(None, backend).matrix,
                &snaps.pairwise_kf_with(None, backend).matrix,
                &ctx,
            );
        }
    }

    fn assert_cells(snaps: &Snapshots, rf: &[u32], wrf: &[f64], kf: &[f64], ctx: &str) {
        let n = snaps.snapshots.len();
        for i in 0..n {
            for j in 0..n {
                let (want_rf, want_wrf, want_kf) =
                    reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);

                assert_eq!(rf[i * n + j] as usize, want_rf, "RF at [{i}][{j}], {ctx}");
                assert!(
                    (wrf[i * n + j] - want_wrf).abs() <= 1e-9 * want_wrf.max(1.0),
                    "WRF {} vs reference {want_wrf} at [{i}][{j}], {ctx}",
                    wrf[i * n + j],
                );
                assert!(
                    (kf[i * n + j] - want_kf).abs() <= 1e-9 * want_kf.max(1.0),
                    "KF {} vs reference {want_kf} at [{i}][{j}], {ctx}",
                    kf[i * n + j],
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
    fn pairwise_wrf_identical_trees_all_zero() {
        // Exactly 0.0 only holds because self and shared are summed in the same
        // column order; reordering either would leave a rounding crumb.
        let t = TREEDIST_TREES[0];
        let snaps = Snapshots::from_newicks(&[t, t, t], false).unwrap();
        let mat = snaps.pairwise_wrf(None);
        assert!(
            mat.iter().all(|&d| d == 0.0),
            "identical trees must have WRF exactly 0 everywhere"
        );
    }

    /// Trees sharing almost nothing, so most splits land in the "held by one
    /// tree only" bucket that [`shared_length_rows`] drops — which the treedist
    /// fixtures barely exercise. Branch lengths are all distinct so a wrongly
    /// dropped column changes the answer visibly.
    const DIVERSE_TREES: [&str; 4] = [
        "(((A:0.11,B:0.12):0.13,(C:0.14,D:0.15):0.16):0.17,((E:0.18,F:0.19):0.21,(G:0.22,H:0.23):0.24):0.25);",
        "(((A:0.31,E:0.32):0.33,(C:0.34,G:0.35):0.36):0.37,((B:0.38,F:0.39):0.41,(D:0.42,H:0.43):0.44):0.45);",
        "(((A:0.51,H:0.52):0.53,(B:0.54,G:0.55):0.56):0.57,((C:0.58,F:0.59):0.61,(D:0.62,E:0.63):0.64):0.65);",
        "(((A:0.71,D:0.72):0.73,(F:0.74,G:0.75):0.76):0.77,((B:0.78,H:0.79):0.81,(C:0.82,E:0.83):0.84):0.85);",
    ];

    #[test]
    fn pairwise_weighted_metrics_match_reference_on_diverse_trees() {
        let snaps = Snapshots::from_newicks(&DIVERSE_TREES, false).unwrap();
        assert_matches_reference(&snaps, "diverse fixtures");
    }

    /// Deterministic LCG, so a failure here is always reproducible.
    fn lcg(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 33
    }

    /// Random binary topology over `n_taxa` leaves, built by repeatedly joining
    /// two randomly chosen subtrees.
    fn random_newick(n_taxa: usize, state: &mut u64) -> String {
        let mut parts: Vec<String> = (0..n_taxa).map(|i| format!("T{i}")).collect();
        while parts.len() > 1 {
            let i = (lcg(state) as usize) % parts.len();
            let a = parts.swap_remove(i);
            let j = (lcg(state) as usize) % parts.len();
            let b = parts.swap_remove(j);
            let la = (lcg(state) % 100 + 1) as f64 / 100.0;
            let lb = (lcg(state) % 100 + 1) as f64 / 100.0;
            parts.push(format!("({a}:{la},{b}:{lb})"));
        }
        format!("{};", parts.pop().unwrap())
    }

    /// Differential test against [`reference_distances`] over a spread of taxon
    /// counts and sharing levels.
    ///
    /// `duplicates` matters because the column filter keys off how many trees
    /// hold a split: distinct trees push splits into the dropped bucket,
    /// duplicates pull them back into shared columns.
    #[test]
    fn pairwise_backends_match_reference_on_random_trees() {
        for &(n_taxa, n_trees, duplicates, seed) in &[
            (8usize, 6usize, 0usize, 1u64),
            (12, 8, 0, 2),
            (20, 10, 0, 3),
            (31, 9, 0, 4),
            (16, 6, 3, 5),   // forces splits into the shared-column bucket
            (24, 5, 5, 6),   // majority duplicates
            (10, 12, 11, 7), // all but one identical
        ] {
            let mut state = seed;
            let mut newicks: Vec<String> = (0..n_trees)
                .map(|_| random_newick(n_taxa, &mut state))
                .collect();
            for k in 0..duplicates {
                newicks.push(newicks[k % n_trees].clone());
            }

            let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
            let snaps = Snapshots::from_newicks(&refs, false).unwrap();
            assert_matches_reference(&snaps, &format!("seed={seed} taxa={n_taxa}"));

            // Duplicated trees must land on exactly zero, not a rounding crumb.
            let n = snaps.snapshots.len();
            let (wrf, kf) = (snaps.pairwise_wrf(None), snaps.pairwise_kf(None));
            for i in 0..n {
                for j in 0..n {
                    if refs[i] == refs[j] {
                        assert_eq!(wrf[i * n + j], 0.0, "WRF identical, seed={seed} [{i}][{j}]");
                        assert_eq!(kf[i * n + j], 0.0, "KF identical, seed={seed} [{i}][{j}]");
                    }
                }
            }
        }
    }

    /// The two kernels must agree cell for cell, not merely to a tolerance —
    /// the selector switches between them silently, so a divergence would show
    /// up as a distance that depends on the dataset's diversity.
    ///
    /// Covers the PHYLIP reference suite, a near-identical set (every column
    /// universal, so the dense side packs nothing), a fully diverse set (every
    /// column unique) and random sets at two taxon counts.
    #[test]
    fn dense_and_sparse_agree_cell_for_cell() {
        let treedist: Vec<&str> = TREEDIST_TREES.to_vec();
        let identical = vec![TREEDIST_TREES[0]; 6];
        let diverse = DIVERSE_TREES.to_vec();

        let mut random: Vec<Vec<String>> = Vec::new();
        for &(n_taxa, n_trees, seed) in &[(9usize, 7usize, 11u64), (40, 6, 12)] {
            let mut state = seed;
            random.push(
                (0..n_trees)
                    .map(|_| random_newick(n_taxa, &mut state))
                    .collect(),
            );
        }
        let random: Vec<Vec<&str>> = random
            .iter()
            .map(|set| set.iter().map(String::as_str).collect())
            .collect();

        let cases = [&treedist, &identical, &diverse]
            .into_iter()
            .chain(random.iter());

        for (case, newicks) in cases.enumerate() {
            let snaps = Snapshots::from_newicks(newicks, false).unwrap();
            assert_eq!(
                snaps.pairwise_rf_with(None, Backend::Dense).matrix,
                snaps.pairwise_rf_with(None, Backend::Sparse).matrix,
                "RF backends disagree, case {case}"
            );
            for (dense, sparse) in [
                (
                    snaps.pairwise_wrf_with(None, Backend::Dense).matrix,
                    snaps.pairwise_wrf_with(None, Backend::Sparse).matrix,
                ),
                (
                    snaps.pairwise_kf_with(None, Backend::Dense).matrix,
                    snaps.pairwise_kf_with(None, Backend::Sparse).matrix,
                ),
            ] {
                for (k, (&d, &sp)) in dense.iter().zip(&sparse).enumerate() {
                    assert!(
                        (d - sp).abs() <= 1e-12 * d.abs().max(1.0),
                        "weighted backends disagree at {k}: {d} vs {sp}, case {case}"
                    );
                }
            }
        }
    }

    /// The memory guard is what stops a diverse collection from allocating a
    /// presence matrix larger than the machine: at 5 000 trees over 500 taxa
    /// with every split unique the dense side wants 1.55 GB.
    #[test]
    fn auto_refuses_a_dense_matrix_over_budget() {
        let (trees, words) = (5_000usize, 38_828usize);
        let occupancy = 497 * trees as u64; // kept splits per tree × trees
        assert_eq!(
            super::choose_kernel(
                Backend::Auto,
                0, // free by the time-based rule; only the budget may refuse it
                super::matrix_bytes(trees, words),
                occupancy,
                trees,
            ),
            Kernel::Sparse
        );
        assert_eq!(
            super::choose_kernel(Backend::Auto, 0, 1 << 20, occupancy, trees),
            Kernel::Dense
        );
    }

    /// Every call reports the kernel that produced its matrix.
    #[test]
    fn auto_reports_the_backend_it_picked() {
        let snaps = Snapshots::from_newicks(&DIVERSE_TREES, false).unwrap();
        for backend in [Backend::Sparse, Backend::Dense] {
            let expected = match backend {
                Backend::Sparse => Kernel::Sparse,
                _ => Kernel::Dense,
            };
            assert_eq!(snaps.pairwise_rf_with(None, backend).kernel, expected);
        }
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
