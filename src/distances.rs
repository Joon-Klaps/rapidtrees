//! Pairwise tree distance backends for [`crate::snapshot::Snapshots`].
//!
//! All three metrics have the form `selfᵢ + selfⱼ − 2·shared`, where `shared` is
//! a popcount (RF), a running minimum (WRF), or a dot product (KF). That shape
//! lets each drop the splits that cannot affect `shared` before the O(n²) sweep.
//! A split held by only one tree is never shared, so every metric drops it; RF
//! also drops the splits held by every tree, which cancel. Both filters are pure
//! optimisations — disabling either changes no distance.
//!
//! Every metric splits its columns by how many trees hold them. A widely held
//! split keeps a dense column: a bit in RF's packed rows, an `f64` in the
//! weighted metrics' rows, swept word by word for every pair. A rarely held one
//! keeps a posting list of the trees that hold it instead, and each row adds
//! those shared terms straight into its output. A posting list of `k` trees
//! costs `k²/2` additions — HashRF's bucket — which is why only the rare splits
//! get one. A dense column costs every pair the same whether or not either
//! tree holds the split, so on a posterior, where most distinct shared splits
//! are held by a handful of trees, the lists take most of the width out of the
//! sweep. RF's bit columns are 64 to a word, so a split has to be much rarer to
//! be worth a list there ([`RF_DENSE_SHARE`]) than in the weighted metrics
//! ([`WEIGHTED_DENSE_SHARE`]).
//!
//! Seven trees, columns in descending holder count as `assign_columns` orders
//! them (`■` = the tree holds the split). The cutoff is `min_dense = 4` here
//! for illustration:
//!
//! ```text
//!             dense   │ posting │ dropped
//!             D1  D2  │ S5  S6  │ U1  U2
//! held by      7   5  │  3   2  │  1   1
//! ────────────────────┼─────────┼────────
//!     T0       ■   ■  │  ·   ·  │  ·   ·
//!     T1       ■   ■  │  ■   ·  │  ·   ·
//!     T2       ■   ■  │  ·   ·  │  ■   ·
//!     T3       ■   ·  │  ·   ■  │  ·   ·
//!     T4       ■   ■  │  ■   ·  │  ·   ·
//!     T5       ■   ·  │  ·   ·  │  ·   ■
//!     T6       ■   ■  │  ■   ■  │  ·   ·
//!                     ▲         ▲
//!                     │         └ held by ≥ 2: can be shared
//!                     └ held by ≥ min_dense
//! ```
//!
//! Dense columns are swept for every pair. `S5` becomes the posting list
//! `[T1, T4, T6]` and costs three pair updates; `S6` becomes `[T3, T6]` and
//! costs one. `U1` and `U2` fold into their tree's `self` term. RF has no
//! posting lists: it keeps `D2`, `S5`, and `S6` as bit columns, and also drops
//! `D1`, since a split held by every tree cancels.
//!
//! There is no per-pair entry point. To compare two trees, build a two-tree
//! `Snapshots` and read the off-diagonal cell.

use crate::par::*;
use crate::snapshot::{InternSnap, Snapshots};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Fill a symmetric `n × n` matrix one row at a time, one rayon task per row.
///
/// `fill_row(i, row)` writes `row[j]` for every `j > i`. `row` arrives holding
/// `T::default()` throughout, so a caller may accumulate into it first. The
/// diagonal stays at `T::default()` and the lower triangle is mirrored from
/// the upper. `progress` is bumped by each row's pair count as that row
/// finishes.
fn fill_symmetric<T, F>(n: usize, progress: Option<&AtomicUsize>, fill_row: F) -> Vec<T>
where
    T: Copy + Default + Send,
    F: Fn(usize, &mut [T]) + Sync,
{
    let mut matrix = vec![T::default(); n * n];

    matrix.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        fill_row(i, row);
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

/// Borrow row `i` of a flat, row-major matrix whose rows are `stride` wide.
#[inline]
fn row_slice<T>(flat: &[T], i: usize, stride: usize) -> &[T] {
    &flat[i * stride..][..stride]
}

/// Replace each slot with the sum of the slots before it and return the total.
fn exclusive_prefix_sum<'a, T>(slots: impl IntoIterator<Item = &'a mut T>) -> T
where
    T: 'a + Copy + Default + std::ops::Add<Output = T>,
{
    let mut total = T::default();
    for slot in slots {
        (*slot, total) = (total, total + *slot);
    }
    total
}

/// Marks a split that [`assign_columns`] gave no column.
const NO_COLUMN: u32 = u32::MAX;

/// The column of split `id`, or `None` if `column_of` gave it none.
#[inline]
fn column(column_of: &[u32], id: u32) -> Option<usize> {
    let col = column_of[id as usize];
    (col != NO_COLUMN).then_some(col as usize)
}

/// Give every split `keep` accepts a packed column index.
///
/// `counts[id]` is how many of the `n_trees` trees hold split `id`. Returns
/// `(column_of, n_columns)`, with [`NO_COLUMN`] for the splits `keep` rejects.
/// Columns run in descending tree count, which clusters the widely-held splits
/// into the low words. A count never exceeds `n_trees`, so this is a counting
/// sort: two passes over the splits and no comparisons. Ties keep ID order.
fn assign_columns(counts: &[u32], n_trees: usize, keep: impl Fn(u32) -> bool) -> (Vec<u32>, usize) {
    // `next[c]` becomes the first column of the splits held by `c` trees.
    let mut next = vec![0u32; n_trees + 1];
    for &count in counts.iter().filter(|&&count| keep(count)) {
        next[count as usize] += 1;
    }
    let kept = exclusive_prefix_sum(next.iter_mut().rev());

    let column_of = counts
        .iter()
        .map(|&count| {
            if keep(count) {
                let col = next[count as usize];
                next[count as usize] += 1;
                col
            } else {
                NO_COLUMN
            }
        })
        .collect();
    (column_of, kept as usize)
}

// ─── posting lists ──────────────────────────────────────────────────────────

/// The trees holding each rarely held split: HashRF's bucket.
///
/// Posting column `c` owns entries `offsets[c]..offsets[c + 1]` of `trees`,
/// ascending because the trees are appended in order, and of `lengths` when
/// they were asked for. A row walks its own posted splits and adds one term
/// per later tree holding each, so a split held by `k` trees costs `k²/2`
/// additions in all, however many trees the collection has.
struct Postings {
    offsets: Vec<usize>,
    trees: Vec<u32>,
    lengths: Vec<f64>,
}

impl Postings {
    /// One list per split that `column_of` gives a column, carrying branch
    /// lengths when `with_lengths`.
    fn build(snaps: &Snapshots, column_of: &[u32], n_columns: usize, with_lengths: bool) -> Self {
        let mut offsets = vec![0usize; n_columns + 1];
        for (&count, &col) in snaps.split_counts.iter().zip(column_of) {
            if col != u32::MAX {
                offsets[col as usize] = count as usize;
            }
        }
        let mut total = 0usize;
        for slot in offsets.iter_mut() {
            (*slot, total) = (total, total + *slot);
        }

        let mut next = offsets.clone();
        let mut trees = vec![0u32; total];
        let mut lengths = if with_lengths {
            vec![0.0f64; total]
        } else {
            Vec::new()
        };
        for (tree, snap) in snaps.snapshots.iter().enumerate() {
            for (k, &id) in snap.split_ids.iter().enumerate() {
                let col = column_of[id as usize];
                if col != u32::MAX {
                    let at = &mut next[col as usize];
                    trees[*at] = tree as u32;
                    if with_lengths {
                        lengths[*at] = snap.lengths[k];
                    }
                    *at += 1;
                }
            }
        }
        Self {
            offsets,
            trees,
            lengths,
        }
    }

    /// The trees after `i` that hold posting column `col`, with their branch
    /// lengths (empty when built without them). Only those are wanted: the
    /// lower triangle is mirrored.
    #[inline]
    fn after(&self, col: u32, i: usize) -> (&[u32], &[f64]) {
        let span = self.offsets[col as usize]..self.offsets[col as usize + 1];
        let trees = &self.trees[span.clone()];
        let from = trees.partition_point(|&tree| tree as usize <= i);
        let lengths = if self.lengths.is_empty() {
            &[][..]
        } else {
            &self.lengths[span][from..]
        };
        (&trees[from..], lengths)
    }
}

// ─── Robinson–Foulds ────────────────────────────────────────────────────────

/// Share of the trees a split must be held by to keep its bit column in RF.
/// Rarer splits get a posting list instead.
///
/// A bit column costs every pair a sixty-fourth of a word, so a split has to be
/// far rarer than in the weighted metrics before its list wins: 3 per cent was
/// the best all-round of 1, 3, 7 and 15 per cent on simulated posteriors from
/// 50 to 20 000 taxa. The boundary changes the speed and nothing else.
const RF_DENSE_SHARE: f64 = 0.03;

/// `RF(i, j) = aᵢ + aⱼ − 2·shared(i, j)`, where `shared` counts the splits
/// both trees hold.
///
/// Splits held by *every* tree add equally to both `a` values and to the shared
/// count, so they cancel exactly and are dropped. Splits held by one tree count
/// towards that tree's `a` but can never be shared, so they get nothing either.
/// On a posterior that is most of the distinct splits. The rest are shared by
/// a popcount over packed bit-rows for the widely held ones and by posting
/// lists for the rare ones.
pub(crate) fn distance_rf(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<u32> {
    let n = snaps.snapshots.len();
    let min_dense = ((n as f64 * RF_DENSE_SHARE).ceil() as u32).max(2);
    distance_rf_split(snaps, progress, min_dense)
}

/// [`distance_rf`] with the bit/posting boundary given: a split held by at
/// least `min_dense` trees (but not all) gets a bit column, one held by fewer
/// (but at least two) gets a posting list.
fn distance_rf_split(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    min_dense: u32,
) -> Vec<u32> {
    let n_trees = snaps.snapshots.len();
    if n_trees == 0 {
        return Vec::new();
    }
    let everyone = n_trees as u32;

    let (bit_slot, kept) = assign_columns(&snaps.split_counts, n_trees, |count| {
        count >= min_dense && count < everyone
    });
    let (posting_col, n_posted) = assign_columns(&snaps.split_counts, n_trees, |count| {
        count >= 2 && count < min_dense.min(everyone)
    });
    let postings = Postings::build(snaps, &posting_col, n_posted, false);
    let everywhere = snaps
        .split_counts
        .iter()
        .filter(|&&count| count == n_trees as u32)
        .count();
    let words = kept.div_ceil(64);

    // One bitmask row per tree: a set bit means "this tree has that split".
    // `spans[i]` is the first and last non-zero word of row `i`, empty as
    // `(1, 0)`. Tracked while the row is built, so it costs nothing extra.
    let mut packed = vec![0u64; n_trees * words];
    let mut spans = vec![(1usize, 0usize); n_trees];
    if words > 0 {
        packed
            .par_chunks_mut(words)
            .zip(spans.par_chunks_mut(1))
            .zip(&snaps.snapshots)
            .for_each(|((row, span), snap)| {
                let (mut lo, mut hi) = (usize::MAX, 0usize);
                for slot in snap
                    .split_ids
                    .iter()
                    .filter_map(|&id| column(&bit_slot, id))
                {
                    let word = slot >> 6;
                    row[word] |= 1u64 << (slot & 63);
                    (lo, hi) = (lo.min(word), hi.max(word));
                }
                if lo != usize::MAX {
                    span[0] = (lo, hi);
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

    fill_symmetric(n_trees, progress, |i, row: &mut [u32]| {
        // Shared splits from the posting lists first, counted in place.
        for &id in &snaps.snapshots[i].split_ids {
            let col = posting_col[id as usize];
            if col != u32::MAX {
                for &j in postings.after(col, i).0 {
                    row[j as usize] += 1;
                }
            }
        }

        for (j, slot) in row.iter_mut().enumerate().skip(i + 1) {
            let (lo, hi) = (spans[i].0.max(spans[j].0), spans[i].1.min(spans[j].1));
            let shared: u32 = if lo > hi {
                0
            } else {
                row_slice(&packed, i, words)[lo..=hi]
                    .iter()
                    .zip(&row_slice(&packed, j, words)[lo..=hi])
                    .map(|(&x, &y)| (x & y).count_ones())
                    .sum()
            };
            *slot = kept_per_tree[i] + kept_per_tree[j] - 2 * (shared + *slot);
        }
    })
}

// ─── weighted metrics (WRF, KF) ─────────────────────────────────────────────

/// Share of the trees a split must be held by to keep a dense column in the
/// weighted metrics. Rarer splits get a posting list instead.
///
/// A dense column costs every pair one element; a posting list of `k` trees
/// costs `k²/2` scattered additions, so the list wins below some share of the
/// trees and loses above it. A quarter was fastest of 5, 25 and 101 per cent on
/// simulated posteriors from 50 to 20 000 taxa. The boundary changes the speed
/// and nothing else: any value gives the same distances.
const WEIGHTED_DENSE_SHARE: f64 = 0.25;

/// `Σ overlap(aₖ, bₖ)` over two rows of equal length, in a fixed order.
///
/// Eight running sums rather than one. A single `f64` sum is a chain the
/// compiler may not reorder, so it can neither vectorise nor pipeline; eight
/// independent ones can do both. The tail shorter than eight lands in the first
/// lanes. The order is set here rather than by the compiler, so `sweep(a, a)`
/// and `sweep(a, b)` add the same terms the same way whenever `a == b`, which
/// is what keeps identical trees at exactly 0.0.
#[inline]
fn sweep(a: &[f64], b: &[f64], overlap: &impl Fn(f64, f64) -> f64) -> f64 {
    const LANES: usize = 8;
    let (a_blocks, a_tail) = a.as_chunks::<LANES>();
    let (b_blocks, b_tail) = b.as_chunks::<LANES>();

    let mut lanes = [0.0f64; LANES];
    let mut add = |xs: &[f64], ys: &[f64]| {
        for ((lane, &x), &y) in lanes.iter_mut().zip(xs).zip(ys) {
            *lane += overlap(x, y);
        }
    };
    for (xs, ys) in a_blocks.iter().zip(b_blocks) {
        add(xs, ys);
    }
    add(a_tail, b_tail);
    lanes.iter().sum()
}

/// A tree's splits that have a column in `column_of`, as `(column, length)`,
/// in the tree's own order.
fn lengths_by_column<'a>(
    snap: &'a InternSnap,
    column_of: &'a [u32],
) -> impl Iterator<Item = (usize, f64)> + 'a {
    snap.split_ids
        .iter()
        .zip(&snap.lengths)
        .filter_map(|(&id, &length)| column(column_of, id).map(|col| (col, length)))
}

/// Every tree's lengths over the dense columns, as one flat row-major
/// `n × stride` matrix with 0.0 where a tree lacks the split.
fn dense_rows(snaps: &Snapshots, dense_col: &[u32], stride: usize) -> Vec<f64> {
    let mut dense = vec![0.0f64; snaps.snapshots.len() * stride];
    if stride > 0 {
        dense
            .par_chunks_mut(stride)
            .zip(&snaps.snapshots)
            .for_each(|(row, snap)| {
                for (col, length) in lengths_by_column(snap, dense_col) {
                    row[col] = length;
                }
            });
    }
    dense
}

/// One tree's posting-list splits as `(column, length)`, in column order.
///
/// Both the tree's own `self` term and every shared term it enters are summed
/// in this order, so two identical trees cancel exactly however their Newick
/// happened to list the children.
fn posted_splits(snap: &InternSnap, posting_col: &[u32]) -> Vec<(usize, f64)> {
    let mut posted: Vec<(usize, f64)> = lengths_by_column(snap, posting_col).collect();
    posted.sort_unstable_by_key(|&(col, _)| col);
    posted
}

/// The posting-list splits, column by column: the trees that hold each one and
/// their lengths.
///
/// Column `c` owns entries `offsets[c]..offsets[c + 1]` of `trees` and
/// `lengths`. The trees run ascending because they are appended in tree order.
struct PostingLists {
    offsets: Vec<usize>,
    trees: Vec<u32>,
    lengths: Vec<f64>,
}

impl PostingLists {
    /// `posting_col` and `n_columns` as [`assign_columns`] returns them.
    fn build(snaps: &Snapshots, posting_col: &[u32], n_columns: usize) -> Self {
        let mut offsets = vec![0usize; n_columns + 1];
        for (&count, &col) in snaps.split_counts.iter().zip(posting_col) {
            if col != NO_COLUMN {
                offsets[col as usize] = count as usize;
            }
        }
        let total = exclusive_prefix_sum(&mut offsets);

        let mut next = offsets.clone();
        let (mut trees, mut lengths) = (vec![0u32; total], vec![0.0f64; total]);
        for (tree, snap) in snaps.snapshots.iter().enumerate() {
            for (col, length) in lengths_by_column(snap, posting_col) {
                let at = &mut next[col];
                (trees[*at], lengths[*at]) = (tree as u32, length);
                *at += 1;
            }
        }
        Self {
            offsets,
            trees,
            lengths,
        }
    }

    /// The trees holding posting column `col`, ascending, and their lengths.
    fn column(&self, col: usize) -> (&[u32], &[f64]) {
        let span = self.offsets[col]..self.offsets[col + 1];
        (&self.trees[span.clone()], &self.lengths[span])
    }
}

/// `finish(selfᵢ + selfⱼ − 2·Σ overlap)` — the shape WRF and KF share.
///
/// `overlap` is the per-split shared term, and a split contributes
/// `overlap(l, l)` to its own tree's `self`. `finish` is applied last.
fn weighted_distances(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Vec<f64> {
    let n = snaps.snapshots.len();
    let min_dense = ((n as f64 * WEIGHTED_DENSE_SHARE).ceil() as u32).max(2);
    weighted_distances_split(snaps, progress, min_dense, term, overlap, finish)
}

/// [`weighted_distances`] with the dense/posting boundary given: a split held
/// by at least `min_dense` trees gets a dense column, one held by fewer (but at
/// least two) gets a posting list, and one held by a single tree gets neither,
/// since it can never be shared. Its length folds into that tree's `self`.
///
/// `self` is summed in the same order as the shared term it has to cancel:
/// dense columns by [`sweep`], posted splits in column order, and the two parts
/// added last in both. Identical trees therefore come out at exactly 0.0. The
/// clamp stops rounding from handing `finish` a negative.
fn weighted_distances_split(
    snaps: &Snapshots,
    progress: Option<&AtomicUsize>,
    min_dense: u32,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Vec<f64> {
    let n = snaps.snapshots.len();
    if n == 0 {
        return Vec::new();
    }
    let counts = &snaps.split_counts;

    // "Everywhere" splits get a dense column; unlike in RF they do not cancel
    // out of a weighted score.
    let (dense_col, stride) = assign_columns(counts, n, |count| count >= min_dense);
    let dense = dense_rows(snaps, &dense_col, stride);
    let (posting_col, n_posted) =
        assign_columns(counts, n, |count| count >= 2 && count < min_dense);
    let postings = PostingLists::build(snaps, &posting_col, n_posted);

    let (self_total, posted): (Vec<f64>, Vec<Vec<(usize, f64)>>) = snaps
        .snapshots
        .par_iter()
        .enumerate()
        .map(|(i, snap)| {
            let row = row_slice(&dense, i, stride);
            let dense_self = sweep(row, row, &overlap);
            let posted = posted_splits(snap, &posting_col);
            let posted_self = posted
                .iter()
                .fold(0.0, |sum, &(_, length)| sum + overlap(length, length));
            let unique_self: f64 = snap
                .split_ids
                .iter()
                .zip(&snap.lengths)
                .filter(|&(&id, _)| counts[id as usize] < 2)
                .map(|(_, &length)| overlap(length, length))
                .sum();
            ((dense_self + posted_self) + unique_self, posted)
        })
        .unzip();

    fill_symmetric(n, progress, |i, row: &mut [f64]| {
        // Shared terms from the posting lists first, accumulated in place.
        // Only trees after `i` are wanted: the lower triangle is mirrored.
        for &(col, length) in &posted[i] {
            let (trees, lengths) = postings.column(col);
            let after = trees.partition_point(|&tree| tree as usize <= i);
            for (&j, &other) in trees[after..].iter().zip(&lengths[after..]) {
                row[j as usize] += overlap(length, other);
            }
        }

        let own = row_slice(&dense, i, stride);
        for (j, slot) in row.iter_mut().enumerate().skip(i + 1) {
            let shared = sweep(own, row_slice(&dense, j, stride), &overlap) + *slot;
            *slot = finish((self_total[i] + self_total[j] - 2.0 * shared).max(0.0));
        }
    })
}

/// `WRF(i, j) = Σ lenᵢ + Σ lenⱼ − 2·Σ min(lenᵢ, lenⱼ)`.
///
/// The `min` form follows from `|a − b| = a + b − 2·min(a, b)`. Assumes
/// non-negative branch lengths; missing lengths parse as 0.0.
pub(crate) fn distance_wrf(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    weighted_distances(snaps, progress, f64::min, |d| d)
}

/// `KF(i, j) = sqrt(Σ lenᵢ² + Σ lenⱼ² − 2·Σ lenᵢ·lenⱼ)` — Euclidean distance in
/// branch-length space.
pub(crate) fn distance_kf(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    weighted_distances(snaps, progress, |a, b| a * b, f64::sqrt)
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
    use super::{TREEDIST_TREES, assign_columns, distance_rf_split, weighted_distances_split};
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

    /// Assert `mat` has an exact-zero diagonal and is symmetric.
    fn assert_symmetric_zero_diagonal(mat: &[f64], n: usize, metric: &str) {
        for (i, row) in mat.chunks(n).enumerate() {
            assert_eq!(row[i], 0.0, "{metric} diagonal [{i}]");
            for (j, v) in row.iter().enumerate() {
                assert!(
                    (v - mat[j * n + i]).abs() < f64::EPSILON,
                    "{metric} symmetry [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn weighted_symmetric_zero_diagonal() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        assert_symmetric_zero_diagonal(&snaps.pairwise_wrf(None), n, "WRF");
        assert_symmetric_zero_diagonal(&snaps.pairwise_kf(None), n, "KF");
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

    /// Assert `d(i, k) ≤ d(i, j) + d(j, k)` for every triple, with a small
    /// epsilon for the float metrics.
    fn assert_triangle_inequality<T: Copy + Into<f64>>(mat: &[T], n: usize, metric: &str) {
        let d = |a: usize, b: usize| -> f64 { mat[a * n + b].into() };
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    assert!(
                        d(i, k) <= d(i, j) + d(j, k) + 1e-10,
                        "{metric} triangle inequality failed for indices {i}, {j}, {k}"
                    );
                }
            }
        }
    }

    #[test]
    fn metrics_satisfy_triangle_inequality() {
        let snaps = three_snapshots();
        let n = snaps.snapshots.len();
        assert_triangle_inequality(&snaps.pairwise_rf(None), n, "RF");
        assert_triangle_inequality(&snaps.pairwise_wrf(None), n, "WRF");
        assert_triangle_inequality(&snaps.pairwise_kf(None), n, "KF");
    }

    #[test]
    fn split_counts_match_a_recount() {
        let snaps = three_snapshots();
        let mut recount = vec![0u32; snaps.n_distinct_splits()];
        for snap in &snaps.snapshots {
            for &id in &snap.split_ids {
                recount[id as usize] += 1;
            }
        }
        assert_eq!(snaps.split_counts, recount);
    }

    #[test]
    fn assign_columns_orders_by_descending_count() {
        // Split 1 is held by one tree, so `count >= 2` drops it.
        let (column_of, kept) = assign_columns(&[3, 1, 2, 3], 3, |count| count >= 2);
        assert_eq!(kept, 3);
        assert_eq!(column_of, vec![0, NO_COLUMN, 2, 1]);
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

    /// Assert `got` is within a relative 1e-9 of the oracle's `want`.
    fn assert_close(got: f64, want: f64, what: &str) {
        assert!(
            (got - want).abs() <= 1e-9 * want.max(1.0),
            "{what}: {got} vs reference {want}"
        );
    }

    /// Assert WRF and KF matrices agree with [`reference_distances`], and that
    /// copies of the same Newick land on exactly 0.0, not a rounding crumb.
    fn assert_weighted_match_reference<S: AsRef<str>>(
        snaps: &Snapshots,
        newicks: &[S],
        (wrf, kf): (&[f64], &[f64]),
        ctx: &str,
    ) {
        let n = snaps.snapshots.len();
        for i in 0..n {
            for j in 0..n {
                let (_, want_wrf, want_kf) =
                    reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);
                let at = i * n + j;
                assert_close(wrf[at], want_wrf, &format!("WRF at [{i}][{j}], {ctx}"));
                assert_close(kf[at], want_kf, &format!("KF at [{i}][{j}], {ctx}"));
                if newicks[i].as_ref() == newicks[j].as_ref() {
                    assert_eq!(wrf[at], 0.0, "WRF identical [{i}][{j}], {ctx}");
                    assert_eq!(kf[at], 0.0, "KF identical [{i}][{j}], {ctx}");
                }
            }
        }
    }

    /// Assert all three metrics agree with [`reference_distances`] on the
    /// snapshots built from `newicks`.
    fn assert_matches_reference<S: AsRef<str>>(snaps: &Snapshots, newicks: &[S], ctx: &str) {
        let n = snaps.snapshots.len();
        let rf = snaps.pairwise_rf(None);
        for i in 0..n {
            for j in 0..n {
                let (want_rf, _, _) = reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);
                assert_eq!(rf[i * n + j] as usize, want_rf, "RF at [{i}][{j}], {ctx}");
            }
        }
        let (wrf, kf) = (snaps.pairwise_wrf(None), snaps.pairwise_kf(None));
        assert_weighted_match_reference(snaps, newicks, (&wrf, &kf), ctx);
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
    fn pairwise_weighted_identical_trees_all_zero() {
        let t = TREEDIST_TREES[0];
        let snaps = Snapshots::from_newicks(&[t, t, t], false).unwrap();

        // Exactly 0.0 only holds because self and shared are summed in the same
        // column order; reordering either would leave a rounding crumb.
        assert!(
            snaps.pairwise_wrf(None).iter().all(|&d| d == 0.0),
            "identical trees must have WRF exactly 0 everywhere"
        );

        // Identical trees make `selfᵢ + selfⱼ − 2·dot` algebraically zero, but
        // rounding can nudge it slightly negative — the clamp must keep the
        // square root from producing NaN.
        assert!(
            snaps.pairwise_kf(None).iter().all(|&d| d == 0.0),
            "identical trees must have KF 0 everywhere (clamp must avoid NaN)"
        );
    }

    /// Trees sharing almost nothing, so most splits land in the "held by one
    /// tree only" bucket, which gets no column and folds into its tree's `self`
    /// term — which the treedist fixtures barely exercise. Branch lengths are
    /// all distinct so a wrongly dropped split changes the answer visibly.
    const DIVERSE_TREES: [&str; 4] = [
        "(((A:0.11,B:0.12):0.13,(C:0.14,D:0.15):0.16):0.17,((E:0.18,F:0.19):0.21,(G:0.22,H:0.23):0.24):0.25);",
        "(((A:0.31,E:0.32):0.33,(C:0.34,G:0.35):0.36):0.37,((B:0.38,F:0.39):0.41,(D:0.42,H:0.43):0.44):0.45);",
        "(((A:0.51,H:0.52):0.53,(B:0.54,G:0.55):0.56):0.57,((C:0.58,F:0.59):0.61,(D:0.62,E:0.63):0.64):0.65);",
        "(((A:0.71,D:0.72):0.73,(F:0.74,G:0.75):0.76):0.77,((B:0.78,H:0.79):0.81,(C:0.82,E:0.83):0.84):0.85);",
    ];

    #[test]
    fn pairwise_weighted_metrics_match_reference_on_diverse_trees() {
        let snaps = Snapshots::from_newicks(&DIVERSE_TREES, false).unwrap();
        assert_matches_reference(&snaps, &DIVERSE_TREES, "diverse fixtures");
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

    /// `n_trees` random trees over `n_taxa` leaves followed by `duplicates`
    /// copies of the first ones, as Newicks and as the snapshots built from them.
    fn random_snapshots(
        n_taxa: usize,
        n_trees: usize,
        duplicates: usize,
        seed: u64,
    ) -> (Vec<String>, Snapshots) {
        let mut state = seed;
        let mut newicks: Vec<String> = (0..n_trees)
            .map(|_| random_newick(n_taxa, &mut state))
            .collect();
        for k in 0..duplicates {
            newicks.push(newicks[k % n_trees].clone());
        }
        let refs: Vec<&str> = newicks.iter().map(String::as_str).collect();
        let snaps = Snapshots::from_newicks(&refs, false).unwrap();
        (newicks, snaps)
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
            let (newicks, snaps) = random_snapshots(n_taxa, n_trees, duplicates, seed);
            assert_matches_reference(&snaps, &newicks, &format!("seed={seed} taxa={n_taxa}"));
        }
    }

    /// WRF and KF with the dense/posting boundary at `min_dense`, mirroring
    /// [`super::distance_wrf`] and [`super::distance_kf`].
    fn weighted_at(snaps: &Snapshots, min_dense: u32) -> (Vec<f64>, Vec<f64>) {
        (
            weighted_distances_split(snaps, None, min_dense, f64::min, |d| d),
            weighted_distances_split(snaps, None, min_dense, |a, b| a * b, f64::sqrt),
        )
    }

    /// The boundary decides which path a split takes and must decide nothing
    /// else: all dense (2), mixed, and all posting lists (`n + 1`) all have to
    /// match the oracle.
    #[test]
    fn weighted_boundary_does_not_change_distances() {
        for &(n_taxa, n_trees, duplicates, seed) in &[
            (12usize, 8usize, 4usize, 11u64),
            (24, 10, 6, 12),
            (40, 6, 9, 13),
        ] {
            let (newicks, snaps) = random_snapshots(n_taxa, n_trees, duplicates, seed);
            let n = snaps.snapshots.len();
            for min_dense in [2, 3, n as u32 / 2, n as u32, n as u32 + 1] {
                let (wrf, kf) = weighted_at(&snaps, min_dense);
                let ctx = format!("seed={seed} min_dense={min_dense}");
                assert_weighted_match_reference(&snaps, &newicks, (&wrf, &kf), &ctx);
            }
        }
    }

    /// Same as [`weighted_boundary_does_not_change_distances`] for RF: all
    /// posting lists (2), mixed, and all bit columns (`n + 1`) must match the
    /// oracle exactly.
    #[test]
    fn rf_boundary_does_not_change_distances() {
        for &(n_taxa, n_trees, duplicates, seed) in &[
            (12usize, 8usize, 4usize, 31u64),
            (24, 10, 6, 32),
            (40, 6, 9, 33),
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
            let n = snaps.snapshots.len();

            for min_dense in [2, 3, n as u32 / 2, n as u32, n as u32 + 1] {
                let rf = distance_rf_split(&snaps, None, min_dense);
                for i in 0..n {
                    for j in 0..n {
                        let (want, _, _) =
                            reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);
                        assert_eq!(
                            rf[i * n + j] as usize,
                            want,
                            "RF seed={seed} min_dense={min_dense} [{i}][{j}]"
                        );
                    }
                }
            }
        }
    }

    /// The same tree written with its children in a different order parses its
    /// splits in a different order. Posted splits are summed in column order,
    /// not parse order, so the copies must still cancel to exactly 0.0 on every
    /// path.
    #[test]
    fn identical_trees_in_any_child_order_cancel_exactly() {
        let trees = [
            "(((A:0.31,B:0.17):0.23,(C:0.41,D:0.13):0.29):0.11,((E:0.37,F:0.19):0.07,(G:0.43,H:0.03):0.47):0.53);",
            "(((H:0.03,G:0.43):0.47,(F:0.19,E:0.37):0.07):0.53,((D:0.13,C:0.41):0.29,(B:0.17,A:0.31):0.23):0.11);",
            "(((A:0.31,E:0.12):0.33,(C:0.34,G:0.35):0.36):0.37,((B:0.38,F:0.39):0.41,(D:0.42,H:0.43):0.44):0.45);",
            "(((C:0.41,D:0.13):0.29,(B:0.17,A:0.31):0.23):0.11,((G:0.43,H:0.03):0.47,(E:0.37,F:0.19):0.07):0.53);",
        ];
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();
        let n = snaps.snapshots.len();
        // 2: all dense. 4: pendants dense, the shared internal splits posted.
        // n + 1: everything posted.
        for min_dense in [2, 4, n as u32 + 1] {
            let (wrf, kf) = weighted_at(&snaps, min_dense);
            for (i, j) in [(0, 1), (0, 3), (1, 3)] {
                assert_eq!(wrf[i * n + j], 0.0, "WRF [{i}][{j}] min_dense={min_dense}");
                assert_eq!(kf[i * n + j], 0.0, "KF [{i}][{j}] min_dense={min_dense}");
            }
            assert!(wrf[2] > 0.0, "tree 2 differs, min_dense={min_dense}");
        }
    }
}
