//! The weighted metrics: weighted Robinson–Foulds and Kuhner–Felsenstein.

use super::boundary::{min_dense, weighted_dense_share};
use super::layout::{Layout, Trees, block_size, column};
use super::matrix::fill_symmetric_banded;
use super::postings::Postings;
use crate::par::*;
use crate::snapshot::Snapshots;
use std::sync::atomic::AtomicUsize;

/// `Σ overlap(aₖ, bₖ)` over two rows of equal length, in a fixed order.
///
/// Eight running sums rather than one. A single `f64` sum is a chain the
/// compiler may not reorder, so it can neither vectorise nor pipeline; eight
/// independent ones can do both. The tail shorter than eight lands in the first
/// lanes. The order is set here rather than by the compiler, so `sweep(a, a)`
/// and `sweep(a, b)` add the same terms the same way whenever `a == b`, which
/// is what keeps identical trees at exactly 0.0.
#[inline]
pub(super) fn sweep(a: &[f64], b: &[f64], overlap: &impl Fn(f64, f64) -> f64) -> f64 {
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

/// Rows per weighted sweep task, rows of the other side held against them,
/// and columns per pass: a tile of `TILE_I × TILE_J` pairs reads each row's
/// panel once instead of once per pair, which turns the sweep from bound by
/// memory bandwidth into bound by arithmetic.
const TILE_I: usize = 8;
const TILE_J: usize = 4;
pub(super) const PANEL: usize = 512;

/// Add `overlap(a[k], b[k])` into `lanes` exactly as [`sweep`] does over the
/// same columns: whole blocks of eight first, then the tail into the first
/// lanes. Called panel by panel, with every panel but the last a multiple of
/// eight wide, the lanes end up holding exactly what [`sweep`] adds up.
#[inline(always)]
fn accumulate(lanes: &mut [f64; 8], a: &[f64], b: &[f64], overlap: &impl Fn(f64, f64) -> f64) {
    let (a_blocks, a_tail) = a.as_chunks::<8>();
    let (b_blocks, b_tail) = b.as_chunks::<8>();
    for (xs, ys) in a_blocks.iter().zip(b_blocks) {
        for ((lane, &x), &y) in lanes.iter_mut().zip(xs).zip(ys) {
            *lane += overlap(x, y);
        }
    }
    for ((lane, &x), &y) in lanes.iter_mut().zip(a_tail).zip(b_tail) {
        *lane += overlap(x, y);
    }
}

/// Every tree's lengths over the dense columns, 0.0 where it lacks the split.
/// Stored a block of trees per buffer, so that they can be allocated as the
/// trees are consumed.
pub(super) struct DenseRows {
    size: usize,
    pub(super) stride: usize,
    blocks: Vec<Vec<f64>>,
}

impl DenseRows {
    /// Tree `i`'s row.
    #[inline]
    pub(super) fn row(&self, i: usize) -> &[f64] {
        let (block, k) = (i / self.size, i % self.size);
        &self.blocks[block][k * self.stride..][..self.stride]
    }
}

/// The weighted kernels' dense rows and posting lists, laid out a block of
/// trees at a time, with every tree's `self` term.
///
/// `self` is summed in the same order as the shared term it has to cancel:
/// dense columns by [`sweep`], posted splits in column order, and the two parts
/// added last in both; a split held by one tree alone folds in after them.
/// Identical trees therefore come out at exactly 0.0.
pub(super) fn weighted_layout(
    trees: Trees,
    counts: &[u32],
    layout: &Layout,
    overlap: &(impl Fn(f64, f64) -> f64 + Sync),
) -> (DenseRows, Postings, Vec<f64>) {
    let n = trees.len();
    let (size, stride) = (block_size(), layout.n_dense);
    // A zero-width row still gets one scratch slot, so that every tree is
    // handed a row; only its first `stride` (none) are read.
    let width = stride.max(1);
    let mut dense = DenseRows {
        size,
        stride,
        blocks: Vec::with_capacity(n.div_ceil(size)),
    };
    let mut postings = Postings::sized(counts, layout, size, true);
    let mut self_total = Vec::with_capacity(n);
    trees.by_blocks(size, |block| {
        let mut rows = vec![0.0f64; block.len() * width];
        let (posted, selfs): (Vec<Vec<(u32, f64)>>, Vec<f64>) = rows
            .par_chunks_mut(width)
            .zip(block.par_iter())
            .map(|(row, snap)| {
                let row = &mut row[..stride];
                let mut own = Vec::new();
                let mut unique_self = 0.0;
                for (id, &length) in snap.ids().zip(&snap.lengths) {
                    match column(&layout.column_of, id) {
                        Some(col) if col < stride => row[col] = length,
                        Some(col) => own.push(((col - stride) as u32, length)),
                        None => unique_self += overlap(length, length),
                    }
                }
                own.sort_unstable_by_key(|&(col, _)| col);
                let dense_self = sweep(row, row, overlap);
                let posted_self = own
                    .iter()
                    .fold(0.0, |sum, &(_, length)| sum + overlap(length, length));
                (own, (dense_self + posted_self) + unique_self)
            })
            .unzip();
        postings.push_block(&posted);
        self_total.extend(selfs);
        dense.blocks.push(rows);
    });
    postings.fill = Vec::new();
    (dense, postings, self_total)
}

/// `finish(selfᵢ + selfⱼ − 2·Σ overlap)` — the shape WRF and KF share.
///
/// `overlap` is the per-split shared term, and a split contributes
/// `overlap(l, l)` to its own tree's `self`. `finish` is applied last.
fn weighted_distances(
    trees: Trees,
    counts: &[u32],
    progress: Option<&AtomicUsize>,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Vec<f64> {
    let min_dense = min_dense(trees.len(), weighted_dense_share());
    weighted_distances_split(trees, counts, progress, min_dense, overlap, finish)
}

/// [`weighted_distances`] with the dense/posting boundary given: a split held
/// by at least `min_dense` trees gets a dense column, one held by fewer (but at
/// least two) gets a posting list, and one held by a single tree gets neither,
/// since it can never be shared. Its length folds into that tree's `self`.
/// The clamp stops rounding from handing `finish` a negative.
pub(super) fn weighted_distances_split(
    trees: Trees,
    counts: &[u32],
    progress: Option<&AtomicUsize>,
    min_dense: u32,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Vec<f64> {
    let n = trees.len();
    if n == 0 {
        return Vec::new();
    }

    // "Everywhere" splits get a dense column; unlike in RF they do not cancel
    // out of a weighted score.
    let layout = Layout::new(counts, n, min_dense, |count| count >= 2);
    let (dense, postings, self_total) = weighted_layout(trees, counts, &layout, &overlap);
    drop(layout);

    fill_symmetric_banded(n, TILE_I, progress, |i0, rows: &mut [f64]| {
        let band = rows.len() / n;
        // Shared terms from the posting lists first, accumulated in place.
        for (r, row) in rows.chunks_mut(n).enumerate() {
            let i = i0 + r;
            let (cols, lengths) = postings.of_tree(i);
            for (&col, &length) in cols.iter().zip(lengths) {
                let (trees, others) = postings.after(col, i);
                for (&j, &other) in trees.iter().zip(others) {
                    row[j as usize] += overlap(length, other);
                }
            }
        }

        // Then the dense columns, a tile of pairs at a time, a panel of
        // columns per pass, each pair keeping its own eight lanes throughout.
        for j0 in (i0 + 1..n).step_by(TILE_J) {
            let width = (n - j0).min(TILE_J);
            let mut lanes = [[[0.0f64; 8]; TILE_J]; TILE_I];
            for p0 in (0..dense.stride).step_by(PANEL) {
                let panel = p0..(p0 + PANEL).min(dense.stride);
                for (r, lanes_r) in lanes.iter_mut().enumerate().take(band) {
                    let a = &dense.row(i0 + r)[panel.clone()];
                    for (c, lane) in lanes_r.iter_mut().enumerate().take(width) {
                        if j0 + c > i0 + r {
                            accumulate(lane, a, &dense.row(j0 + c)[panel.clone()], &overlap);
                        }
                    }
                }
            }
            for (r, lanes_r) in lanes.iter().enumerate().take(band) {
                let i = i0 + r;
                for (c, lane) in lanes_r.iter().enumerate().take(width) {
                    let j = j0 + c;
                    if j > i {
                        let slot = &mut rows[r * n + j];
                        let shared = lane.iter().sum::<f64>() + *slot;
                        *slot = finish((self_total[i] + self_total[j] - 2.0 * shared).max(0.0));
                    }
                }
            }
        }
    })
}

/// The weighted metrics' trees and split counts, taken out of a collection the
/// caller is done with.
fn owned_trees(snaps: Snapshots) -> (Trees<'static>, Vec<u32>) {
    let Snapshots {
        snapshots,
        split_counts,
        ..
    } = snaps;
    (Trees::Owned(snapshots), split_counts)
}

/// `WRF(i, j) = Σ lenᵢ + Σ lenⱼ − 2·Σ min(lenᵢ, lenⱼ)`.
///
/// The `min` form follows from `|a − b| = a + b − 2·min(a, b)`. Assumes
/// non-negative branch lengths; missing lengths parse as 0.0.
pub(crate) fn distance_wrf(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    let trees = Trees::Borrowed(&snaps.snapshots);
    weighted_distances(trees, &snaps.split_counts, progress, f64::min, |d| d)
}

/// [`distance_wrf`] for a caller that is done with `snaps`: each tree's split
/// IDs and branch lengths are dropped, a block at a time, as the kernel lays
/// them out.
pub(crate) fn distance_wrf_owned(snaps: Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    let (trees, counts) = owned_trees(snaps);
    weighted_distances(trees, &counts, progress, f64::min, |d| d)
}

/// `KF(i, j) = sqrt(Σ lenᵢ² + Σ lenⱼ² − 2·Σ lenᵢ·lenⱼ)` — Euclidean distance in
/// branch-length space.
pub(crate) fn distance_kf(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    let trees = Trees::Borrowed(&snaps.snapshots);
    weighted_distances(
        trees,
        &snaps.split_counts,
        progress,
        |a, b| a * b,
        f64::sqrt,
    )
}

/// [`distance_kf`] for a caller that is done with `snaps`; see
/// [`distance_wrf_owned`].
pub(crate) fn distance_kf_owned(snaps: Snapshots, progress: Option<&AtomicUsize>) -> Vec<f64> {
    let (trees, counts) = owned_trees(snaps);
    weighted_distances(trees, &counts, progress, |a, b| a * b, f64::sqrt)
}
