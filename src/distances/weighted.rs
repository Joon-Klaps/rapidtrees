//! The weighted metrics: weighted Robinson–Foulds and Kuhner–Felsenstein.

use super::boundary::{min_dense, weighted_dense_share};
use super::layout::{Column, Input, Laid, Layout, lay_out};
use super::matrix::fill_symmetric;
use std::sync::atomic::AtomicUsize;

/// Running sums per weighted sweep: see [`sweep`].
const LANES: usize = 8;

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
    accumulate([0.0; LANES], a, b, overlap).iter().sum()
}

/// Rows per weighted sweep task, rows of the other side held against them,
/// and columns per pass: a tile of `TILE_I × TILE_J` pairs reads each row's
/// panel once instead of once per pair, which turns the sweep from bound by
/// memory bandwidth into bound by arithmetic.
const TILE_I: usize = 8;
const TILE_J: usize = 4;
pub(super) const PANEL: usize = 512;
// Every panel but the last must end on a block boundary for `accumulate` to
// add in `sweep`'s order.
const _: () = assert!(PANEL.is_multiple_of(LANES));

/// `lanes` with `overlap(a[k], b[k])` added in [`sweep`]'s order: whole
/// blocks of eight first, then the tail into the first lanes. Called panel by
/// panel, with every panel but the last a multiple of eight wide, the lanes end
/// up holding exactly what one call over the whole rows would.
///
/// The lanes go in and come out by value. Behind a `&mut` the compiler cannot
/// always prove that they do not alias the rows, and then stores every lane
/// back on every block and gives up on vectorising: twice the instructions.
#[inline(always)]
fn accumulate(
    mut lanes: [f64; LANES],
    a: &[f64],
    b: &[f64],
    overlap: &impl Fn(f64, f64) -> f64,
) -> [f64; LANES] {
    let (a_blocks, a_tail) = a.as_chunks::<LANES>();
    let (b_blocks, b_tail) = b.as_chunks::<LANES>();
    for (xs, ys) in a_blocks.iter().zip(b_blocks) {
        for ((lane, &x), &y) in lanes.iter_mut().zip(xs).zip(ys) {
            *lane += overlap(x, y);
        }
    }
    for ((lane, &x), &y) in lanes.iter_mut().zip(a_tail).zip(b_tail) {
        *lane += overlap(x, y);
    }
    lanes
}

/// The weighted kernels' trees laid out with the dense/posting boundary at
/// `min_dense`: each tree's lengths over the dense columns, 0.0 where it lacks
/// the split, its posted splits with their lengths, and its `self` term.
///
/// `self` is summed in the same order as the shared term it has to cancel:
/// dense columns by [`sweep`], posted splits in column order, and the two parts
/// added last in both; a split held by one tree alone folds in after them.
/// Identical trees therefore come out at exactly 0.0.
pub(super) fn weighted_lay_out(
    input: Input<'_>,
    min_dense: u32,
    overlap: &(impl Fn(f64, f64) -> f64 + Sync),
) -> Laid<f64, f64> {
    // "Everywhere" splits get a dense column; unlike in RF they do not cancel
    // out of a weighted score.
    let layout = Layout::new(&input.counts, input.trees.len(), min_dense, |count| {
        count >= 2
    });
    let stride = layout.n_dense;
    lay_out(
        input,
        layout,
        stride,
        true,
        |layout, snap, row: &mut [f64]| {
            let mut posted = Vec::new();
            let mut unique_self = 0.0;
            for (id, &length) in snap.ids().zip(&snap.lengths) {
                match layout.place(id) {
                    Some(Column::Dense(col)) => row[col] = length,
                    Some(Column::Posted(col)) => posted.push((col, length)),
                    None => unique_self += overlap(length, length),
                }
            }
            posted.sort_unstable_by_key(|&(col, _)| col);
            let dense_self = sweep(row, row, overlap);
            let posted_self = posted
                .iter()
                .fold(0.0, |sum, &(_, length)| sum + overlap(length, length));
            (posted, (dense_self + posted_self) + unique_self)
        },
    )
}

/// `finish(selfᵢ + selfⱼ − 2·Σ overlap)` — the shape WRF and KF share.
///
/// `overlap` is the per-split shared term, and a split contributes
/// `overlap(l, l)` to its own tree's `self`. `finish` is applied last.
fn weighted_distances(
    input: Input<'_>,
    progress: Option<&AtomicUsize>,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Vec<f64> {
    let min_dense = min_dense(input.trees.len(), weighted_dense_share());
    weighted_distances_split(input, progress, min_dense, overlap, finish)
}

/// [`weighted_distances`] with the dense/posting boundary given: a split held
/// by at least `min_dense` trees gets a dense column, one held by fewer (but at
/// least two) gets a posting list, and one held by a single tree gets neither,
/// since it can never be shared. Its length folds into that tree's `self`.
/// The clamp stops rounding from handing `finish` a negative.
pub(super) fn weighted_distances_split(
    input: Input<'_>,
    progress: Option<&AtomicUsize>,
    min_dense: u32,
    overlap: impl Fn(f64, f64) -> f64 + Sync,
    finish: impl Fn(f64) -> f64 + Sync,
) -> Vec<f64> {
    let n = input.trees.len();
    if n == 0 {
        return Vec::new();
    }
    let Laid {
        rows: dense,
        postings,
        selfs: self_total,
    } = weighted_lay_out(input, min_dense, &overlap);

    fill_symmetric(n, TILE_I, progress, |i0, rows: &mut [f64]| {
        let i1 = i0 + rows.len() / n;
        // Shared terms from the posting lists first, accumulated in place.
        for (i, row) in (i0..i1).zip(rows.chunks_mut(n)) {
            postings.add_shared(i, row, &overlap);
        }

        // Then the dense columns, a tile of pairs at a time, a panel of
        // columns per pass, each pair keeping its own eight lanes throughout.
        for j0 in (i0 + 1..n).step_by(TILE_J) {
            let j1 = (j0 + TILE_J).min(n);
            // Row `i` skips the tile's trees up to and including itself.
            let at_or_before = |i: usize| (i + 1).saturating_sub(j0);
            let mut lanes = [[[0.0f64; LANES]; TILE_J]; TILE_I];
            for p0 in (0..dense.stride).step_by(PANEL) {
                let panel = p0..(p0 + PANEL).min(dense.stride);
                for (i, lanes_i) in (i0..i1).zip(&mut lanes) {
                    let a = &dense.row(i)[panel.clone()];
                    for (j, lane) in (j0..j1).zip(lanes_i).skip(at_or_before(i)) {
                        *lane = accumulate(*lane, a, &dense.row(j)[panel.clone()], &overlap);
                    }
                }
            }
            for ((i, row), lanes_i) in (i0..i1).zip(rows.chunks_mut(n)).zip(&lanes) {
                let tile = (j0..j1).zip(lanes_i).zip(&mut row[j0..j1]);
                for ((j, lane), slot) in tile.skip(at_or_before(i)) {
                    let shared = lane.iter().sum::<f64>() + *slot;
                    *slot = finish((self_total[i] + self_total[j] - 2.0 * shared).max(0.0));
                }
            }
        }
    })
}

/// `WRF(i, j) = Σ lenᵢ + Σ lenⱼ − 2·Σ min(lenᵢ, lenⱼ)`.
///
/// The `min` form follows from `|a − b| = a + b − 2·min(a, b)`. Assumes
/// non-negative branch lengths; missing lengths parse as 0.0.
pub(crate) fn distance_wrf(input: Input<'_>, progress: Option<&AtomicUsize>) -> Vec<f64> {
    weighted_distances(input, progress, f64::min, |d| d)
}

/// `KF(i, j) = sqrt(Σ lenᵢ² + Σ lenⱼ² − 2·Σ lenᵢ·lenⱼ)` — Euclidean distance in
/// branch-length space.
pub(crate) fn distance_kf(input: Input<'_>, progress: Option<&AtomicUsize>) -> Vec<f64> {
    weighted_distances(input, progress, |a, b| a * b, f64::sqrt)
}
