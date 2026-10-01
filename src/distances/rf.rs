//! Robinson–Foulds.

use super::boundary::{min_dense, rf_dense_share};
use super::layout::{Column, Input, Laid, Layout, Rows, lay_out};
use super::matrix::{fill_symmetric, row_slice};
use std::sync::atomic::AtomicUsize;

/// Fewest trees at which RF gives any split a posting list.
const RF_MIN_TREES_FOR_POSTINGS: usize = 256;

/// One bitmask row per tree over RF's bit columns: bit `c` of row `i` is set
/// when tree `i` holds the split in column `c`.
///
/// `spans[i]` is the half-open range of words where row `i` has bits set,
/// `(0, 0)` when it has none. Columns run in descending tree count, so a row's
/// bits cluster and a pair only sweeps the words where both spans overlap.
///
/// The rows are kept flat rather than in [`Rows`]' blocks: a pair reads two
/// short rows and does little else, so a block lookup per row would show in
/// the sweep, and bit rows are small enough to copy once.
struct BitRows {
    pub(super) words: usize,
    pub(super) packed: Vec<u64>,
    pub(super) spans: Vec<(usize, usize)>,
}

impl BitRows {
    /// The `n` trees' `rows`, flattened, with the span of each.
    fn new(rows: Rows<u64>, n: usize) -> Self {
        let words = rows.stride;
        let packed = rows.into_flat();
        let spans = (0..n)
            .map(|i| {
                let row = row_slice(&packed, i, words);
                let hi = row
                    .iter()
                    .rposition(|&word| word != 0)
                    .map_or(0, |last| last + 1);
                let lo = row.iter().position(|&word| word != 0).unwrap_or(hi);
                (lo, hi)
            })
            .collect();
        Self {
            words,
            packed,
            spans,
        }
    }

    /// How many bit columns rows `i` and `j` both have set.
    #[inline]
    fn shared(&self, i: usize, j: usize) -> u32 {
        let lo = self.spans[i].0.max(self.spans[j].0);
        let hi = self.spans[i].1.min(self.spans[j].1);
        if lo >= hi {
            return 0;
        }
        let a = &row_slice(&self.packed, i, self.words)[lo..hi];
        let b = &row_slice(&self.packed, j, self.words)[lo..hi];
        a.iter().zip(b).map(|(&x, &y)| (x & y).count_ones()).sum()
    }
}

/// `RF(i, j) = selfᵢ + selfⱼ − 2·shared(i, j)`, where `self` counts a tree's
/// splits and `shared` the splits both trees hold.
///
/// Splits held by *every* tree add equally to both `self` values and to the
/// shared count, so they cancel exactly and are dropped. Splits held by one tree
/// count towards that tree's `self` but can never be shared, so they get nothing
/// either. On a posterior that is most of the distinct splits. The rest are
/// shared by a popcount over packed bit-rows for the widely held ones and by
/// posting lists for the rare ones.
pub(crate) fn distance_rf(input: Input<'_>, progress: Option<&AtomicUsize>) -> Vec<u32> {
    let n = input.trees.len();
    let min_dense = if n < RF_MIN_TREES_FOR_POSTINGS {
        2 // every shareable split gets a bit column
    } else {
        min_dense(n, rf_dense_share())
    };
    distance_rf_split(input, progress, min_dense)
}

/// [`distance_rf`] with the bit/posting boundary given: a split held by at
/// least `min_dense` trees (but not all) gets a bit column, one held by fewer
/// (but at least two) gets a posting list.
pub(super) fn distance_rf_split(
    input: Input<'_>,
    progress: Option<&AtomicUsize>,
    min_dense: u32,
) -> Vec<u32> {
    let n = input.trees.len();
    if n == 0 {
        return Vec::new();
    }

    // A split held by all `n` trees cancels, so it gets no column of either kind.
    let all = n as u32;
    let layout = Layout::new(&input.counts, n, min_dense.min(all), |count| {
        (2..all).contains(&count)
    });
    // Every tree holds each of the `n_universal` splits, so its `self` is its
    // split count minus that fixed number.
    let n_universal = input.counts.iter().filter(|&&count| count == all).count();
    let words = layout.n_dense.div_ceil(64);
    let Laid {
        rows,
        postings,
        selfs: self_count,
    } = lay_out(
        input,
        layout,
        words,
        false,
        |layout, snap, row: &mut [u64]| {
            let mut posted = Vec::new();
            for id in snap.ids() {
                match layout.place(id) {
                    Some(Column::Dense(col)) => row[col / 64] |= 1u64 << (col % 64),
                    Some(Column::Posted(col)) => posted.push((col, 0.0)),
                    None => {}
                }
            }
            posted.sort_unstable_by_key(|&(col, _)| col);
            (posted, (snap.n_splits() - n_universal) as u32)
        },
    );
    let bits = BitRows::new(rows, n);

    fill_symmetric(n, 1, progress, |i, row: &mut [u32]| {
        // Shared splits from the posting lists first, counted in place.
        for &col in postings.of_tree(i).0 {
            for &j in postings.after(col, i).0 {
                row[j as usize] += 1;
            }
        }

        for (j, slot) in row.iter_mut().enumerate().skip(i + 1) {
            let shared = bits.shared(i, j) + *slot;
            *slot = self_count[i] + self_count[j] - 2 * shared;
        }
    })
}
