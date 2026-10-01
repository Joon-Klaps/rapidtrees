//! Robinson–Foulds.

use super::boundary::{min_dense, rf_dense_share};
use super::layout::{Layout, Trees, block_size, column};
use super::matrix::{fill_symmetric, row_slice};
use super::postings::Postings;
use crate::par::*;
use crate::snapshot::Snapshots;
use std::sync::atomic::AtomicUsize;

/// Fewest trees at which RF gives any split a posting list.
const RF_MIN_TREES_FOR_POSTINGS: usize = 256;

/// One bitmask row per tree over RF's bit columns: bit `c` of row `i` is set
/// when tree `i` holds the split in column `c`.
///
/// `spans[i]` is the half-open range of words where row `i` has bits set,
/// `(0, 0)` when it has none. Columns run in descending tree count, so a row's
/// bits cluster and a pair only sweeps the words where both spans overlap.
struct BitRows {
    words: usize,
    packed: Vec<u64>,
    spans: Vec<(usize, usize)>,
}

impl BitRows {
    /// The `n` rows packed `words` to a row in `packed`.
    fn from_packed(packed: Vec<u64>, words: usize, n: usize) -> Self {
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

/// RF's bit rows and posting lists, laid out a block of trees at a time, with
/// every tree's `self` count: its splits less the `n_universal` that every
/// tree holds.
fn rf_layout(
    trees: Trees,
    counts: &[u32],
    layout: &Layout,
    n_universal: usize,
) -> (BitRows, Postings, Vec<u32>) {
    let n = trees.len();
    let words = layout.n_dense.div_ceil(64);
    // A zero-width row still gets one scratch word, so that every tree is
    // handed a row; no bit is ever set in it.
    let width = words.max(1);
    let mut packed = vec![0u64; n * width];
    let size = block_size();
    let mut postings = Postings::sized(counts, layout, size, false);
    let mut self_count = Vec::with_capacity(n);
    let mut first = 0;
    trees.by_blocks(size, |block| {
        let rows = &mut packed[first * width..(first + block.len()) * width];
        let posted: Vec<Vec<(u32, f64)>> = rows
            .par_chunks_mut(width)
            .zip(block.par_iter())
            .map(|(row, snap)| {
                let mut own = Vec::new();
                for &id in &snap.split_ids {
                    match column(&layout.column_of, id) {
                        Some(col) if col < layout.n_dense => row[col / 64] |= 1u64 << (col % 64),
                        Some(col) => own.push(((col - layout.n_dense) as u32, 0.0)),
                        None => {}
                    }
                }
                own.sort_unstable_by_key(|&(col, _)| col);
                own
            })
            .collect();
        postings.push_block(&posted);
        self_count.extend(
            block
                .iter()
                .map(|snap| (snap.n_splits() - n_universal) as u32),
        );
        first += block.len();
    });
    postings.fill = Vec::new();
    if words == 0 {
        packed = Vec::new();
    }
    (BitRows::from_packed(packed, words, n), postings, self_count)
}

/// The bit/posting boundary RF uses for `n` trees: a split held by fewer than
/// [`rf_dense_share`] of them gets a posting list, once there are enough trees
/// for lists to pay.
fn rf_min_dense(n: usize) -> u32 {
    if n < RF_MIN_TREES_FOR_POSTINGS {
        2 // every shareable split gets a bit column
    } else {
        min_dense(n, rf_dense_share())
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
pub(crate) fn distance_rf(snaps: &Snapshots, progress: Option<&AtomicUsize>) -> Vec<u32> {
    let trees = Trees::Borrowed(&snaps.snapshots);
    let min_dense = rf_min_dense(trees.len());
    distance_rf_split(trees, &snaps.split_counts, progress, min_dense)
}

/// [`distance_rf`] for a caller that is done with `snaps`: each block of trees
/// is dropped once its bit rows and posting entries are written, before the
/// pairwise sweep.
pub(crate) fn distance_rf_owned(snaps: Snapshots, progress: Option<&AtomicUsize>) -> Vec<u32> {
    let Snapshots {
        snapshots,
        split_counts,
        ..
    } = snaps;
    let min_dense = rf_min_dense(snapshots.len());
    distance_rf_split(Trees::Owned(snapshots), &split_counts, progress, min_dense)
}

/// [`distance_rf`] with the bit/posting boundary given: a split held by at
/// least `min_dense` trees (but not all) gets a bit column, one held by fewer
/// (but at least two) gets a posting list.
pub(super) fn distance_rf_split(
    trees: Trees,
    counts: &[u32],
    progress: Option<&AtomicUsize>,
    min_dense: u32,
) -> Vec<u32> {
    let n = trees.len();
    if n == 0 {
        return Vec::new();
    }

    // A split held by all `n` trees cancels, so it gets no column of either kind.
    let all = n as u32;
    let layout = Layout::new(counts, n, min_dense.min(all), |count| {
        (2..all).contains(&count)
    });
    // Every tree holds each of the `n_universal` splits, so its `self` is its
    // split count minus that fixed number.
    let n_universal = counts.iter().filter(|&&count| count == all).count();
    let (bits, postings, self_count) = rf_layout(trees, counts, &layout, n_universal);
    drop(layout);

    fill_symmetric(n, progress, |i, row: &mut [u32]| {
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
