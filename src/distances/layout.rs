//! Column layout of the splits, and the trees a kernel reads.

use super::matrix::exclusive_prefix_sum;
use crate::par::current_num_threads;
use crate::snapshot::InternSnap;

/// Marks a split that [`assign_columns`] gave no column.
pub(super) const NO_COLUMN: u32 = u32::MAX;

/// The column of split `id`, or `None` if `column_of` gave it none.
#[inline]
pub(super) fn column(column_of: &[u32], id: u32) -> Option<usize> {
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
pub(super) fn assign_columns(
    counts: &[u32],
    n_trees: usize,
    keep: impl Fn(u32) -> bool,
) -> (Vec<u32>, usize) {
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

/// One column map for both halves of a metric's layout.
///
/// [`assign_columns`] runs in descending holder count, so when it is given
/// every split the metric can share, the widely held ones take the low
/// columns: `0..n_dense` are dense and `n_dense..n_columns` are posted. One
/// `u32` per distinct split describes both, rather than one map each.
pub(super) struct Layout {
    pub(super) column_of: Vec<u32>,
    pub(super) n_dense: usize,
    pub(super) n_columns: usize,
}

impl Layout {
    /// Columns for the splits `shareable` accepts, given how many of the
    /// `n_trees` trees hold each; those held by at least `min_dense` are dense.
    pub(super) fn new(
        counts: &[u32],
        n_trees: usize,
        min_dense: u32,
        shareable: impl Fn(u32) -> bool,
    ) -> Self {
        let (column_of, n_columns) = assign_columns(counts, n_trees, &shareable);
        let n_dense = counts
            .iter()
            .filter(|&&count| shareable(count) && count >= min_dense)
            .count();
        Self {
            column_of,
            n_dense,
            n_columns,
        }
    }
}

/// Trees laid out per block: a few rounds of the pool, so a block is laid out
/// in parallel and its scratch stays small.
pub(super) fn block_size() -> usize {
    4 * current_num_threads().max(1)
}

/// The trees a kernel reads: borrowed from a collection the caller keeps, or
/// owned, in which case each block of trees is dropped as soon as the kernel
/// has laid it out, so the trees and the kernel's layout of them are never both
/// whole.
pub(super) enum Trees<'a> {
    Borrowed(&'a [InternSnap]),
    Owned(Vec<InternSnap>),
}

impl Trees<'_> {
    pub(super) fn len(&self) -> usize {
        match self {
            Trees::Borrowed(trees) => trees.len(),
            Trees::Owned(trees) => trees.len(),
        }
    }

    /// Hand the trees to `lay_out` in order, a block of `size` at a time.
    pub(super) fn by_blocks(self, size: usize, mut lay_out: impl FnMut(&[InternSnap])) {
        match self {
            Trees::Borrowed(trees) => trees.chunks(size).for_each(lay_out),
            Trees::Owned(trees) => {
                let mut trees = trees.into_iter();
                loop {
                    let block: Vec<InternSnap> = trees.by_ref().take(size).collect();
                    if block.is_empty() {
                        break;
                    }
                    lay_out(&block);
                }
            }
        }
    }
}
