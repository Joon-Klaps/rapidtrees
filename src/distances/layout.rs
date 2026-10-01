//! Column layout of the splits, and the trees a kernel reads.

use super::matrix::{exclusive_prefix_sum, row_slice};
use super::postings::Postings;
use crate::par::*;
use crate::snapshot::{InternSnap, Snapshots};
use std::borrow::Cow;

/// Marks a split that [`assign_columns`] gave no column.
pub(super) const NO_COLUMN: u32 = u32::MAX;

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

/// Where a [`Layout`] puts a split.
pub(super) enum Column {
    /// A dense column: a bit in RF's rows, an `f64` in the weighted metrics'.
    Dense(usize),
    /// A posting list, numbered from 0 after the dense columns.
    Posted(u32),
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

    /// Where split `id` goes, or `None` if it has no column.
    #[inline]
    pub(super) fn place(&self, id: u32) -> Option<Column> {
        match self.column_of[id as usize] {
            NO_COLUMN => None,
            col if (col as usize) < self.n_dense => Some(Column::Dense(col as usize)),
            col => Some(Column::Posted(col - self.n_dense as u32)),
        }
    }
}

/// What a kernel reads: the trees, and how many of them hold each split.
///
/// Borrowed from a collection the caller keeps, or owned when the caller is
/// done with it. Owned trees are dropped a block at a time as the kernel lays
/// them out, so the trees and the kernel's layout of them are never both whole.
pub(crate) struct Input<'a> {
    pub(super) trees: Cow<'a, [InternSnap]>,
    pub(super) counts: Cow<'a, [u32]>,
}

impl<'a> From<&'a Snapshots> for Input<'a> {
    fn from(snaps: &'a Snapshots) -> Self {
        Self {
            trees: Cow::Borrowed(&snaps.snapshots),
            counts: Cow::Borrowed(&snaps.split_counts),
        }
    }
}

impl From<Snapshots> for Input<'static> {
    fn from(snaps: Snapshots) -> Self {
        Self {
            trees: Cow::Owned(snaps.snapshots),
            counts: Cow::Owned(snaps.split_counts),
        }
    }
}

/// Trees per layout block, as a power of two: two to four rounds of the pool,
/// so a block is laid out in parallel and its scratch stays small, and finding
/// a tree's block is a shift rather than a division.
pub(super) fn block_shift() -> u32 {
    (4 * current_num_threads().max(1)).ilog2()
}

/// Hand `trees` to `lay_out` in order, `size` at a time. Owned trees are
/// dropped a block at a time, as soon as `lay_out` has read them.
pub(super) fn by_blocks(
    trees: Cow<'_, [InternSnap]>,
    size: usize,
    mut lay_out: impl FnMut(&[InternSnap]),
) {
    match trees {
        Cow::Borrowed(trees) => trees.chunks(size).for_each(lay_out),
        Cow::Owned(trees) => {
            let mut trees = trees.into_iter().peekable();
            while trees.peek().is_some() {
                let block: Vec<InternSnap> = trees.by_ref().take(size).collect();
                lay_out(&block);
            }
        }
    }
}

/// Per-tree data stored a block of trees per item. Each block is sized as it
/// is laid out, so nothing is reserved for all the trees while they are still
/// held. Every block but the last holds `1 << shift` trees.
pub(super) struct Blocks<T> {
    pub(super) shift: u32,
    pub(super) blocks: Vec<T>,
}

impl<T> Blocks<T> {
    /// No blocks yet, with room for those of `n` trees.
    pub(super) fn with_capacity(shift: u32, n: usize) -> Self {
        Self {
            shift,
            blocks: Vec::with_capacity(n.div_ceil(1 << shift)),
        }
    }

    /// Add the next block.
    pub(super) fn push(&mut self, block: T) {
        self.blocks.push(block);
    }

    /// The block holding tree `i`, and where in it `i` is.
    #[inline]
    pub(super) fn locate(&self, i: usize) -> (&T, usize) {
        (&self.blocks[i >> self.shift], i & ((1 << self.shift) - 1))
    }
}

/// One row of `stride` values per tree over a layout's dense columns.
pub(super) struct Rows<T> {
    pub(super) stride: usize,
    pub(super) blocks: Blocks<Vec<T>>,
}

impl<T> Rows<T> {
    /// Tree `i`'s row.
    #[inline]
    pub(super) fn row(&self, i: usize) -> &[T] {
        let (block, k) = self.blocks.locate(i);
        row_slice(block, k, self.stride)
    }

    /// Every row in one flat, row-major buffer.
    pub(super) fn into_flat(self) -> Vec<T>
    where
        T: Clone,
    {
        self.blocks.blocks.concat()
    }
}

/// The trees as a kernel sweeps them: each tree's dense row, the posting
/// lists, and each tree's `self` term.
pub(super) struct Laid<T, S> {
    pub(super) rows: Rows<T>,
    pub(super) postings: Postings,
    pub(super) selfs: Vec<S>,
}

/// Lay the trees out over `layout` a block at a time, `width` values to a row.
///
/// `tree(layout, snap, row)` writes `snap`'s dense columns into its zeroed
/// `row` and returns its posted splits in column order, with their branch
/// lengths when `with_lengths`, and its `self` term. A block's trees are laid
/// out in parallel and then added to the posting lists in tree order, which
/// keeps every list ascending without a sort or atomics.
pub(super) fn lay_out<T, S>(
    input: Input<'_>,
    layout: Layout,
    width: usize,
    with_lengths: bool,
    tree: impl Fn(&Layout, &InternSnap, &mut [T]) -> (Vec<(u32, f64)>, S) + Sync,
) -> Laid<T, S>
where
    T: Copy + Default + Send,
    S: Send,
{
    let Input { trees, counts } = input;
    let n = trees.len();
    let shift = block_shift();
    let mut rows = Rows {
        stride: width,
        blocks: Blocks::with_capacity(shift, n),
    };
    let mut postings = Postings::builder(&counts, &layout, n, shift, with_lengths);
    let mut selfs = Vec::with_capacity(n);
    // `par_chunks_mut` takes no zero width, so a zero-width row is handed one
    // scratch value, which is dropped once its block is laid out.
    let scratch = width.max(1);
    by_blocks(trees, 1 << shift, |block| {
        let mut buffer = vec![T::default(); block.len() * scratch];
        let (posted, block_selfs): (Vec<_>, Vec<_>) = buffer
            .par_chunks_mut(scratch)
            .zip(block.par_iter())
            .map(|(row, snap)| tree(&layout, snap, &mut row[..width]))
            .unzip();
        postings.push_block(&posted);
        selfs.extend(block_selfs);
        rows.blocks
            .push(if width == 0 { Vec::new() } else { buffer });
    });
    Laid {
        rows,
        postings: postings.finish(),
        selfs,
    }
}
