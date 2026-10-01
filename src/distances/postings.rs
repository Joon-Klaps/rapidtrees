//! Posting lists: the trees holding each rarely held split.

use super::layout::{Layout, column};
use super::matrix::exclusive_prefix_sum;

/// The trees holding each rarely held split: HashRF's bucket.
///
/// Posted column `c` (numbered from 0 after the dense ones) owns entries
/// `offsets[c]..offsets[c + 1]` of `trees`, ascending because the trees are
/// added in order, and of `lengths` when they were asked for. A row walks its
/// own posted splits and adds one term per later tree holding each, so a split
/// held by `k` trees costs `k²/2` additions in all, however many trees the
/// collection has.
///
/// Each tree's own posted splits are recorded once, in column order, so a row
/// walks only those rather than scanning all of its tree's splits. They are
/// stored a block of trees at a time, sized as each block is laid out, so they
/// are never reserved while the trees they come from are still held.
pub(super) struct Postings {
    offsets: Vec<usize>,
    trees: Vec<u32>,
    lengths: Vec<f64>,
    /// Trees per block of `own`: tree `i` is tree `i % block` of `own[i / block]`.
    block: usize,
    own: Vec<OwnPosted>,
    /// While the lists are filled: the next free entry of each.
    pub(super) fill: Vec<usize>,
}

/// One block of trees' own posted splits: tree `k` of the block owns entries
/// `starts[k]..starts[k + 1]` of `cols` and, with lengths, of `lengths`.
struct OwnPosted {
    starts: Vec<usize>,
    cols: Vec<u32>,
    lengths: Vec<f64>,
}

impl Postings {
    /// Empty lists for `layout`'s posted columns, carrying branch lengths when
    /// `with_lengths`, for trees that arrive `block` at a time. A posted
    /// split's list holds exactly the trees that hold it, so `counts` sizes
    /// every list before any tree is read. When there are none, nothing is
    /// scanned or allocated beyond an empty index.
    pub(super) fn sized(counts: &[u32], layout: &Layout, block: usize, with_lengths: bool) -> Self {
        let n_posted = layout.n_columns - layout.n_dense;
        let mut offsets = vec![0usize; n_posted + 1];
        if n_posted > 0 {
            for (id, &count) in counts.iter().enumerate() {
                let posted = column(&layout.column_of, id as u32)
                    .and_then(|col| col.checked_sub(layout.n_dense));
                if let Some(col) = posted {
                    offsets[col] = count as usize;
                }
            }
        }
        let total = exclusive_prefix_sum(&mut offsets);
        Self {
            fill: offsets[..n_posted].to_vec(),
            offsets,
            trees: vec![0; total],
            lengths: vec![0.0; if with_lengths { total } else { 0 }],
            block,
            own: Vec::new(),
        }
    }

    /// Add the next block of trees' posted splits, each tree's in column
    /// order. Blocks come in order, so every list stays ascending.
    pub(super) fn push_block(&mut self, posted: &[Vec<(u32, f64)>]) {
        let first = self.own.len() * self.block;
        let total = posted.iter().map(Vec::len).sum();
        let with_lengths = !self.lengths.is_empty();
        let mut own = OwnPosted {
            starts: Vec::with_capacity(posted.len() + 1),
            cols: Vec::with_capacity(total),
            lengths: Vec::with_capacity(if with_lengths { total } else { 0 }),
        };
        own.starts.push(0);
        for (k, splits) in posted.iter().enumerate() {
            for &(col, length) in splits {
                let at = &mut self.fill[col as usize];
                self.trees[*at] = (first + k) as u32;
                if let Some(slot) = self.lengths.get_mut(*at) {
                    *slot = length;
                }
                *at += 1;
            }
            own.cols.extend(splits.iter().map(|&(col, _)| col));
            if with_lengths {
                own.lengths.extend(splits.iter().map(|&(_, length)| length));
            }
            own.starts.push(own.cols.len());
        }
        self.own.push(own);
    }

    /// Tree `i`'s posted columns in column order, with their branch lengths
    /// (empty when built without them).
    #[inline]
    pub(super) fn of_tree(&self, i: usize) -> (&[u32], &[f64]) {
        let own = &self.own[i / self.block];
        let k = i % self.block;
        let span = own.starts[k]..own.starts[k + 1];
        let lengths = if own.lengths.is_empty() {
            &[][..]
        } else {
            &own.lengths[span.clone()]
        };
        (&own.cols[span], lengths)
    }

    /// The trees after `i` that hold posted column `col`, with their branch
    /// lengths (empty when built without them). Only those are wanted: the
    /// lower triangle is mirrored.
    #[inline]
    pub(super) fn after(&self, col: u32, i: usize) -> (&[u32], &[f64]) {
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
