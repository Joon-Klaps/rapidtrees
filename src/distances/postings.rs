//! Posting lists: the trees holding each rarely held split.

use super::layout::{Blocks, Column, Layout};
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
/// walks only those rather than scanning all of its tree's splits. The order
/// is set by [`PostingsBuilder::finish`], not by how a tree lists its splits.
pub(super) struct Postings {
    pub(super) offsets: Vec<usize>,
    pub(super) trees: Vec<u32>,
    pub(super) lengths: Vec<f64>,
    pub(super) own: Blocks<OwnPosted>,
}

/// One block of trees' own posted splits: tree `k` of the block owns entries
/// `starts[k]..starts[k + 1]` of `cols` and, with lengths, of `lengths`.
pub(super) struct OwnPosted {
    pub(super) starts: Vec<usize>,
    pub(super) cols: Vec<u32>,
    pub(super) lengths: Vec<f64>,
}

/// [`Postings`] while the trees are added to them, a block at a time.
pub(super) struct PostingsBuilder {
    pub(super) postings: Postings,
    pub(super) with_lengths: bool,
    /// The next free entry of each list.
    pub(super) fill: Vec<usize>,
    /// How many trees have been added so far.
    pub(super) added: usize,
}

impl Postings {
    /// Empty lists for `layout`'s posted columns, carrying branch lengths when
    /// `with_lengths`, to be filled with `n` trees in blocks of `1 << shift`. A
    /// posted split's list holds exactly the trees that hold it, so `counts`
    /// sizes every list before any tree is read.
    pub(super) fn builder(
        counts: &[u32],
        layout: &Layout,
        n: usize,
        shift: u32,
        with_lengths: bool,
    ) -> PostingsBuilder {
        let n_posted = layout.n_columns - layout.n_dense;
        let mut offsets = vec![0usize; n_posted + 1];
        if n_posted > 0 {
            for (id, &count) in counts.iter().enumerate() {
                if let Some(Column::Posted(col)) = layout.place(id as u32) {
                    offsets[col as usize] = count as usize;
                }
            }
        }
        let total = exclusive_prefix_sum(&mut offsets);
        PostingsBuilder {
            fill: offsets[..n_posted].to_vec(),
            postings: Self {
                offsets,
                trees: vec![0; total],
                lengths: vec![0.0; if with_lengths { total } else { 0 }],
                own: Blocks::with_capacity(shift, n),
            },
            with_lengths,
            added: 0,
        }
    }

    /// Tree `i`'s posted columns in column order, with their branch lengths
    /// (empty when built without them).
    #[inline]
    pub(super) fn of_tree(&self, i: usize) -> (&[u32], &[f64]) {
        let (own, k) = self.own.locate(i);
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

    /// Add `overlap(lengthᵢ, lengthⱼ)` into `row[j]` for each posted split tree
    /// `i` holds and each later tree `j` holding it too, in tree `i`'s column
    /// order. Reads the branch lengths, so the lists must be built with them.
    #[inline]
    pub(super) fn add_shared(&self, i: usize, row: &mut [f64], overlap: &impl Fn(f64, f64) -> f64) {
        let (cols, lengths) = self.of_tree(i);
        for (&col, &length) in cols.iter().zip(lengths) {
            let (trees, others) = self.after(col, i);
            for (&j, &other) in trees.iter().zip(others) {
                row[j as usize] += overlap(length, other);
            }
        }
    }
}

impl PostingsBuilder {
    /// Add the next block of trees, each tree's posted splits in any order.
    /// Blocks come in order, so every list stays ascending. A tree's own
    /// posted splits are only counted here; [`Self::finish`] fills them in.
    pub(super) fn push_block(&mut self, posted: &[Vec<(u32, f64)>]) {
        let mut starts = Vec::with_capacity(posted.len() + 1);
        starts.push(0);
        for (tree, splits) in (self.added..).zip(posted) {
            for &(col, length) in splits {
                let at = &mut self.fill[col as usize];
                self.postings.trees[*at] = tree as u32;
                if self.with_lengths {
                    self.postings.lengths[*at] = length;
                }
                *at += 1;
            }
            starts.push(starts[starts.len() - 1] + splits.len());
        }
        let total = starts[posted.len()];
        self.added += posted.len();
        self.postings.own.push(OwnPosted {
            starts,
            cols: vec![0; total],
            lengths: vec![0.0; if self.with_lengths { total } else { 0 }],
        });
    }

    /// The filled lists, with each tree's own posted splits in column order.
    ///
    /// Every list holds its trees in ascending order, so walking the lists in
    /// column order and handing each entry to its tree is a transpose: each
    /// tree receives its columns in ascending order, in time linear in the
    /// entries and without a sort.
    pub(super) fn finish(self) -> Postings {
        let mut postings = self.postings;
        let shift = postings.own.shift;
        // Where each tree's next own posted split goes in its block.
        let mut next: Vec<usize> = (postings.own.blocks.iter())
            .flat_map(|own| own.starts[..own.starts.len() - 1].iter().copied())
            .collect();
        for (col, span) in postings.offsets.windows(2).enumerate() {
            for at in span[0]..span[1] {
                let tree = postings.trees[at] as usize;
                let own = &mut postings.own.blocks[tree >> shift];
                let slot = &mut next[tree];
                own.cols[*slot] = col as u32;
                if self.with_lengths {
                    own.lengths[*slot] = postings.lengths[at];
                }
                *slot += 1;
            }
        }
        postings
    }
}
