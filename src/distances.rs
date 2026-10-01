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
//! be worth a list there (`RF_DENSE_SHARE`) than in the weighted metrics
//! (`WEIGHTED_DENSE_SHARE`).
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
//! costs one. `U1` and `U2` fold into their tree's `self` term. RF lays its
//! columns out the same way with its own cutoff, and also drops `D1`, since a
//! split held by every tree cancels.
//!
//! There is no per-pair entry point. To compare two trees, build a two-tree
//! `Snapshots` and read the off-diagonal cell.

use crate::par::*;
use crate::snapshot::{InternSnap, Snapshots};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Share of the trees a split must be held by to keep its bit column in RF.
/// Rarer splits get a posting list instead.
const RF_DENSE_SHARE: f64 = 0.03;
const WEIGHTED_DENSE_SHARE: f64 = 0.25;

/// Fewest trees at which RF gives any split a posting list.
const RF_MIN_TREES_FOR_POSTINGS: usize = 256;

/// Marks a split that [`assign_columns`] gave no column.
const NO_COLUMN: u32 = u32::MAX;

// ─── dense/posting boundary ─────────────────────────────────────────────────

/// `default`, unless `value` is a finite, non-negative number.
///
fn parse_share(value: Option<&str>, default: f64) -> f64 {
    value
        .and_then(|v| v.trim().parse::<f64>().ok())
        .filter(|share| share.is_finite() && *share >= 0.0)
        .unwrap_or(default)
}

/// [`RF_DENSE_SHARE`], or `RAPIDTREES_RF_DENSE_SHARE` when that is set. Read
/// once per process. The variable exists to re-tune the boundary on other data
/// or hardware; it changes the speed and never a distance.
fn rf_dense_share() -> f64 {
    static SHARE: OnceLock<f64> = OnceLock::new();
    *SHARE.get_or_init(|| {
        parse_share(
            std::env::var("RAPIDTREES_RF_DENSE_SHARE").ok().as_deref(),
            RF_DENSE_SHARE,
        )
    })
}

/// [`WEIGHTED_DENSE_SHARE`], or `RAPIDTREES_WEIGHTED_DENSE_SHARE` when that is
/// set. See [`rf_dense_share`].
fn weighted_dense_share() -> f64 {
    static SHARE: OnceLock<f64> = OnceLock::new();
    *SHARE.get_or_init(|| {
        parse_share(
            std::env::var("RAPIDTREES_WEIGHTED_DENSE_SHARE")
                .ok()
                .as_deref(),
            WEIGHTED_DENSE_SHARE,
        )
    })
}

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
    fill_symmetric_banded(n, 1, progress, fill_row)
}

/// [`fill_symmetric`] with one rayon task per band of `band` rows:
/// `fill_band(i0, rows)` writes `rows[r * n + j]` for every row `i0 + r` of the
/// band and every `j > i0 + r`.
fn fill_symmetric_banded<T, F>(
    n: usize,
    band: usize,
    progress: Option<&AtomicUsize>,
    fill_band: F,
) -> Vec<T>
where
    T: Copy + Default + Send,
    F: Fn(usize, &mut [T]) + Sync,
{
    let mut matrix = vec![T::default(); n * n];

    matrix
        .par_chunks_mut((n * band).max(1))
        .enumerate()
        .for_each(|(b, rows)| {
            let i0 = b * band;
            fill_band(i0, rows);
            if let Some(counter) = progress {
                let pairs = (i0..i0 + rows.len() / n).map(|i| n - i - 1).sum();
                counter.fetch_add(pairs, Ordering::Relaxed);
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

/// The fewest trees that must hold a split for it to keep a dense column:
/// `share` of the `n_trees`, and never fewer than two, since a split held by
/// one tree is never shared.
fn min_dense(n_trees: usize, share: f64) -> u32 {
    ((n_trees as f64 * share).ceil() as u32).max(2)
}

/// One column map for both halves of a metric's layout.
///
/// [`assign_columns`] runs in descending holder count, so when it is given
/// every split the metric can share, the widely held ones take the low
/// columns: `0..n_dense` are dense and `n_dense..n_columns` are posted. One
/// `u32` per distinct split describes both, rather than one map each.
struct Layout {
    column_of: Vec<u32>,
    n_dense: usize,
    n_columns: usize,
}

impl Layout {
    /// Columns for the splits `shareable` accepts, given how many of the
    /// `n_trees` trees hold each; those held by at least `min_dense` are dense.
    fn new(
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
// ─── the trees a kernel reads ───────────────────────────────────────────────

/// Trees laid out per block: a few rounds of the pool, so a block is laid out
/// in parallel and its scratch stays small.
fn block_size() -> usize {
    4 * current_num_threads().max(1)
}

/// The trees a kernel reads: borrowed from a collection the caller keeps, or
/// owned, in which case each block of trees is dropped as soon as the kernel
/// has laid it out, so the trees and the kernel's layout of them are never both
/// whole.
enum Trees<'a> {
    Borrowed(&'a [InternSnap]),
    Owned(Vec<InternSnap>),
}

impl Trees<'_> {
    fn len(&self) -> usize {
        match self {
            Trees::Borrowed(trees) => trees.len(),
            Trees::Owned(trees) => trees.len(),
        }
    }

    /// Hand the trees to `lay_out` in order, a block of `size` at a time.
    fn by_blocks(self, size: usize, mut lay_out: impl FnMut(&[InternSnap])) {
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

// ─── posting lists ──────────────────────────────────────────────────────────

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
struct Postings {
    offsets: Vec<usize>,
    trees: Vec<u32>,
    lengths: Vec<f64>,
    /// Trees per block of `own`: tree `i` is tree `i % block` of `own[i / block]`.
    block: usize,
    own: Vec<OwnPosted>,
    /// While the lists are filled: the next free entry of each.
    fill: Vec<usize>,
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
    fn sized(counts: &[u32], layout: &Layout, block: usize, with_lengths: bool) -> Self {
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
    fn push_block(&mut self, posted: &[Vec<(u32, f64)>]) {
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
    fn of_tree(&self, i: usize) -> (&[u32], &[f64]) {
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
fn distance_rf_split(
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

// ─── weighted metrics (WRF, KF) ─────────────────────────────────────────────

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

/// Rows per weighted sweep task, rows of the other side held against them,
/// and columns per pass: a tile of `TILE_I × TILE_J` pairs reads each row's
/// panel once instead of once per pair, which turns the sweep from bound by
/// memory bandwidth into bound by arithmetic.
const TILE_I: usize = 8;
const TILE_J: usize = 4;
const PANEL: usize = 512;

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
struct DenseRows {
    size: usize,
    stride: usize,
    blocks: Vec<Vec<f64>>,
}

impl DenseRows {
    /// Tree `i`'s row.
    #[inline]
    fn row(&self, i: usize) -> &[f64] {
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
fn weighted_layout(
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
fn weighted_distances_split(
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
    use super::{
        Layout, NO_COLUMN, PANEL, TREEDIST_TREES, Trees, assign_columns, distance_rf_split,
        fill_symmetric, min_dense, parse_share, sweep, weighted_distances_split, weighted_layout,
    };
    use crate::{
        distances::{RF_DENSE_SHARE, WEIGHTED_DENSE_SHARE, rf_dense_share, weighted_dense_share},
        snapshot::{InternSnap, Snapshots},
    };
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
            for id in snap.ids() {
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
            assert_eq!(snap.n_splits(), snap.lengths.len());
        }
    }

    /// Test oracle: all three metrics for one pair, straight from the
    /// definitions over the union of both trees' splits (absent = length 0).
    ///
    /// Deliberately naive — no filtering, no `selfᵢ + selfⱼ − 2·shared` rewrite —
    /// so it cannot reproduce a backend bug by construction.
    fn reference_distances(a: &InternSnap, b: &InternSnap) -> (usize, f64, f64) {
        let lengths_by_split = |s: &InternSnap| -> BTreeMap<u32, f64> {
            s.ids().zip(s.lengths.iter().copied()).collect()
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

    /// Assert an RF matrix agrees exactly with [`reference_distances`].
    fn assert_rf_matches_reference(snaps: &Snapshots, rf: &[u32], ctx: &str) {
        let n = snaps.snapshots.len();
        for i in 0..n {
            for j in 0..n {
                let (want, _, _) = reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);
                assert_eq!(rf[i * n + j] as usize, want, "RF at [{i}][{j}], {ctx}");
            }
        }
    }

    /// Assert all three metrics agree with [`reference_distances`] on the
    /// snapshots built from `newicks`.
    fn assert_matches_reference<S: AsRef<str>>(snaps: &Snapshots, newicks: &[S], ctx: &str) {
        assert_rf_matches_reference(snaps, &snaps.pairwise_rf(None), ctx);
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
    /// [`super::distance_wrf`] and [`super::distance_kf`]. The trees are read
    /// borrowed and owned, and the two have to agree to the bit.
    fn weighted_at(snaps: &Snapshots, min_dense: u32) -> (Vec<f64>, Vec<f64>) {
        let counts = &snaps.split_counts;
        let both = |trees: fn(&Snapshots) -> Trees| {
            (
                weighted_distances_split(trees(snaps), counts, None, min_dense, f64::min, |d| d),
                weighted_distances_split(
                    trees(snaps),
                    counts,
                    None,
                    min_dense,
                    |a, b| a * b,
                    f64::sqrt,
                ),
            )
        };
        let borrowed = both(|snaps| Trees::Borrowed(&snaps.snapshots));
        let owned = both(|snaps| Trees::Owned(snaps.snapshots.clone()));
        assert_eq!(
            borrowed, owned,
            "borrowed and owned trees disagree at min_dense={min_dense}"
        );
        borrowed
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

    #[test]
    fn share_override_accepts_only_finite_non_negative_numbers() {
        assert_eq!(parse_share(None, 0.25), 0.25);
        assert_eq!(parse_share(Some(" 0.1 "), 0.25), 0.1);
        assert_eq!(parse_share(Some("0"), 0.25), 0.0);
        assert_eq!(parse_share(Some("2"), 0.25), 2.0);
        for bad in ["", "quarter", "-0.1", "NaN", "inf"] {
            assert_eq!(parse_share(Some(bad), 0.25), 0.25, "{bad:?}");
        }
    }

    /// Each share is its default unless the environment overrides it. The test
    /// sets no variable: the shares are read once per process, and tests share
    /// that process, so a test that set one would leak into the others.
    #[test]
    fn dense_shares_follow_the_environment() {
        for (share, var, default) in [
            (
                rf_dense_share(),
                "RAPIDTREES_RF_DENSE_SHARE",
                RF_DENSE_SHARE,
            ),
            (
                weighted_dense_share(),
                "RAPIDTREES_WEIGHTED_DENSE_SHARE",
                WEIGHTED_DENSE_SHARE,
            ),
        ] {
            let expected = parse_share(std::env::var(var).ok().as_deref(), default);
            assert_eq!(share, expected, "{var}");
        }
    }

    /// The tiled sweep adds each pair's terms in the order [`sweep`] does, so it
    /// gives the untiled formula to the bit. Seventeen trees leave partial
    /// bands and tiles, and with 603 taxa the dense columns cross a panel
    /// boundary and end in a tail.
    #[test]
    fn tiled_sweep_matches_the_untiled_formula_to_the_bit() {
        let (_, snaps) = random_snapshots(603, 13, 4, 51);
        let (counts, n) = (&snaps.split_counts, snaps.snapshots.len());
        let min_dense = min_dense(n, weighted_dense_share());
        type Overlap = fn(f64, f64) -> f64;
        type Finish = fn(f64) -> f64;
        let metrics: [(Overlap, Finish); 2] = [(f64::min, |d| d), (|a, b| a * b, f64::sqrt)];
        for (overlap, finish) in metrics {
            let trees = || Trees::Borrowed(&snaps.snapshots);
            let tiled = weighted_distances_split(trees(), counts, None, min_dense, overlap, finish);

            let layout = Layout::new(counts, n, min_dense, |count| count >= 2);
            let (dense, postings, self_total) = weighted_layout(trees(), counts, &layout, &overlap);
            assert!(
                dense.stride > PANEL && dense.stride % 8 != 0,
                "{}",
                dense.stride
            );
            let untiled = fill_symmetric(n, None, |i, row: &mut [f64]| {
                let (cols, lengths) = postings.of_tree(i);
                for (&col, &length) in cols.iter().zip(lengths) {
                    let (trees, others) = postings.after(col, i);
                    for (&j, &other) in trees.iter().zip(others) {
                        row[j as usize] += overlap(length, other);
                    }
                }
                for (j, slot) in row.iter_mut().enumerate().skip(i + 1) {
                    let shared = sweep(dense.row(i), dense.row(j), &overlap) + *slot;
                    *slot = finish((self_total[i] + self_total[j] - 2.0 * shared).max(0.0));
                }
            });
            assert_eq!(tiled, untiled);
        }
    }

    #[test]
    fn no_trees_give_an_empty_matrix() {
        let snaps = Snapshots::from_newicks(&[], false).unwrap();
        assert!(snaps.pairwise_rf(None).is_empty());
        assert!(snaps.pairwise_wrf(None).is_empty());
        assert!(snaps.pairwise_kf(None).is_empty());
    }

    /// Handing a collection to a kernel gives the matrix that borrowing it
    /// gives, for every metric at the default boundaries: over several layout
    /// blocks, and with enough trees for RF to give rare splits posting lists.
    #[test]
    fn consuming_kernels_match_borrowing_ones() {
        for &(n_taxa, n_trees, duplicates, seed) in
            &[(20usize, 12usize, 4usize, 41u64), (30, 300, 40, 42)]
        {
            let build = || random_snapshots(n_taxa, n_trees, duplicates, seed).1;
            let snaps = build();
            assert_eq!(
                snaps.pairwise_rf(None),
                build().into_pairwise_rf(None),
                "RF, seed={seed}"
            );
            assert_eq!(
                snaps.pairwise_wrf(None),
                build().into_pairwise_wrf(None),
                "WRF, seed={seed}"
            );
            assert_eq!(
                snaps.pairwise_kf(None),
                build().into_pairwise_kf(None),
                "KF, seed={seed}"
            );
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
            let (_, snaps) = random_snapshots(n_taxa, n_trees, duplicates, seed);
            let n = snaps.snapshots.len();
            for min_dense in [2, 3, n as u32 / 2, n as u32, n as u32 + 1] {
                let counts = &snaps.split_counts;
                let rf =
                    distance_rf_split(Trees::Borrowed(&snaps.snapshots), counts, None, min_dense);
                let owned = distance_rf_split(
                    Trees::Owned(snaps.snapshots.clone()),
                    counts,
                    None,
                    min_dense,
                );
                assert_eq!(
                    rf, owned,
                    "borrowed and owned trees disagree at min_dense={min_dense}"
                );
                assert_rf_matches_reference(
                    &snaps,
                    &rf,
                    &format!("seed={seed} min_dense={min_dense}"),
                );
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
