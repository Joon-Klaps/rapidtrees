//! GPU buffer layouts and WGSL kernel sources — with no GPU dependency.
//!
//! This crate never talks to a GPU. `wgpu` on `wasm32-unknown-unknown` links
//! `wasm-bindgen`, and a browser WebGPU path needs an async device request and
//! an async buffer readback, neither of which fits behind the synchronous
//! `Snapshots` API the CLI and Python bindings use. So the device, the
//! pipelines and the dispatch live in the consumer — `treetracer-core` for the
//! browser — and what lives here is the half that is pure data: how to lay the
//! snapshots out in memory so a shader can read them, and the shader text
//! itself.
//!
//! Keeping the layout here rather than in the consumer is what stops the bit
//! packing from being reimplemented against the public API and drifting out of
//! agreement with [`crate::distances`]. The tests at the bottom of this file
//! check exactly that: they replay each shader's arithmetic on the CPU over
//! these buffers and require the result to match `pairwise_*` — so a layout bug
//! fails `cargo test` without needing a GPU.
//!
//! # Sizing
//!
//! WebGPU's default `maxStorageBufferBindingSize` is 128 MiB, and every buffer
//! below has to clear it. The `storage_bytes` methods report the input side;
//! the output matrix is a separate `n² × 4` bytes, which is what binds first —
//! it passes 128 MiB at about 5 800 trees, above which the consumer has to
//! dispatch in row blocks.

use crate::distances::{assign_columns, shared_length_rows, split_tree_counts};
use crate::par::*;
use crate::snapshot::Snapshots;

/// WGSL source for the Robinson-Foulds kernel. Reads [`SplitBitRows`].
pub const RF_WGSL: &str = include_str!("../shaders/rf.wgsl");

/// WGSL source for the weighted Robinson-Foulds kernel. Reads [`WeightedRows`]
/// built with [`WeightedMetric::Wrf`].
pub const WRF_WGSL: &str = include_str!("../shaders/wrf.wgsl");

/// WGSL source for the Kuhner-Felsenstein kernel. Reads [`WeightedRows`] built
/// with [`WeightedMetric::Kf`].
pub const KF_WGSL: &str = include_str!("../shaders/kf.wgsl");

/// Bit-packed split presence, one row per tree, laid out for [`RF_WGSL`].
///
/// A set bit means "this tree holds that split". Only splits that *some but not
/// all* trees hold get a bit: a split held by every tree adds equally to both
/// `kept_per_tree` terms and to the shared popcount, so it cancels out of RF
/// exactly and is dropped before packing.
///
/// Words are `u32` because WGSL has no 64-bit integer type. The CPU path in
/// [`crate::distances`] packs the same bits into `u64` instead.
#[derive(Debug, Clone)]
pub struct SplitBitRows {
    /// Flat `n × words` bit matrix, row-major.
    pub packed: Vec<u32>,
    /// `u32` words per row.
    pub words: usize,
    /// Per tree, how many *kept* splits it holds — the `a` term in
    /// `RF(i, j) = aᵢ + aⱼ − 2·shared`.
    pub kept_per_tree: Vec<u32>,
    /// Number of trees.
    pub n: usize,
}

impl SplitBitRows {
    /// Bytes the `packed` buffer occupies on the device.
    #[inline]
    pub fn storage_bytes(&self) -> usize {
        self.packed.len() * size_of::<u32>()
    }
}

/// Which weighted metric a [`WeightedRows`] was built for.
///
/// The variants differ only in the `term` applied to the branch lengths that
/// land in `unique_self`, but that term has to match the shader, so building
/// rows for one metric and running the other's kernel over them is wrong.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WeightedMetric {
    /// Weighted Robinson-Foulds: `term(l) = l`, paired with [`WRF_WGSL`].
    Wrf,
    /// Kuhner-Felsenstein: `term(l) = l²`, paired with [`KF_WGSL`].
    Kf,
}

/// Dense `f32` branch lengths over shared-candidate splits, for [`WRF_WGSL`]
/// and [`KF_WGSL`].
///
/// Only splits at least two trees hold get a column. A split held by one tree
/// alone can never overlap, so instead of costing a near-empty column it is
/// folded into that tree's `unique_self`. On diverse posteriors this cuts
/// ~48 000 columns to ~2 600; on near-identical trees it is close to a no-op.
///
/// Lengths are downcast from the `f64` the CPU path carries, because WGSL has
/// no `f64`. Expect ~1e-5 relative error against [`Snapshots::pairwise_wrf`]
/// and [`Snapshots::pairwise_kf`] — GPU results are not bit-identical to CPU
/// ones and should not be compared for exact equality.
#[derive(Debug, Clone)]
pub struct WeightedRows {
    /// Flat `n × stride` length matrix, row-major; `0.0` where a tree lacks the
    /// split.
    pub rows: Vec<f32>,
    /// Row stride — how many split columns survived the filter.
    pub stride: usize,
    /// Per tree, `Σ term(length)` over the splits no other tree holds. The
    /// kernel adds `unique_self[i] + unique_self[j]` to every cell.
    pub unique_self: Vec<f32>,
    /// Number of trees.
    pub n: usize,
    /// Which metric these rows were built for.
    pub metric: WeightedMetric,
}

impl WeightedRows {
    /// Bytes the `rows` buffer occupies on the device.
    #[inline]
    pub fn storage_bytes(&self) -> usize {
        self.rows.len() * size_of::<f32>()
    }
}

/// Pack split presence into `u32` bit rows for the RF kernel.
///
/// Mirrors the packing inside [`crate::distances`]'s RF backend, but 32 bits to
/// a word instead of 64. See [`SplitBitRows`] for the layout and for why
/// universal splits are dropped.
pub fn split_bit_rows(snaps: &Snapshots) -> SplitBitRows {
    let n = snaps.snapshots.len();
    if n == 0 {
        return SplitBitRows {
            packed: Vec::new(),
            words: 0,
            kept_per_tree: Vec::new(),
            n: 0,
        };
    }

    let counts = split_tree_counts(snaps);
    let n_splits = counts.len();
    let (bit_slot, kept) = assign_columns(&counts, |count| count < n as u32);
    let everywhere = n_splits - kept;
    let words = kept.div_ceil(32);

    let mut packed = vec![0u32; n * words];
    if words > 0 {
        packed
            .par_chunks_mut(words)
            .zip(&snaps.snapshots)
            .for_each(|(row, snap)| {
                for &id in &snap.split_ids {
                    let slot = bit_slot[id as usize];
                    if slot != u32::MAX {
                        let slot = slot as usize;
                        row[slot >> 5] |= 1u32 << (slot & 31);
                    }
                }
            });
    }

    // Every tree holds all the universal splits, so its kept count is just its
    // total minus that fixed number.
    let kept_per_tree = snaps
        .snapshots
        .iter()
        .map(|snap| (snap.split_ids.len() - everywhere) as u32)
        .collect();

    SplitBitRows {
        packed,
        words,
        kept_per_tree,
        n,
    }
}

/// Build dense `f32` length rows for the WRF or KF kernel.
///
/// Delegates the column filtering and the `unique_self` fold to the same
/// [`shared_length_rows`] the CPU path uses, then downcasts to `f32`, so the
/// two backends cannot disagree about which splits earn a column.
pub fn weighted_rows(snaps: &Snapshots, metric: WeightedMetric) -> WeightedRows {
    let n = snaps.snapshots.len();
    if n == 0 {
        return WeightedRows {
            rows: Vec::new(),
            stride: 0,
            unique_self: Vec::new(),
            n: 0,
            metric,
        };
    }

    // Must match `distance_wrf` / `distance_kf` in `distances.rs`.
    let shared = match metric {
        WeightedMetric::Wrf => shared_length_rows(snaps, |l| l),
        WeightedMetric::Kf => shared_length_rows(snaps, |l| l * l),
    };

    WeightedRows {
        rows: shared.rows.iter().map(|&l| l as f32).collect(),
        stride: shared.stride,
        unique_self: shared.unique_self.iter().map(|&s| s as f32).collect(),
        n,
        metric,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distances::TREEDIST_TREES;
    use crate::distances::tests::{DIVERSE_TREES, random_newick};

    /// Replay [`RF_WGSL`] on the CPU: one cell per (i, j), exactly as the
    /// kernel computes it. Integer throughout, so this must match bit for bit.
    fn replay_rf(l: &SplitBitRows) -> Vec<u32> {
        let n = l.n;
        let mut matrix = vec![0u32; n * n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let shared: u32 = (0..l.words)
                    .map(|k| (l.packed[i * l.words + k] & l.packed[j * l.words + k]).count_ones())
                    .sum();
                matrix[i * n + j] = l.kept_per_tree[i] + l.kept_per_tree[j] - 2 * shared;
            }
        }
        matrix
    }

    /// Replay [`WRF_WGSL`] / [`KF_WGSL`] on the CPU, including the
    /// `unique[i] + unique[j]` correction the shared-column layout depends on.
    fn replay_weighted(l: &WeightedRows) -> Vec<f32> {
        let n = l.n;
        let mut matrix = vec![0f32; n * n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let (row_i, row_j) = (i * l.stride, j * l.stride);
                let acc: f32 = (0..l.stride)
                    .map(|k| {
                        let d = l.rows[row_i + k] - l.rows[row_j + k];
                        match l.metric {
                            WeightedMetric::Wrf => d.abs(),
                            WeightedMetric::Kf => d * d,
                        }
                    })
                    .sum();
                let total = acc + l.unique_self[i] + l.unique_self[j];
                matrix[i * n + j] = match l.metric {
                    WeightedMetric::Wrf => total,
                    WeightedMetric::Kf => total.sqrt(),
                };
            }
        }
        matrix
    }

    /// The GPU buffers must reproduce every `pairwise_*` value on `snaps`.
    ///
    /// This is the test that makes the layout trustworthy without a GPU: if a
    /// column is dropped that should not have been, or `unique_self` is folded
    /// with the wrong `term`, the replayed cell diverges from the CPU answer.
    fn assert_layout_matches_cpu(snaps: &Snapshots, ctx: &str) {
        let n = snaps.snapshots.len();

        let rf_cpu = snaps.pairwise_rf(None);
        let rf_gpu = replay_rf(&split_bit_rows(snaps));
        for (idx, (&want, &got)) in rf_cpu.iter().zip(&rf_gpu).enumerate() {
            assert_eq!(
                got as usize,
                want,
                "RF at [{}][{}], {ctx}",
                idx / n,
                idx % n
            );
        }

        for (metric, cpu) in [
            (WeightedMetric::Wrf, snaps.pairwise_wrf(None)),
            (WeightedMetric::Kf, snaps.pairwise_kf(None)),
        ] {
            let gpu = replay_weighted(&weighted_rows(snaps, metric));
            for (idx, (&want, &got)) in cpu.iter().zip(&gpu).enumerate() {
                // f32 downcast, so this is a tolerance check, not equality.
                assert!(
                    (got as f64 - want).abs() <= 1e-4 * want.max(1.0),
                    "{metric:?} {got} vs CPU {want} at [{}][{}], {ctx}",
                    idx / n,
                    idx % n,
                );
            }
        }
    }

    #[test]
    fn layout_matches_cpu_on_treedist_fixtures() {
        let snaps = Snapshots::from_newicks(&TREEDIST_TREES, false).unwrap();
        assert_layout_matches_cpu(&snaps, "treedist fixtures");
    }

    #[test]
    fn layout_matches_cpu_on_diverse_trees() {
        // Most splits are held by a single tree here, so this is the case that
        // exercises the `unique_self` fold rather than the shared columns.
        let snaps = Snapshots::from_newicks(&DIVERSE_TREES, false).unwrap();
        assert_layout_matches_cpu(&snaps, "diverse fixtures");
    }

    #[test]
    fn layout_matches_cpu_on_random_trees() {
        let mut state = 0x5EED_1234_u64;
        for n_taxa in [4usize, 8, 16, 32] {
            for duplicates in [1usize, 3] {
                let mut newicks = Vec::new();
                for _ in 0..6 {
                    let t = random_newick(n_taxa, &mut state);
                    for _ in 0..duplicates {
                        newicks.push(t.clone());
                    }
                }
                let refs: Vec<&str> = newicks.iter().map(String::as_str).collect();
                let snaps = Snapshots::from_newicks(&refs, false).unwrap();
                assert_layout_matches_cpu(&snaps, &format!("{n_taxa} taxa, ×{duplicates}"));
            }
        }
    }

    #[test]
    fn identical_trees_give_zero_width_rows_and_zero_distance() {
        // Every split is universal, so RF drops all of them and `words` is 0.
        // The kernel must still produce 0 rather than indexing an empty buffer.
        let t = TREEDIST_TREES[0];
        let snaps = Snapshots::from_newicks(&[t, t, t], false).unwrap();

        let bits = split_bit_rows(&snaps);
        assert_eq!(
            bits.words, 0,
            "universal splits should leave no bit columns"
        );
        assert!(bits.packed.is_empty());
        assert!(replay_rf(&bits).iter().all(|&d| d == 0));

        for metric in [WeightedMetric::Wrf, WeightedMetric::Kf] {
            let rows = weighted_rows(&snaps, metric);
            assert!(
                rows.unique_self.iter().all(|&s| s == 0.0),
                "{metric:?}: no split is unique when all trees are identical"
            );
            assert!(replay_weighted(&rows).iter().all(|&d| d == 0.0));
        }
    }

    #[test]
    fn empty_snapshots_produce_empty_buffers() {
        let snaps = Snapshots::from_newicks(&[] as &[&str], false).unwrap();
        let bits = split_bit_rows(&snaps);
        assert_eq!((bits.n, bits.words, bits.storage_bytes()), (0, 0, 0));

        for metric in [WeightedMetric::Wrf, WeightedMetric::Kf] {
            let rows = weighted_rows(&snaps, metric);
            assert_eq!((rows.n, rows.stride, rows.storage_bytes()), (0, 0, 0));
        }
    }

    #[test]
    fn weighted_metric_is_recorded_on_the_rows() {
        // The rows differ only in the `term` applied to `unique_self`, so the
        // tag is the only thing stopping a caller running the wrong kernel.
        let snaps = Snapshots::from_newicks(&DIVERSE_TREES, false).unwrap();
        for metric in [WeightedMetric::Wrf, WeightedMetric::Kf] {
            assert_eq!(weighted_rows(&snaps, metric).metric, metric);
        }
    }

    #[test]
    fn shader_sources_are_present_and_declare_an_entry_point() {
        // Guards the `include_str!` paths: a typo there is a compile error, but
        // pointing at the wrong file is not.
        for (name, src) in [("rf", RF_WGSL), ("wrf", WRF_WGSL), ("kf", KF_WGSL)] {
            assert!(src.contains("fn main("), "{name}.wgsl has no entry point");
            assert!(
                src.contains("@compute @workgroup_size(16, 16, 1)"),
                "{name}.wgsl is not the 2-D compute kernel the dispatch assumes"
            );
        }
    }
}
