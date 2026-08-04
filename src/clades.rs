//! Clade (bipartition) frequencies and consensus-tree selection.
//!
//! Ports `treetracer/clade_freq/__init__.py` and
//! `consensus_tree.compute_consensus_tree_index`.
//!
//! Both work off the presence matrix that already falls out of building
//! `Snapshots`: `presence[i * n_bip + j] == 1` iff tree `i` contains
//! bipartition `j`. Nothing here re-parses a newick.

use crate::par::*;

/// Per-run frequency of every bipartition, plus what the bipartition contains.
pub struct CladeFrequencies {
    /// `freq[run * n_bip + j]` — share of run `run`'s trees containing split `j`.
    pub freq: Vec<f64>,
    pub n_runs: usize,
    pub n_bipartitions: usize,
    /// Number of taxa on the canonical side of each split — a cheap proxy for
    /// how deep in the tree it sits.
    pub clade_sizes: Vec<u32>,
}

/// Compute per-run bipartition frequencies from a presence matrix.
///
/// `presence` is row-major `n_trees x n_bip`; `run_ids[i]` names the run tree
/// `i` came from.
///
/// Frequencies are per run rather than pooled because that is the comparison
/// that matters: a clade at 0.98 in one chain and 0.11 in another is the
/// signature of two chains that have not converged on the same posterior,
/// and pooling would average that away into an unremarkable 0.55.
pub fn clade_frequencies(
    presence: &[u8],
    n_trees: usize,
    n_bip: usize,
    run_ids: &[u32],
    n_runs: usize,
) -> CladeFrequencies {
    let mut counts = vec![0u32; n_runs * n_bip];
    let mut run_totals = vec![0u32; n_runs];

    for i in 0..n_trees {
        let run = *run_ids.get(i).unwrap_or(&0) as usize;
        if run >= n_runs {
            continue;
        }
        run_totals[run] += 1;
        let row = &presence[i * n_bip..(i + 1) * n_bip];
        let target = &mut counts[run * n_bip..(run + 1) * n_bip];
        for (slot, &p) in target.iter_mut().zip(row) {
            *slot += p as u32;
        }
    }

    let mut freq = vec![0.0f64; n_runs * n_bip];
    for run in 0..n_runs {
        let total = run_totals[run].max(1) as f64;
        for j in 0..n_bip {
            freq[run * n_bip + j] = counts[run * n_bip + j] as f64 / total;
        }
    }

    CladeFrequencies {
        freq,
        n_runs,
        n_bipartitions: n_bip,
        clade_sizes: Vec::new(),
    }
}

/// How often each bipartition occurs within a subset of trees.
///
/// This is what a *saved* summary tree carries: the counts vector plus the
/// number of trees behind it. Two saved trees drawn from the same analysis
/// index into the same bipartition basis, so comparing them is column
/// alignment rather than key matching — which is why the desktop app requires
/// both sides of a comparison to come from one distance matrix.
pub struct CladeCounts {
    pub counts: Vec<u32>,
    pub n_trees: usize,
}

/// Count bipartition occurrences over an arbitrary subset of trees.
pub fn clade_counts_for(
    presence: &[u8],
    n_bip: usize,
    indices: &[usize],
    n_trees_total: usize,
) -> CladeCounts {
    let mut counts = vec![0u32; n_bip];
    let mut used = 0usize;
    for &i in indices {
        if i >= n_trees_total {
            continue;
        }
        used += 1;
        let row = &presence[i * n_bip..(i + 1) * n_bip];
        for (slot, &p) in counts.iter_mut().zip(row) {
            *slot += p as u32;
        }
    }
    CladeCounts {
        counts,
        n_trees: used,
    }
}

/// The maximum-clade-credibility tree of a set.
pub struct MccResult {
    /// Index of the winning tree, into the pooled tree order.
    pub index: usize,
    /// Its score: the sum over its splits of `log(P(split))`.
    pub log_clade_credibility: f64,
}

/// Pick the tree whose clades are, jointly, the most probable.
///
/// Score is `sum_j presence[i][j] * log(freq_j)` — the log product of the
/// posterior frequencies of the clades tree `i` contains. This is a *credible*
/// tree rather than a constructed consensus: it is one of the sampled trees,
/// so it is guaranteed to be a real, resolved topology with consistent branch
/// lengths, which a majority-rule consensus is not.
///
/// Splits absent from every tree are pinned to frequency 1 so their `log` is
/// 0; they are zero in every presence row anyway, so they cannot contribute.
pub fn max_clade_credibility(presence: &[u8], n_trees: usize, n_bip: usize) -> Option<MccResult> {
    if n_trees == 0 || n_bip == 0 || presence.len() != n_trees * n_bip {
        return None;
    }

    let mut counts = vec![0u32; n_bip];
    for i in 0..n_trees {
        let row = &presence[i * n_bip..(i + 1) * n_bip];
        for (slot, &p) in counts.iter_mut().zip(row) {
            *slot += p as u32;
        }
    }

    let log_freq: Vec<f64> = counts
        .iter()
        .map(|&c| {
            if c > 0 {
                (c as f64 / n_trees as f64).ln()
            } else {
                0.0
            }
        })
        .collect();

    let mut scores = vec![0.0f64; n_trees];
    scores.par_chunks_mut(1).enumerate().for_each(|(i, slot)| {
        let row = &presence[i * n_bip..(i + 1) * n_bip];
        let mut acc = 0.0;
        for (j, &p) in row.iter().enumerate() {
            if p != 0 {
                acc += log_freq[j];
            }
        }
        slot[0] = acc;
    });

    let mut best = 0usize;
    for i in 1..n_trees {
        if scores[i] > scores[best] {
            best = i;
        }
    }

    Some(MccResult {
        index: best,
        log_clade_credibility: scores[best],
    })
}

/// Maximum clade credibility restricted to a subset of trees.
///
/// Clade frequencies are computed *within the subset*, not inherited from the
/// full sample. That is the point of summarising a selection: a cluster in
/// tree space is its own little posterior, and the tree that best represents
/// it is the one whose clades are most probable among its neighbours — not the
/// one that scores best globally.
///
/// The returned index is into the original tree order, not into `indices`.
pub fn max_clade_credibility_for(
    presence: &[u8],
    n_trees_total: usize,
    n_bip: usize,
    indices: &[usize],
) -> Option<MccResult> {
    if indices.is_empty() || n_bip == 0 || presence.len() != n_trees_total * n_bip {
        return None;
    }

    let counts = clade_counts_for(presence, n_bip, indices, n_trees_total);
    if counts.n_trees == 0 {
        return None;
    }

    let log_freq: Vec<f64> = counts
        .counts
        .iter()
        .map(|&c| {
            if c > 0 {
                (c as f64 / counts.n_trees as f64).ln()
            } else {
                0.0
            }
        })
        .collect();

    let mut best = None::<(usize, f64)>;
    for &i in indices {
        if i >= n_trees_total {
            continue;
        }
        let row = &presence[i * n_bip..(i + 1) * n_bip];
        let mut acc = 0.0;
        for (j, &p) in row.iter().enumerate() {
            if p != 0 {
                acc += log_freq[j];
            }
        }
        if best.is_none_or(|(_, s)| acc > s) {
            best = Some((i, acc));
        }
    }

    best.map(|(index, score)| MccResult {
        index,
        log_clade_credibility: score,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frequencies_are_per_run() {
        // 4 trees, 3 splits. Run 0 = trees 0,1; run 1 = trees 2,3.
        let presence: Vec<u8> = vec![
            1, 0, 0, //
            1, 1, 0, //
            0, 1, 1, //
            0, 1, 1, //
        ];
        let run_ids = vec![0, 0, 1, 1];
        let f = clade_frequencies(&presence, 4, 3, &run_ids, 2);

        // Run 0: split 0 in both trees, split 1 in one, split 2 in none.
        assert_eq!(f.freq[0], 1.0);
        assert_eq!(f.freq[1], 0.5);
        assert_eq!(f.freq[2], 0.0);
        // Run 1: split 0 absent, splits 1 and 2 in both.
        assert_eq!(f.freq[3], 0.0);
        assert_eq!(f.freq[4], 1.0);
        assert_eq!(f.freq[5], 1.0);
    }

    #[test]
    fn mcc_picks_the_most_credible_tree() {
        // Split 0 is in 3 of 4 trees, split 1 in 1 of 4.
        // Tree 0 has the common split only — it should win.
        let presence: Vec<u8> = vec![
            1, 0, //
            1, 0, //
            1, 1, //
            0, 1, //
        ];
        let r = max_clade_credibility(&presence, 4, 2).expect("mcc");
        assert!(r.index == 0 || r.index == 1, "picked {}", r.index);
        // log(3/4) for the single common split.
        assert!((r.log_clade_credibility - (0.75f64).ln()).abs() < 1e-12);
    }

    #[test]
    fn subset_mcc_uses_subset_frequencies() {
        // Trees 0,1 share split 0. Trees 2,3 share split 1. Globally each
        // split sits at 2/4, so the full-set MCC is a coin toss; restricted to
        // {2,3} the answer must come from that pair.
        let presence: Vec<u8> = vec![
            1, 0, //
            1, 0, //
            0, 1, //
            0, 1, //
        ];
        let r = max_clade_credibility_for(&presence, 4, 2, &[2, 3]).expect("subset mcc");
        assert!(r.index == 2 || r.index == 3, "picked {}", r.index);
        // Split 1 is in both selected trees, so its frequency is 1 and the
        // score is log(1) = 0 — not log(2/4) as the global view would give.
        assert!(
            r.log_clade_credibility.abs() < 1e-12,
            "{}",
            r.log_clade_credibility
        );
    }

    #[test]
    fn subset_counts_only_count_the_subset() {
        let presence: Vec<u8> = vec![1, 0, 1, 1, 0, 1, 0, 0];
        let c = clade_counts_for(&presence, 2, &[0, 1], 4);
        assert_eq!(c.n_trees, 2);
        assert_eq!(c.counts, vec![2, 1]);
        // Out-of-range indices are skipped rather than panicking.
        let c2 = clade_counts_for(&presence, 2, &[0, 99], 4);
        assert_eq!(c2.n_trees, 1);
    }

    #[test]
    fn subset_mcc_rejects_empty_selection() {
        let presence: Vec<u8> = vec![1, 0, 0, 1];
        assert!(max_clade_credibility_for(&presence, 2, 2, &[]).is_none());
    }

    #[test]
    fn mcc_rejects_empty() {
        assert!(max_clade_credibility(&[], 0, 0).is_none());
        assert!(max_clade_credibility(&[1, 0], 1, 3).is_none());
    }
}
