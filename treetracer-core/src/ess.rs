//! Effective sample size for MCMC traces, and pseudo-ESS for tree chains.
//!
//! Ports `treetracer/ess/ess.py` and `treetracer/ess/pseudo_ess.py`.
//!
//! # ESS
//!
//! The Gelman variogram estimator on a split chain (BDA3 §11.5). A single
//! trace is cut in half and the halves treated as two pseudo-chains, which
//! gives a Gelman–Rubin marginal variance that notices drift between the early
//! and late parts of the chain. Autocorrelation then comes from the variogram
//!
//! ```text
//! V_t = <(x_{j,i+t} - x_{j,i})^2>
//! rho_t = 1 - V_t / (2 sigma^2)
//! ESS = m n / (1 + 2 sum_t rho_t)
//! ```
//!
//! summed until the first *even* lag where `rho_{t-1} + rho_t` goes negative —
//! Geyer's initial-positive-sequence rule written on consecutive lags.
//!
//! # Pseudo-ESS
//!
//! Tree topology has no univariate trace to run ESS on, but the RF distance
//! from every tree in the chain to some reference tree is a perfectly good
//! scalar (Lanfear et al. 2016). Different references give different answers,
//! so we take a sample of references, compute one ESS per reference, and
//! report the distribution. The conservative reading of "the chain's ESS" is
//! the minimum.
//!
//! Column `j` of the cached RF matrix *is* the trace for reference `j`, so
//! this costs no recomputation.

use statrs::distribution::{ContinuousCDF, Normal};

/// Rank-normalise a trace, as Stan does before computing ESS.
///
/// Ranks are averaged over ties, then mapped through the normal quantile
/// function at `(r - 3/8) / (n + 1/4)` — the Blom plotting position scipy's
/// `norm.ppf` is fed in the Python implementation.
///
/// This matters for RF traces specifically: they are bounded, discrete and
/// heavily tied, which is exactly the regime where a raw variogram ESS
/// misbehaves.
pub fn rank_normalize(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| x[a].partial_cmp(&x[b]).unwrap_or(std::cmp::Ordering::Equal));

    // Average ranks within tied groups, matching `scipy.stats.rankdata`'s
    // "average" method.
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && x[order[j + 1]] == x[order[i]] {
            j += 1;
        }
        // Ranks are 1-based; the group spans positions i..=j.
        let avg = ((i + j) as f64) / 2.0 + 1.0;
        for &idx in &order[i..=j] {
            ranks[idx] = avg;
        }
        i = j + 1;
    }

    let normal = Normal::new(0.0, 1.0).expect("standard normal");
    ranks
        .iter()
        .map(|&r| normal.inverse_cdf((r - 0.375) / (n as f64 + 0.25)))
        .collect()
}

/// Gelman–Rubin marginal posterior variance for `m` chains of `n` samples.
///
/// `s^2 = ((n-1)/n) W + B/n`, combining within- and between-chain variance.
fn split_chain_variance(chains: &[&[f64]]) -> f64 {
    let m = chains.len();
    let n = chains[0].len();
    if m < 2 || n < 2 {
        return 0.0;
    }

    let chain_means: Vec<f64> = chains
        .iter()
        .map(|c| c.iter().sum::<f64>() / n as f64)
        .collect();
    let grand_mean = chain_means.iter().sum::<f64>() / m as f64;

    let b_over_n = chain_means
        .iter()
        .map(|cm| (cm - grand_mean).powi(2))
        .sum::<f64>()
        / (m as f64 - 1.0);

    let w = chains
        .iter()
        .zip(&chain_means)
        .map(|(c, &cm)| c.iter().map(|v| (v - cm).powi(2)).sum::<f64>())
        .sum::<f64>()
        / (m as f64 * (n as f64 - 1.0));

    w * (n as f64 - 1.0) / n as f64 + b_over_n
}

/// Effective sample size of a 1-D trace. `NaN` when the trace is too short.
pub fn effective_sample_size(x: &[f64]) -> f64 {
    let n_total = x.len();
    if n_total < 4 {
        return f64::NAN;
    }

    // Even-length split; an odd final sample is dropped, the standard
    // split-Rhat convention.
    let half = n_total / 2;
    let first = &x[0..half];
    let second = &x[half..2 * half];
    let chains: [&[f64]; 2] = [first, second];
    let m = 2usize;
    let n = half;

    let post_var = split_chain_variance(&chains);
    if post_var == 0.0 {
        return (m * n) as f64;
    }

    let mut rho = vec![1.0f64; n];
    let mut t = 1usize;
    let mut negative = false;
    while !negative && t < n {
        let mut acc = 0.0;
        for chain in chains.iter() {
            for i in 0..(n - t) {
                let d = chain[i + t] - chain[i];
                acc += d * d;
            }
        }
        let v_t = acc / (m as f64 * (n - t) as f64);
        rho[t] = 1.0 - v_t / (2.0 * post_var);
        if t.is_multiple_of(2) && rho[t - 1] + rho[t] < 0.0 {
            negative = true;
        }
        t += 1;
    }

    let tau = 1.0 + 2.0 * rho[1..t].iter().sum::<f64>();
    if tau <= 0.0 {
        return (m * n) as f64;
    }
    (m * n) as f64 / tau
}

/// A pseudo-ESS result: one ESS estimate per reference tree, plus summaries.
pub struct PseudoEss {
    pub ess_values: Vec<f64>,
    pub ref_indices: Vec<usize>,
    pub min: f64,
    pub median: f64,
    pub max: f64,
}

/// Deterministic reference sampling.
///
/// Evenly spaced over the chain rather than pseudo-random. Two reasons: the
/// same input always yields the same diagnostic, which matters for a number
/// users will quote; and even spacing actually covers the chain better than a
/// random draw of the same size. Lanfear et al. sample at random, but nothing
/// in the method depends on the references being random — only on their being
/// a spread of distinct trees.
pub fn evenly_spaced_refs(n_trees: usize, n_refs: usize) -> Vec<usize> {
    let k = n_refs.min(n_trees);
    if k == 0 {
        return Vec::new();
    }
    (0..k).map(|i| (i * n_trees) / k).collect()
}

/// Pseudo-ESS over a set of reference trees.
///
/// `distmat` is the row-major `n x n` RF matrix. `refs` selects the reference
/// columns; pass [`evenly_spaced_refs`] unless you have a reason not to.
pub fn pseudo_ess(distmat: &[u32], n: usize, refs: &[usize]) -> PseudoEss {
    if n < 4 || refs.is_empty() || distmat.len() != n * n {
        return PseudoEss {
            ess_values: Vec::new(),
            ref_indices: Vec::new(),
            min: f64::NAN,
            median: f64::NAN,
            max: f64::NAN,
        };
    }

    let mut ess_values = Vec::with_capacity(refs.len());
    let mut trace = vec![0.0f64; n];
    for &r in refs {
        if r >= n {
            ess_values.push(f64::NAN);
            continue;
        }
        // Column r — the RF distance from every tree to reference r. The
        // matrix is symmetric so a row read would do equally well, but the
        // column is what the Python indexes and keeping them aligned makes
        // the cross-check exact.
        for (i, slot) in trace.iter_mut().enumerate() {
            *slot = distmat[i * n + r] as f64;
        }
        let normalized = rank_normalize(&trace);
        ess_values.push(effective_sample_size(&normalized));
    }

    let mut valid: Vec<f64> = ess_values.iter().copied().filter(|v| !v.is_nan()).collect();
    valid.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let (min, median, max) = if valid.is_empty() {
        (f64::NAN, f64::NAN, f64::NAN)
    } else {
        let mid = valid.len() / 2;
        let med = if valid.len().is_multiple_of(2) {
            (valid[mid - 1] + valid[mid]) / 2.0
        } else {
            valid[mid]
        };
        (valid[0], med, valid[valid.len() - 1])
    };

    PseudoEss {
        ess_values,
        ref_indices: refs.to_vec(),
        min,
        median,
        max,
    }
}

/// Pseudo-ESS restricted to a subset of trees.
///
/// The trace for a reference must be built from the subset alone. Taking the
/// full column and filtering afterwards would leave the reference's own zero
/// in place for a tree outside the subset, and — more importantly — an ESS is
/// a property of one chain's mixing. Pooling four runs into a single trace
/// measures how far apart the runs sit, not how well any of them mixed.
pub fn pseudo_ess_for(distmat: &[u32], n: usize, indices: &[usize], n_refs: usize) -> PseudoEss {
    let m = indices.len();
    if m < 4 || distmat.len() != n * n {
        return PseudoEss {
            ess_values: Vec::new(),
            ref_indices: Vec::new(),
            min: f64::NAN,
            median: f64::NAN,
            max: f64::NAN,
        };
    }

    // Compact the subset into its own dense matrix so the existing estimator
    // applies unchanged.
    let mut sub = vec![0u32; m * m];
    for (a, &ia) in indices.iter().enumerate() {
        for (b, &ib) in indices.iter().enumerate() {
            if ia < n && ib < n {
                sub[a * m + b] = distmat[ia * n + ib];
            }
        }
    }

    let refs = evenly_spaced_refs(m, n_refs.max(1));
    pseudo_ess(&sub, m, &refs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ess_of_independent_samples_is_near_n() {
        // Deterministic pseudo-random, uncorrelated.
        let mut state = 12345u64;
        let x: Vec<f64> = (0..2000)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((state >> 11) as f64) / ((1u64 << 53) as f64)
            })
            .collect();
        let ess = effective_sample_size(&x);
        // Uncorrelated: ESS should land in the same ballpark as N.
        assert!(ess > 1000.0 && ess <= 2200.0, "ess = {ess}");
    }

    #[test]
    fn ess_of_correlated_samples_is_much_less_than_n() {
        // An AR(1) chain with strong autocorrelation.
        let mut state = 999u64;
        let mut prev = 0.0f64;
        let x: Vec<f64> = (0..2000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let u = ((state >> 11) as f64) / ((1u64 << 53) as f64) - 0.5;
                prev = 0.95 * prev + u;
                prev
            })
            .collect();
        let ess = effective_sample_size(&x);
        assert!(ess < 400.0, "strongly correlated chain gave ess = {ess}");
        assert!(ess > 0.0);
    }

    #[test]
    fn subset_ess_uses_only_the_subset() {
        // Two tight blocks far apart. Pooled, the trace jumps between blocks
        // and looks badly mixed; within a block it is flat.
        let n = 40;
        let mut d = vec![0u32; n * n];
        for i in 0..n {
            for j in 0..n {
                let same_block = (i < 20) == (j < 20);
                d[i * n + j] = if i == j {
                    0
                } else if same_block {
                    1
                } else {
                    100
                };
            }
        }
        let block: Vec<usize> = (0..20).collect();
        let within = pseudo_ess_for(&d, n, &block, 10);
        let pooled = pseudo_ess(&d, n, &evenly_spaced_refs(n, 10));
        assert!(!within.median.is_nan());
        assert!(!pooled.median.is_nan());
        // The within-block chain is the better-behaved one.
        assert!(
            within.median >= pooled.median,
            "within {} vs pooled {}",
            within.median,
            pooled.median
        );
    }

    #[test]
    fn subset_ess_needs_enough_trees() {
        let d = vec![0u32; 16];
        assert!(pseudo_ess_for(&d, 4, &[0, 1], 5).median.is_nan());
    }

    #[test]
    fn short_traces_are_nan() {
        assert!(effective_sample_size(&[1.0, 2.0, 3.0]).is_nan());
    }

    #[test]
    fn constant_trace_returns_n() {
        let x = vec![7.0; 100];
        assert_eq!(effective_sample_size(&x), 100.0);
    }

    #[test]
    fn rank_normalize_is_monotone_and_centred() {
        let x = vec![5.0, 1.0, 3.0, 3.0, 9.0];
        let r = rank_normalize(&x);
        // Ties get equal values.
        assert!((r[2] - r[3]).abs() < 1e-12);
        // Order is preserved.
        assert!(r[1] < r[2] && r[2] < r[0] && r[0] < r[4]);

        // Symmetry only holds without ties: a tied pair collapses two
        // plotting positions onto their midpoint, which shifts the mean.
        let untied = vec![10.0, 20.0, 30.0, 40.0, 50.0];
        let ru = rank_normalize(&untied);
        let mean: f64 = ru.iter().sum::<f64>() / ru.len() as f64;
        assert!(mean.abs() < 1e-12, "mean = {mean}");
        // ...and it is antisymmetric end to end.
        assert!((ru[0] + ru[4]).abs() < 1e-12);
    }

    #[test]
    fn evenly_spaced_refs_spans_and_dedups() {
        let refs = evenly_spaced_refs(100, 10);
        assert_eq!(refs.len(), 10);
        assert_eq!(refs[0], 0);
        assert!(refs.windows(2).all(|w| w[0] < w[1]));
        assert!(*refs.last().unwrap() >= 90);
        // Asking for more references than trees is capped.
        assert_eq!(evenly_spaced_refs(5, 100).len(), 5);
    }
}
