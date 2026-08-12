//! Small statistics helpers for convergence diagnostics.
//!
//! Ports the `scipy.stats.gaussian_kde` calls in TreeTracer's
//! `callbacks/diagnostics.py`, which draw the marginal density beside each
//! trace plot. Keeping them here rather than reimplementing in TypeScript
//! means the desktop and browser builds share one definition of the curve.

/// Gaussian kernel density estimate on a regular grid.
///
/// Bandwidth follows Scott's rule exactly as `scipy.stats.gaussian_kde` does by
/// default: for 1-D data the factor is `n^(-1/5)`, applied to the sample
/// standard deviation. Matching scipy matters because these curves are
/// compared against the desktop app's output.
///
/// Returns `(grid, density)`, both of length `n_grid`. `None` if there is too
/// little data or no spread to estimate from.
///
/// The kernel is summed over every sample with no tail cutoff. Truncating at
/// 6 sigma looks free — `exp(-18)` is about 1.5e-8 — but it showed up as a
/// 1e-10 disagreement with scipy on a real RF trace, and the inner loop is
/// `n_grid × n` exponentials, so a few milliseconds at the sizes this sees.
/// Exactness is worth more than that here.
pub fn gaussian_kde(values: &[f64], n_grid: usize) -> Option<(Vec<f64>, Vec<f64>)> {
    let n = values.len();
    if n < 2 || n_grid == 0 {
        return None;
    }

    let mean = values.iter().sum::<f64>() / n as f64;
    // scipy's gaussian_kde uses the unbiased covariance (ddof = 1).
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
    let std = var.sqrt();
    if !std.is_finite() || std <= 0.0 {
        // Every value identical — a density here is a spike, not a curve.
        return None;
    }

    // Scott's factor for d = 1.
    let factor = (n as f64).powf(-1.0 / 5.0);
    let bandwidth = std * factor;

    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !min.is_finite() || !max.is_finite() {
        return None;
    }

    let step = if n_grid > 1 {
        (max - min) / (n_grid as f64 - 1.0)
    } else {
        0.0
    };

    let norm = 1.0 / (n as f64 * bandwidth * (2.0 * std::f64::consts::PI).sqrt());
    let inv_bw = 1.0 / bandwidth;

    let mut grid = Vec::with_capacity(n_grid);
    let mut density = Vec::with_capacity(n_grid);
    for i in 0..n_grid {
        let x = min + step * i as f64;
        let mut acc = 0.0;
        for &v in values {
            let z = (x - v) * inv_bw;
            acc += (-0.5 * z * z).exp();
        }
        grid.push(x);
        density.push(acc * norm);
    }

    Some((grid, density))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A standard normal sample should integrate to ~1 and peak near zero.
    #[test]
    fn integrates_to_one() {
        // Deterministic pseudo-normal sample via Box-Muller on a fixed LCG.
        let mut state = 42u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((state >> 11) as f64) / ((1u64 << 53) as f64)
        };
        let values: Vec<f64> = (0..2000)
            .map(|_| {
                let u1: f64 = next().max(1e-12);
                let u2: f64 = next();
                (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
            })
            .collect();

        let (grid, density) = gaussian_kde(&values, 512).expect("kde");
        let step = grid[1] - grid[0];
        let area: f64 = density.iter().sum::<f64>() * step;
        assert!((area - 1.0).abs() < 0.02, "area = {area}");

        // Peak should sit near the mean of a standard normal.
        let (peak_idx, _) = density
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap();
        assert!(grid[peak_idx].abs() < 0.35, "peak at {}", grid[peak_idx]);
    }

    #[test]
    fn rejects_degenerate_input() {
        assert!(gaussian_kde(&[], 100).is_none());
        assert!(gaussian_kde(&[1.0], 100).is_none());
        // No spread: not estimable.
        assert!(gaussian_kde(&[3.0, 3.0, 3.0, 3.0], 100).is_none());
    }

    #[test]
    fn density_is_non_negative_and_finite() {
        let values: Vec<f64> = (0..200).map(|i| (i % 17) as f64 * 1.5).collect();
        let (_, density) = gaussian_kde(&values, 200).expect("kde");
        assert!(density.iter().all(|d| d.is_finite() && *d >= 0.0));
    }
}
