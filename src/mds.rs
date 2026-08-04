//! Classical multidimensional scaling (PCoA) on a precomputed distance matrix.
//!
//! This is a Rust port of the `_pcoa_fast` path in TreeTracer's `rf/mds.py`,
//! which used vectorised double-centering plus `scipy.sparse.linalg.eigsh`.
//! Having it here rather than in either frontend means the desktop app and the
//! browser build share one implementation, and there is a single place to
//! validate against SciPy.
//!
//! Two properties carried over from the Python version:
//!
//! 1. **No materialised centering matrix.** The textbook `B = -½ H D² H` costs
//!    two O(n³) matrix products plus an n×n scratch buffer. The identity
//!    `B_ij = -½(D²_ij − rowmean_i − rowmean_j + grandmean)` is valid because
//!    D² is symmetric, and is pure O(n²) done in place.
//!
//! 2. **Top-k only.** A full eigendecomposition computes n eigenpairs so that
//!    n−k of them can be thrown away. Lanczos iteration finds the largest k
//!    directly from matrix-vector products.
//!
//! Everything stays in `f64`. At n≈8000 the row and grand-mean accumulators
//! over D² reach ~10¹⁴ for 1000-taxon RF distances; in `f32` that would lose
//! roughly seven significant digits and corrupt the centering subtraction.

use crate::par::*;

/// Result of a PCoA embedding.
pub struct Pcoa {
    /// Row-major `n × k` coordinates.
    pub coords: Vec<f64>,
    /// The k retained eigenvalues, descending, clamped at zero.
    pub eigenvalues: Vec<f64>,
    /// Fraction of total positive eigenvalue mass on each retained axis.
    pub explained: Vec<f64>,
    pub n: usize,
    pub k: usize,
}

/// Double-centre a squared-distance matrix in place, producing the Gram matrix.
fn gram_from_distances(distances: &[u32], n: usize) -> Vec<f64> {
    let mut b = vec![0.0f64; n * n];
    for (i, row) in b.chunks_mut(n).enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            let d = distances[i * n + j] as f64;
            *slot = d * d;
        }
    }

    let mut row_mean = vec![0.0f64; n];
    for (i, mean) in row_mean.iter_mut().enumerate() {
        let row = &b[i * n..(i + 1) * n];
        *mean = row.iter().sum::<f64>() / n as f64;
    }
    let grand_mean = row_mean.iter().sum::<f64>() / n as f64;

    let means = &row_mean;
    b.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        let ri = means[i];
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = -0.5 * (*slot - ri - means[j] + grand_mean);
        }
    });
    b
}

/// `out = a * x` for a dense symmetric `n × n` matrix in row-major order.
fn matvec(a: &[f64], x: &[f64], out: &mut [f64], n: usize) {
    out.par_chunks_mut(1).enumerate().for_each(|(i, slot)| {
        let row = &a[i * n..(i + 1) * n];
        let mut acc = 0.0;
        for (j, &v) in row.iter().enumerate() {
            acc += v * x[j];
        }
        slot[0] = acc;
    });
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Find a unit vector orthogonal to every vector in `basis`.
///
/// Used to restart Lanczos after it hits an invariant subspace. Trying the
/// standard unit vectors in turn is enough: if the basis does not already span
/// R^n, at least one `e_i` must have a component outside it.
fn fresh_orthogonal(basis: &[Vec<f64>], n: usize, tol: f64) -> Option<Vec<f64>> {
    if basis.len() >= n {
        return None;
    }
    for seed in 0..n {
        let mut candidate = vec![0.0f64; n];
        candidate[seed] = 1.0;
        for _ in 0..2 {
            for qi in basis {
                let proj = dot(&candidate, qi);
                for (cv, qv) in candidate.iter_mut().zip(qi) {
                    *cv -= proj * qv;
                }
            }
        }
        let nb = norm(&candidate);
        if nb > tol {
            let inv = 1.0 / nb;
            for x in candidate.iter_mut() {
                *x *= inv;
            }
            return Some(candidate);
        }
    }
    None
}

/// One Lanczos run: reduce `b` to a symmetric tridiagonal of size ≤ `m`.
///
/// Returns `(alpha, beta, basis, steps)` where `alpha[..steps]` is the
/// diagonal, `beta[..steps]` the subdiagonal (with `beta[steps-1]` the norm of
/// the final residual vector, which drives the convergence test), and `basis`
/// the orthonormal Krylov vectors.
fn lanczos(b: &[f64], n: usize, m: usize) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, usize) {
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut alpha = vec![0.0f64; m];
    let mut beta = vec![0.0f64; m];

    // Deterministic start vector: a fixed pattern rather than random, so the
    // same distance matrix always yields the same embedding. (The sign of each
    // axis is still arbitrary, as in any eigendecomposition.)
    let mut v: Vec<f64> = (0..n)
        .map(|i| 1.0 + ((i % 7) as f64) * 0.1 - ((i % 3) as f64) * 0.05)
        .collect();
    let v_norm = norm(&v);
    for x in v.iter_mut() {
        *x /= v_norm;
    }
    q.push(v);

    let mut w = vec![0.0f64; n];
    let tol = 1e-10;
    let mut steps = 0usize;

    for j in 0..m {
        if j >= q.len() {
            break; // basis could not be extended; Krylov space exhausted
        }
        matvec(b, &q[j], &mut w, n);
        let a = dot(&q[j], &w);
        alpha[j] = a;
        steps = j + 1;

        for (idx, qi) in q.iter().enumerate() {
            let coeff = if idx == j {
                a
            } else if idx + 1 == j {
                beta[j - 1]
            } else {
                0.0
            };
            if coeff != 0.0 {
                for (wv, qv) in w.iter_mut().zip(qi) {
                    *wv -= coeff * qv;
                }
            }
        }

        // Two passes of Gram-Schmidt against the whole basis. One pass is not
        // enough in floating point — "twice is enough" is the classic result.
        for _ in 0..2 {
            for qi in q.iter() {
                let proj = dot(&w, qi);
                for (wv, qv) in w.iter_mut().zip(qi) {
                    *wv -= proj * qv;
                }
            }
        }

        let nb = norm(&w);
        if j + 1 == m {
            beta[j] = nb;
            break;
        }

        if nb > tol {
            beta[j] = nb;
            let inv = 1.0 / nb;
            q.push(w.iter().map(|x| x * inv).collect());
        } else {
            // The Krylov space closed on an invariant subspace. Left here,
            // plain Lanczos would stop and report one eigenvector per
            // *distinct* eigenvalue — silently missing every repeated one,
            // because a single start vector cannot span a degenerate
            // eigenspace. Restarting on a fresh direction orthogonal to
            // everything found so far opens a new block of the tridiagonal
            // matrix; `beta[j] = 0` is what keeps the blocks independent.
            beta[j] = 0.0;
            match fresh_orthogonal(&q, n, tol) {
                Some(next) => q.push(next),
                None => break, // basis already spans the whole space
            }
        }
    }

    (alpha, beta, q, steps)
}

/// Compute a PCoA embedding from a dense square distance matrix.
///
/// `distances` is row-major `n × n`; `k` is the number of axes wanted. Returns
/// `None` for a degenerate input (n < 2 or k == 0).
///
/// # Accuracy
///
/// The Krylov subspace is grown until every retained eigenpair meets a
/// residual tolerance, so results track a full eigendecomposition to near
/// machine precision rather than to whatever a single fixed-size pass happens
/// to reach.
///
/// A single Lanczos pass sized like ARPACK's default `ncv` (`2k+20`) is *not*
/// enough on this workload. Tree-space Gram matrices have a flat spectrum —
/// on a 75-tree sample the top four eigenvalues sat within 12% of each other
/// — and closely spaced eigenvalues are exactly the case where an
/// unrestarted pass stalls, leaving coordinates about 1e-5 off. ARPACK hides
/// this behind implicit restarting; growing the subspace and rechecking is
/// the simpler equivalent, and terminates because at `m == n` the reduction
/// is exact.
pub fn pcoa(distances: &[u32], n: usize, k: usize) -> Option<Pcoa> {
    if n < 2 || k == 0 || distances.len() != n * n {
        return None;
    }
    let k = k.min(n - 1);
    let b = gram_from_distances(distances, n);

    // Residual tolerance, relative to the dominant eigenvalue.
    const RESID_TOL: f64 = 1e-12;

    let mut m = (2 * k + 20).min(n);
    loop {
        let (alpha, beta, q, steps) = lanczos(&b, n, m);

        // Lanczos has reduced the n×n problem to a `steps`×`steps` symmetric
        // tridiagonal one — a few dozen rows. At that size a *full* dense
        // eigendecomposition is free, so hand it to nalgebra rather than
        // carrying a hand-written QL implementation here.
        //
        // nalgebra rather than `ndarray-linalg` or `nalgebra-lapack`: those
        // bind to Fortran LAPACK and do not compile for wasm32 at all.
        // nalgebra's `SymmetricEigen` is pure Rust.
        let s = steps;
        let mut t = nalgebra::DMatrix::<f64>::zeros(s, s);
        for i in 0..s {
            t[(i, i)] = alpha[i];
            if i + 1 < s {
                t[(i, i + 1)] = beta[i];
                t[(i + 1, i)] = beta[i];
            }
        }
        let eig = nalgebra::SymmetricEigen::new(t);
        let diag: Vec<f64> = eig.eigenvalues.iter().copied().collect();
        let z = eig.eigenvectors;

        // Rank the Ritz values, largest first.
        let mut order: Vec<usize> = (0..s).collect();
        order.sort_by(|&a, &b2| {
            diag[b2]
                .partial_cmp(&diag[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let take = k.min(s);

        // Standard Lanczos residual bound: for Ritz pair (theta_j, Q s_j),
        // ||B y_j - theta_j y_j|| == |beta_last| * |s[last, j]|. Free — it
        // needs no extra matrix-vector product.
        let beta_last = beta[s - 1];
        let scale = order
            .first()
            .map(|&c| diag[c].abs())
            .unwrap_or(0.0)
            .max(1e-30);
        let max_resid = order
            .iter()
            .take(take)
            .map(|&col| (beta_last * z[(s - 1, col)]).abs())
            .fold(0.0f64, f64::max);

        let converged = max_resid <= RESID_TOL * scale;
        if !converged && m < n {
            m = (m * 2).min(n);
            continue;
        }

        // Explained variance is measured against the trace of the Gram matrix,
        // which equals the sum of *all* n eigenvalues and costs O(n) to read
        // off the diagonal.
        //
        // The obvious alternative — dividing by the eigenvalues this run
        // happened to find — makes the answer depend on the internal Krylov
        // size, so growing the subspace for accuracy would silently change
        // every reported percentage. Using the trace keeps the figure a
        // property of the data.
        //
        // Note RF distances are not Euclidean, so B has some negative
        // eigenvalues and the trace is smaller than the positive mass alone.
        // Fractions are therefore mildly optimistic; they are a guide to the
        // relative weight of each axis, not an exact variance decomposition.
        let total: f64 = (0..n).map(|i| b[i * n + i]).sum();
        let mut coords = vec![0.0f64; n * take];
        let mut eigenvalues = Vec::with_capacity(take);
        let mut explained = Vec::with_capacity(take);

        for (axis, &col) in order.iter().take(take).enumerate() {
            let lambda = diag[col].max(0.0);
            let scale_axis = lambda.sqrt();
            eigenvalues.push(lambda);
            explained.push(if total > 0.0 { lambda / total } else { 0.0 });

            // Ritz vector: the Lanczos basis combined by this eigenvector of T.
            for (i, coord_row) in coords.chunks_mut(take).enumerate() {
                let mut acc = 0.0;
                for (j, qj) in q.iter().enumerate().take(s) {
                    acc += qj[i] * z[(j, col)];
                }
                coord_row[axis] = acc * scale_axis;
            }
        }

        return Some(Pcoa {
            coords,
            eigenvalues,
            explained,
            n,
            k: take,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Four points on a line at 0, 1, 2, 3. PCoA should recover one dominant
    /// axis whose coordinates are the original spacing up to sign and offset.
    #[test]
    fn recovers_a_line() {
        let pts = [0u32, 1, 2, 3];
        let n = 4;
        let mut d = vec![0u32; n * n];
        for i in 0..n {
            for j in 0..n {
                d[i * n + j] = pts[i].abs_diff(pts[j]);
            }
        }
        let out = pcoa(&d, n, 2).expect("pcoa");
        let axis: Vec<f64> = (0..n).map(|i| out.coords[i * out.k]).collect();

        // Monotone along the line, and the second axis carries ~no variance.
        let ascending = axis.windows(2).all(|w| w[0] < w[1]);
        let descending = axis.windows(2).all(|w| w[0] > w[1]);
        assert!(ascending || descending, "axis 1 not monotone: {axis:?}");
        assert!(
            out.eigenvalues[1] < out.eigenvalues[0] * 1e-6,
            "second axis should be negligible: {:?}",
            out.eigenvalues
        );

        // Equal spacing preserved.
        let d01 = (axis[1] - axis[0]).abs();
        let d12 = (axis[2] - axis[1]).abs();
        assert!((d01 - d12).abs() < 1e-9, "spacing not uniform: {axis:?}");
    }

    /// Distances between the corners of a unit square, taken along the edges
    /// only (a 2-D configuration), should need exactly two axes.
    #[test]
    fn recovers_two_dimensions() {
        // Squared Euclidean layout of (0,0) (2,0) (0,2) (2,2), rounded to
        // integers so it fits the u32 input: edges 2, diagonals ~3.
        let n = 4;
        let d: Vec<u32> = vec![
            0, 2, 2, 3, //
            2, 0, 3, 2, //
            2, 3, 0, 2, //
            3, 2, 2, 0,
        ];
        let out = pcoa(&d, n, 3).expect("pcoa");
        assert!(out.eigenvalues[0] > 0.0);
        assert!(out.eigenvalues[1] > 0.0);
        // Two axes should carry essentially all the positive mass.
        assert!(
            out.explained[0] + out.explained[1] > 0.99,
            "explained: {:?}",
            out.explained
        );
    }

    #[test]
    fn rejects_degenerate_input() {
        assert!(pcoa(&[0], 1, 2).is_none());
        assert!(pcoa(&[0, 1, 1, 0], 2, 0).is_none());
    }
}
