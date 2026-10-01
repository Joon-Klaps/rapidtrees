//! Weighted Robinson–Foulds and Kuhner–Felsenstein entry points.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator};
use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;

use super::input::{IterInput, collect_snapshots_from_iter, validate_iter_args};
use super::progress::{ProgressCounter, with_counter};
use super::{PyRfSnapshotResult, n_pairs};
use crate::snapshot::{Retain, Snapshots};

/// A pairwise weighted metric, as a plain function pointer so WRF and KF can
/// share one body.
type WeightedMetric = fn(&Snapshots, Option<&AtomicUsize>) -> Vec<f64>;

/// A weighted metric that consumes the collection: for entry points that
/// return the matrix alone.
type ConsumingMetric = fn(Snapshots, Option<&AtomicUsize>) -> Vec<f64>;

/// Shared body of the plain WRF and KF entry points.
///
/// These paths export nothing and the dense kernels read only `split_ids` +
/// `lengths`, so the bipartition table is never built in the first place, and
/// the kernel is handed the collection to drop as it goes.
fn weighted_pairwise(
    py: Python<'_>,
    input: IterInput<'_, '_>,
    progress: Option<Py<ProgressCounter>>,
    metric: ConsumingMetric,
) -> PyResult<Vec<f64>> {
    let snaps = collect_snapshots_from_iter(
        input,
        Retain {
            lengths: true,
            bipartitions: false,
        },
    )?;

    let total = n_pairs(snaps.len());
    with_counter(py, progress, total, move |counter| {
        metric(snaps, Some(counter))
    })
}

/// Shared body of the WRF and KF snapshot-exporting entry points.
///
/// The branch-length matrix is metric-agnostic, so the two differ only in which
/// distance fills the first buffer.
fn weighted_pairwise_with_snapshots(
    py: Python<'_>,
    names: Vec<String>,
    input: IterInput<'_, '_>,
    progress: Option<Py<ProgressCounter>>,
    metric: WeightedMetric,
) -> PyResult<PyRfSnapshotResult> {
    let snaps = collect_snapshots_from_iter(input, Retain::everything())?;

    let total = n_pairs(snaps.len());
    let matrix = with_counter(py, progress, total, |counter| metric(&snaps, Some(counter)))?;
    let matrix_bytes: Vec<u8> = matrix.iter().flat_map(|&v| v.to_ne_bytes()).collect();

    let n_bip = snaps.n_distinct_splits();
    let leaf_names = snaps.leaf_names.clone();
    let (bl_bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
    let bip_bytes = snaps.build_bipartition_bytes(&col_to_bip_id);

    Ok((
        names,
        PyBytes::new(py, &matrix_bytes).into(),
        leaf_names,
        n_bip,
        PyBytes::new(py, &bl_bytes).into(),
        PyBytes::new(py, &bip_bytes).into(),
    ))
}

/// Compute pairwise WRF distances and export a branch-length matrix in a single pass.
///
/// Parses each newick once, building both the pairwise WRF distance matrix and a
/// per-edge branch-length matrix encoding the branch length of every edge in each tree.
///
/// # Branch-length matrix format
///
/// Shape `(n_trees, n_bip)`, encoded as a flat row-major `float64` byte buffer.
/// `branch_length_bytes[i, j]` is the branch length of edge `j` in tree `i`, or `0.0`
/// if that edge is absent. Pendant (leaf) edges are always present in every tree.
/// Column order matches `bipartition_clade_bytes`, and is stable across calls.
///
/// Reconstruct and compute Fréchet traces on the Python side:
/// ```python
/// bl = np.frombuffer(branch_length_bytes, dtype=np.float64).reshape(n_trees, n_bip)
/// wrf_trace = np.sum(np.abs(bl[ref_idx, :] - bl), axis=1)   # wRF to one reference
/// kf_trace  = np.sqrt(np.sum((bl[ref_idx, :] - bl)**2, axis=1))  # KF to one reference
/// ```
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///     progress: Optional ``rapidtrees.ProgressCounter`` whose ``.value()`` /
///         ``.total()`` / ``.fraction()`` reflect live progress. Typically
///         polled from another Python thread. See
///         ``pairwise_rf_from_newick_iter`` for details.
///
/// Returns:
///     6-tuple (tree_names, wrf_matrix_bytes, leaf_names, n_bip,
///     branch_length_bytes, bipartition_clade_bytes).
///     wrf_matrix_bytes is flat float64 bytes (row-major n×n).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
pub(super) fn pairwise_wrf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;
    let input = IterInput {
        newick_iter,
        translate_maps: &translate_maps,
        map_indices: &map_indices,
        rooted,
    };
    weighted_pairwise_with_snapshots(py, names, input, progress, Snapshots::pairwise_wrf)
}

/// Compute pairwise KF distances and export a branch-length matrix in a single pass.
///
/// Identical to `pairwise_wrf_with_snapshots_from_newick_iter` except the distance
/// matrix uses the Kuhner–Felsenstein (Branch Score) metric instead of Weighted RF.
/// The branch-length matrix is metric-agnostic and identical in both functions.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///     progress: Optional ``rapidtrees.ProgressCounter`` whose ``.value()`` /
///         ``.total()`` / ``.fraction()`` reflect live progress. Typically
///         polled from another Python thread. See
///         ``pairwise_rf_from_newick_iter`` for details.
///
/// Returns:
///     6-tuple (tree_names, kf_matrix_bytes, leaf_names, n_bip,
///     branch_length_bytes, bipartition_clade_bytes).
///     kf_matrix_bytes is flat float64 bytes (row-major n×n).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
pub(super) fn pairwise_kf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;
    let input = IterInput {
        newick_iter,
        translate_maps: &translate_maps,
        map_indices: &map_indices,
        rooted,
    };
    weighted_pairwise_with_snapshots(py, names, input, progress, Snapshots::pairwise_kf)
}

/// Compute pairwise Weighted Robinson-Foulds distances from a lazy Python iterator.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///     progress: Optional ``rapidtrees.ProgressCounter`` whose ``.value()`` /
///         ``.total()`` / ``.fraction()`` reflect live progress. Typically
///         polled from another Python thread. See
///         ``pairwise_rf_from_newick_iter`` for details.
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of Weighted RF distances.
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
pub(super) fn pairwise_wrf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
) -> PyResult<(Vec<String>, Vec<f64>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;
    let input = IterInput {
        newick_iter,
        translate_maps: &translate_maps,
        map_indices: &map_indices,
        rooted,
    };
    let matrix = weighted_pairwise(py, input, progress, Snapshots::into_pairwise_wrf)?;
    Ok((names, matrix))
}

/// Compute pairwise Kuhner-Felsenstein (Branch Score) distances from a lazy Python iterator.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///     progress: Optional ``rapidtrees.ProgressCounter`` whose ``.value()`` /
///         ``.total()`` / ``.fraction()`` reflect live progress. Typically
///         polled from another Python thread. See
///         ``pairwise_rf_from_newick_iter`` for details.
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of KF distances.
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
pub(super) fn pairwise_kf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
) -> PyResult<(Vec<String>, Vec<f64>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;
    let input = IterInput {
        newick_iter,
        translate_maps: &translate_maps,
        map_indices: &map_indices,
        rooted,
    };
    let matrix = weighted_pairwise(py, input, progress, Snapshots::into_pairwise_kf)?;
    Ok((names, matrix))
}
