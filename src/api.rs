//! Python binding layer for tree distance calculations.
//!
//! All computation lives in [`crate::distances`] and [`crate::snapshot`]; this module
//! contains only the PyO3-specific glue: iterating Python iterators, wrapping
//! byte buffers into `PyBytes`, and registering functions in the Python module.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator};
use std::collections::HashMap;

use crate::progress::{ProgressCounter, with_counter};
use crate::snapshot::Snapshots;

/// Number of upper-triangle pairs for `n` trees — what a [`ProgressCounter`]
/// counts up to as "100% done".
#[inline]
fn n_pairs(n: usize) -> usize {
    n.saturating_mul(n.saturating_sub(1)) / 2
}

type PyRfSnapshotResult = (
    Vec<String>,
    Py<PyAny>,
    Vec<String>,
    usize,
    Py<PyAny>,
    Py<PyAny>,
);

/// Collect tree snapshots from a lazy Python iterator of newick strings.
///
/// `store_lengths` controls whether branch lengths are retained: weighted metrics
/// (WRF/KF) need them, RF-only paths pass `false` to skip the per-tree `lengths`
/// allocation.
fn collect_snapshots_from_iter(
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: &[HashMap<String, String>],
    map_indices: &[usize],
    rooted: bool,
    store_lengths: bool,
) -> PyResult<Snapshots> {
    let newicks: Vec<String> = newick_iter
        .map(|item| item?.extract::<String>())
        .collect::<PyResult<_>>()?;

    if newicks.len() < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }

    let entries = newicks
        .iter()
        .zip(map_indices.iter())
        .map(|(n, &idx)| (n.as_str(), &translate_maps[idx]));
    Snapshots::from_newick_iter_opts(entries, rooted, store_lengths).map_err(PyValueError::new_err)
}

/// Validate argument consistency for iterator-based functions.
fn validate_iter_args(
    names: &[String],
    map_indices: &[usize],
    translate_maps: &[HashMap<String, String>],
) -> PyResult<()> {
    let n_names = names.len();

    if n_names < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }

    if n_names != map_indices.len() {
        return Err(PyValueError::new_err(format!(
            "names length ({}) must equal map_indices length ({})",
            n_names,
            map_indices.len()
        )));
    }

    for (i, &idx) in map_indices.iter().enumerate() {
        if idx >= translate_maps.len() {
            return Err(PyValueError::new_err(format!(
                "map_indices[{}] = {} is out of bounds (only {} translate maps provided)",
                i,
                idx,
                translate_maps.len()
            )));
        }
    }
    Ok(())
}

/// Compute pairwise Robinson-Foulds distances from a lazy Python iterator of newick strings.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///     progress: Optional ``rapidtrees.ProgressCounter``. Instantiate one,
///         hand it in, then read ``.value()`` / ``.total()`` / ``.fraction()``
///         from any other Python thread to observe live progress. While the
///         counter is in use the GIL is released for the duration of the
///         pairwise loop, so a polling thread runs unimpeded. When ``None``
///         (the default) the call has zero progress-related overhead and the
///         GIL stays held throughout. Note: when ``use_gpu=True`` and a GPU is
///         used, progress stays at 0 until the GPU finishes (no per-row updates).
///     use_gpu: If ``True``, attempt to accelerate computation on a GPU
///         (autodetect). RF is exact on GPU; WRF/KF run in f32 (~1e-5 relative
///         precision). Falls back to CPU silently when no suitable GPU is found
///         or the problem is below the minimum size threshold. Defaults to
///         ``False`` (CPU, full f64 precision).
///
/// Returns:
///     (tree_names, rf_matrix_bytes) — rf_matrix_bytes is flat u32 bytes (row-major n×n).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None, use_gpu=false))]
fn pairwise_rf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
    use_gpu: bool,
) -> PyResult<(Vec<String>, Py<PyAny>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let mut snaps =
        collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted, false)?;

    // This path never exports snapshots, and the dense RF computation reads
    // only `split_ids` — release the bipartition bitsets before the O(n²) loop.
    snaps.bipartitions = Vec::new();

    let n = snaps.snapshots.len();
    let rf_matrix = with_counter(py, progress, n_pairs(n), |counter| {
        crate::distances::dispatch_rf(&snaps, Some(counter), use_gpu)
    })?;
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(n)
        .flat_map(|row| row.iter().flat_map(|&v| (v as u32).to_ne_bytes()))
        .collect();

    let py_rf = PyBytes::new(py, &rf_bytes);

    Ok((names, py_rf.into()))
}

/// Compute pairwise RF distances and export binary tree snapshots in a single pass.
///
/// Parses each newick once, building both the RF distance matrix and a binary
/// presence matrix encoding which bipartitions appear in each tree.
///
/// # Presence matrix format
///
/// Shape `(n_trees, n_bipartitions)`, encoded as a flat row-major `uint8` byte buffer.
/// Column ordering is deterministic (ascending `Bitset` order) and stable across
/// calls on the same tree set. Reconstruct on the Python side:
/// ```python
/// presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n_trees, n_bip).copy()
/// ```
///
/// # Bipartition clade bytes
///
/// `bipartition_clade_bytes` is a flat `bytes` buffer encoding the canonical leaf
/// membership of every bipartition. Shape: `(n_bipartitions, ceil(n_leaves / 8))`,
/// row-major. Bit `i` of row `j` — **little-endian bit order within each byte** — is `1`
/// if `leaf_names[i]` is on the canonical side of bipartition `j`. Column order matches
/// the presence matrix (ascending `Bitset` order, stable across calls).
///
/// For unrooted trees the canonical side is the half that does **not** contain the first
/// leaf alphabetically, so bit 0 of every row is always `0`.
///
/// Decode with NumPy:
/// ```python
/// import math
/// bytes_per_bip = math.ceil(len(leaf_names) / 8)
/// bip_arr  = np.frombuffer(bipartition_clade_bytes, dtype=np.uint8)
///                .reshape(n_bip, bytes_per_bip)
/// bip_bool = np.unpackbits(bip_arr, axis=1, bitorder='little')[:, :len(leaf_names)]
/// # bip_bool[j, i] == 1  →  leaf_names[i] is in bipartition j
///
/// # Human-readable column labels
/// col_labels = [
///     "|".join(n for i, n in enumerate(leaf_names) if bip_bool[j, i])
///     for j in range(n_bip)
/// ]
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
///     6-tuple (tree_names, rf_matrix_bytes, leaf_names, n_bipartitions, presence_bytes,
///     bipartition_clade_bytes).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None, use_gpu=false))]
fn pairwise_rf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
    use_gpu: bool,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps =
        collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted, false)?;

    let n = snaps.snapshots.len();
    let rf_matrix = with_counter(py, progress, n_pairs(n), |counter| {
        crate::distances::dispatch_rf(&snaps, Some(counter), use_gpu)
    })?;
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(n)
        .flat_map(|row| row.iter().flat_map(|&v| (v as u32).to_ne_bytes()))
        .collect();

    let n_bipartitions = snaps.bipartitions.len();
    let (presence_vec, col_to_bip_id) = snaps.build_presence_matrix();
    let leaf_names = snaps.leaf_names.clone();
    let bip_clade_bytes = snaps.build_bipartition_bytes(&col_to_bip_id);

    let py_rf = PyBytes::new(py, &rf_bytes);
    let py_pres = PyBytes::new(py, &presence_vec);
    let py_bip = PyBytes::new(py, &bip_clade_bytes);
    Ok((
        names,
        py_rf.into(),
        leaf_names,
        n_bipartitions,
        py_pres.into(),
        py_bip.into(),
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
/// Column order matches `bipartition_clade_bytes` (ascending `Bitset` order, stable).
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
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None, use_gpu=false))]
fn pairwise_wrf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
    use_gpu: bool,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps =
        collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted, true)?;

    let n = snaps.snapshots.len();
    let wrf_matrix = with_counter(py, progress, n_pairs(n), |counter| {
        crate::distances::dispatch_wrf(&snaps, Some(counter), use_gpu)
    })?;
    let wrf_bytes: Vec<u8> = wrf_matrix.iter().flat_map(|&v| v.to_ne_bytes()).collect();

    let n_bip = snaps.bipartitions.len();
    let leaf_names = snaps.leaf_names.clone();
    let (bl_bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
    let bip_bytes = snaps.build_bipartition_bytes(&col_to_bip_id);

    let py_wrf = PyBytes::new(py, &wrf_bytes);
    let py_bl = PyBytes::new(py, &bl_bytes);
    let py_bip = PyBytes::new(py, &bip_bytes);
    Ok((
        names,
        py_wrf.into(),
        leaf_names,
        n_bip,
        py_bl.into(),
        py_bip.into(),
    ))
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
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None, use_gpu=false))]
fn pairwise_kf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
    use_gpu: bool,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps =
        collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted, true)?;

    let n = snaps.snapshots.len();
    let kf_matrix = with_counter(py, progress, n_pairs(n), |counter| {
        crate::distances::dispatch_kf(&snaps, Some(counter), use_gpu)
    })?;
    let kf_bytes: Vec<u8> = kf_matrix.iter().flat_map(|&v| v.to_ne_bytes()).collect();

    let n_bip = snaps.bipartitions.len();
    let leaf_names = snaps.leaf_names.clone();
    let (bl_bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
    let bip_bytes = snaps.build_bipartition_bytes(&col_to_bip_id);

    let py_kf = PyBytes::new(py, &kf_bytes);
    let py_bl = PyBytes::new(py, &bl_bytes);
    let py_bip = PyBytes::new(py, &bip_bytes);
    Ok((
        names,
        py_kf.into(),
        leaf_names,
        n_bip,
        py_bl.into(),
        py_bip.into(),
    ))
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
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None, use_gpu=false))]
fn pairwise_wrf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
    use_gpu: bool,
) -> PyResult<(Vec<String>, Vec<f64>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let mut snaps =
        collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted, true)?;

    // This path never exports snapshots, and the dense WRF computation reads
    // only `split_ids` + `lengths` — release the bipartition bitsets first.
    snaps.bipartitions = Vec::new();

    let n = snaps.snapshots.len();
    let matrix = with_counter(py, progress, n_pairs(n), |counter| {
        crate::distances::dispatch_wrf(&snaps, Some(counter), use_gpu)
    })?;
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
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None, use_gpu=false))]
fn pairwise_kf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
    use_gpu: bool,
) -> PyResult<(Vec<String>, Vec<f64>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let mut snaps =
        collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted, true)?;

    // This path never exports snapshots, and the dense KF computation reads
    // only `split_ids` + `lengths` — release the bipartition bitsets first.
    snaps.bipartitions = Vec::new();

    let n = snaps.snapshots.len();
    let matrix = with_counter(py, progress, n_pairs(n), |counter| {
        crate::distances::dispatch_kf(&snaps, Some(counter), use_gpu)
    })?;
    Ok((names, matrix))
}

/// Python module definition
#[pymodule]
fn rapidtrees(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pairwise_rf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(
        pairwise_rf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        pairwise_wrf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        pairwise_kf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(pairwise_wrf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_kf_from_newick_iter, m)?)?;
    m.add_class::<ProgressCounter>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::snapshot::Snapshots;
    use std::collections::HashMap;

    fn empty_map() -> HashMap<String, String> {
        HashMap::new()
    }

    #[test]
    fn test_snapshots_from_newicks_basic() {
        let snaps = Snapshots::from_newicks(
            &["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"],
            false,
        )
        .unwrap();
        assert_eq!(snaps.len(), 2);
        assert_eq!(snaps.leaf_names, vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn test_snapshots_from_newicks_mismatched_leaves_errors() {
        let result = Snapshots::from_newicks(
            &["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,B:1):1,(C:1,E:1):1);"],
            false,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_snapshots_from_newicks_leaf_names_sorted() {
        let snaps = Snapshots::from_newicks(&["((D:1,C:1):1,(B:1,A:1):1);"], false).unwrap();
        assert_eq!(
            snaps.leaf_names,
            vec!["A", "B", "C", "D"],
            "leaf names must be sorted"
        );
    }

    #[test]
    fn test_snapshots_from_newick_iter_uses_translate() {
        let translate: HashMap<String, String> = [
            ("1".to_string(), "A".to_string()),
            ("2".to_string(), "B".to_string()),
            ("3".to_string(), "C".to_string()),
            ("4".to_string(), "D".to_string()),
        ]
        .into();
        let newick = "((1:1,2:1):1,(3:1,4:1):1);";
        let entries = [(newick, &translate), (newick, &translate)];
        let snaps = Snapshots::from_newick_iter(entries, false).unwrap();
        assert_eq!(snaps.len(), 2);
        assert_eq!(snaps.leaf_names, vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn test_load_beast_raw_nonexistent_returns_empty() {
        let _ = empty_map(); // suppress unused warning
        let (_, pairs) = crate::io::load_beast_raw("nonexistent.trees", 0, 0, false);
        assert!(pairs.is_empty());
    }
}

/// Embedded-Python integration tests that exercise every PyO3 entry point
/// from Rust.
///
/// **Why** — `cargo test --lib` would normally never enter any `#[pyfunction]`
/// body or `ProgressCounter`'s `#[pymethods]`: they are only reachable through
/// the Python interpreter that loads the wheel. That made codecov report all
/// PyO3 surface as `FNDA:0` (untested) even though pytest exercised it.
///
/// This module pulls the Python interpreter into the test binary instead.
/// `append_to_inittab!(rapidtrees)` registers the `#[pymodule] rapidtrees`
/// function into the interpreter's built-in module table, and
/// `prepare_freethreaded_python()` initialises CPython. From then on we can
/// `py.import("rapidtrees")` and call the real `pairwise_*` / `ProgressCounter`
/// functions exactly like Python would. cargo test now measures coverage of
/// the whole PyO3 surface without needing pytest to be run under
/// `cargo llvm-cov`.
///
/// Requires `cargo test --features python` (no `extension-module`, so the
/// test binary links libpython) and `DYLD_LIBRARY_PATH` / `LD_LIBRARY_PATH`
/// pointing at the python lib dir — pixi tasks set those for you.
#[cfg(test)]
mod py_integration_tests {
    use pyo3::prelude::*;
    use pyo3::types::{PyDict, PyList};
    use std::sync::Once;

    // Bring the `#[pymodule] fn rapidtrees` from the parent module into scope
    // so the `append_to_inittab!` macro can find it by name.
    use super::rapidtrees;

    /// Initialise the embedded interpreter exactly once across all tests
    /// running in this process.
    fn ensure_python() {
        static INIT: Once = Once::new();
        INIT.call_once(|| {
            pyo3::append_to_inittab!(rapidtrees);
            Python::initialize();
        });
    }

    /// Three small newicks with a known RF profile, ready to hand to a
    /// pairwise function from Python.
    fn fixture<'py>(py: Python<'py>) -> (Bound<'py, PyList>, Bound<'py, PyList>) {
        let trees = PyList::new(
            py,
            [
                "(A:0.1,(B:0.1,C:0.1):0.1);",
                "(A:0.1,(C:0.1,B:0.1):0.1);",
                "((A:0.1,B:0.1):0.1,C:0.1);",
            ],
        )
        .unwrap();
        let names = PyList::new(py, ["t0", "t1", "t2"]).unwrap();
        (names, trees)
    }

    fn call_pairwise<'py>(
        py: Python<'py>,
        func_name: &str,
        progress: Option<&Bound<'py, PyAny>>,
    ) -> Bound<'py, PyAny> {
        let rapidtrees = py.import("rapidtrees").expect("import rapidtrees");
        let func = rapidtrees.getattr(func_name).expect(func_name);
        let (names, trees) = fixture(py);
        let translate_maps = PyList::new(py, [PyDict::new(py)]).unwrap();
        let map_indices = PyList::new(py, [0i64, 0, 0]).unwrap();
        let kwargs = PyDict::new(py);
        if let Some(pc) = progress {
            kwargs.set_item("progress", pc).unwrap();
        }
        func.call(
            (
                names,
                trees.try_iter().unwrap(),
                translate_maps,
                map_indices,
            ),
            Some(&kwargs),
        )
        .unwrap_or_else(|e| panic!("{func_name} call failed: {e}"))
    }

    #[test]
    fn rapidtrees_module_exports_expected_symbols() {
        ensure_python();
        Python::attach(|py| {
            let m = py.import("rapidtrees").unwrap();
            for name in [
                "pairwise_rf_from_newick_iter",
                "pairwise_wrf_from_newick_iter",
                "pairwise_kf_from_newick_iter",
                "pairwise_rf_with_snapshots_from_newick_iter",
                "pairwise_wrf_with_snapshots_from_newick_iter",
                "pairwise_kf_with_snapshots_from_newick_iter",
                "ProgressCounter",
            ] {
                assert!(
                    m.hasattr(name).unwrap(),
                    "rapidtrees module is missing `{name}`"
                );
            }
        });
    }

    /// Build a fresh ProgressCounter for tests that want to drive the
    /// counter-bearing branch of `with_counter`. Without this the
    /// counter-handling closure for that metric+variant monomorphization
    /// stays at FNDA:0 in lcov.
    fn fresh_counter<'py>(py: Python<'py>) -> Bound<'py, PyAny> {
        py.import("rapidtrees")
            .unwrap()
            .getattr("ProgressCounter")
            .unwrap()
            .call0()
            .unwrap()
    }

    #[test]
    fn pairwise_rf_via_python_returns_bytes() {
        ensure_python();
        Python::attach(|py| {
            let pc = fresh_counter(py);
            let result = call_pairwise(py, "pairwise_rf_from_newick_iter", Some(&pc));
            let tuple: (Vec<String>, Vec<u8>) = result.extract().unwrap();
            assert_eq!(tuple.0.len(), 3);
            assert_eq!(tuple.1.len(), 3 * 3 * 4); // n*n*sizeof(u32)
            // Also exercises with_counter's `Some` branch + final `store` line.
            let frac: f64 = pc.call_method0("fraction").unwrap().extract().unwrap();
            assert_eq!(frac, 1.0);
        });
    }

    #[test]
    fn pairwise_wrf_via_python() {
        ensure_python();
        Python::attach(|py| {
            let pc = fresh_counter(py);
            let result = call_pairwise(py, "pairwise_wrf_from_newick_iter", Some(&pc));
            let (_names, matrix): (Vec<String>, Vec<f64>) = result.extract().unwrap();
            assert_eq!(matrix.len(), 9);
        });
    }

    #[test]
    fn pairwise_kf_via_python() {
        ensure_python();
        Python::attach(|py| {
            let pc = fresh_counter(py);
            let result = call_pairwise(py, "pairwise_kf_from_newick_iter", Some(&pc));
            let (_names, matrix): (Vec<String>, Vec<f64>) = result.extract().unwrap();
            assert_eq!(matrix.len(), 9);
        });
    }

    #[test]
    fn pairwise_rf_with_snapshots_via_python() {
        ensure_python();
        Python::attach(|py| {
            let pc = fresh_counter(py);
            let _ = call_pairwise(py, "pairwise_rf_with_snapshots_from_newick_iter", Some(&pc));
        });
    }

    #[test]
    fn pairwise_wrf_with_snapshots_via_python() {
        ensure_python();
        Python::attach(|py| {
            let pc = fresh_counter(py);
            let _ = call_pairwise(
                py,
                "pairwise_wrf_with_snapshots_from_newick_iter",
                Some(&pc),
            );
        });
    }

    #[test]
    fn pairwise_kf_with_snapshots_via_python() {
        ensure_python();
        Python::attach(|py| {
            let pc = fresh_counter(py);
            let _ = call_pairwise(py, "pairwise_kf_with_snapshots_from_newick_iter", Some(&pc));
        });
    }

    /// One call without a counter, exercising the fast path inside
    /// `with_counter` (the `None` arm — different code path from the
    /// counter-bearing branch the rest of the tests cover).
    #[test]
    fn pairwise_without_counter_takes_fast_path() {
        ensure_python();
        Python::attach(|py| {
            let result = call_pairwise(py, "pairwise_rf_from_newick_iter", None);
            let tuple: (Vec<String>, Vec<u8>) = result.extract().unwrap();
            assert_eq!(tuple.0.len(), 3);
        });
    }

    #[test]
    fn progress_counter_lifecycle() {
        ensure_python();
        Python::attach(|py| {
            let m = py.import("rapidtrees").unwrap();
            let pc_class = m.getattr("ProgressCounter").unwrap();
            let pc = pc_class.call0().unwrap();

            // Initial state — exercises ProgressCounter::value/total/fraction.
            let zero: usize = pc.call_method0("value").unwrap().extract().unwrap();
            assert_eq!(zero, 0);
            let zero_total: usize = pc.call_method0("total").unwrap().extract().unwrap();
            assert_eq!(zero_total, 0);
            let zero_frac: f64 = pc.call_method0("fraction").unwrap().extract().unwrap();
            assert_eq!(zero_frac, 0.0);

            // __repr__ — codecov was flagging this as uncovered.
            let repr = pc.repr().unwrap().to_string();
            assert!(repr.contains("ProgressCounter"));
            assert!(repr.contains("value=0"));
            assert!(repr.contains("total=0"));

            // Drive the counter through a real pairwise call. After the call
            // value == total and fraction == 1.0.
            let _ = call_pairwise(py, "pairwise_rf_from_newick_iter", Some(&pc));
            let total: usize = pc.call_method0("total").unwrap().extract().unwrap();
            let value: usize = pc.call_method0("value").unwrap().extract().unwrap();
            let frac: f64 = pc.call_method0("fraction").unwrap().extract().unwrap();
            assert_eq!(total, 3); // n*(n-1)/2 for n=3
            assert_eq!(value, total);
            assert_eq!(frac, 1.0);

            // reset() — only line in __pymethod_reset__ that the unit tests
            // alone never hit.
            pc.call_method0("reset").unwrap();
            let after: usize = pc.call_method0("value").unwrap().extract().unwrap();
            let after_total: usize = pc.call_method0("total").unwrap().extract().unwrap();
            assert_eq!(after, 0);
            assert_eq!(after_total, 0);
        });
    }

    #[test]
    fn invalid_args_raise_value_error() {
        ensure_python();
        Python::attach(|py| {
            // names < 2 → ValueError from validate_iter_args.
            let m = py.import("rapidtrees").unwrap();
            let func = m.getattr("pairwise_rf_from_newick_iter").unwrap();
            let names = PyList::new(py, ["only-one"]).unwrap();
            let trees = PyList::new(py, ["(A:0.1,B:0.1);"]).unwrap();
            let translate_maps = PyList::new(py, [PyDict::new(py)]).unwrap();
            let map_indices = PyList::new(py, [0i64]).unwrap();
            let err = func
                .call1((
                    names,
                    trees.try_iter().unwrap(),
                    translate_maps,
                    map_indices,
                ))
                .expect_err("expected ValueError for n<2");
            assert!(err.is_instance_of::<pyo3::exceptions::PyValueError>(py));
        });
    }
}
