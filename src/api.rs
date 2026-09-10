//! Python binding layer for tree distance calculations.
//!
//! All computation lives in [`crate::distances`] and [`crate::snapshot`]; this module
//! contains only the PyO3-specific glue: iterating Python iterators, wrapping
//! byte buffers into `PyBytes`, and registering functions in the Python module.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator};
use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;

use crate::progress::{ProgressCounter, with_counter};
use crate::snapshot::{Retain, Snapshots};

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

/// The tree source as handed in from Python: the newicks, how to rename their
/// taxa, and whether to compare clades or bipartitions.
///
/// Every entry point takes these four together and does nothing with them but
/// forward them, so they travel as one value.
struct IterInput<'py, 'a> {
    newick_iter: Bound<'py, PyIterator>,
    translate_maps: &'a [HashMap<String, String>],
    map_indices: &'a [usize],
    rooted: bool,
}

/// Collect tree snapshots from a lazy Python iterator of newick strings.
///
/// `retain` says what to build beyond the split IDs — see [`Retain`]. Each entry
/// point asks only for what it actually returns.
fn collect_snapshots_from_iter(input: IterInput<'_, '_>, retain: Retain) -> PyResult<Snapshots> {
    let newicks: Vec<String> = input
        .newick_iter
        .map(|item| item?.extract::<String>())
        .collect::<PyResult<_>>()?;

    if newicks.len() < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }

    let entries = newicks
        .iter()
        .zip(input.map_indices.iter())
        .map(|(n, &idx)| (n.as_str(), &input.translate_maps[idx]));
    Snapshots::from_newick_iter_opts(entries, input.rooted, retain).map_err(PyValueError::new_err)
}

/// A pairwise weighted metric, as a plain function pointer so WRF and KF can
/// share one body.
type WeightedMetric = fn(&Snapshots, Option<&AtomicUsize>) -> Vec<f64>;

/// Shared body of the plain WRF and KF entry points.
///
/// These paths export nothing and the dense kernels read only `split_ids` +
/// `lengths`, so the bipartition table is never built in the first place.
fn weighted_pairwise(
    py: Python<'_>,
    input: IterInput<'_, '_>,
    progress: Option<Py<ProgressCounter>>,
    metric: WeightedMetric,
) -> PyResult<Vec<f64>> {
    let snaps = collect_snapshots_from_iter(
        input,
        Retain {
            lengths: true,
            bipartitions: false,
        },
    )?;

    let total = n_pairs(snaps.len());
    with_counter(py, progress, total, |counter| metric(&snaps, Some(counter)))
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
///         GIL stays held throughout.
///
/// Returns:
///     (tree_names, rf_matrix_bytes) — rf_matrix_bytes is flat u32 bytes (row-major n×n).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
fn pairwise_rf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
) -> PyResult<(Vec<String>, Py<PyAny>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let input = IterInput {
        newick_iter,
        translate_maps: &translate_maps,
        map_indices: &map_indices,
        rooted,
    };

    // This path exports nothing and the dense RF kernel reads only `split_ids`,
    // so neither lengths nor bipartitions are built.
    let snaps = collect_snapshots_from_iter(
        input,
        Retain {
            lengths: false,
            bipartitions: false,
        },
    )?;

    let n = snaps.len();
    let rf_matrix = with_counter(py, progress, n_pairs(n), |counter| {
        snaps.pairwise_rf(Some(counter))
    })?;
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(n)
        .flat_map(|row| row.iter().flat_map(|&v: &u32| v.to_ne_bytes()))
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
/// Column ordering is deterministic and stable across calls on the same tree set. Reconstruct on the Python side:
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
/// the presence matrix, and is stable across calls.
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
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
fn pairwise_rf_with_snapshots_from_newick_iter(
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

    let snaps = collect_snapshots_from_iter(
        input,
        Retain {
            lengths: false,
            bipartitions: true,
        },
    )?;

    let n = snaps.snapshots.len();
    let rf_matrix = with_counter(py, progress, n_pairs(n), |counter| {
        snaps.pairwise_rf(Some(counter))
    })?;
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(n)
        .flat_map(|row| row.iter().flat_map(|&v: &u32| v.to_ne_bytes()))
        .collect();

    let n_bipartitions = snaps.n_distinct_splits();
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
fn pairwise_wrf_with_snapshots_from_newick_iter(
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
fn pairwise_kf_with_snapshots_from_newick_iter(
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
fn pairwise_wrf_from_newick_iter(
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
    let matrix = weighted_pairwise(py, input, progress, Snapshots::pairwise_wrf)?;
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
fn pairwise_kf_from_newick_iter(
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
    let matrix = weighted_pairwise(py, input, progress, Snapshots::pairwise_kf)?;
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
