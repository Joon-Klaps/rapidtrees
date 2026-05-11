//! Python binding layer for tree distance calculations.
//!
//! All computation lives in [`crate::distances`] and [`crate::snapshot`]; this module
//! contains only the PyO3-specific glue: iterating Python iterators, wrapping
//! byte buffers into `PyBytes`, and registering functions in the Python module.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator};
use std::collections::HashMap;

use crate::snapshot::Snapshots;

type PyRfSnapshotResult = (
    Vec<String>,
    Py<PyAny>,
    Vec<String>,
    usize,
    Py<PyAny>,
    Py<PyAny>,
);

/// Collect tree snapshots from a lazy Python iterator of newick strings.
fn collect_snapshots_from_iter(
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: &[HashMap<String, String>],
    map_indices: &[usize],
    rooted: bool,
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
    Snapshots::from_newick_iter(entries, rooted).map_err(PyValueError::new_err)
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
///
/// Returns:
///     (tree_names, rf_matrix_bytes) — rf_matrix_bytes is flat u32 bytes (row-major n×n).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false))]
fn pairwise_rf_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
) -> PyResult<(Vec<String>, Py<PyAny>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps = collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted)?;

    let rf_matrix = snaps.pairwise_rf();
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(snaps.snapshots.len())
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
///
/// Returns:
///     6-tuple (tree_names, rf_matrix_bytes, leaf_names, n_bipartitions, presence_bytes,
///     bipartition_clade_bytes).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false))]
fn pairwise_rf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps = collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted)?;

    let rf_matrix = snaps.pairwise_rf();
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(snaps.snapshots.len())
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

/// Compute pairwise Weighted Robinson-Foulds distances from a lazy Python iterator.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of Weighted RF distances.
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false))]
fn pairwise_wrf_from_newick_iter(
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<f64>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps = collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted)?;

    Ok((names, snaps.pairwise_wrf()))
}

/// Compute pairwise Kuhner-Felsenstein (Branch Score) distances from a lazy Python iterator.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of KF distances.
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false))]
fn pairwise_kf_from_newick_iter(
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<f64>)> {
    validate_iter_args(&names, &map_indices, &translate_maps)?;

    let snaps = collect_snapshots_from_iter(newick_iter, &translate_maps, &map_indices, rooted)?;

    Ok((names, snaps.pairwise_kf()))
}

/// Python module definition
#[pymodule]
fn rapidtrees(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pairwise_rf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(
        pairwise_rf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(pairwise_wrf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_kf_from_newick_iter, m)?)?;
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
