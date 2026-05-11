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
    Vec<Vec<u32>>,
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
/// # Bipartition leaf indices
///
/// `bipartition_leaf_indices` is a list of length `n_bipartitions`. Each entry is a list
/// of integer indices into `leaf_names` identifying the leaves on the canonical side of
/// that bipartition. For unrooted trees the canonical side is defined as the half that
/// does not contain the first leaf alphabetically.
///
/// Use it to build human-readable column labels for the presence matrix:
/// ```python
/// col_labels = [
///     "|".join(leaf_names[i] for i in indices)
///     for indices in bipartition_leaf_indices
/// ]
/// ```
///
/// For post-hoc analysis with pandas, construct a named DataFrame directly:
/// ```python
/// import pandas as pd
/// import numpy as np
///
/// tree_names, rf_bytes, leaf_names, n_bip, pres_bytes, bip_leaf_indices = (
///     pairwise_rf_with_snapshots_from_newick_iter(...)
/// )
/// presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(len(tree_names), n_bip).copy()
/// col_labels = [
///     "|".join(leaf_names[i] for i in indices)
///     for indices in bip_leaf_indices
/// ]
/// df = pd.DataFrame(presence, index=tree_names, columns=col_labels)
/// # Each column is now a named bipartition, e.g. "C|D|E", and each row is a tree.
/// # Columns that vary across trees (df.std() > 0) are the informative splits.
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
///     bipartition_leaf_indices).
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
    let (presence_vec, sorted_bipartitions) = snaps.build_presence_matrix();
    let leaf_names = snaps.leaf_names.clone();
    let n_leaves = leaf_names.len();
    let bipartition_leaf_indices: Vec<Vec<u32>> = sorted_bipartitions
        .iter()
        .map(|bitset| {
            (0..n_leaves)
                .filter(|&i| bitset.get(i))
                .map(|i| i as u32)
                .collect()
        })
        .collect();

    let py_rf = PyBytes::new(py, &rf_bytes);
    let py_pres = PyBytes::new(py, &presence_vec);
    Ok((
        names,
        py_rf.into(),
        leaf_names,
        n_bipartitions,
        py_pres.into(),
        bipartition_leaf_indices,
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
