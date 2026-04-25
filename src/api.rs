//! Python binding layer for tree distance calculations.
//!
//! All computation lives in [`crate::distances`] and [`crate::io`]; this module
//! contains only the PyO3-specific glue: iterating Python iterators, wrapping
//! byte buffers into `PyBytes`, and registering functions in the Python module.

use phylotree::tree::Tree as PhyloTree;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator};
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::bitset::Bitset;
use crate::distances::{pairwise_kf_matrix, pairwise_rf_matrix, pairwise_wrf_matrix};
use crate::io::{
    read_and_snapshot_trees, rename_leaf_nodes, snapshot_newicks, strip_beast_annotations,
};
use crate::snapshot::TreeSnapshot;

type PyRfSnapshotResult = (Vec<String>, Py<PyAny>, Vec<String>, usize, Py<PyAny>);

/// Compute pairwise Robinson-Foulds distances from multiple tree files.
///
/// Args:
///     paths: List of file paths to BEAST/NEXUS tree files.
///     burnin_trees: Number of trees to skip at the start of each file (default: 0).
///     burnin_states: Minimum STATE value to keep trees (default: 0).
///     use_real_taxa: Use TRANSLATE block for taxon names when available (default: True).
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of RF distances.
///
/// Raises:
///     ValueError: If no trees are found, trees have different leaf sets, or fewer than 2 trees.
#[pyfunction]
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true, rooted=false))]
fn pairwise_rf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<usize>>)> {
    let (names, snapshots) =
        read_and_snapshot_trees(&paths, burnin_trees, burnin_states, use_real_taxa, rooted)
            .map_err(PyValueError::new_err)?;
    Ok((names, pairwise_rf_matrix(&snapshots)))
}

/// Compute pairwise Weighted Robinson-Foulds distances from multiple tree files.
///
/// This metric considers branch lengths when comparing trees.
///
/// Args:
///     paths: List of file paths to BEAST/NEXUS tree files.
///     burnin_trees: Number of trees to skip at the start of each file (default: 0).
///     burnin_states: Minimum STATE value to keep trees (default: 0).
///     use_real_taxa: Use TRANSLATE block for taxon names when available (default: True).
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of Weighted RF distances.
///
/// Raises:
///     ValueError: If no trees are found, trees have different leaf sets, or fewer than 2 trees.
#[pyfunction]
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true, rooted=false))]
fn pairwise_weighted_rf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<f64>>)> {
    let (names, snapshots) =
        read_and_snapshot_trees(&paths, burnin_trees, burnin_states, use_real_taxa, rooted)
            .map_err(PyValueError::new_err)?;
    Ok((names, pairwise_wrf_matrix(&snapshots)))
}

/// Compute pairwise Kuhner-Felsenstein (Branch Score) distances from multiple tree files.
///
/// This metric uses squared differences of branch lengths: sqrt(Σ(length_a - length_b)²).
///
/// Args:
///     paths: List of file paths to BEAST/NEXUS tree files.
///     burnin_trees: Number of trees to skip at the start of each file (default: 0).
///     burnin_states: Minimum STATE value to keep trees (default: 0).
///     use_real_taxa: Use TRANSLATE block for taxon names when available (default: True).
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, distance_matrix) — a 2D list of KF distances.
///
/// Raises:
///     ValueError: If no trees are found, trees have different leaf sets, or fewer than 2 trees.
#[pyfunction]
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true, rooted=false))]
fn pairwise_kf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<f64>>)> {
    let (names, snapshots) =
        read_and_snapshot_trees(&paths, burnin_trees, burnin_states, use_real_taxa, rooted)
            .map_err(PyValueError::new_err)?;
    Ok((names, pairwise_kf_matrix(&snapshots)))
}

/// Compute pairwise Robinson-Foulds distances from newick strings with translate maps.
///
/// Args:
///     names: List of tree identifiers (one per newick).
///     newicks: List of newick strings (may contain BEAST annotations).
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, distance_matrix) — names are the input names, matrix is a 2D list.
///
/// Raises:
///     ValueError: If lengths mismatch, indices are out of bounds, or fewer than 2 trees.
#[pyfunction]
#[pyo3(signature = (names, newicks, translate_maps, map_indices, rooted=false))]
fn pairwise_rf_from_newicks(
    names: Vec<String>,
    newicks: Vec<String>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<usize>>)> {
    if names.len() != newicks.len() {
        return Err(PyValueError::new_err(format!(
            "names length ({}) must equal newicks length ({})",
            names.len(),
            newicks.len()
        )));
    }
    let (snapshots, _leaf_names) =
        snapshot_newicks(&newicks, &translate_maps, &map_indices, rooted)
            .map_err(PyValueError::new_err)?;
    Ok((names, pairwise_rf_matrix(&snapshots)))
}

/// Collect tree snapshots from a lazy Python iterator of newick strings.
///
/// This mirrors `snapshot_newicks`, but it must validate incrementally because a
/// lazy Python iterator cannot be batch-validated through the shared `io.rs` helper.
///
/// This function must live in `api.rs` because it consumes a `Bound<'_, PyIterator>`.
fn collect_snapshots_from_iter(
    n: usize,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: &[HashMap<String, String>],
    map_indices: &[usize],
    rooted: bool,
) -> PyResult<(Vec<TreeSnapshot>, Vec<String>)> {
    let mut snapshots: Vec<TreeSnapshot> = Vec::with_capacity(n);
    let mut sorted_leaf_names: Vec<String> = Vec::new();
    let mut reference_leaves: Option<HashSet<String>> = None;
    let mut tree_count: usize = 0;

    for item in newick_iter {
        let item = item?;
        let newick: String = item.extract()?;

        if tree_count >= n {
            return Err(PyValueError::new_err(format!(
                "Iterator yielded more than {} newick strings (expected {})",
                n, n
            )));
        }

        let map_idx = map_indices[tree_count];
        let clean = strip_beast_annotations(&newick);
        let mut tree = PhyloTree::from_newick(&clean).map_err(|e| {
            PyValueError::new_err(format!(
                "Failed to parse newick at index {}: {}",
                tree_count, e
            ))
        })?;
        rename_leaf_nodes(&mut tree, &translate_maps[map_idx]);

        let leaves: HashSet<String> = tree
            .get_leaves()
            .iter()
            .filter_map(|&id| tree.get(&id).ok()?.name.clone())
            .collect();

        match &reference_leaves {
            None => {
                let mut sorted: Vec<String> = leaves.iter().cloned().collect();
                sorted.sort_unstable();
                sorted_leaf_names = sorted;
                reference_leaves = Some(leaves);
            }
            Some(ref_leaves) => {
                if leaves.len() != ref_leaves.len() {
                    return Err(PyValueError::new_err(format!(
                        "Tree {} has {} leaves, but tree 0 has {} leaves. All trees must have the same number of leaves.",
                        tree_count,
                        leaves.len(),
                        ref_leaves.len()
                    )));
                }
                if leaves != *ref_leaves {
                    return Err(PyValueError::new_err(format!(
                        "Tree {} has different leaf set than tree 0. All trees must have the same taxa.",
                        tree_count
                    )));
                }
            }
        }

        let snapshot = TreeSnapshot::from_tree(&tree, rooted).map_err(|e| {
            PyValueError::new_err(format!(
                "Failed to create tree snapshot at index {}: {}",
                tree_count, e
            ))
        })?;
        snapshots.push(snapshot);
        tree_count += 1;
        // `newick` and `tree` are dropped here — only the snapshot is retained
    }

    if tree_count != n {
        return Err(PyValueError::new_err(format!(
            "Iterator yielded {} newick strings, but names has {} entries",
            tree_count, n
        )));
    }

    Ok((snapshots, sorted_leaf_names))
}

/// Validate argument consistency for iterator-based functions.
fn validate_iter_args(
    n: usize,
    map_indices: &[usize],
    translate_maps: &[HashMap<String, String>],
) -> PyResult<()> {
    if n != map_indices.len() {
        return Err(PyValueError::new_err(format!(
            "names length ({}) must equal map_indices length ({})",
            n,
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

/// Compute pairwise RF distances from a lazy Python iterator of newick strings.
///
/// Thin wrapper around [`pairwise_rf_with_snapshots_from_newick_iter`] that
/// discards the snapshot outputs, ensuring both functions share a single
/// implementation.
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings. Must yield exactly
///         len(names) strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     (tree_names, matrix_bytes) — matrix_bytes is flat u32 bytes (row-major n×n).
///
/// Raises:
///     ValueError: If lengths mismatch, indices are out of bounds, trees have
///                 different leaf sets, or fewer than 2 trees.
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
    let (names, rf_bytes, _leaf_names, _n_bip, _pres_bytes) =
        pairwise_rf_with_snapshots_from_newick_iter(
            py,
            names,
            newick_iter,
            translate_maps,
            map_indices,
            rooted,
        )?;
    Ok((names, rf_bytes))
}

/// Compute pairwise RF distances and export binary tree snapshots in a single pass.
///
/// Parses each newick once, building both the RF distance matrix and a binary
/// presence matrix encoding which bipartitions appear in each tree. The presence
/// matrix enables all topology-based convergence diagnostics (Pseudo ESS, Fréchet
/// ESS, ASDSF) without re-parsing the original tree files.
///
/// # Presence matrix format
///
/// Shape `(n_trees, n_bipartitions)`, encoded as a flat row-major `uint8` byte buffer.
/// Column ordering is deterministic (ascending `Bitset` order via `BTreeSet`)
/// and stable across calls on the same tree set. Reconstruct on the Python side:
/// ```python
/// presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n_trees, n_bip).copy()
/// ```
///
/// # RF trace from presence matrix
///
/// ```python
/// rf_trace = np.sum(presence ^ presence[r:r+1], axis=1)
/// # Equals rf_matrix[r] exactly.
/// ```
///
/// Args:
///     names: Tree identifiers (one per newick).
///     newick_iter: Python iterator yielding newick strings. Must yield exactly
///         len(names) strings.
///     translate_maps: List of translate maps (number → taxon name).
///     map_indices: Per-tree index into translate_maps. Defaults to all-zero.
///     rooted: If True compare clades; if False compare bipartitions (default: False).
///
/// Returns:
///     5-tuple (tree_names, rf_matrix_bytes, leaf_names, n_bipartitions, presence_bytes):
///     - tree_names: input names in order.
///     - rf_matrix_bytes: flat bytes of n×n u32 RF distances (row-major).
///     - leaf_names: taxa sorted alphabetically (defines the bipartition columns).
///     - n_bipartitions: number of unique bipartitions across all trees.
///     - presence_bytes: flat bytes of uint8, shape (n_trees, n_bipartitions).
///
/// Raises:
///     ValueError: If lengths mismatch, indices are out of bounds, trees have
///                 different leaf sets, or fewer than 2 trees.
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
    let n = names.len();
    validate_iter_args(n, &map_indices, &translate_maps)?;

    let (snapshots, leaf_names) =
        collect_snapshots_from_iter(n, newick_iter, &translate_maps, &map_indices, rooted)?;

    if snapshots.len() < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }

    let rf_matrix = pairwise_rf_matrix(&snapshots);
    let rf_bytes: Vec<u8> = rf_matrix
        .iter()
        .flat_map(|row| row.iter().flat_map(|&v| (v as u32).to_ne_bytes()))
        .collect();
    let (n_bipartitions, presence_vec) = {
        let sorted_bips: Vec<Bitset> = snapshots
            .iter()
            .flat_map(|s| s.parts.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let n_bip = sorted_bips.len();
        let mut bytes = vec![0u8; snapshots.len() * n_bip];
        for (i, snap) in snapshots.iter().enumerate() {
            for part in &snap.parts {
                if let Ok(j) = sorted_bips.binary_search(part) {
                    bytes[i * n_bip + j] = 1;
                }
            }
        }
        (n_bip, bytes)
    };

    let py_rf = PyBytes::new(py, &rf_bytes);
    let py_pres = PyBytes::new(py, &presence_vec);
    Ok((
        names,
        py_rf.into(),
        leaf_names,
        n_bipartitions,
        py_pres.into(),
    ))
}

/// Python module definition
#[pymodule]
fn rapidtrees(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pairwise_rf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_weighted_rf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_kf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_rf_from_newicks, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_rf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(
        pairwise_rf_with_snapshots_from_newick_iter,
        m
    )?)?;
    Ok(())
}
