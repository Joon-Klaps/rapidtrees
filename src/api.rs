//! Python binding layer for tree distance calculations.
//!
//! Provides Python functions for computing pairwise tree distances
//! from BEAST/NEXUS tree files.

use phylotree::tree::Tree as PhyloTree;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashSet;

use crate::distances::{kf_from_snapshots, rf_from_snapshots, weighted_rf_from_snapshots};
use crate::io::read_beast_trees;
use crate::snapshot::TreeSnapshot;

/// Compute pairwise Robinson-Foulds distances from multiple tree files.
///
/// Args:
///     paths: List of file paths to BEAST/NEXUS tree files
///     burnin_trees: Number of trees to skip at the beginning of each file (default: 0)
///     burnin_states: Minimum STATE value to keep trees (default: 0)
///     use_real_taxa: Use TRANSLATE block for taxon names when available (default: True)
///
/// Returns:
///     A tuple of (tree_names, distance_matrix) where:
///     - tree_names is a list of tree identifiers
///     - distance_matrix is a 2D list of RF distances
///
/// Raises:
///     ValueError: If no trees are found, trees have different leaf sets, or sanity checks fail
#[pyfunction]
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true))]
fn pairwise_rf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> PyResult<(Vec<String>, Vec<Vec<usize>>)> {
    // Read all trees from all files
    let (tree_names, trees) = read_all_trees(&paths, burnin_trees, burnin_states, use_real_taxa)?;

    // Perform sanity checks
    sanity_check_trees(&trees)?;

    // Build snapshots
    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(TreeSnapshot::from_tree)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| PyValueError::new_err(format!("Failed to create tree snapshot: {}", e)))?;

    // Compute pairwise distances
    let n = snapshots.len();
    let mut matrix = vec![vec![0usize; n]; n];

    // Parallel computation across all pairs
    let pairs: Vec<(usize, usize, usize)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| {
            let dist = rf_from_snapshots(&snapshots[i], &snapshots[j]);
            (i, j, dist)
        })
        .collect();

    // Fill matrix (symmetric)
    for (i, j, dist) in pairs {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    Ok((tree_names, matrix))
}

/// Compute pairwise Weighted Robinson-Foulds distances from multiple tree files.
///
/// This metric considers branch lengths when comparing trees.
///
/// Args:
///     paths: List of file paths to BEAST/NEXUS tree files
///     burnin_trees: Number of trees to skip at the beginning of each file (default: 0)
///     burnin_states: Minimum STATE value to keep trees (default: 0)
///     use_real_taxa: Use TRANSLATE block for taxon names when available (default: True)
///
/// Returns:
///     A tuple of (tree_names, distance_matrix) where:
///     - tree_names is a list of tree identifiers
///     - distance_matrix is a 2D list of weighted RF distances
///
/// Raises:
///     ValueError: If no trees are found, trees have different leaf sets, or sanity checks fail
#[pyfunction]
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true))]
fn pairwise_weighted_rf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> PyResult<(Vec<String>, Vec<Vec<f64>>)> {
    let (tree_names, trees) = read_all_trees(&paths, burnin_trees, burnin_states, use_real_taxa)?;
    sanity_check_trees(&trees)?;

    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(TreeSnapshot::from_tree)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| PyValueError::new_err(format!("Failed to create tree snapshot: {}", e)))?;

    let n = snapshots.len();
    let mut matrix = vec![vec![0.0f64; n]; n];

    let pairs: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| {
            let dist = weighted_rf_from_snapshots(&snapshots[i], &snapshots[j]);
            (i, j, dist)
        })
        .collect();

    for (i, j, dist) in pairs {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    Ok((tree_names, matrix))
}

/// Compute pairwise Kuhner-Felsenstein (Branch Score) distances from multiple tree files.
///
/// This metric uses squared differences of branch lengths: sqrt(Σ(length_a - length_b)²)
///
/// Args:
///     paths: List of file paths to BEAST/NEXUS tree files
///     burnin_trees: Number of trees to skip at the beginning of each file (default: 0)
///     burnin_states: Minimum STATE value to keep trees (default: 0)
///     use_real_taxa: Use TRANSLATE block for taxon names when available (default: True)
///
/// Returns:
///     A tuple of (tree_names, distance_matrix) where:
///     - tree_names is a list of tree identifiers
///     - distance_matrix is a 2D list of KF distances
///
/// Raises:
///     ValueError: If no trees are found, trees have different leaf sets, or sanity checks fail
#[pyfunction]
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true))]
fn pairwise_kf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> PyResult<(Vec<String>, Vec<Vec<f64>>)> {
    let (tree_names, trees) = read_all_trees(&paths, burnin_trees, burnin_states, use_real_taxa)?;
    sanity_check_trees(&trees)?;

    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(TreeSnapshot::from_tree)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| PyValueError::new_err(format!("Failed to create tree snapshot: {}", e)))?;

    let n = snapshots.len();
    let mut matrix = vec![vec![0.0f64; n]; n];

    let pairs: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| {
            let dist = kf_from_snapshots(&snapshots[i], &snapshots[j]);
            (i, j, dist)
        })
        .collect();

    for (i, j, dist) in pairs {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    Ok((tree_names, matrix))
}

/// Helper function to read trees from multiple files
fn read_all_trees(
    paths: &[String],
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> PyResult<(Vec<String>, Vec<PhyloTree>)> {
    let mut all_tree_names = Vec::new();
    let mut all_trees = Vec::new();

    for (file_idx, path) in paths.iter().enumerate() {
        let (_taxons, named_trees) = read_beast_trees(
            std::path::PathBuf::from(path),
            burnin_trees,
            burnin_states,
            use_real_taxa,
        );

        if named_trees.is_empty() {
            return Err(PyValueError::new_err(format!(
                "No trees found in file '{}' after burnin removal",
                path
            )));
        }

        // Add trees with file prefix in name
        for (name, tree) in named_trees {
            let full_name = format!("file{}_{}", file_idx, name);
            all_tree_names.push(full_name);
            all_trees.push(tree);
        }
    }

    if all_trees.is_empty() {
        return Err(PyValueError::new_err(
            "No trees found in any of the provided files",
        ));
    }

    Ok((all_tree_names, all_trees))
}

/// Perform sanity checks on trees
fn sanity_check_trees(trees: &[PhyloTree]) -> PyResult<()> {
    if trees.is_empty() {
        return Err(PyValueError::new_err("No trees to compare"));
    }

    if trees.len() < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }

    // Check that all trees have the same leaf set
    let first_leaves: HashSet<String> = trees[0]
        .get_leaves()
        .iter()
        .filter_map(|&id| trees[0].get(&id).ok()?.name.clone())
        .collect();

    let first_leaf_count = first_leaves.len();

    for (idx, tree) in trees.iter().enumerate().skip(1) {
        let leaves: HashSet<String> = tree
            .get_leaves()
            .iter()
            .filter_map(|&id| tree.get(&id).ok()?.name.clone())
            .collect();

        if leaves.len() != first_leaf_count {
            return Err(PyValueError::new_err(format!(
                "Tree {} has {} leaves, but tree 0 has {} leaves. All trees must have the same number of leaves.",
                idx,
                leaves.len(),
                first_leaf_count
            )));
        }

        if leaves != first_leaves {
            return Err(PyValueError::new_err(format!(
                "Tree {} has different leaf set than tree 0. All trees must have the same taxa.",
                idx
            )));
        }
    }

    Ok(())
}

/// Python module definition
#[pymodule]
fn rapidtrees(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pairwise_rf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_weighted_rf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_kf, m)?)?;
    Ok(())
}
