//! Python binding layer for tree distance calculations.
//!
//! Provides Python functions for computing pairwise tree distances
//! from BEAST/NEXUS tree files.

use phylotree::tree::Tree as PhyloTree;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::distances::{kf_from_snapshots, rf_from_snapshots, weighted_rf_from_snapshots};
use crate::io::{read_beast_trees, rename_leaf_nodes, strip_beast_annotations};
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
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true, rooted=false))]
fn pairwise_rf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<usize>>)> {
    // Read all trees from all files
    let (tree_names, trees) = read_all_trees(&paths, burnin_trees, burnin_states, use_real_taxa)?;

    // Perform sanity checks
    sanity_check_trees(&trees)?;

    // Build snapshots
    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(|t| TreeSnapshot::from_tree(t, rooted))
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
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true, rooted=false))]
fn pairwise_weighted_rf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<f64>>)> {
    let (tree_names, trees) = read_all_trees(&paths, burnin_trees, burnin_states, use_real_taxa)?;
    sanity_check_trees(&trees)?;

    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(|t| TreeSnapshot::from_tree(t, rooted))
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
#[pyo3(signature = (paths, burnin_trees=0, burnin_states=0, use_real_taxa=true, rooted=false))]
fn pairwise_kf(
    paths: Vec<String>,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<f64>>)> {
    let (tree_names, trees) = read_all_trees(&paths, burnin_trees, burnin_states, use_real_taxa)?;
    sanity_check_trees(&trees)?;

    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(|t| TreeSnapshot::from_tree(t, rooted))
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

/// Parse newick strings and apply translate maps to produce PhyloTrees.
///
/// Each newick is stripped of BEAST annotations, parsed, then leaf nodes are
/// renamed using the translate map selected by `map_indices[i]`.
pub(crate) fn parse_and_translate(
    newicks: &[String],
    translate_maps: &[HashMap<String, String>],
    map_indices: &[usize],
) -> Result<Vec<PhyloTree>, String> {
    if newicks.len() != map_indices.len() {
        return Err(format!(
            "newicks length ({}) must equal map_indices length ({})",
            newicks.len(),
            map_indices.len()
        ));
    }

    for (i, &idx) in map_indices.iter().enumerate() {
        if idx >= translate_maps.len() {
            return Err(format!(
                "map_indices[{}] = {} is out of bounds (only {} translate maps provided)",
                i,
                idx,
                translate_maps.len()
            ));
        }
    }

    let mut trees = Vec::with_capacity(newicks.len());
    for (i, newick) in newicks.iter().enumerate() {
        let clean = strip_beast_annotations(newick);
        let mut tree = PhyloTree::from_newick(&clean)
            .map_err(|e| format!("Failed to parse newick at index {}: {}", i, e))?;
        rename_leaf_nodes(&mut tree, &translate_maps[map_indices[i]]);
        trees.push(tree);
    }

    Ok(trees)
}

/// Compute pairwise Robinson-Foulds distances from newick strings with translate maps.
///
/// Args:
///     names: List of tree identifiers (one per newick)
///     newicks: List of newick strings (may contain BEAST annotations)
///     translate_maps: List of translate maps (number → taxon name)
///     map_indices: Per-tree index into translate_maps
///
/// Returns:
///     A tuple of (tree_names, distance_matrix) where:
///     - tree_names is the input names list
///     - distance_matrix is a 2D list of RF distances
///
/// Raises:
///     ValueError: If lengths mismatch, indices are out of bounds, or fewer than 2 trees
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

    let trees = parse_and_translate(&newicks, &translate_maps, &map_indices)
        .map_err(PyValueError::new_err)?;

    sanity_check_trees(&trees)?;

    let snapshots: Vec<TreeSnapshot> = trees
        .iter()
        .map(|t| TreeSnapshot::from_tree(t, rooted))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| PyValueError::new_err(format!("Failed to create tree snapshot: {}", e)))?;

    let n = snapshots.len();
    let mut matrix = vec![vec![0usize; n]; n];

    let pairs: Vec<(usize, usize, usize)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| {
            let dist = rf_from_snapshots(&snapshots[i], &snapshots[j]);
            (i, j, dist)
        })
        .collect();

    for (i, j, dist) in pairs {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    Ok((names, matrix))
}

/// Compute pairwise Robinson-Foulds distances from a Python iterator of newick strings.
///
/// Unlike `pairwise_rf_from_newicks`, this accepts a lazy Python iterator so that
/// only one newick string needs to be in memory at a time. Each newick is parsed
/// into a compact `TreeSnapshot` immediately and the raw string is discarded.
///
/// Args:
///     names: List of tree identifiers (one per newick)
///     newick_iter: Python iterator yielding newick strings
///     translate_maps: List of translate maps (number → taxon name)
///     map_indices: Per-tree index into translate_maps
///     rooted: If True compare clades; if False compare bipartitions (default: False)
///
/// Returns:
///     A tuple of (tree_names, distance_matrix)
///
/// Raises:
///     ValueError: If lengths mismatch, indices are out of bounds, trees have
///                 different leaf sets, or fewer than 2 trees
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false))]
fn pairwise_rf_from_newick_iter(
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
) -> PyResult<(Vec<String>, Vec<Vec<usize>>)> {
    let n = names.len();

    if n != map_indices.len() {
        return Err(PyValueError::new_err(format!(
            "names length ({}) must equal map_indices length ({})",
            n,
            map_indices.len()
        )));
    }

    // Validate map_indices upfront
    for (i, &idx) in map_indices.iter().enumerate() {
        if idx >= translate_maps.len() {
            return Err(PyValueError::new_err(format!(
                "map_indices[{}] = {} is out of bounds (only {} translate maps provided)",
                i, idx, translate_maps.len()
            )));
        }
    }

    let mut snapshots: Vec<TreeSnapshot> = Vec::with_capacity(n);
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

        // Incremental sanity check: verify leaf set matches the first tree
        let leaves: HashSet<String> = tree
            .get_leaves()
            .iter()
            .filter_map(|&id| tree.get(&id).ok()?.name.clone())
            .collect();

        match &reference_leaves {
            None => {
                reference_leaves = Some(leaves);
            }
            Some(ref_leaves) => {
                if leaves.len() != ref_leaves.len() {
                    return Err(PyValueError::new_err(format!(
                        "Tree {} has {} leaves, but tree 0 has {} leaves. All trees must have the same number of leaves.",
                        tree_count, leaves.len(), ref_leaves.len()
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
        // `newick` and `tree` are dropped here — only snapshot remains
    }

    if tree_count != n {
        return Err(PyValueError::new_err(format!(
            "Iterator yielded {} newick strings, but names has {} entries",
            tree_count, n
        )));
    }

    if tree_count < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }

    // Parallel pairwise RF computation (same as pairwise_rf_from_newicks)
    let mut matrix = vec![vec![0usize; tree_count]; tree_count];

    let pairs: Vec<(usize, usize, usize)> = (0..tree_count)
        .into_par_iter()
        .flat_map_iter(|i| (i + 1..tree_count).map(move |j| (i, j)))
        .map(|(i, j)| {
            let dist = rf_from_snapshots(&snapshots[i], &snapshots[j]);
            (i, j, dist)
        })
        .collect();

    for (i, j, dist) in pairs {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    Ok((names, matrix))
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
fn rust_python_tree_distances(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pairwise_rf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_weighted_rf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_kf, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_rf_from_newicks, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_rf_from_newick_iter, m)?)?;
    Ok(())
}
