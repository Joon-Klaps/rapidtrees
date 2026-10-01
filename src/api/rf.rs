//! Robinson–Foulds entry points.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyIterator};
use std::collections::HashMap;

use super::input::{IterInput, collect_snapshots_from_iter, validate_iter_args};
use super::progress::{ProgressCounter, with_counter};
use super::{PyRfSnapshotResult, n_pairs};
use crate::snapshot::{Retain, Snapshots};

/// Shared body of the RF snapshot-exporting entry points.
///
/// They differ only in how they encode which clades each tree holds: `membership` builds that payload and returns it with the column order the clade bitmasks must follow. The tuple comes back in the dense endpoint's order, payload fifth and clade bitmasks sixth.
fn rf_with_snapshots(
    py: Python<'_>,
    names: Vec<String>,
    input: IterInput<'_, '_>,
    progress: Option<Py<ProgressCounter>>,
    membership: impl FnOnce(&Snapshots) -> PyResult<(Py<PyAny>, Vec<usize>)>,
) -> PyResult<PyRfSnapshotResult> {
    validate_iter_args(&names, input.map_indices, input.translate_maps)?;
    let snaps = collect_snapshots_from_iter(
        input,
        Retain {
            lengths: false,
            bipartitions: true,
        },
    )?;

    let rf_matrix = with_counter(py, progress, n_pairs(snaps.len()), |counter| {
        snaps.pairwise_rf(Some(counter))
    })?;
    let rf_bytes: Vec<u8> = rf_matrix.iter().flat_map(|v| v.to_ne_bytes()).collect();
    let (payload, col_to_bip_id) = membership(&snaps)?;
    let clade_bytes = snaps.build_bipartition_bytes(&col_to_bip_id);
    let n_clades = snaps.n_distinct_splits();

    Ok((
        names,
        PyBytes::new(py, &rf_bytes).into(),
        snaps.leaf_names,
        n_clades,
        payload,
        PyBytes::new(py, &clade_bytes).into(),
    ))
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
pub(super) fn pairwise_rf_from_newick_iter(
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
    let rf_matrix = with_counter(py, progress, n_pairs(n), move |counter| {
        snaps.into_pairwise_rf(Some(counter))
    })?;
    let rf_bytes: Vec<u8> = rf_matrix
        .chunks(n)
        .flat_map(|row| row.iter().flat_map(|&v: &u32| v.to_ne_bytes()))
        .collect();

    let py_rf = PyBytes::new(py, &rf_bytes);

    Ok((names, py_rf.into()))
}

/// Compute pairwise RF distances and export tree snapshots in a single pass.
///
/// Parses each newick once, building both the RF distance matrix and the tree-by-bipartition
/// presence matrix, stored as compressed sparse rows (CSR) so a tree costs one entry per
/// bipartition it holds rather than one byte per bipartition in the collection.
///
/// # Presence matrix format
///
/// The fifth element is a dict with these keys:
///
/// - `format_version`: `1`.
/// - `encoding`: `"csr"`.
/// - `n_entries`: total number of (tree, bipartition) entries.
/// - `row_offsets`: native-endian `uint64` bytes, `n_trees + 1` values.
/// - `column_indices`: native-endian `uint32` bytes, `n_entries` values.
///
/// Tree `i` holds the columns `column_indices[row_offsets[i]:row_offsets[i + 1]]`, ascending and in the column order of `bipartition_clade_bytes`, which is stable across calls on the same tree set. Reconstruct the dense matrix on the Python side:
/// ```python
/// offsets = np.frombuffer(sparse["row_offsets"], dtype=np.uint64)
/// columns = np.frombuffer(sparse["column_indices"], dtype=np.uint32)
/// presence = np.zeros((n_trees, n_bip), dtype=np.uint8)
/// for i in range(n_trees):
///     presence[i, columns[offsets[i]:offsets[i + 1]]] = 1
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
///     6-tuple (tree_names, rf_matrix_bytes, leaf_names, n_bipartitions, sparse_presence,
///     bipartition_clade_bytes).
///
/// Raises:
///     ValueError: If fewer than 2 trees, leaf sets differ, or argument lengths mismatch.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, translate_maps, map_indices, rooted=false, progress=None))]
pub(super) fn pairwise_rf_with_snapshots_from_newick_iter(
    py: Python<'_>,
    names: Vec<String>,
    newick_iter: Bound<'_, PyIterator>,
    translate_maps: Vec<HashMap<String, String>>,
    map_indices: Vec<usize>,
    rooted: bool,
    progress: Option<Py<ProgressCounter>>,
) -> PyResult<PyRfSnapshotResult> {
    let input = IterInput {
        newick_iter,
        translate_maps: &translate_maps,
        map_indices: &map_indices,
        rooted,
    };
    rf_with_snapshots(py, names, input, progress, |snaps| {
        let (row_offsets, column_indices, col_to_bip_id) = snaps.build_sparse_presence_matrix();
        let sparse = PyDict::new(py);
        sparse.set_item("format_version", 1)?;
        sparse.set_item("encoding", "csr")?;
        sparse.set_item("n_entries", column_indices.len() / size_of::<u32>())?;
        sparse.set_item("row_offsets", PyBytes::new(py, &row_offsets))?;
        sparse.set_item("column_indices", PyBytes::new(py, &column_indices))?;
        Ok((sparse.into_any().unbind(), col_to_bip_id))
    })
}
