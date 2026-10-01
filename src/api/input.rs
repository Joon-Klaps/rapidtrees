//! The tree source as it arrives from Python, checked and turned into [`Snapshots`].

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use std::collections::HashMap;

use crate::snapshot::{Retain, Snapshots};

/// The tree source as handed in from Python: the newicks, how to rename their
/// taxa, and whether to compare clades or bipartitions.
///
/// Every entry point takes these four together and does nothing with them but
/// forward them, so they travel as one value.
pub(super) struct IterInput<'py, 'a> {
    pub(super) newick_iter: Bound<'py, PyIterator>,
    pub(super) translate_maps: &'a [HashMap<String, String>],
    pub(super) map_indices: &'a [usize],
    pub(super) rooted: bool,
}

/// Collect tree snapshots from a lazy Python iterator of newick strings.
///
/// The iterator is pulled a chunk of trees at a time and each chunk's strings
/// are dropped once parsed, so a generator reading a file streams straight
/// through rather than being collected first. The first item that is not a
/// string, or that the iterator fails to produce, raises its own error.
///
/// `retain` says what to build beyond the split IDs — see [`Retain`]. Each entry
/// point asks only for what it actually returns.
pub(super) fn collect_snapshots_from_iter(
    input: IterInput<'_, '_>,
    retain: Retain,
) -> PyResult<Snapshots> {
    let translate_maps = input.translate_maps;
    let mut failed = None;
    let entries = input
        .newick_iter
        .zip(input.map_indices)
        .map_while(
            |(item, &idx)| match item.and_then(|obj| obj.extract::<String>()) {
                Ok(newick) => Some((newick, &translate_maps[idx])),
                Err(e) => {
                    failed = Some(e);
                    None
                }
            },
        );
    let built = Snapshots::from_newick_iter_opts(entries, input.rooted, retain);
    if let Some(e) = failed {
        return Err(e);
    }
    let snaps = built.map_err(PyValueError::new_err)?;
    if snaps.len() < 2 {
        return Err(PyValueError::new_err(
            "Need at least 2 trees to compute pairwise distances",
        ));
    }
    Ok(snaps)
}

/// Validate argument consistency for iterator-based functions.
pub(super) fn validate_iter_args(
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
