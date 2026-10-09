//! Python binding layer for tree distance calculations.
//!
//! All computation lives in [`crate::distances`] and `crate::snapshot`; this module
//! contains only the PyO3-specific glue: iterating Python iterators, wrapping
//! byte buffers into `PyBytes`, and registering functions in the Python module.
//!
//! - `input`: argument validation and pulling snapshots out of a Python iterator.
//! - `rf`: the Robinson–Foulds entry points.
//! - `weighted`: the weighted RF and Kuhner–Felsenstein entry points.
//! - `progress`: `ProgressCounter`, the live progress handle Python polls.

mod input;
mod progress;
mod rf;
mod weighted;

#[cfg(test)]
mod py_tests;

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};

use crate::snapshot::Snapshots;
use progress::ProgressCounter;

/// Number of upper-triangle pairs for `n` trees — what a [`ProgressCounter`]
/// counts up to as "100% done".
#[inline]
pub(super) fn n_pairs(n: usize) -> usize {
    n.saturating_mul(n.saturating_sub(1)) / 2
}

pub(super) type PyRfSnapshotResult = (
    Vec<String>,
    Py<PyAny>,
    Vec<String>,
    usize,
    Py<PyAny>,
    Py<PyAny>,
);

/// The CSR dict a snapshot entry point returns fifth, with the column order its clade bitmasks must follow. `lengths` is a key only when `with_lengths` is set.
pub(super) fn csr_dict(
    py: Python<'_>,
    snaps: &Snapshots,
    with_lengths: bool,
) -> PyResult<(Py<PyAny>, Vec<usize>)> {
    let (row_offsets, column_indices, lengths, col_to_bip_id) = snaps.build_csr(with_lengths);
    let sparse = PyDict::new(py);
    sparse.set_item("format_version", 1)?;
    sparse.set_item("encoding", "csr")?;
    sparse.set_item("n_entries", column_indices.len() / size_of::<u32>())?;
    sparse.set_item("row_offsets", PyBytes::new(py, &row_offsets))?;
    sparse.set_item("column_indices", PyBytes::new(py, &column_indices))?;
    if with_lengths {
        sparse.set_item("lengths", PyBytes::new(py, &lengths))?;
    }
    Ok((sparse.into_any().unbind(), col_to_bip_id))
}

/// Python module definition
#[pymodule]
fn rapidtrees(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rf::pairwise_rf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(
        rf::pairwise_rf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        weighted::pairwise_wrf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        weighted::pairwise_kf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        weighted::pairwise_wrf_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(weighted::pairwise_kf_from_newick_iter, m)?)?;
    m.add_class::<ProgressCounter>()?;
    Ok(())
}
