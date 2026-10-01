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

/// Python module definition
#[pymodule]
fn rapidtrees(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rf::pairwise_rf_from_newick_iter, m)?)?;
    m.add_function(wrap_pyfunction!(
        rf::pairwise_rf_with_snapshots_from_newick_iter,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        rf::pairwise_rf_with_sparse_snapshots_from_newick_iter,
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
