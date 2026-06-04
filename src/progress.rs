//! Progress reporting for long-running pairwise computations.
//!
//! A single Python-facing handle, [`ProgressCounter`], that the caller
//! instantiates and passes to a pairwise function. Internally it wraps an
//! `Arc<AtomicUsize>` that rayon workers bump as they finish rows; Python
//! reads `.value()` / `.total()` / `.fraction()` on its own cadence
//! (typically from another thread while the pairwise call blocks).
//!
//! Reading is a single atomic load (~1 ns). There is no callback, no
//! monitor thread, no GIL acquisition from Rust — the caller decides when
//! and how often to look at the counter.
//!
//! When no counter is supplied the GIL stays held and the call has zero
//! progress-related overhead.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use pyo3::prelude::*;

/// Thread-safe pair-completion counter that Python code can read on its own
/// schedule while a Rust pairwise computation runs in another thread.
///
/// Construct with `rapidtrees.ProgressCounter()`, pass to any of the
/// `pairwise_*_from_newick_iter` functions as `progress=...`. The function
/// (a) resets the counter and sets its `total` to `n*(n-1)/2` at entry,
/// (b) atomically increments it as rayon workers finish rows, and
/// (c) releases the GIL while the rayon loop runs so a polling Python thread
/// can keep reading the counter.
///
/// Reading is lock-free: `.value()` is a single atomic load (~1 ns).
#[pyclass(module = "rapidtrees")]
#[derive(Clone)]
pub struct ProgressCounter {
    /// Pairs completed so far. Bumped by rayon workers via `fetch_add`.
    value: Arc<AtomicUsize>,
    /// `n*(n-1)/2` for the current call (0 before any call is in flight).
    total: Arc<AtomicUsize>,
}

#[pymethods]
impl ProgressCounter {
    /// Construct a fresh counter at `value=0`, `total=0`.
    #[new]
    fn new() -> Self {
        Self {
            value: Arc::new(AtomicUsize::new(0)),
            total: Arc::new(AtomicUsize::new(0)),
        }
    }

    /// Pairs completed so far (lock-free atomic load).
    fn value(&self) -> usize {
        self.value.load(Ordering::Relaxed)
    }

    /// Total pairs the in-flight call will compute (``n*(n-1)/2``), or ``0``
    /// before any call has started.
    fn total(&self) -> usize {
        self.total.load(Ordering::Relaxed)
    }

    /// ``value() / total()`` clamped to ``[0.0, 1.0]``. Returns ``0.0`` if
    /// no call has started yet (``total == 0``).
    fn fraction(&self) -> f64 {
        frac_done(self.value(), self.total())
    }

    /// Reset `value` and `total` to ``0``. Called automatically at the start
    /// of every pairwise call that takes this counter; exposed so consumers
    /// can also reset between manual runs.
    fn reset(&self) {
        self.value.store(0, Ordering::Relaxed);
        self.total.store(0, Ordering::Relaxed);
    }

    fn __repr__(&self) -> String {
        format!(
            "ProgressCounter(value={}, total={})",
            self.value(),
            self.total()
        )
    }
}

impl ProgressCounter {
    /// Take a snapshot of the internal atomics so they can be shared with
    /// rayon worker threads without holding a `Py<…>` reference.
    pub(crate) fn arc_value(&self) -> Arc<AtomicUsize> {
        Arc::clone(&self.value)
    }

    /// Initialise the counter for a new pairwise call.
    pub(crate) fn begin(&self, total: usize) {
        self.value.store(0, Ordering::Relaxed);
        self.total.store(total, Ordering::Relaxed);
    }
}

/// Convert a (done, total) pair into a fraction in `[0.0, 1.0]`.
///
/// Returns `0.0` when `total == 0` (avoids a divide-by-zero) and clamps the
/// result so a counter that briefly overshoots `total` does not yield a
/// fraction above `1.0`.
#[inline]
fn frac_done(done: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        (done as f64 / total as f64).clamp(0.0, 1.0)
    }
}

/// Run `f` while exposing live progress through an optional
/// [`ProgressCounter`].
///
/// `f` receives a reference to the [`AtomicUsize`] that it must increment as
/// work units complete. When `counter` is `Some`, that atomic is shared with
/// the user's `ProgressCounter` (so `.value()` reads the live count) and the
/// GIL is released for the duration of `f` — without this a Python polling
/// thread couldn't run.
///
/// When `counter` is `None` we take the fast path: no Arc, no GIL release,
/// same cost as a plain function call.
pub(crate) fn with_counter<R, F>(
    py: Python<'_>,
    counter: Option<Py<ProgressCounter>>,
    total: usize,
    f: F,
) -> PyResult<R>
where
    R: Send,
    F: Send + FnOnce(&AtomicUsize) -> R,
{
    let Some(pc) = counter else {
        // Fast path: GIL stays held, single stack-allocated counter.
        let local = AtomicUsize::new(0);
        return Ok(f(&local));
    };

    // Share the user's atomic so .value() reflects what the rayon workers see.
    let work_counter = {
        let pc_ref = pc.borrow(py);
        pc_ref.begin(total);
        pc_ref.arc_value()
    };

    let result = py.detach(|| f(&work_counter));

    // Pin the final count at exactly `total` so callers polling for
    // `value == total` (or `fraction == 1.0`) always observe completion.
    work_counter.store(total, Ordering::Relaxed);

    Ok(result)
}

#[cfg(test)]
mod tests {
    //! Unit tests for the pure-Rust portions of the progress helper.
    //!
    //! `with_progress` itself takes a `Python<'_>` and a Python callback, so
    //! its end-to-end behaviour is exercised by the Python test suite
    //! (`tests/test_python_api.py::TestProgressCallback`). These tests cover
    //! the GIL-free pieces — the fraction calculation and the row-level
    //! counter increments the dense `pairwise_*` paths perform — which are
    //! the only places the progress code can silently go wrong without
    //! Python ever noticing.

    use super::frac_done;
    use crate::snapshot::Snapshots;
    use std::collections::HashMap;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn frac_done_zero_total_is_zero() {
        // Empty work set → reporting 0% is the only sensible answer (and avoids
        // the divide-by-zero a naive implementation would produce).
        assert_eq!(frac_done(0, 0), 0.0);
        assert_eq!(frac_done(7, 0), 0.0);
    }

    #[test]
    fn frac_done_endpoints() {
        assert_eq!(frac_done(0, 10), 0.0);
        assert_eq!(frac_done(10, 10), 1.0);
    }

    #[test]
    fn frac_done_midpoint() {
        assert!((frac_done(5, 10) - 0.5).abs() < 1e-12);
        assert!((frac_done(3, 4) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn frac_done_clamps_overshoot() {
        // The monitor briefly seeing `done > total` (a race against the rayon
        // workers finalising the last row) must never produce a fraction
        // above 1.0 — UIs commonly assert this and crash otherwise.
        assert_eq!(frac_done(11, 10), 1.0);
        assert_eq!(frac_done(usize::MAX, 10), 1.0);
    }

    fn small_snapshot_set() -> Snapshots {
        // Six trees over the same five taxa; differences are immaterial — we
        // only need the pairwise computation to actually visit every upper-
        // triangle cell so we can check the counter.
        let newicks: Vec<&str> = vec![
            "((A:0.1,B:0.1):0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1):0.1);",
            "((A:0.1,C:0.1):0.1,(B:0.1,(D:0.1,E:0.1):0.1):0.1);",
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,(D:0.1,E:0.1):0.1);",
            "(A:0.1,((B:0.1,D:0.1):0.1,(C:0.1,E:0.1):0.1):0.1);",
            "((A:0.1,(B:0.1,C:0.1):0.1):0.1,(D:0.1,E:0.1):0.1);",
        ];
        let empty_map: HashMap<String, String> = HashMap::new();
        let entries = newicks.iter().map(|n| (*n, &empty_map)).collect::<Vec<_>>();
        Snapshots::from_newick_iter(entries, false).expect("fixture snapshots should parse")
    }

    #[test]
    fn counter_reaches_total_pairs() {
        // Six trees → 15 upper-triangle pairs. The counter must equal that
        // exact value once the parallel loop returns, otherwise a progress
        // callback would never see frac == 1.0 before the final manual fire.
        // Exercised on the KF path, whose dense sweep shares the same
        // row-level counter increment as RF and WRF.
        let snaps = small_snapshot_set();
        let n = 6;
        let expected = n * (n - 1) / 2;
        let counter = AtomicUsize::new(0);
        let _ = snaps.pairwise_kf(Some(&counter));
        assert_eq!(counter.load(Ordering::Relaxed), expected);
    }

    #[test]
    fn counter_none_matches_default_path() {
        // Regression guard: passing a counter must not change the matrix the
        // uncounted path produces. If someone "optimises" the counted variant
        // by skipping a row when the counter is absent, this catches it.
        let snaps = small_snapshot_set();
        let counter = AtomicUsize::new(0);
        let with_counter = snaps.pairwise_kf(Some(&counter));
        let without = snaps.pairwise_kf(None);
        assert_eq!(with_counter, without);
    }

    #[test]
    fn pairwise_rf_counted_matches_uncounted_and_increments() {
        // Exercises `pairwise_rf(Some(counter))` end-to-end: the matrix must be
        // identical to the uncounted `pairwise_rf(None)`, and the counter must
        // land on exactly n*(n-1)/2.
        let snaps = small_snapshot_set();
        let n = 6;
        let counter = AtomicUsize::new(0);
        let counted = snaps.pairwise_rf(Some(&counter));
        let uncounted = snaps.pairwise_rf(None);
        assert_eq!(
            counted, uncounted,
            "pairwise_rf(Some(..)) must match pairwise_rf(None)"
        );
        assert_eq!(counter.load(Ordering::Relaxed), n * (n - 1) / 2);
    }

    #[test]
    fn pairwise_wrf_counted_matches_uncounted_and_increments() {
        // Same contract for WRF — confirms the f64 metric path also hits the
        // `counter.fetch_add` line in the dense `fill_symmetric` sweep.
        let snaps = small_snapshot_set();
        let n = 6;
        let counter = AtomicUsize::new(0);
        let counted = snaps.pairwise_wrf(Some(&counter));
        let uncounted = snaps.pairwise_wrf(None);
        assert_eq!(counted, uncounted);
        assert_eq!(counter.load(Ordering::Relaxed), n * (n - 1) / 2);
    }

    #[test]
    fn pairwise_kf_counted_matches_uncounted_and_increments() {
        // Same contract for KF.
        let snaps = small_snapshot_set();
        let n = 6;
        let counter = AtomicUsize::new(0);
        let counted = snaps.pairwise_kf(Some(&counter));
        let uncounted = snaps.pairwise_kf(None);
        assert_eq!(counted, uncounted);
        assert_eq!(counter.load(Ordering::Relaxed), n * (n - 1) / 2);
    }

    // -- ProgressCounter (pure-Rust portions) ---------------------------------

    use super::ProgressCounter;

    #[test]
    fn progress_counter_new_is_zero() {
        let pc = ProgressCounter::new();
        assert_eq!(pc.value(), 0);
        assert_eq!(pc.total(), 0);
        assert_eq!(pc.fraction(), 0.0);
    }

    #[test]
    fn progress_counter_begin_sets_total_and_resets_value() {
        let pc = ProgressCounter::new();
        // Simulate a prior call leaving the value non-zero.
        pc.arc_value().store(123, Ordering::Relaxed);
        pc.begin(500);
        assert_eq!(pc.value(), 0);
        assert_eq!(pc.total(), 500);
        assert_eq!(pc.fraction(), 0.0);
    }

    #[test]
    fn progress_counter_fraction_endpoints_and_overshoot() {
        let pc = ProgressCounter::new();
        pc.begin(10);
        assert_eq!(pc.fraction(), 0.0);
        pc.arc_value().store(5, Ordering::Relaxed);
        assert!((pc.fraction() - 0.5).abs() < 1e-12);
        pc.arc_value().store(10, Ordering::Relaxed);
        assert_eq!(pc.fraction(), 1.0);
        // Overshoot must still report 1.0 (matches the callback monitor's
        // clamping behaviour — UIs assume the fraction never exceeds 1.0).
        pc.arc_value().store(99, Ordering::Relaxed);
        assert_eq!(pc.fraction(), 1.0);
    }

    #[test]
    fn progress_counter_reset_clears_both_fields() {
        let pc = ProgressCounter::new();
        pc.begin(42);
        pc.arc_value().store(7, Ordering::Relaxed);
        pc.reset();
        assert_eq!(pc.value(), 0);
        assert_eq!(pc.total(), 0);
    }

    #[test]
    fn progress_counter_arc_value_is_shared() {
        // Cloning the inner Arc must give us a handle to the SAME atomic that
        // ProgressCounter's `.value()` reads. This is the invariant
        // `with_counter` relies on to make rayon writes visible to Python.
        let pc = ProgressCounter::new();
        pc.begin(100);
        let shared = pc.arc_value();
        shared.fetch_add(33, Ordering::Relaxed);
        assert_eq!(pc.value(), 33);
        assert!((pc.fraction() - 0.33).abs() < 1e-12);
    }
}
