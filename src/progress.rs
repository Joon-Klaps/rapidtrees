//! Progress reporting for long-running pairwise computations.
//!
//! Bridges a thread-safe `AtomicUsize` work counter (incremented from rayon
//! worker threads) and an optional Python callback. A dedicated monitor
//! thread polls the counter on a fixed interval and invokes the callback
//! with a fraction in `[0.0, 1.0]`. The main computation thread releases the
//! GIL with [`pyo3::Python::detach`] so the monitor can re-acquire it
//! independently.
//!
//! When no callback is supplied, [`with_progress`] short-circuits and runs
//! the work with a throwaway counter and no GIL release — preserving the
//! original zero-overhead path.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::thread;
use std::time::Duration;

use pyo3::prelude::*;

/// How often the monitor thread wakes up to poll the counter.
const POLL_INTERVAL: Duration = Duration::from_millis(100);

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

/// Run `f` while reporting progress to an optional Python `callback`.
///
/// `f` receives a reference to an [`AtomicUsize`] that it must increment as
/// work units complete. `total` is the count that corresponds to "100% done"
/// — the callback receives `done / total` clamped to `[0.0, 1.0]`.
///
/// Behaviour:
/// - **`callback == None`**: `f` is invoked directly with a dummy counter.
///   No monitor thread, no GIL release, no callbacks fired. Zero overhead.
/// - **`callback == Some(cb)`**:
///   1. A monitor thread is spawned that polls every `POLL_INTERVAL`.
///   2. `f` runs inside `py.allow_threads(|| ...)` so the monitor can grab
///      the GIL and invoke `cb(frac)` periodically.
///   3. The monitor also calls `Python::check_signals()`, so a Python
///      `KeyboardInterrupt` propagates through the callback as a `PyErr`.
///   4. After `f` returns, a final `cb(1.0)` is fired and any error from
///      the callback (or `check_signals`) is propagated as the function's
///      `Err` return.
pub(crate) fn with_progress<R, F>(
    py: Python<'_>,
    callback: Option<Py<PyAny>>,
    total: usize,
    f: F,
) -> PyResult<R>
where
    R: Send,
    F: Send + FnOnce(&AtomicUsize) -> R,
{
    // Fast path: no callback supplied → run inline, no monitor, no GIL release.
    let Some(cb) = callback else {
        let counter = AtomicUsize::new(0);
        return Ok(f(&counter));
    };

    let counter = Arc::new(AtomicUsize::new(0));
    let terminate = Arc::new(AtomicBool::new(false));
    // Shared slot for an error raised by the Python callback (or check_signals).
    let cb_error: Arc<std::sync::Mutex<Option<PyErr>>> = Arc::new(std::sync::Mutex::new(None));

    let monitor_counter = Arc::clone(&counter);
    let monitor_terminate = Arc::clone(&terminate);
    let monitor_error = Arc::clone(&cb_error);
    let monitor_cb = cb.clone_ref(py);

    let monitor = thread::spawn(move || {
        let mut last_frac: f64 = -1.0;
        loop {
            thread::sleep(POLL_INTERVAL);
            if monitor_terminate.load(Ordering::Relaxed) {
                return;
            }
            let done = monitor_counter.load(Ordering::Relaxed);
            let frac = frac_done(done, total);
            if frac == last_frac {
                continue;
            }
            last_frac = frac;

            let outcome: PyResult<()> = Python::attach(|py| {
                monitor_cb.call1(py, (frac,))?;
                py.check_signals()?;
                Ok(())
            });
            if let Err(err) = outcome {
                *monitor_error.lock().expect("progress mutex poisoned") = Some(err);
                return;
            }
        }
    });

    let result = py.detach(|| f(&counter));

    terminate.store(true, Ordering::Relaxed);
    let _ = monitor.join();

    // Surface any error raised inside the monitor (callback or KeyboardInterrupt).
    if let Some(err) = cb_error.lock().expect("progress mutex poisoned").take() {
        return Err(err);
    }

    // Final callback at exactly 1.0 so consumers can finalise their progress UI.
    cb.call1(py, (1.0_f64,))?;
    py.check_signals()?;

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
    //! counter increments used by `pairwise_symmetric_counted` — which are
    //! the only places the progress code can silently go wrong without
    //! Python ever noticing.

    use super::frac_done;
    use crate::distances::pairwise_symmetric_counted;
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
        let snaps = small_snapshot_set();
        let n = 6;
        let expected = n * (n - 1) / 2;
        let counter = AtomicUsize::new(0);
        let _ =
            pairwise_symmetric_counted(&snaps, crate::distances::rf_distance_fast, Some(&counter));
        assert_eq!(counter.load(Ordering::Relaxed), expected);
    }

    #[test]
    fn counter_none_matches_default_path() {
        // Regression guard: passing `None` for the counter must produce the
        // same matrix as the original (uncounted) `pairwise_symmetric`. If
        // someone "optimises" the counted variant by skipping a row when the
        // counter is absent, this catches it.
        let snaps = small_snapshot_set();
        let with_none: Vec<usize> =
            pairwise_symmetric_counted(&snaps, crate::distances::rf_distance_fast, None);
        let with_counter: Vec<usize> = {
            let counter = AtomicUsize::new(0);
            pairwise_symmetric_counted(&snaps, crate::distances::rf_distance_fast, Some(&counter))
        };
        assert_eq!(with_none, with_counter);
    }

    #[test]
    fn pairwise_rf_counted_matches_uncounted_and_increments() {
        // Exercises `Snapshots::pairwise_rf_counted` end-to-end: the matrix
        // must be identical to the uncounted `pairwise_rf`, and the counter
        // must land on exactly n*(n-1)/2.
        let snaps = small_snapshot_set();
        let n = 6;
        let counter = AtomicUsize::new(0);
        let counted = snaps.pairwise_rf_counted(&counter);
        let uncounted = snaps.pairwise_rf();
        assert_eq!(
            counted, uncounted,
            "pairwise_rf_counted must produce the same matrix as pairwise_rf"
        );
        assert_eq!(counter.load(Ordering::Relaxed), n * (n - 1) / 2);
    }

    #[test]
    fn pairwise_wrf_counted_matches_uncounted_and_increments() {
        // Same contract for WRF — confirms the f64 metric path also hits the
        // `counter.fetch_add` line in `pairwise_symmetric_counted`.
        let snaps = small_snapshot_set();
        let n = 6;
        let counter = AtomicUsize::new(0);
        let counted = snaps.pairwise_wrf_counted(&counter);
        let uncounted = snaps.pairwise_wrf();
        assert_eq!(counted, uncounted);
        assert_eq!(counter.load(Ordering::Relaxed), n * (n - 1) / 2);
    }

    #[test]
    fn pairwise_kf_counted_matches_uncounted_and_increments() {
        // Same contract for KF.
        let snaps = small_snapshot_set();
        let n = 6;
        let counter = AtomicUsize::new(0);
        let counted = snaps.pairwise_kf_counted(&counter);
        let uncounted = snaps.pairwise_kf();
        assert_eq!(counted, uncounted);
        assert_eq!(counter.load(Ordering::Relaxed), n * (n - 1) / 2);
    }
}
