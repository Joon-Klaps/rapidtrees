//! Embedded-Python integration tests that exercise every PyO3 entry point
//! from Rust.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::sync::Once;

// Bring the `#[pymodule] fn rapidtrees` from the parent module into scope
// so the `append_to_inittab!` macro can find it by name.
use super::rapidtrees;

/// Initialise the embedded interpreter exactly once across all tests
/// running in this process.
fn ensure_python() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        pyo3::append_to_inittab!(rapidtrees);
        Python::initialize();
    });
}

/// Three small newicks with a known RF profile, ready to hand to a
/// pairwise function from Python.
fn fixture<'py>(py: Python<'py>) -> (Bound<'py, PyList>, Bound<'py, PyList>) {
    let trees = PyList::new(
        py,
        [
            "(A:0.1,(B:0.1,C:0.1):0.1);",
            "(A:0.1,(C:0.1,B:0.1):0.1);",
            "((A:0.1,B:0.1):0.1,C:0.1);",
        ],
    )
    .unwrap();
    let names = PyList::new(py, ["t0", "t1", "t2"]).unwrap();
    (names, trees)
}

fn call_pairwise<'py>(
    py: Python<'py>,
    func_name: &str,
    progress: Option<&Bound<'py, PyAny>>,
) -> Bound<'py, PyAny> {
    let rapidtrees = py.import("rapidtrees").expect("import rapidtrees");
    let func = rapidtrees.getattr(func_name).expect(func_name);
    let (names, trees) = fixture(py);
    let translate_maps = PyList::new(py, [PyDict::new(py)]).unwrap();
    let map_indices = PyList::new(py, [0i64, 0, 0]).unwrap();
    let kwargs = PyDict::new(py);
    if let Some(pc) = progress {
        kwargs.set_item("progress", pc).unwrap();
    }
    func.call(
        (
            names,
            trees.try_iter().unwrap(),
            translate_maps,
            map_indices,
        ),
        Some(&kwargs),
    )
    .unwrap_or_else(|e| panic!("{func_name} call failed: {e}"))
}

#[test]
fn rapidtrees_module_exports_expected_symbols() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        for name in [
            "pairwise_rf_from_newick_iter",
            "pairwise_wrf_from_newick_iter",
            "pairwise_kf_from_newick_iter",
            "pairwise_rf_with_snapshots_from_newick_iter",
            "pairwise_wrf_with_snapshots_from_newick_iter",
            "pairwise_kf_with_snapshots_from_newick_iter",
            "ProgressCounter",
        ] {
            assert!(
                m.hasattr(name).unwrap(),
                "rapidtrees module is missing `{name}`"
            );
        }
    });
}

/// Build a fresh ProgressCounter for tests that want to drive the
/// counter-bearing branch of `with_counter`. Without this the
/// counter-handling closure for that metric+variant monomorphization
/// stays at FNDA:0 in lcov.
fn fresh_counter<'py>(py: Python<'py>) -> Bound<'py, PyAny> {
    py.import("rapidtrees")
        .unwrap()
        .getattr("ProgressCounter")
        .unwrap()
        .call0()
        .unwrap()
}

#[test]
fn pairwise_rf_via_python_returns_bytes() {
    ensure_python();
    Python::attach(|py| {
        let pc = fresh_counter(py);
        let result = call_pairwise(py, "pairwise_rf_from_newick_iter", Some(&pc));
        let tuple: (Vec<String>, Vec<u8>) = result.extract().unwrap();
        assert_eq!(tuple.0.len(), 3);
        assert_eq!(tuple.1.len(), 3 * 3 * 4); // n*n*sizeof(u32)
        // Also exercises with_counter's `Some` branch + final `store` line.
        let frac: f64 = pc.call_method0("fraction").unwrap().extract().unwrap();
        assert_eq!(frac, 1.0);
    });
}

#[test]
fn pairwise_wrf_via_python() {
    ensure_python();
    Python::attach(|py| {
        let pc = fresh_counter(py);
        let result = call_pairwise(py, "pairwise_wrf_from_newick_iter", Some(&pc));
        let (_names, matrix): (Vec<String>, Vec<f64>) = result.extract().unwrap();
        assert_eq!(matrix.len(), 9);
    });
}

#[test]
fn pairwise_kf_via_python() {
    ensure_python();
    Python::attach(|py| {
        let pc = fresh_counter(py);
        let result = call_pairwise(py, "pairwise_kf_from_newick_iter", Some(&pc));
        let (_names, matrix): (Vec<String>, Vec<f64>) = result.extract().unwrap();
        assert_eq!(matrix.len(), 9);
    });
}

#[test]
fn pairwise_rf_with_snapshots_via_python() {
    ensure_python();
    Python::attach(|py| {
        let pc = fresh_counter(py);
        let _ = call_pairwise(py, "pairwise_rf_with_snapshots_from_newick_iter", Some(&pc));
    });
}

#[test]
fn pairwise_wrf_with_snapshots_via_python() {
    ensure_python();
    Python::attach(|py| {
        let pc = fresh_counter(py);
        let _ = call_pairwise(
            py,
            "pairwise_wrf_with_snapshots_from_newick_iter",
            Some(&pc),
        );
    });
}

#[test]
fn pairwise_kf_with_snapshots_via_python() {
    ensure_python();
    Python::attach(|py| {
        let pc = fresh_counter(py);
        let _ = call_pairwise(py, "pairwise_kf_with_snapshots_from_newick_iter", Some(&pc));
    });
}

/// One call without a counter, exercising the fast path inside
/// `with_counter` (the `None` arm — different code path from the
/// counter-bearing branch the rest of the tests cover).
#[test]
fn pairwise_without_counter_takes_fast_path() {
    ensure_python();
    Python::attach(|py| {
        let result = call_pairwise(py, "pairwise_rf_from_newick_iter", None);
        let tuple: (Vec<String>, Vec<u8>) = result.extract().unwrap();
        assert_eq!(tuple.0.len(), 3);
    });
}

#[test]
fn progress_counter_lifecycle() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        let pc_class = m.getattr("ProgressCounter").unwrap();
        let pc = pc_class.call0().unwrap();

        // Initial state — exercises ProgressCounter::value/total/fraction.
        let zero: usize = pc.call_method0("value").unwrap().extract().unwrap();
        assert_eq!(zero, 0);
        let zero_total: usize = pc.call_method0("total").unwrap().extract().unwrap();
        assert_eq!(zero_total, 0);
        let zero_frac: f64 = pc.call_method0("fraction").unwrap().extract().unwrap();
        assert_eq!(zero_frac, 0.0);

        // __repr__ — codecov was flagging this as uncovered.
        let repr = pc.repr().unwrap().to_string();
        assert!(repr.contains("ProgressCounter"));
        assert!(repr.contains("value=0"));
        assert!(repr.contains("total=0"));

        // Drive the counter through a real pairwise call. After the call
        // value == total and fraction == 1.0.
        let _ = call_pairwise(py, "pairwise_rf_from_newick_iter", Some(&pc));
        let total: usize = pc.call_method0("total").unwrap().extract().unwrap();
        let value: usize = pc.call_method0("value").unwrap().extract().unwrap();
        let frac: f64 = pc.call_method0("fraction").unwrap().extract().unwrap();
        assert_eq!(total, 3); // n*(n-1)/2 for n=3
        assert_eq!(value, total);
        assert_eq!(frac, 1.0);

        // reset() — only line in __pymethod_reset__ that the unit tests
        // alone never hit.
        pc.call_method0("reset").unwrap();
        let after: usize = pc.call_method0("value").unwrap().extract().unwrap();
        let after_total: usize = pc.call_method0("total").unwrap().extract().unwrap();
        assert_eq!(after, 0);
        assert_eq!(after_total, 0);
    });
}

/// A parse/leaf-set failure has to surface as `ValueError` from *every*
/// entry point, not just the RF one — each reaches
/// `collect_snapshots_from_iter` through a different helper.
#[test]
fn bad_trees_raise_value_error_from_every_entry_point() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        // Same leaf count, different taxa: tree 1 has `Z` where tree 0 has `C`.
        let trees = PyList::new(py, ["(A:1,(B:1,C:1):1);", "(A:1,(B:1,Z:1):1);"]).unwrap();
        let names = PyList::new(py, ["t0", "t1"]).unwrap();

        for func_name in [
            "pairwise_rf_from_newick_iter",
            "pairwise_wrf_from_newick_iter",
            "pairwise_kf_from_newick_iter",
            "pairwise_rf_with_snapshots_from_newick_iter",
            "pairwise_wrf_with_snapshots_from_newick_iter",
            "pairwise_kf_with_snapshots_from_newick_iter",
        ] {
            let err = m
                .getattr(func_name)
                .unwrap()
                .call1((
                    names.clone(),
                    trees.clone().try_iter().unwrap(),
                    PyList::new(py, [PyDict::new(py)]).unwrap(),
                    PyList::new(py, [0i64, 0]).unwrap(),
                ))
                .expect_err(&format!("{func_name} accepted mismatched leaf sets"));
            assert!(
                err.is_instance_of::<pyo3::exceptions::PyValueError>(py),
                "{func_name} raised something other than ValueError"
            );
        }
    });
}

/// `validate_iter_args` counts `names`, but the newick iterator is lazy and
/// may yield fewer — so the tree count is re-checked after draining it.
#[test]
fn short_iterator_raises_value_error() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        let func = m.getattr("pairwise_rf_from_newick_iter").unwrap();
        // Two names promised, one tree delivered.
        let err = func
            .call1((
                PyList::new(py, ["t0", "t1"]).unwrap(),
                PyList::new(py, ["(A:1,B:1);"]).unwrap().try_iter().unwrap(),
                PyList::new(py, [PyDict::new(py)]).unwrap(),
                PyList::new(py, [0i64, 0]).unwrap(),
            ))
            .expect_err("expected ValueError when the iterator is shorter than names");
        assert!(err.is_instance_of::<pyo3::exceptions::PyValueError>(py));
    });
}

/// The newick iterator is pulled lazily: a generator streams straight
/// through and gives the matrix a list gives.
#[test]
fn generator_input_matches_list_input() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        let func = m.getattr("pairwise_rf_from_newick_iter").unwrap();
        let (names, trees) = fixture(py);
        let call = |newicks| {
            func.call1((
                names.clone(),
                newicks,
                PyList::new(py, [PyDict::new(py)]).unwrap(),
                PyList::new(py, [0i64, 0, 0]).unwrap(),
            ))
            .unwrap()
        };
        let from_list = call(trees.try_iter().unwrap().into_any());
        let generator = py
            .eval(
                c"(t for t in ['(A:0.1,(B:0.1,C:0.1):0.1);', '(A:0.1,(C:0.1,B:0.1):0.1);', '((A:0.1,B:0.1):0.1,C:0.1);'])",
                None,
                None,
            )
            .unwrap();
        assert!(from_list.eq(call(generator)).unwrap());
    });
}

/// An item that is not a string, and an exception the iterator raises
/// itself, surface as they are rather than as a parse error.
#[test]
fn iterator_errors_propagate_unchanged() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        let func = m.getattr("pairwise_rf_from_newick_iter").unwrap();
        let call = |code: &std::ffi::CStr| {
            func.call1((
                PyList::new(py, ["t0", "t1", "t2"]).unwrap(),
                py.eval(code, None, None).unwrap(),
                PyList::new(py, [PyDict::new(py)]).unwrap(),
                PyList::new(py, [0i64, 0, 0]).unwrap(),
            ))
            .expect_err("a bad iterator was accepted")
        };
        let err = call(c"iter(['(A:1,B:1,C:1);', 7, '(A:1,B:1,C:1);'])");
        assert!(err.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
        let err = call(c"(t if i < 2 else 1 // 0 for i, t in enumerate(['(A:1,B:1,C:1);'] * 3))");
        assert!(err.is_instance_of::<pyo3::exceptions::PyZeroDivisionError>(py));
    });
}

/// Argument validation rejects a `map_indices` entry pointing past the end
/// of `translate_maps`, and a `names` length that disagrees with it.
#[test]
fn argument_shape_mismatches_raise_value_error() {
    ensure_python();
    Python::attach(|py| {
        let m = py.import("rapidtrees").unwrap();
        let func = m.getattr("pairwise_rf_from_newick_iter").unwrap();
        let trees = PyList::new(py, ["(A:1,B:1);", "(A:1,B:1);"]).unwrap();
        let names = PyList::new(py, ["t0", "t1"]).unwrap();

        // map_indices shorter than names.
        let err = func
            .call1((
                names.clone(),
                trees.clone().try_iter().unwrap(),
                PyList::new(py, [PyDict::new(py)]).unwrap(),
                PyList::new(py, [0i64]).unwrap(),
            ))
            .expect_err("expected ValueError for a names/map_indices length mismatch");
        assert!(err.is_instance_of::<pyo3::exceptions::PyValueError>(py));

        // map_indices pointing past the end of translate_maps.
        let err = func
            .call1((
                names,
                trees.try_iter().unwrap(),
                PyList::new(py, [PyDict::new(py)]).unwrap(),
                PyList::new(py, [0i64, 7]).unwrap(),
            ))
            .expect_err("expected ValueError for an out-of-bounds map index");
        assert!(err.is_instance_of::<pyo3::exceptions::PyValueError>(py));
    });
}

#[test]
fn invalid_args_raise_value_error() {
    ensure_python();
    Python::attach(|py| {
        // names < 2 → ValueError from validate_iter_args.
        let m = py.import("rapidtrees").unwrap();
        let func = m.getattr("pairwise_rf_from_newick_iter").unwrap();
        let names = PyList::new(py, ["only-one"]).unwrap();
        let trees = PyList::new(py, ["(A:0.1,B:0.1);"]).unwrap();
        let translate_maps = PyList::new(py, [PyDict::new(py)]).unwrap();
        let map_indices = PyList::new(py, [0i64]).unwrap();
        let err = func
            .call1((
                names,
                trees.try_iter().unwrap(),
                translate_maps,
                map_indices,
            ))
            .expect_err("expected ValueError for n<2");
        assert!(err.is_instance_of::<pyo3::exceptions::PyValueError>(py));
    });
}
