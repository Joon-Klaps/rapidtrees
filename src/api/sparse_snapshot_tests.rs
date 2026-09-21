//! Embedded-Python smoke tests for the sparse-snapshot endpoint.

use super::py_integration_tests::{call_pairwise, ensure_python, fresh_counter};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

fn dict_item<'py>(dict: &Bound<'py, PyDict>, key: &str) -> Bound<'py, PyAny> {
    dict.get_item(key)
        .unwrap_or_else(|error| panic!("failed to read sparse-snapshot key {key:?}: {error}"))
        .unwrap_or_else(|| panic!("sparse-snapshot dictionary is missing {key:?}"))
}

fn decode_u32(bytes: &[u8]) -> Vec<u32> {
    let (values, remainder) = bytes.as_chunks::<4>();
    assert!(remainder.is_empty());
    values.iter().copied().map(u32::from_ne_bytes).collect()
}

fn decode_u64(bytes: &[u8]) -> Vec<u64> {
    let (values, remainder) = bytes.as_chunks::<8>();
    assert!(remainder.is_empty());
    values.iter().copied().map(u64::from_ne_bytes).collect()
}

#[test]
fn sparse_snapshot_endpoint_is_registered() {
    ensure_python();
    Python::attach(|py| {
        let module = py.import("rapidtrees").unwrap();
        assert!(
            module
                .hasattr("pairwise_rf_with_sparse_snapshots_from_newick_iter")
                .unwrap()
        );
    });
}

#[test]
fn sparse_snapshot_endpoint_returns_versioned_csr_buffers() {
    ensure_python();
    Python::attach(|py| {
        let progress = fresh_counter(py);
        let result = call_pairwise(
            py,
            "pairwise_rf_with_sparse_snapshots_from_newick_iter",
            Some(&progress),
        );
        let output = result.cast::<PyTuple>().expect("six-element result tuple");
        assert_eq!(output.len(), 6);

        let tree_names: Vec<String> = output.get_item(0).unwrap().extract().unwrap();
        let rf_bytes: Vec<u8> = output.get_item(1).unwrap().extract().unwrap();
        let leaf_names: Vec<String> = output.get_item(2).unwrap().extract().unwrap();
        let n_clades: usize = output.get_item(3).unwrap().extract().unwrap();
        let clade_bytes: Vec<u8> = output.get_item(4).unwrap().extract().unwrap();
        let sparse_object = output.get_item(5).unwrap();
        let sparse = sparse_object
            .cast::<PyDict>()
            .expect("sparse snapshot must be a dictionary");

        assert_eq!(tree_names, ["t0", "t1", "t2"]);
        assert_eq!(leaf_names, ["A", "B", "C"]);
        assert_eq!(rf_bytes.len(), 3 * 3 * size_of::<u32>());
        assert_eq!(clade_bytes.len(), n_clades * leaf_names.len().div_ceil(8));
        assert_eq!(sparse.len(), 5);

        let version: u8 = dict_item(sparse, "format_version").extract().unwrap();
        let encoding: String = dict_item(sparse, "encoding").extract().unwrap();
        let n_entries: usize = dict_item(sparse, "n_entries").extract().unwrap();
        let offset_bytes: Vec<u8> = dict_item(sparse, "row_offsets").extract().unwrap();
        let column_bytes: Vec<u8> = dict_item(sparse, "column_indices").extract().unwrap();
        let offsets = decode_u64(&offset_bytes);
        let columns = decode_u32(&column_bytes);

        assert_eq!(version, 1);
        assert_eq!(encoding, "csr");
        assert_eq!(offsets.len(), tree_names.len() + 1);
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets.last().copied(), Some(n_entries as u64));
        assert_eq!(columns.len(), n_entries);
        assert!(columns.iter().all(|&column| column < n_clades as u32));
        for row in offsets.windows(2) {
            let values = &columns[row[0] as usize..row[1] as usize];
            assert!(values.windows(2).all(|pair| pair[0] < pair[1]));
        }

        let fraction: f64 = progress
            .call_method0("fraction")
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(fraction, 1.0);
    });
}
