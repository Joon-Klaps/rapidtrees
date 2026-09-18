//! Embedded-Python smoke tests for the rooted-facts endpoint.

use super::py_integration_tests::{call_pairwise, ensure_python, fresh_counter};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

fn dict_item<'py>(dict: &Bound<'py, PyDict>, key: &str) -> Bound<'py, PyAny> {
    dict.get_item(key)
        .unwrap_or_else(|error| panic!("failed to read rooted-facts key {key:?}: {error}"))
        .unwrap_or_else(|| panic!("rooted-facts dictionary is missing {key:?}"))
}

fn decode_u32(bytes: &[u8]) -> Vec<u32> {
    let (values, remainder) = bytes.as_chunks::<4>();
    assert!(remainder.is_empty());
    values.iter().copied().map(u32::from_ne_bytes).collect()
}

fn decode_f64(bytes: &[u8]) -> Vec<f64> {
    let (values, remainder) = bytes.as_chunks::<8>();
    assert!(remainder.is_empty());
    values.iter().copied().map(f64::from_ne_bytes).collect()
}

#[test]
fn rooted_facts_endpoint_is_registered() {
    ensure_python();
    Python::attach(|py| {
        let module = py.import("rapidtrees").unwrap();
        assert!(
            module
                .hasattr("pairwise_rf_with_rooted_facts_from_newick_iter")
                .unwrap()
        );
    });
}

#[test]
fn rooted_facts_endpoint_returns_versioned_fixed_width_buffers() {
    ensure_python();
    Python::attach(|py| {
        let progress = fresh_counter(py);
        let result = call_pairwise(
            py,
            "pairwise_rf_with_rooted_facts_from_newick_iter",
            Some(&progress),
        );
        let output = result
            .cast::<PyTuple>()
            .expect("seven-element result tuple");
        assert_eq!(output.len(), 7);

        let tree_names: Vec<String> = output.get_item(0).unwrap().extract().unwrap();
        let rf_bytes: Vec<u8> = output.get_item(1).unwrap().extract().unwrap();
        let leaf_names: Vec<String> = output.get_item(2).unwrap().extract().unwrap();
        let n_clades: usize = output.get_item(3).unwrap().extract().unwrap();
        let presence: Vec<u8> = output.get_item(4).unwrap().extract().unwrap();
        let clade_bytes: Vec<u8> = output.get_item(5).unwrap().extract().unwrap();
        let facts_object = output.get_item(6).unwrap();
        let facts = facts_object
            .cast::<PyDict>()
            .expect("rooted facts must be a dictionary");

        assert_eq!(tree_names, ["t0", "t1", "t2"]);
        assert_eq!(leaf_names, ["A", "B", "C"]);
        assert_eq!(n_clades, 5);
        assert_eq!(rf_bytes.len(), 3 * 3 * size_of::<u32>());
        assert_eq!(presence.len(), 3 * n_clades);
        assert_eq!(clade_bytes.len(), n_clades * leaf_names.len().div_ceil(8));
        assert!(presence.chunks_exact(n_clades).all(|row| {
            row.iter()
                .map(|&present| usize::from(present))
                .sum::<usize>()
                == 4
        }));

        assert_eq!(facts.len(), 8);
        for key in [
            "format_version",
            "root_column",
            "nodes_per_tree",
            "splits_per_tree",
            "node_columns",
            "node_heights",
            "root_heights",
            "split_columns",
        ] {
            assert!(facts.contains(key).unwrap(), "missing {key:?}");
        }

        let format_version: u8 = dict_item(facts, "format_version").extract().unwrap();
        let root_column: u32 = dict_item(facts, "root_column").extract().unwrap();
        let nodes_per_tree: usize = dict_item(facts, "nodes_per_tree").extract().unwrap();
        let splits_per_tree: usize = dict_item(facts, "splits_per_tree").extract().unwrap();
        let node_column_bytes: Vec<u8> = dict_item(facts, "node_columns").extract().unwrap();
        let node_height_bytes: Vec<u8> = dict_item(facts, "node_heights").extract().unwrap();
        let root_height_bytes: Vec<u8> = dict_item(facts, "root_heights").extract().unwrap();
        let split_column_bytes: Vec<u8> = dict_item(facts, "split_columns").extract().unwrap();

        assert_eq!(format_version, 1);
        assert_eq!(root_column, n_clades as u32);
        assert_eq!(nodes_per_tree, 4);
        assert_eq!(splits_per_tree, 2);
        assert_eq!(
            node_column_bytes.len(),
            3 * nodes_per_tree * size_of::<u32>()
        );
        assert_eq!(
            node_height_bytes.len(),
            3 * nodes_per_tree * size_of::<f64>()
        );
        assert_eq!(root_height_bytes.len(), 3 * size_of::<f64>());
        assert_eq!(
            split_column_bytes.len(),
            3 * splits_per_tree * 3 * size_of::<u32>()
        );

        let node_columns = decode_u32(&node_column_bytes);
        assert!(node_columns.iter().all(|&column| column < root_column));

        let node_heights = decode_f64(&node_height_bytes);
        assert!(node_heights.iter().all(|height| height.is_finite()));
        let root_heights = decode_f64(&root_height_bytes);
        assert!(
            root_heights
                .iter()
                .all(|height| (height - 0.2).abs() < 1e-12)
        );

        let split_columns = decode_u32(&split_column_bytes);
        let triples = split_columns.chunks_exact(3).collect::<Vec<_>>();
        assert_eq!(
            triples
                .iter()
                .filter(|triple| triple[0] == root_column)
                .count(),
            3
        );
        assert!(triples.iter().all(|triple| {
            triple[0] <= root_column && triple[1] < triple[2] && triple[2] < root_column
        }));

        let fraction: f64 = progress
            .call_method0("fraction")
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(fraction, 1.0);
    });
}
