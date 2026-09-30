//! Tests for compressed sparse-row snapshot export.

use super::Snapshots;

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

fn expand_sparse(n_trees: usize, n_clades: usize, offsets: &[u64], columns: &[u32]) -> Vec<u8> {
    let mut dense = vec![0u8; n_trees * n_clades];
    for tree in 0..n_trees {
        for &column in &columns[offsets[tree] as usize..offsets[tree + 1] as usize] {
            dense[tree * n_clades + column as usize] = 1;
        }
    }
    dense
}

#[test]
fn sparse_presence_round_trips_dense_in_both_rooting_modes() {
    let trees = [
        "((A:1,B:1):1,(C:1,D:1):1);",
        "((A:1,C:1):1,(B:1,D:1):1);",
        "((A:1,D:1):1,(B:1,C:1):1);",
    ];

    for rooted in [false, true] {
        let snapshots = Snapshots::from_newicks(&trees, rooted).unwrap();
        let (dense, dense_columns) = snapshots.build_presence_matrix();
        let (sparse, sparse_columns) = snapshots.build_sparse_presence_matrix().unwrap();
        let offsets = decode_u64(&sparse.row_offsets);
        let columns = decode_u32(&sparse.column_indices);

        assert_eq!(sparse_columns, dense_columns);
        assert_eq!(offsets.len(), trees.len() + 1);
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets.last().copied(), Some(sparse.n_entries as u64));
        assert_eq!(columns.len(), sparse.n_entries);
        assert!(offsets.windows(2).all(|pair| pair[0] <= pair[1]));
        for bounds in offsets.windows(2) {
            let row = &columns[bounds[0] as usize..bounds[1] as usize];
            assert!(row.windows(2).all(|pair| pair[0] < pair[1]));
        }
        assert_eq!(
            expand_sparse(
                trees.len(),
                snapshots.n_distinct_splits(),
                &offsets,
                &columns
            ),
            dense
        );
    }
}

#[test]
fn sparse_presence_supports_variable_width_polytomies() {
    let trees = ["(A:1,B:1,C:1,D:1);", "((A:1,B:1):1,(C:1,D:1):1);"];
    let snapshots = Snapshots::from_newicks(&trees, true).unwrap();
    let (dense, _) = snapshots.build_presence_matrix();
    let (sparse, _) = snapshots.build_sparse_presence_matrix().unwrap();
    let offsets = decode_u64(&sparse.row_offsets);
    let columns = decode_u32(&sparse.column_indices);

    assert_eq!(offsets, [0, 4, 10]);
    assert_eq!(columns.len(), 10);
    assert_eq!(
        expand_sparse(
            trees.len(),
            snapshots.n_distinct_splits(),
            &offsets,
            &columns
        ),
        dense
    );
}

#[test]
fn empty_sparse_presence_has_one_zero_offset() {
    let snapshots = Snapshots::from_newicks(&[], false).unwrap();
    let (sparse, columns) = snapshots.build_sparse_presence_matrix().unwrap();

    assert!(columns.is_empty());
    assert_eq!(sparse.n_entries, 0);
    assert_eq!(decode_u64(&sparse.row_offsets), [0]);
    assert!(sparse.column_indices.is_empty());
}
