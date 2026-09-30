"""Public contract tests for compressed sparse-row RF snapshots."""

import numpy as np
import pytest

try:
    import rapidtrees as rtd

    RUST_MODULE_AVAILABLE = True
except ImportError:
    RUST_MODULE_AVAILABLE = False


pytestmark = pytest.mark.skipif(
    not RUST_MODULE_AVAILABLE,
    reason="rapidtrees module not installed. Run: maturin develop --release --features python",
)


SOURCE_TREES = (
    "((A:1,B:2):3,(C:4,D:5):6);",
    "((A:2,C:3):4,(B:5,D:6):7);",
    "((A:1,D:4):2,(B:3,C:2):5);",
)
SOURCE_NAMES = ("t0", "t1", "t2")
SPARSE_KEYS = {
    "format_version",
    "encoding",
    "n_entries",
    "row_offsets",
    "column_indices",
}


def _call(trees=SOURCE_TREES, names=SOURCE_NAMES, *, rooted=False):
    trees = tuple(trees)
    return rtd.pairwise_rf_with_sparse_snapshots_from_newick_iter(
        list(names),
        iter(trees),
        [{}],
        [0] * len(trees),
        rooted=rooted,
    )


def _decode(result):
    assert isinstance(result, tuple)
    assert len(result) == 6
    tree_names, rf_bytes, leaf_names, n_clades, clade_bytes, sparse = result

    assert isinstance(tree_names, list)
    assert isinstance(rf_bytes, bytes)
    assert isinstance(leaf_names, list)
    assert isinstance(n_clades, int)
    assert isinstance(clade_bytes, bytes)
    assert isinstance(sparse, dict)
    assert set(sparse) == SPARSE_KEYS
    assert sparse["format_version"] == 1
    assert sparse["encoding"] == "csr"
    assert isinstance(sparse["row_offsets"], bytes)
    assert isinstance(sparse["column_indices"], bytes)

    n_trees = len(tree_names)
    rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n_trees, n_trees)
    offsets = np.frombuffer(sparse["row_offsets"], dtype=np.uint64)
    columns = np.frombuffer(sparse["column_indices"], dtype=np.uint32)
    presence = np.zeros((n_trees, n_clades), dtype=np.uint8)

    assert offsets.shape == (n_trees + 1,)
    assert offsets[0] == 0
    assert np.all(offsets[:-1] <= offsets[1:])
    assert int(offsets[-1]) == sparse["n_entries"] == columns.size
    assert np.all(columns < n_clades)
    for tree_index, (start, stop) in enumerate(
        zip(offsets[:-1], offsets[1:], strict=True)
    ):
        row = columns[int(start) : int(stop)]
        assert np.all(row[:-1] < row[1:])
        presence[tree_index, row.astype(np.intp, copy=False)] = 1

    bytes_per_clade = (len(leaf_names) + 7) // 8
    assert len(clade_bytes) == n_clades * bytes_per_clade
    return {
        "tree_names": tree_names,
        "rf": rf,
        "leaf_names": leaf_names,
        "n_clades": n_clades,
        "clade_bytes": clade_bytes,
        "offsets": offsets,
        "columns": columns,
        "presence": presence,
    }


@pytest.mark.parametrize("rooted", [False, True])
def test_sparse_snapshots_reconstruct_established_dense_bytes(rooted):
    sparse_result = _call(rooted=rooted)
    dense_result = rtd.pairwise_rf_with_snapshots_from_newick_iter(
        list(SOURCE_NAMES),
        iter(SOURCE_TREES),
        [{}],
        [0] * len(SOURCE_TREES),
        rooted=rooted,
    )
    decoded = _decode(sparse_result)

    assert sparse_result[:4] == dense_result[:4]
    assert sparse_result[4] == dense_result[5]
    assert decoded["presence"].tobytes() == dense_result[4]


def test_sparse_snapshots_support_variable_width_polytomies():
    trees = (
        "(A:1,B:1,C:1,D:1);",
        "((A:1,B:1):1,(C:1,D:1):1);",
    )
    decoded = _decode(_call(trees, ("star", "binary"), rooted=True))

    np.testing.assert_array_equal(decoded["offsets"], [0, 4, 10])
    np.testing.assert_array_equal(decoded["presence"].sum(axis=1), [4, 6])


def test_rooted_binary_columns_match_rooted_facts_exactly():
    sparse_result = _call(rooted=True)
    rooted_facts_result = rtd.pairwise_rf_with_rooted_facts_from_newick_iter(
        list(SOURCE_NAMES),
        iter(SOURCE_TREES),
        [{}],
        [0] * len(SOURCE_TREES),
    )
    sparse = sparse_result[5]
    facts = rooted_facts_result[5]

    assert sparse_result[:5] == rooted_facts_result[:5]
    assert sparse["column_indices"] == facts["clade_columns"]
    np.testing.assert_array_equal(
        np.frombuffer(sparse["row_offsets"], dtype=np.uint64),
        np.arange(len(SOURCE_TREES) + 1, dtype=np.uint64)
        * facts["nodes_per_tree"],
    )


def test_sparse_snapshot_output_is_byte_deterministic():
    assert _call(rooted=False) == _call(rooted=False)
    assert _call(rooted=True) == _call(rooted=True)
