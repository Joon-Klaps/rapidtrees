"""Public contract tests for rooted RF snapshots with MrHIPSTR tree facts."""

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
FACT_KEYS = {
    "format_version",
    "root_column",
    "nodes_per_tree",
    "splits_per_tree",
    "node_columns",
    "node_heights",
    "root_heights",
    "split_columns",
}


def _call(
    trees=SOURCE_TREES,
    names=SOURCE_NAMES,
    translate_maps=None,
    map_indices=None,
):
    """Call the rooted-facts endpoint with a fresh lazy iterator."""
    trees = tuple(trees)
    names = list(names)
    translate_maps = [{}] if translate_maps is None else translate_maps
    map_indices = [0] * len(trees) if map_indices is None else map_indices
    return rtd.pairwise_rf_with_rooted_facts_from_newick_iter(
        names,
        iter(trees),
        translate_maps,
        map_indices,
    )


def _decode(result):
    """Decode the complete public wire format and assert its scalar contract."""
    assert isinstance(result, tuple)
    assert len(result) == 7
    tree_names, rf_bytes, leaf_names, n_clades, presence_bytes, clade_bytes, facts = result

    assert isinstance(tree_names, list)
    assert isinstance(rf_bytes, bytes)
    assert isinstance(leaf_names, list)
    assert isinstance(n_clades, int)
    assert isinstance(presence_bytes, bytes)
    assert isinstance(clade_bytes, bytes)
    assert isinstance(facts, dict)
    assert set(facts) == FACT_KEYS
    assert facts["format_version"] == 1
    assert facts["root_column"] == n_clades

    n_trees = len(tree_names)
    n_leaves = len(leaf_names)
    nodes_per_tree = 2 * n_leaves - 2
    splits_per_tree = n_leaves - 1
    assert facts["nodes_per_tree"] == nodes_per_tree
    assert facts["splits_per_tree"] == splits_per_tree

    for key in ("node_columns", "node_heights", "root_heights", "split_columns"):
        assert isinstance(facts[key], bytes)

    rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n_trees, n_trees)
    presence = np.frombuffer(presence_bytes, dtype=np.uint8).reshape(n_trees, n_clades)
    node_columns = np.frombuffer(facts["node_columns"], dtype=np.uint32).reshape(
        n_trees, nodes_per_tree
    )
    node_heights = np.frombuffer(facts["node_heights"], dtype=np.float64).reshape(
        n_trees, nodes_per_tree
    )
    root_heights = np.frombuffer(facts["root_heights"], dtype=np.float64).reshape(n_trees)
    split_columns = np.frombuffer(facts["split_columns"], dtype=np.uint32).reshape(
        n_trees, splits_per_tree, 3
    )

    bytes_per_clade = (n_leaves + 7) // 8
    packed_clades = np.frombuffer(clade_bytes, dtype=np.uint8).reshape(
        n_clades, bytes_per_clade
    )
    clades = np.unpackbits(packed_clades, axis=1, bitorder="little")[:, :n_leaves].astype(
        np.bool_, copy=False
    )

    assert rf.dtype == np.dtype(np.uint32)
    assert presence.dtype == np.dtype(np.uint8)
    assert node_columns.dtype == np.dtype(np.uint32)
    assert node_heights.dtype == np.dtype(np.float64)
    assert root_heights.dtype == np.dtype(np.float64)
    assert split_columns.dtype == np.dtype(np.uint32)
    assert clades.dtype == np.dtype(np.bool_)
    assert rf.shape == (n_trees, n_trees)
    assert presence.shape == (n_trees, n_clades)
    assert node_columns.shape == (n_trees, nodes_per_tree)
    assert node_heights.shape == (n_trees, nodes_per_tree)
    assert root_heights.shape == (n_trees,)
    assert split_columns.shape == (n_trees, splits_per_tree, 3)
    assert clades.shape == (n_clades, n_leaves)

    return {
        "tree_names": tree_names,
        "leaf_names": leaf_names,
        "n_clades": n_clades,
        "root_column": facts["root_column"],
        "nodes_per_tree": nodes_per_tree,
        "splits_per_tree": splits_per_tree,
        "rf": rf,
        "presence": presence,
        "clades": clades,
        "node_columns": node_columns,
        "node_heights": node_heights,
        "root_heights": root_heights,
        "split_columns": split_columns,
    }


def _height_map(decoded, tree_index):
    """Index one exported height row by public taxon-name clade."""
    mapping = {}
    for column, height in zip(
        decoded["node_columns"][tree_index],
        decoded["node_heights"][tree_index],
        strict=True,
    ):
        members = frozenset(
            name
            for name, present in zip(
                decoded["leaf_names"], decoded["clades"][column], strict=True
            )
            if present
        )
        mapping[members] = float(height)
    return mapping


def test_exact_types_dtypes_and_shapes():
    decoded = _decode(_call())

    assert decoded["tree_names"] == list(SOURCE_NAMES)
    assert decoded["leaf_names"] == ["A", "B", "C", "D"]
    assert decoded["n_clades"] == 10
    assert decoded["root_column"] == 10
    assert decoded["nodes_per_tree"] == 6
    assert decoded["splits_per_tree"] == 3


def test_first_six_outputs_match_existing_rooted_endpoint_byte_for_byte():
    rooted_facts_result = _call()
    established_result = rtd.pairwise_rf_with_snapshots_from_newick_iter(
        list(SOURCE_NAMES),
        iter(SOURCE_TREES),
        [{}],
        [0, 0, 0],
        rooted=True,
    )

    assert rooted_facts_result[:6] == established_result


def test_nodes_and_observed_splits_reference_public_clade_columns():
    decoded = _decode(_call())
    presence = decoded["presence"]
    clades = decoded["clades"]
    root_column = decoded["root_column"]
    all_taxa = np.ones(len(decoded["leaf_names"]), dtype=np.bool_)

    for tree_index in range(len(decoded["tree_names"])):
        node_columns = decoded["node_columns"][tree_index]
        present_columns = np.flatnonzero(presence[tree_index])
        assert set(node_columns.tolist()) == set(present_columns.tolist())
        assert len(set(node_columns.tolist())) == decoded["nodes_per_tree"]

        root_splits = 0
        for parent, left, right in decoded["split_columns"][tree_index]:
            parent = int(parent)
            left = int(left)
            right = int(right)
            assert left < right < root_column
            assert presence[tree_index, left] == 1
            assert presence[tree_index, right] == 1
            assert not np.any(clades[left] & clades[right])

            if parent == root_column:
                root_splits += 1
                parent_clade = all_taxa
            else:
                assert 0 <= parent < root_column
                assert presence[tree_index, parent] == 1
                parent_clade = clades[parent]

            np.testing.assert_array_equal(clades[left] | clades[right], parent_clade)

        assert root_splits == 1


def test_non_ultrametric_heights_follow_max_tip_distance_definition():
    decoded = _decode(_call())
    np.testing.assert_allclose(decoded["root_heights"], [11.0, 13.0, 8.0])

    assert _height_map(decoded, 0) == {
        frozenset({"A"}): 7.0,
        frozenset({"B"}): 6.0,
        frozenset({"A", "B"}): 8.0,
        frozenset({"C"}): 1.0,
        frozenset({"D"}): 0.0,
        frozenset({"C", "D"}): 5.0,
    }
    assert _height_map(decoded, 2) == {
        frozenset({"A"}): 5.0,
        frozenset({"D"}): 2.0,
        frozenset({"A", "D"}): 6.0,
        frozenset({"B"}): 0.0,
        frozenset({"C"}): 1.0,
        frozenset({"B", "C"}): 3.0,
    }


def test_translate_maps_and_beast_annotations_preserve_facts():
    translate_a = {"1": "A", "2": "B", "3": "C", "4": "D"}
    translate_b = {"1": "C", "2": "A", "3": "D", "4": "B"}
    trees = (
        "((1:[&rate=0.1]1,2:2):3,(3:4,4:[&rate=0.2]5):6);",
        "((2:1,4:[&rate=0.3]2):3,(1:4,3:5):6);",
    )
    decoded = _decode(
        _call(
            trees=trees,
            names=("translated-a", "translated-b"),
            translate_maps=[translate_a, translate_b],
            map_indices=[0, 1],
        )
    )

    assert decoded["leaf_names"] == ["A", "B", "C", "D"]
    assert decoded["rf"][0, 1] == 0
    np.testing.assert_array_equal(decoded["presence"][0], decoded["presence"][1])
    np.testing.assert_allclose(decoded["root_heights"], [11.0, 11.0])
    assert _height_map(decoded, 0) == _height_map(decoded, 1)


def test_repeated_calls_are_byte_for_byte_deterministic():
    first = _call()
    second = _call()

    assert first == second


@pytest.mark.parametrize(
    ("trees", "message"),
    [
        (("(A:1,B:1,C:1);",) * 2, "exactly two children"),
        (("(A:1,B);",) * 2, "missing an explicit branch length"),
        (
            ("((A:1e308,B:1e308):1e308,C:1);",) * 2,
            "non-finite cumulative root distance",
        ),
        (("(A:1,B:1);", "(A:1,C:1);"), "leaf set"),
    ],
)
def test_descriptive_value_errors(trees, message):
    with pytest.raises(ValueError, match=message):
        _call(trees=trees, names=("bad-0", "bad-1"))


def test_endpoint_has_no_ambiguous_rooted_mode_flag():
    with pytest.raises(TypeError, match="rooted"):
        rtd.pairwise_rf_with_rooted_facts_from_newick_iter(
            list(SOURCE_NAMES),
            iter(SOURCE_TREES),
            [{}],
            [0, 0, 0],
            rooted=True,
        )
