"""
Tests for pairwise_rf_from_newicks.

Uses the same 12 PHYLIP reference trees (10 taxa A-J) from
https://evolution.genetics.washington.edu/phylip/doc/treedist.html
that are also tested in the Rust unit tests (distances.rs).
"""

import re

import pytest

try:
    import rust_python_tree_distances as rtd
    RUST_MODULE_AVAILABLE = True
except ImportError:
    RUST_MODULE_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not RUST_MODULE_AVAILABLE,
    reason="rust_python_tree_distances module not installed. Run: maturin develop --release --features python",
)

# 12 PHYLIP reference trees with real taxon names (A-J)
REFERENCE_NEWICKS = [
    "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
]

# Expected RF distances (symmetric 12x12 matrix) from PHYLIP treedist
EXPECTED_RF = [
    [0, 4, 2, 10, 10, 10, 10, 10, 10, 10, 2, 10],
    [4, 0, 2, 10, 8, 10, 8, 10, 8, 10, 2, 10],
    [2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
    [10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
    [10, 8, 10, 2, 0, 4, 2, 4, 2, 2, 10, 4],
    [10, 10, 10, 2, 4, 0, 2, 2, 4, 2, 10, 2],
    [10, 8, 10, 4, 2, 2, 0, 4, 2, 4, 10, 4],
    [10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
    [10, 8, 10, 4, 2, 4, 2, 2, 0, 4, 10, 2],
    [10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
    [2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
    [10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
]


def natural_map():
    """Translate map: 1->A, 2->B, ..., 10->J (alphabetical order)."""
    taxa = list("ABCDEFGHIJ")
    return {str(i + 1): taxa[i] for i in range(10)}


def shuffled_map():
    """Translate map with a permuted assignment: 1->E, 2->H, 3->C, ..."""
    permutation = list("EHCJAGBFDI")
    return {str(i + 1): permutation[i] for i in range(10)}


def numberise_newick(newick, translate_map):
    """Replace taxon names with their numeric keys using a translate map.

    The translate_map maps numbers to names (e.g. {"1": "A", "2": "B"}).
    This function inverts the map and replaces names with numbers.
    """
    inv = {v: k for k, v in translate_map.items()}
    # Sort by name length descending to avoid partial replacements (e.g. "I" matching inside "BI")
    for name in sorted(inv.keys(), key=len, reverse=True):
        newick = re.sub(rf"(?<=[,(]){re.escape(name)}(?=[:)])", inv[name], newick)
    return newick


class TestSameTranslateMap:
    """All trees use the same translate map (natural ordering)."""

    def test_all_pairs_match_reference(self):
        tmap = natural_map()
        numbered = [numberise_newick(n, tmap) for n in REFERENCE_NEWICKS]
        names = [f"tree_{i}" for i in range(12)]

        tree_names, matrix = rtd.pairwise_rf_from_newicks(
            names, numbered, [tmap], [0] * 12
        )

        assert tree_names == names
        assert len(matrix) == 12
        for i in range(12):
            for j in range(12):
                assert matrix[i][j] == EXPECTED_RF[i][j], (
                    f"RF[{i}][{j}] = {matrix[i][j]}, expected {EXPECTED_RF[i][j]}"
                )


class TestCrossFileDifferentMaps:
    """Trees from two 'files' with different translate maps."""

    def test_representative_pairs(self):
        map_a = natural_map()
        map_b = shuffled_map()

        numbered_a = [numberise_newick(n, map_a) for n in REFERENCE_NEWICKS]
        numbered_b = [numberise_newick(n, map_b) for n in REFERENCE_NEWICKS]

        # Mix: first 6 use map_a, last 6 use map_b
        newicks = numbered_a[:6] + numbered_b[6:]
        names = [f"tree_{i}" for i in range(12)]
        map_indices = [0] * 6 + [1] * 6

        tree_names, matrix = rtd.pairwise_rf_from_newicks(
            names, newicks, [map_a, map_b], map_indices
        )

        assert tree_names == names
        # Check representative pairs across the two map groups
        pairs = [(0, 1), (0, 6), (3, 9), (5, 11), (2, 10)]
        for i, j in pairs:
            assert matrix[i][j] == EXPECTED_RF[i][j], (
                f"RF[{i}][{j}] = {matrix[i][j]}, expected {EXPECTED_RF[i][j]}"
            )


class TestIdenticalTopologyDifferentNumbering:
    """Same tree topology numbered with two different maps should yield RF=0."""

    def test_all_trees_rf_zero(self):
        map_a = natural_map()
        map_b = shuffled_map()

        for idx in range(12):
            newick_a = numberise_newick(REFERENCE_NEWICKS[idx], map_a)
            newick_b = numberise_newick(REFERENCE_NEWICKS[idx], map_b)

            names = ["a", "b"]
            tree_names, matrix = rtd.pairwise_rf_from_newicks(
                names, [newick_a, newick_b], [map_a, map_b], [0, 1]
            )
            assert matrix[0][1] == 0, (
                f"Tree {idx}: same topology with different numbering should have RF=0, got {matrix[0][1]}"
            )


class TestBeastAnnotations:
    """BEAST-style [&...] annotations are stripped before parsing."""

    def _inject_annotations(self, newick):
        """Insert fake BEAST annotations after every colon."""
        return newick.replace(":", ":[&rate=0.42]")

    def test_annotated_vs_clean_different_tree(self):
        tmap = natural_map()
        clean_0 = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        annotated_3 = self._inject_annotations(numberise_newick(REFERENCE_NEWICKS[3], tmap))

        names = ["clean_0", "annotated_3"]
        _, matrix = rtd.pairwise_rf_from_newicks(
            names, [clean_0, annotated_3], [tmap], [0, 0]
        )
        assert matrix[0][1] == EXPECTED_RF[0][3]

    def test_annotated_vs_clean_same_tree(self):
        tmap = natural_map()
        clean = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        annotated = self._inject_annotations(clean)

        names = ["clean", "annotated"]
        _, matrix = rtd.pairwise_rf_from_newicks(
            names, [clean, annotated], [tmap], [0, 0]
        )
        assert matrix[0][1] == 0


class TestValidation:
    """Input validation errors."""

    def test_names_newicks_length_mismatch(self):
        tmap = natural_map()
        newick = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        with pytest.raises(ValueError, match="names length"):
            rtd.pairwise_rf_from_newicks(
                ["a", "b"], [newick], [tmap], [0]
            )

    def test_map_index_out_of_bounds(self):
        tmap = natural_map()
        newicks = [
            numberise_newick(REFERENCE_NEWICKS[0], tmap),
            numberise_newick(REFERENCE_NEWICKS[1], tmap),
        ]
        with pytest.raises(ValueError, match="out of bounds"):
            rtd.pairwise_rf_from_newicks(
                ["a", "b"], newicks, [tmap], [0, 5]
            )

    def test_fewer_than_two_trees(self):
        tmap = natural_map()
        newick = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        with pytest.raises(ValueError, match="at least 2"):
            rtd.pairwise_rf_from_newicks(
                ["a"], [newick], [tmap], [0]
            )
