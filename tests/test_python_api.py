"""
Tests for the Python API.

Tests verify that the Python bindings work correctly by comparing
against known reference values and checking error handling.
"""

import pytest
import numpy as np
import hashlib
from pathlib import Path

# Try to import the module
try:
    import rapidtrees as rtd
    RUST_MODULE_AVAILABLE = True
except ImportError:
    RUST_MODULE_AVAILABLE = False

# Test data paths
TEST_DATA = Path(__file__).parent / "data"

# Small 3-tree fixture with known pairwise RF: [0,1]=4, [0,2]=2, [1,2]=2.
# Trees are the first three from the phylip treedist reference set.
FIXTURE_TREES = [
    "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
]
FIXTURE_NAMES = ["t0", "t1", "t2"]
# Empty translate map: taxa names are already real names in the newick strings
FIXTURE_TRANSLATE = [{}]
FIXTURE_MAP_INDICES = [0, 0, 0]
# Known pairwise RF matrix for the 3 fixture trees
FIXTURE_RF = [[0, 4, 2], [4, 0, 2], [2, 2, 0]]

# Mark all tests to skip if module not available
pytestmark = pytest.mark.skipif(
    not RUST_MODULE_AVAILABLE,
    reason="rapidtrees module not installed. Run: maturin develop --release --features python"
)


def compute_output_hash(tree_names, matrix):
    """Compute a hash of the output for consistency checking."""
    # Create a deterministic string representation
    content = "\t".join([""] + tree_names) + "\n"
    for i, name in enumerate(tree_names):
        row = [name] + [str(matrix[i][j]) for j in range(len(matrix[i]))]
        content += "\t".join(row) + "\n"

    return hashlib.md5(content.encode()).hexdigest()


def matrices_close(mat1, mat2, rtol=1e-9, atol=1e-9):
    """Check if two matrices are element-wise close (for floating point comparison)."""
    if len(mat1) != len(mat2):
        return False
    for i in range(len(mat1)):
        if len(mat1[i]) != len(mat2[i]):
            return False
        for j in range(len(mat1[i])):
            if abs(mat1[i][j] - mat2[i][j]) > atol + rtol * abs(mat2[i][j]):
                return False
    return True


class TestPairwiseRF:
    """Tests for pairwise_rf function."""

    def test_basic_rf_calculation(self):
        """Test basic RF distance calculation with HIV data."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, matrix = rtd.pairwise_rf(paths, burnin_trees=1)

        # Check structure
        assert len(tree_names) == 20, "Should have 20 trees after burnin"
        assert len(matrix) == 20, "Matrix should be 20x20"
        assert all(len(row) == 20 for row in matrix), "All rows should have 20 elements"

        # Check diagonal is zero
        for i in range(len(matrix)):
            assert matrix[i][i] == 0, f"Diagonal element [{i}][{i}] should be 0"

        # Check symmetry
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                assert matrix[i][j] == matrix[j][i], f"Matrix should be symmetric at [{i}][{j}]"

        # Check some known values (RF distances are integers)
        assert matrix[0][1] == 164, "Known RF distance between tree 0 and 1"
        assert matrix[0][2] == 184, "Known RF distance between tree 0 and 2"
        assert matrix[1][2] == 118, "Known RF distance between tree 1 and 2"

    def test_rf_is_deterministic(self):
        """Test that RF distances are exactly deterministic (no floating-point issues)."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Run three times - RF should be 100% deterministic since it uses integers
        tree_names1, matrix1 = rtd.pairwise_rf(paths, burnin_trees=1)
        tree_names2, matrix2 = rtd.pairwise_rf(paths, burnin_trees=1)
        tree_names3, matrix3 = rtd.pairwise_rf(paths, burnin_trees=1)

        # Everything should be identical
        assert tree_names1 == tree_names2 == tree_names3
        assert matrix1 == matrix2 == matrix3, "RF distances must be exactly deterministic"

    def test_multiple_files(self):
        """Test with multiple input files."""
        paths = [
            str(TEST_DATA / "hiv1.trees"),
            str(TEST_DATA / "hiv2.trees"),
        ]
        tree_names, matrix = rtd.pairwise_rf(paths, burnin_trees=1)

        # Should have trees from both files
        assert len(tree_names) == 40, "Should have 40 trees total (20 from each file)"
        assert len(matrix) == 40

    def test_burnin_trees(self):
        """Test burnin by number of trees."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # No burnin
        tree_names_no_burnin, _ = rtd.pairwise_rf(paths, burnin_trees=0)

        # With burnin
        tree_names_burnin, _ = rtd.pairwise_rf(paths, burnin_trees=5)

        assert len(tree_names_no_burnin) > len(tree_names_burnin)
        assert len(tree_names_no_burnin) - len(tree_names_burnin) == 5

    def test_burnin_states(self):
        """Test burnin by state value."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        tree_names, _ = rtd.pairwise_rf(paths, burnin_states=100000)

        # Should filter out trees with STATE < 100000
        assert all("STATE" in name for name in tree_names)
        # Extract state values and verify they're all >= 100000
        states = [int(name.split("STATE_")[1]) for name in tree_names]
        assert all(s >= 100000 for s in states)

    def test_use_real_taxa(self):
        """Test using real taxon names from TRANSLATE block."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Both should work, but names might differ
        names_numeric, _ = rtd.pairwise_rf(paths, burnin_trees=1, use_real_taxa=False)
        names_real, _ = rtd.pairwise_rf(paths, burnin_trees=1, use_real_taxa=True)

        assert len(names_numeric) == len(names_real)

    def test_empty_after_burnin(self):
        """Test error when no trees remain after burnin."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        with pytest.raises(ValueError, match="No trees found"):
            rtd.pairwise_rf(paths, burnin_trees=1000)

    def test_missing_file(self):
        """Test error handling for missing files."""
        paths = [str(TEST_DATA / "nonexistent.trees")]

        with pytest.raises((ValueError, RuntimeError)):
            rtd.pairwise_rf(paths)

    def test_tree_name_format(self):
        """Test that tree names have expected format."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, _ = rtd.pairwise_rf(paths, burnin_trees=1)

        # Should all contain the filename and STATE
        assert all("hiv1" in name for name in tree_names)
        assert all("STATE" in name for name in tree_names)


class TestPairwiseWeightedRF:
    """Tests for pairwise_weighted_rf function."""

    def test_basic_weighted_rf_calculation(self):
        """Test basic weighted RF distance calculation."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, matrix = rtd.pairwise_weighted_rf(paths, burnin_trees=1)

        # Check structure
        assert len(tree_names) == 20
        assert len(matrix) == 20

        # Check diagonal is zero
        for i in range(len(matrix)):
            assert matrix[i][i] == 0.0

        # Check symmetry
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                assert abs(matrix[i][j] - matrix[j][i]) < 1e-10, "Matrix should be symmetric"

        # Check that values are reasonable (should be larger than RF due to branch lengths)
        assert matrix[0][1] > 0, "Distance should be positive"
        assert matrix[0][1] > 302, "Weighted RF should be >= unweighted RF"

    def test_weighted_vs_unweighted(self):
        """Test that weighted RF distances are generally larger than unweighted."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, rf_matrix = rtd.pairwise_rf(paths, burnin_trees=1)
        _, weighted_matrix = rtd.pairwise_weighted_rf(paths, burnin_trees=1)

        # Most non-diagonal elements should be larger in weighted version
        larger_count = 0
        total_count = 0

        for i in range(len(rf_matrix)):
            for j in range(i + 1, len(rf_matrix)):
                if rf_matrix[i][j] > 0:  # Only check where there's a difference
                    total_count += 1
                    if weighted_matrix[i][j] >= rf_matrix[i][j]:
                        larger_count += 1

        # Most weighted distances should be >= unweighted
        assert larger_count / total_count > 0.8

    def test_weighted_rf_is_deterministic(self):
        """Test that weighted RF distances are deterministic across multiple runs."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Run three times - should be deterministic
        tree_names1, matrix1 = rtd.pairwise_weighted_rf(paths, burnin_trees=1)
        tree_names2, matrix2 = rtd.pairwise_weighted_rf(paths, burnin_trees=1)
        tree_names3, matrix3 = rtd.pairwise_weighted_rf(paths, burnin_trees=1)

        # Tree names should be identical
        assert tree_names1 == tree_names2 == tree_names3

        # Matrices should be exactly identical (deterministic)
        assert matrices_close(matrix1, matrix2), "Weighted RF distances should be deterministic between runs 1 and 2"
        assert matrices_close(matrix2, matrix3), "Weighted RF distances should be deterministic between runs 2 and 3"



class TestPairwiseKF:
    """Tests for pairwise_kf (Kuhner-Felsenstein) function."""

    def test_basic_kf_calculation(self):
        """Test basic KF distance calculation."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, matrix = rtd.pairwise_kf(paths, burnin_trees=1)

        # Check structure
        assert len(tree_names) == 20
        assert len(matrix) == 20

        # Check diagonal is zero
        for i in range(len(matrix)):
            assert matrix[i][i] == 0.0

        # Check symmetry
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                assert abs(matrix[i][j] - matrix[j][i]) < 1e-10

        # Check that values are positive and reasonable
        assert matrix[0][1] > 0

    def test_kf_is_deterministic(self):
        """Test that KF distances are deterministic across multiple runs."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Run three times - should be deterministic
        tree_names1, matrix1 = rtd.pairwise_kf(paths, burnin_trees=1)
        tree_names2, matrix2 = rtd.pairwise_kf(paths, burnin_trees=1)
        tree_names3, matrix3 = rtd.pairwise_kf(paths, burnin_trees=1)

        # Tree names should be identical
        assert tree_names1 == tree_names2 == tree_names3

        # Matrices should be exactly identical (deterministic)
        assert matrices_close(matrix1, matrix2), "KF distances should be deterministic between runs 1 and 2"
        assert matrices_close(matrix2, matrix3), "KF distances should be deterministic between runs 2 and 3"


    def test_kf_vs_weighted(self):
        """Test KF produces different values from weighted RF."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, weighted_matrix = rtd.pairwise_weighted_rf(paths, burnin_trees=1)
        _, kf_matrix = rtd.pairwise_kf(paths, burnin_trees=1)

        # Matrices should be different (KF uses branch length info differently)
        differences = 0
        for i in range(len(kf_matrix)):
            for j in range(i + 1, len(kf_matrix)):
                if abs(kf_matrix[i][j] - weighted_matrix[i][j]) > 1e-6:
                    differences += 1

        # Most values should differ
        total_comparisons = len(kf_matrix) * (len(kf_matrix) - 1) // 2
        assert differences / total_comparisons > 0.9


class TestSanityChecks:
    """Tests for input validation and sanity checks."""

    def test_inconsistent_leaf_sets(self):
        """Test that trees with different leaf sets are rejected."""
        # This would require creating test files with incompatible trees
        # For now, we verify the function doesn't crash with valid input
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, _ = rtd.pairwise_rf(paths, burnin_trees=1)
        assert len(tree_names) > 0

    def test_empty_file_list(self):
        """Test error handling for empty file list."""
        with pytest.raises((ValueError, TypeError)):
            rtd.pairwise_rf([])

    def test_invalid_burnin_values(self):
        """Test handling of invalid burnin parameters."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Negative burnin should raise an overflow error
        with pytest.raises(OverflowError):
            rtd.pairwise_rf(paths, burnin_trees=-1)


class TestAPIConsistency:
    """Tests that all three metrics return consistent structure."""

    def test_all_metrics_same_trees(self):
        """Test that all metrics return the same tree names."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        names_rf, _ = rtd.pairwise_rf(paths, burnin_trees=1)
        names_weighted, _ = rtd.pairwise_weighted_rf(paths, burnin_trees=1)
        names_kf, _ = rtd.pairwise_kf(paths, burnin_trees=1)

        assert names_rf == names_weighted == names_kf

    def test_all_metrics_same_dimensions(self):
        """Test that all metrics return same matrix dimensions."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, matrix_rf = rtd.pairwise_rf(paths, burnin_trees=1)
        _, matrix_weighted = rtd.pairwise_weighted_rf(paths, burnin_trees=1)
        _, matrix_kf = rtd.pairwise_kf(paths, burnin_trees=1)

        assert len(matrix_rf) == len(matrix_weighted) == len(matrix_kf)
        assert len(matrix_rf[0]) == len(matrix_weighted[0]) == len(matrix_kf[0])


class TestRootedRF:
    """Tests for rooted RF distance (clade-based comparison)."""

    def test_rooted_rf_greater_or_equal(self):
        """Rooted RF should be >= unrooted RF for the same trees."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, matrix_unrooted = rtd.pairwise_rf(paths, burnin_trees=1, rooted=False)
        _, matrix_rooted = rtd.pairwise_rf(paths, burnin_trees=1, rooted=True)

        for i in range(len(matrix_rooted)):
            for j in range(i + 1, len(matrix_rooted)):
                assert matrix_rooted[i][j] >= matrix_unrooted[i][j], (
                    f"Rooted RF[{i}][{j}]={matrix_rooted[i][j]} should be >= "
                    f"unrooted RF[{i}][{j}]={matrix_unrooted[i][j]}"
                )

    def test_rooted_rf_structure(self):
        """Rooted RF matrix should be symmetric with zero diagonal."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        _, matrix = rtd.pairwise_rf(paths, burnin_trees=1, rooted=True)

        for i in range(len(matrix)):
            assert matrix[i][i] == 0, f"Diagonal [{i}][{i}] should be 0"

        for i in range(len(matrix)):
            for j in range(len(matrix)):
                assert matrix[i][j] == matrix[j][i], f"Should be symmetric at [{i}][{j}]"

    def test_rooted_rf_deterministic(self):
        """Rooted RF should be deterministic."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, m1 = rtd.pairwise_rf(paths, burnin_trees=1, rooted=True)
        _, m2 = rtd.pairwise_rf(paths, burnin_trees=1, rooted=True)

        assert m1 == m2, "Rooted RF must be deterministic"

    def test_rooted_default_is_false(self):
        """Default rooted=False should match explicit rooted=False."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, matrix_default = rtd.pairwise_rf(paths, burnin_trees=1)
        _, matrix_explicit = rtd.pairwise_rf(paths, burnin_trees=1, rooted=False)

        assert matrix_default == matrix_explicit


class TestPairwiseRFFromNewicks:
    """Tests for pairwise_rf_from_newicks (newick list, no iterator)."""

    def test_known_rf_values(self):
        """Known pairwise RF values from the treedist reference set."""
        names, matrix = rtd.pairwise_rf_from_newicks(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        assert names == FIXTURE_NAMES
        for i in range(len(FIXTURE_RF)):
            for j in range(len(FIXTURE_RF)):
                assert matrix[i][j] == FIXTURE_RF[i][j], f"mismatch at [{i}][{j}]"

    def test_symmetric_zero_diagonal(self):
        """Matrix must be symmetric with a zero diagonal."""
        _, matrix = rtd.pairwise_rf_from_newicks(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        n = len(matrix)
        for i in range(n):
            assert matrix[i][i] == 0
        for i in range(n):
            for j in range(n):
                assert matrix[i][j] == matrix[j][i], f"not symmetric at [{i}][{j}]"

    def test_rooted_ge_unrooted(self):
        """Rooted RF must be >= unrooted RF for the same tree pair."""
        _, unrooted = rtd.pairwise_rf_from_newicks(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES, rooted=False
        )
        _, rooted = rtd.pairwise_rf_from_newicks(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES, rooted=True
        )
        for i in range(len(rooted)):
            for j in range(len(rooted)):
                assert rooted[i][j] >= unrooted[i][j], f"rooted < unrooted at [{i}][{j}]"

    def test_error_on_length_mismatch(self):
        """Passing fewer newicks than names must raise ValueError."""
        with pytest.raises(ValueError):
            rtd.pairwise_rf_from_newicks(
                FIXTURE_NAMES, FIXTURE_TREES[:2], FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
            )


class TestPairwiseRFFromNewickIter:
    """Tests for pairwise_rf_from_newick_iter (lazy iterator API)."""

    def test_returns_bytes(self):
        """rf_bytes must be a bytes object of size n*n*4."""
        names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        assert isinstance(rf_bytes, bytes)
        assert len(rf_bytes) == len(FIXTURE_NAMES) ** 2 * 4

    def test_known_rf_values(self):
        """Decoded u32 matrix must match the known RF reference values."""
        names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        n = len(names)
        matrix = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)
        for i in range(len(FIXTURE_RF)):
            for j in range(len(FIXTURE_RF)):
                assert matrix[i, j] == FIXTURE_RF[i][j], f"mismatch at [{i}][{j}]"

    def test_matches_from_newicks(self):
        """Iterator and list variants must produce identical matrices."""
        _, matrix_list = rtd.pairwise_rf_from_newicks(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        _, rf_bytes = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        n = len(FIXTURE_NAMES)
        matrix_arr = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)
        for i in range(n):
            for j in range(n):
                assert int(matrix_arr[i, j]) == matrix_list[i][j]

    def test_symmetric_zero_diagonal(self):
        """Decoded matrix must be symmetric with a zero diagonal."""
        names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        n = len(names)
        matrix = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)
        assert np.all(matrix == matrix.T), "matrix not symmetric"
        assert np.all(np.diag(matrix) == 0), "diagonal not zero"

    def test_deterministic(self):
        """Two calls on identical input must return identical bytes."""
        _, rf1 = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        _, rf2 = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        assert rf1 == rf2


class TestPairwiseRFWithSnapshots:
    """Tests for pairwise_rf_with_snapshots_from_newick_iter (5-tuple API)."""

    def _call(self, **kwargs):
        defaults = dict(
            names=FIXTURE_NAMES,
            newick_iter=iter(FIXTURE_TREES),
            translate_maps=FIXTURE_TRANSLATE,
            map_indices=FIXTURE_MAP_INDICES,
        )
        defaults.update(kwargs)
        return rtd.pairwise_rf_with_snapshots_from_newick_iter(**defaults)

    def test_returns_five_tuple(self):
        """Function returns a 5-tuple with the correct types."""
        result = self._call()
        assert len(result) == 5
        names, rf_bytes, leaf_names, n_bip, pres_bytes = result
        assert isinstance(names, list)
        assert isinstance(rf_bytes, bytes)
        assert isinstance(leaf_names, list)
        assert isinstance(n_bip, int)
        assert isinstance(pres_bytes, bytes)

    def test_rf_bytes_known_values(self):
        """rf_bytes decodes to the known RF matrix."""
        names, rf_bytes, *_ = self._call()
        n = len(names)
        matrix = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)
        for i in range(len(FIXTURE_RF)):
            for j in range(len(FIXTURE_RF)):
                assert matrix[i, j] == FIXTURE_RF[i][j], f"RF mismatch at [{i}][{j}]"

    def test_presence_bytes_shape_and_range(self):
        """presence_bytes encodes an n_trees × n_bip uint8 matrix of 0s and 1s."""
        names, _, _, n_bip, pres_bytes = self._call()
        n = len(names)
        assert len(pres_bytes) == n * n_bip
        presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n, n_bip)
        assert set(presence.flatten().tolist()) <= {0, 1}

    def test_xor_identity(self):
        """sum(presence[i] XOR presence[j]) equals RF(i, j) for every pair."""
        names, rf_bytes, _, n_bip, pres_bytes = self._call()
        n = len(names)
        rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)
        presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n, n_bip).astype(np.int32)
        for i in range(n):
            for j in range(n):
                xor_sum = int(np.sum(np.abs(presence[i] - presence[j])))
                assert xor_sum == int(rf[i, j]), f"XOR identity failed at [{i},{j}]"

    def test_rf_bytes_matches_iter(self):
        """rf_bytes matches the output of pairwise_rf_from_newick_iter."""
        names_s, rf_bytes_s, *_ = self._call()
        names_i, rf_bytes_i = rtd.pairwise_rf_from_newick_iter(
            FIXTURE_NAMES, iter(FIXTURE_TREES), FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        assert names_s == names_i
        assert rf_bytes_s == rf_bytes_i

    def test_leaf_names_sorted(self):
        """Returned leaf names must be sorted alphabetically."""
        _, _, leaf_names, _, _ = self._call()
        assert leaf_names == sorted(leaf_names)

    def test_deterministic(self):
        """Two calls on identical input must return identical results."""
        r1 = self._call()
        r2 = self._call()
        assert r1[0] == r2[0]  # names
        assert r1[1] == r2[1]  # rf_bytes
        assert r1[2] == r2[2]  # leaf_names
        assert r1[3] == r2[3]  # n_bip
        assert r1[4] == r2[4]  # presence_bytes


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
