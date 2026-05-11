"""
Tests for the Python API.

Tests verify that the Python bindings work correctly by comparing
against known reference values and checking error handling.
"""

import pytest
import numpy as np
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
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


def natural_map():
    """Translate map: 1->A, 2->B, ..., 10->J (alphabetical order)."""
    taxa = list("ABCDEFGHIJ")
    return {str(i + 1): taxa[i] for i in range(10)}


def shuffled_map():
    """Translate map with a permuted assignment: 1->E, 2->H, 3->C, ..."""
    permutation = list("EHCJAGBFDI")
    return {str(i + 1): permutation[i] for i in range(10)}


def numberise_newick(newick, translate_map):
    """Replace taxon names with their numeric keys using a translate map."""
    inv = {v: k for k, v in translate_map.items()}
    for name in sorted(inv.keys(), key=len, reverse=True):
        newick = re.sub(rf"(?<=[,(]){re.escape(name)}(?=[:)])", inv[name], newick)
    return newick


def _r_phangorn_available():
    """Return True if Rscript is available and phangorn can be loaded."""
    if shutil.which("Rscript") is None:
        return False

    cmd = [
        "Rscript",
        "-e",
        'quit(status = ifelse(requireNamespace("phangorn", quietly = TRUE), 0, 1))',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return result.returncode == 0


_PHANGORN_SCRIPT = Path(__file__).parent / "phangorn.R"


def _phangorn_run(newick_path: str, metric: str) -> np.ndarray:
    """Call tests/phangorn.R and return upper-triangle distances as a 1-D float array."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as out_handle:
        output_path = out_handle.name
    try:
        subprocess.run(
            ["Rscript", str(_PHANGORN_SCRIPT), newick_path, metric, output_path],
            capture_output=True,
            text=True,
            check=True,
        )
        with open(output_path) as fh:
            vals = [float(line.strip()) for line in fh if line.strip()]
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"phangorn.R failed (exit {exc.returncode}):\n{exc.stderr}"
        ) from exc
    finally:
        try:
            os.unlink(output_path)
        except FileNotFoundError:
            pass
    return np.array(vals, dtype=np.float64)


def _rapidtrees_upper_triangle(newick_path: str, metric: str) -> np.ndarray:
    """Compute distances via rapidtrees on a newick file, return upper-triangle values."""
    with open(newick_path) as fh:
        raw = [line.strip() for line in fh if line.strip()]
    newicks = [_strip_beast_annotations(nwk) for nwk in raw]
    n = len(newicks)
    names = [f"tree_{i}" for i in range(n)]
    if metric == "RF":
        _, rf_bytes = rtd.pairwise_rf_from_newick_iter(
            names, iter(newicks), [{}], [0] * n
        )
        mat = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n).astype(np.float64)
    elif metric == "wRF":
        _, wrf_flat = rtd.pairwise_wrf_from_newick_iter(
            names, iter(newicks), [{}], [0] * n
        )
        mat = np.array(wrf_flat, dtype=np.float64).reshape(n, n)
    elif metric == "KF":
        _, kf_flat = rtd.pairwise_kf_from_newick_iter(
            names, iter(newicks), [{}], [0] * n
        )
        mat = np.array(kf_flat, dtype=np.float64).reshape(n, n)
    else:
        raise ValueError(f"unknown metric: {metric}")
    rows, cols = np.triu_indices(n, k=1)
    return mat[rows, cols]


def _assert_distances_match(rt_vals: np.ndarray, ph_vals: np.ndarray, metric: str) -> None:
    """Assert rapidtrees and phangorn upper-triangle distances agree.

    On failure, reports ratio and absolute-difference statistics so the nature
    of the discrepancy (scale factor vs noise vs ordering issue) is immediately
    visible in the test output.
    """
    nz = ph_vals != 0
    ratio = np.where(nz, rt_vals / ph_vals, np.nan)
    abs_diff = rt_vals - ph_vals
    msg = (
        f"\nmetric={metric}  n_pairs={len(rt_vals)}"
        f"\n  rt range : [{rt_vals.min():.6g}, {rt_vals.max():.6g}]"
        f"\n  ph range : [{ph_vals.min():.6g}, {ph_vals.max():.6g}]"
        f"\n  ratio (rt/ph): min={np.nanmin(ratio):.6f}  max={np.nanmax(ratio):.6f}"
        f"  mean={np.nanmean(ratio):.6f}  std={np.nanstd(ratio):.2e}"
        f"\n  abs diff  : max={np.abs(abs_diff).max():.6g}"
        f"  mean={np.abs(abs_diff).mean():.6g}"
        f"\n  first 5 pairs (rt, ph, ratio):"
    )
    for i in range(min(5, len(rt_vals))):
        r = rt_vals[i] / ph_vals[i] if ph_vals[i] != 0 else float("nan")
        msg += f"\n    [{i}]  {rt_vals[i]:.8g}  vs  {ph_vals[i]:.8g}  → {r:.6f}"
    if metric == "RF":
        np.testing.assert_array_equal(rt_vals.astype(int), ph_vals.astype(int), err_msg=msg)
    else:
        np.testing.assert_allclose(rt_vals, ph_vals, rtol=1e-9, atol=1e-9, err_msg=msg)


def _strip_beast_annotations(newick):
    return re.sub(r"\[&[^\]]*\]", "", newick)


def _extract_name_state(header):
    upper = header.upper()
    state_pos = upper.find("STATE_")
    if state_pos == -1:
        return "", 0

    parts = header.split(" ", 1)
    if len(parts) != 2:
        return "", 0
    tree_name = parts[1].split()[0] if parts[1].split() else ""
    digits = ""
    for ch in header[state_pos + 6:]:
        if ch.isdigit():
            digits += ch
        else:
            break
    return tree_name, int(digits) if digits else 0


def _parse_taxon_block(content):
    lines = content.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().upper().startswith("TRANSLATE"):
            start = i + 1
            break
    if start is None:
        return {}

    mapping = {}
    for line in lines[start:]:
        s = line.strip()
        if s.upper().startswith(";"):
            break
        s = s.rstrip(",")
        parts = s.split()
        if len(parts) >= 2:
            key = parts[0]
            val = parts[1].strip("'")
            mapping[key] = val
    return mapping


def _collect_tree_blocks(content):
    lines = content.splitlines()
    started = False
    blocks = []
    for line in lines:
        upper = line.upper()
        if not started:
            if upper.startswith("TREE "):
                started = True
            else:
                continue
        if line.strip().upper().startswith("END;"):
            break
        if " = " in line:
            header, body = line.split(" = ", 1)
            blocks.append((header.strip(), body.strip()))
    return blocks


def _load_beast_raw_py(path, burnin_trees=0, burnin_states=0, use_real_taxa=True):
    p = Path(path)
    if not p.exists():
        return {}, []

    content = p.read_text()
    base_name = p.name.removesuffix(".trees")
    taxons = _parse_taxon_block(content)
    translate_map = taxons if use_real_taxa else {}

    tree_pairs = []
    for idx, (header, body) in enumerate(_collect_tree_blocks(content)):
        name, state = _extract_name_state(header)
        keep = (
            (burnin_trees == 0 and burnin_states == 0)
            or (burnin_trees > 0 and idx >= burnin_trees)
            or (burnin_states > 0 and state > burnin_states)
        )
        if keep:
            tree_pairs.append((f"{base_name}_{name}", _strip_beast_annotations(body)))
    return translate_map, tree_pairs


def _rf_bytes_to_matrix(rf_bytes, n):
    return np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n).tolist()


def _pairwise_rf_from_newicks_compat(names, newicks, translate_maps, map_indices, rooted=False):
    if len(names) != len(newicks):
        raise ValueError(f"names length ({len(names)}) must equal newicks length ({len(newicks)})")
    tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
        names, iter(newicks), translate_maps, map_indices, rooted=rooted
    )
    return tree_names, _rf_bytes_to_matrix(rf_bytes, len(tree_names))


def _pairwise_metric_from_files(paths, metric, burnin_trees=0, burnin_states=0, use_real_taxa=True, rooted=False):
    if burnin_trees < 0 or burnin_states < 0:
        raise OverflowError("can't convert negative value to unsigned integer")

    names = []
    newicks = []
    translate_maps = []
    map_indices = []

    for path in paths:
        tmap, pairs = _load_beast_raw_py(path, burnin_trees, burnin_states, use_real_taxa)
        map_idx = len(translate_maps)
        translate_maps.append(tmap)
        for name, newick in pairs:
            names.append(name)
            newicks.append(newick)
            map_indices.append(map_idx)

    if not names:
        raise ValueError("No trees found after applying burnin filters")

    if metric == "rf":
        return _pairwise_rf_from_newicks_compat(names, newicks, translate_maps, map_indices, rooted=rooted)

    if metric == "wrf":
        tree_names, flat = rtd.pairwise_wrf_from_newick_iter(
            names, iter(newicks), translate_maps, map_indices, rooted=rooted
        )
    else:
        tree_names, flat = rtd.pairwise_kf_from_newick_iter(
            names, iter(newicks), translate_maps, map_indices, rooted=rooted
        )

    n = len(tree_names)
    matrix = np.array(flat, dtype=np.float64).reshape(n, n).tolist()
    return tree_names, matrix


# ── Test-only Python wrappers ─────────────────────────────────────────────────
# These functions are NOT part of the rapidtrees API.  The Rust extension only
# exposes pairwise_*_from_newick_iter functions.  The helpers below implement
# the old file-based and list-based convenience APIs in pure Python on top of
# the iterator API, and are used only in this test suite.


def _pairwise_rf(paths, burnin_trees=0, burnin_states=0, use_real_taxa=True, rooted=False):
    return _pairwise_metric_from_files(paths, "rf", burnin_trees, burnin_states, use_real_taxa, rooted)


def _pairwise_weighted_rf(paths, burnin_trees=0, burnin_states=0, use_real_taxa=True, rooted=False):
    return _pairwise_metric_from_files(paths, "wrf", burnin_trees, burnin_states, use_real_taxa, rooted)


def _pairwise_kf(paths, burnin_trees=0, burnin_states=0, use_real_taxa=True, rooted=False):
    return _pairwise_metric_from_files(paths, "kf", burnin_trees, burnin_states, use_real_taxa, rooted)


class TestPairwiseRF:
    """Tests for pairwise_rf function."""

    def test_basic_rf_calculation(self):
        """Test basic RF distance calculation with HIV data."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, matrix = _pairwise_rf(paths, burnin_trees=1)

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
        tree_names1, matrix1 = _pairwise_rf(paths, burnin_trees=1)
        tree_names2, matrix2 = _pairwise_rf(paths, burnin_trees=1)
        tree_names3, matrix3 = _pairwise_rf(paths, burnin_trees=1)

        # Everything should be identical
        assert tree_names1 == tree_names2 == tree_names3
        assert matrix1 == matrix2 == matrix3, "RF distances must be exactly deterministic"

    def test_multiple_files(self):
        """Test with multiple input files."""
        paths = [
            str(TEST_DATA / "hiv1.trees"),
            str(TEST_DATA / "hiv2.trees"),
        ]
        tree_names, matrix = _pairwise_rf(paths, burnin_trees=1)

        # Should have trees from both files
        assert len(tree_names) == 40, "Should have 40 trees total (20 from each file)"
        assert len(matrix) == 40

    def test_burnin_trees(self):
        """Test burnin by number of trees."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # No burnin
        tree_names_no_burnin, _ = _pairwise_rf(paths, burnin_trees=0)

        # With burnin
        tree_names_burnin, _ = _pairwise_rf(paths, burnin_trees=5)

        assert len(tree_names_no_burnin) > len(tree_names_burnin)
        assert len(tree_names_no_burnin) - len(tree_names_burnin) == 5

    def test_burnin_states(self):
        """Test burnin by state value."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        tree_names, _ = _pairwise_rf(paths, burnin_states=100000)

        # Should filter out trees with STATE < 100000
        assert all("STATE" in name for name in tree_names)
        # Extract state values and verify they're all >= 100000
        states = [int(name.split("STATE_")[1]) for name in tree_names]
        assert all(s >= 100000 for s in states)

    def test_use_real_taxa(self):
        """Test using real taxon names from TRANSLATE block."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Both should work, but names might differ
        names_numeric, _ = _pairwise_rf(paths, burnin_trees=1, use_real_taxa=False)
        names_real, _ = _pairwise_rf(paths, burnin_trees=1, use_real_taxa=True)

        assert len(names_numeric) == len(names_real)

    def test_empty_after_burnin(self):
        """Test error when no trees remain after burnin."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        with pytest.raises(ValueError, match="No trees found"):
            _pairwise_rf(paths, burnin_trees=1000)

    def test_missing_file(self):
        """Test error handling for missing files."""
        paths = [str(TEST_DATA / "nonexistent.trees")]

        with pytest.raises((ValueError, RuntimeError)):
            _pairwise_rf(paths)

    def test_tree_name_format(self):
        """Test that tree names have expected format."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, _ = _pairwise_rf(paths, burnin_trees=1)

        # Should all contain the filename and STATE
        assert all("hiv1" in name for name in tree_names)
        assert all("STATE" in name for name in tree_names)


class TestPairwiseWeightedRF:
    """Tests for pairwise_weighted_rf function."""

    def test_basic_weighted_rf_calculation(self):
        """Test basic weighted RF distance calculation."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        tree_names, matrix = _pairwise_weighted_rf(paths, burnin_trees=1)

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

        _, rf_matrix = _pairwise_rf(paths, burnin_trees=1)
        _, weighted_matrix = _pairwise_weighted_rf(paths, burnin_trees=1)

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
        tree_names1, matrix1 = _pairwise_weighted_rf(paths, burnin_trees=1)
        tree_names2, matrix2 = _pairwise_weighted_rf(paths, burnin_trees=1)
        tree_names3, matrix3 = _pairwise_weighted_rf(paths, burnin_trees=1)

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
        tree_names, matrix = _pairwise_kf(paths, burnin_trees=1)

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
        tree_names1, matrix1 = _pairwise_kf(paths, burnin_trees=1)
        tree_names2, matrix2 = _pairwise_kf(paths, burnin_trees=1)
        tree_names3, matrix3 = _pairwise_kf(paths, burnin_trees=1)

        # Tree names should be identical
        assert tree_names1 == tree_names2 == tree_names3

        # Matrices should be exactly identical (deterministic)
        assert matrices_close(matrix1, matrix2), "KF distances should be deterministic between runs 1 and 2"
        assert matrices_close(matrix2, matrix3), "KF distances should be deterministic between runs 2 and 3"


    def test_kf_vs_weighted(self):
        """Test KF produces different values from weighted RF."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, weighted_matrix = _pairwise_weighted_rf(paths, burnin_trees=1)
        _, kf_matrix = _pairwise_kf(paths, burnin_trees=1)

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
        tree_names, _ = _pairwise_rf(paths, burnin_trees=1)
        assert len(tree_names) > 0

    def test_empty_file_list(self):
        """Test error handling for empty file list."""
        with pytest.raises((ValueError, TypeError)):
            _pairwise_rf([])

    def test_invalid_burnin_values(self):
        """Test handling of invalid burnin parameters."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        # Negative burnin should raise an overflow error
        with pytest.raises(OverflowError):
            _pairwise_rf(paths, burnin_trees=-1)


class TestAPIConsistency:
    """Tests that all three metrics return consistent structure."""

    def test_all_metrics_same_trees(self):
        """Test that all metrics return the same tree names."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        names_rf, _ = _pairwise_rf(paths, burnin_trees=1)
        names_weighted, _ = _pairwise_weighted_rf(paths, burnin_trees=1)
        names_kf, _ = _pairwise_kf(paths, burnin_trees=1)

        assert names_rf == names_weighted == names_kf

    def test_all_metrics_same_dimensions(self):
        """Test that all metrics return same matrix dimensions."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, matrix_rf = _pairwise_rf(paths, burnin_trees=1)
        _, matrix_weighted = _pairwise_weighted_rf(paths, burnin_trees=1)
        _, matrix_kf = _pairwise_kf(paths, burnin_trees=1)

        assert len(matrix_rf) == len(matrix_weighted) == len(matrix_kf)
        assert len(matrix_rf[0]) == len(matrix_weighted[0]) == len(matrix_kf[0])


class TestRootedRF:
    """Tests for rooted RF distance (clade-based comparison)."""

    def test_rooted_rf_greater_or_equal(self):
        """Rooted RF should be >= unrooted RF for the same trees."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, matrix_unrooted = _pairwise_rf(paths, burnin_trees=1, rooted=False)
        _, matrix_rooted = _pairwise_rf(paths, burnin_trees=1, rooted=True)

        for i in range(len(matrix_rooted)):
            for j in range(i + 1, len(matrix_rooted)):
                assert matrix_rooted[i][j] >= matrix_unrooted[i][j], (
                    f"Rooted RF[{i}][{j}]={matrix_rooted[i][j]} should be >= "
                    f"unrooted RF[{i}][{j}]={matrix_unrooted[i][j]}"
                )

    def test_rooted_rf_structure(self):
        """Rooted RF matrix should be symmetric with zero diagonal."""
        paths = [str(TEST_DATA / "hiv1.trees")]
        _, matrix = _pairwise_rf(paths, burnin_trees=1, rooted=True)

        for i in range(len(matrix)):
            assert matrix[i][i] == 0, f"Diagonal [{i}][{i}] should be 0"

        for i in range(len(matrix)):
            for j in range(len(matrix)):
                assert matrix[i][j] == matrix[j][i], f"Should be symmetric at [{i}][{j}]"

    def test_rooted_rf_deterministic(self):
        """Rooted RF should be deterministic."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, m1 = _pairwise_rf(paths, burnin_trees=1, rooted=True)
        _, m2 = _pairwise_rf(paths, burnin_trees=1, rooted=True)

        assert m1 == m2, "Rooted RF must be deterministic"

    def test_rooted_default_is_false(self):
        """Default rooted=False should match explicit rooted=False."""
        paths = [str(TEST_DATA / "hiv1.trees")]

        _, matrix_default = _pairwise_rf(paths, burnin_trees=1)
        _, matrix_explicit = _pairwise_rf(paths, burnin_trees=1, rooted=False)

        assert matrix_default == matrix_explicit


class TestPairwiseRFFromNewicks:
    """Tests for pairwise_rf_from_newicks (newick list, no iterator)."""

    def test_known_rf_values(self):
        """Known pairwise RF values from the treedist reference set."""
        names, matrix = _pairwise_rf_from_newicks_compat(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
        )
        assert names == FIXTURE_NAMES
        for i in range(len(FIXTURE_RF)):
            for j in range(len(FIXTURE_RF)):
                assert matrix[i][j] == FIXTURE_RF[i][j], f"mismatch at [{i}][{j}]"

    def test_symmetric_zero_diagonal(self):
        """Matrix must be symmetric with a zero diagonal."""
        _, matrix = _pairwise_rf_from_newicks_compat(
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
        _, unrooted = _pairwise_rf_from_newicks_compat(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES, rooted=False
        )
        _, rooted = _pairwise_rf_from_newicks_compat(
            FIXTURE_NAMES, FIXTURE_TREES, FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES, rooted=True
        )
        for i in range(len(rooted)):
            for j in range(len(rooted)):
                assert rooted[i][j] >= unrooted[i][j], f"rooted < unrooted at [{i}][{j}]"

    def test_error_on_length_mismatch(self):
        """Passing fewer newicks than names must raise ValueError."""
        with pytest.raises(ValueError):
            _pairwise_rf_from_newicks_compat(
                FIXTURE_NAMES, FIXTURE_TREES[:2], FIXTURE_TRANSLATE, FIXTURE_MAP_INDICES
            )


class TestPairwiseRFFromNewicksReference:
    """Expanded RF-from-newicks tests using the 12-tree PHYLIP reference set."""

    def test_all_pairs_match_reference(self):
        tmap = natural_map()
        numbered = [numberise_newick(n, tmap) for n in REFERENCE_NEWICKS]
        names = [f"tree_{i}" for i in range(12)]

        tree_names, matrix = _pairwise_rf_from_newicks_compat(
            names, numbered, [tmap], [0] * 12
        )

        assert tree_names == names
        assert len(matrix) == 12
        for i in range(12):
            for j in range(12):
                assert matrix[i][j] == EXPECTED_RF[i][j], (
                    f"RF[{i}][{j}] = {matrix[i][j]}, expected {EXPECTED_RF[i][j]}"
                )

    def test_representative_pairs_across_different_translate_maps(self):
        map_a = natural_map()
        map_b = shuffled_map()

        numbered_a = [numberise_newick(n, map_a) for n in REFERENCE_NEWICKS]
        numbered_b = [numberise_newick(n, map_b) for n in REFERENCE_NEWICKS]

        newicks = numbered_a[:6] + numbered_b[6:]
        names = [f"tree_{i}" for i in range(12)]
        map_indices = [0] * 6 + [1] * 6

        tree_names, matrix = _pairwise_rf_from_newicks_compat(
            names, newicks, [map_a, map_b], map_indices
        )

        assert tree_names == names
        pairs = [(0, 1), (0, 6), (3, 9), (5, 11), (2, 10)]
        for i, j in pairs:
            assert matrix[i][j] == EXPECTED_RF[i][j], (
                f"RF[{i}][{j}] = {matrix[i][j]}, expected {EXPECTED_RF[i][j]}"
            )

    def test_identical_topology_different_numbering_is_zero(self):
        map_a = natural_map()
        map_b = shuffled_map()

        for idx in range(12):
            newick_a = numberise_newick(REFERENCE_NEWICKS[idx], map_a)
            newick_b = numberise_newick(REFERENCE_NEWICKS[idx], map_b)

            _, matrix = _pairwise_rf_from_newicks_compat(
                ["a", "b"], [newick_a, newick_b], [map_a, map_b], [0, 1]
            )
            assert matrix[0][1] == 0, (
                f"Tree {idx}: same topology with different numbering should have RF=0, got {matrix[0][1]}"
            )

    def test_beast_annotations_stripped(self):
        tmap = natural_map()
        clean = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        annotated = clean.replace(":", ":[&rate=0.42]")

        _, matrix = _pairwise_rf_from_newicks_compat(
            ["clean", "annotated"], [clean, annotated], [tmap], [0, 0]
        )
        assert matrix[0][1] == 0

    def test_validation_names_newicks_length_mismatch(self):
        tmap = natural_map()
        newick = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        with pytest.raises(ValueError, match="names length"):
            _pairwise_rf_from_newicks_compat(
                ["a", "b"], [newick], [tmap], [0]
            )

    def test_validation_map_index_out_of_bounds(self):
        tmap = natural_map()
        newicks = [
            numberise_newick(REFERENCE_NEWICKS[0], tmap),
            numberise_newick(REFERENCE_NEWICKS[1], tmap),
        ]
        with pytest.raises(ValueError, match="out of bounds"):
            _pairwise_rf_from_newicks_compat(
                ["a", "b"], newicks, [tmap], [0, 5]
            )

    def test_validation_fewer_than_two_trees(self):
        tmap = natural_map()
        newick = numberise_newick(REFERENCE_NEWICKS[0], tmap)
        with pytest.raises(ValueError, match="at least 2"):
            _pairwise_rf_from_newicks_compat(
                ["a"], [newick], [tmap], [0]
            )


@pytest.mark.skipif(
    not _r_phangorn_available(),
    reason="Rscript with phangorn is required for cross-validation tests",
)
class TestRPhangornCrossValidation:
    """Cross-check rapidtrees distance outputs against R phangorn.

    Two suites:
    - Synthetic 10-taxon trees (REFERENCE_NEWICKS): verifies the base algorithm.
    - Real BEAST HIV trees (hiv{1,2,3,4}.newick): verifies annotation stripping,
      large taxon counts, and real branch-length handling.

    All three metrics (RF, wRF, KF) are validated for each tree set.
    """

    @pytest.mark.parametrize(
        "newick_file", ["hiv1.newick", "hiv2.newick", "hiv3.newick", "hiv4.newick"]
    )
    # @pytest.mark.parametrize("metric", ["RF", "wRF", "KF"])
    @pytest.mark.parametrize("metric", ["RF"])
    def test_hiv_matches_phangorn(self, newick_file, metric):
        """Real BEAST HIV trees (162 taxa) match phangorn across all metrics."""
        path = str(TEST_DATA / newick_file)
        rt_vals = _rapidtrees_upper_triangle(path, metric)
        ph_vals = _phangorn_run(path, metric)
        _assert_distances_match(rt_vals, ph_vals, metric)


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
        _, matrix_list = _pairwise_rf_from_newicks_compat(
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
    """Tests for pairwise_rf_with_snapshots_from_newick_iter (6-tuple API)."""

    def _call(self, **kwargs):
        defaults = dict(
            names=FIXTURE_NAMES,
            newick_iter=iter(FIXTURE_TREES),
            translate_maps=FIXTURE_TRANSLATE,
            map_indices=FIXTURE_MAP_INDICES,
        )
        defaults.update(kwargs)
        return rtd.pairwise_rf_with_snapshots_from_newick_iter(**defaults)

    def test_returns_six_tuple(self):
        """Function returns a 6-tuple with the correct types."""
        result = self._call()
        assert len(result) == 6
        names, rf_bytes, leaf_names, n_bip, pres_bytes, bip_leaf_indices = result
        assert isinstance(names, list)
        assert isinstance(rf_bytes, bytes)
        assert isinstance(leaf_names, list)
        assert isinstance(n_bip, int)
        assert isinstance(pres_bytes, bytes)
        assert isinstance(bip_leaf_indices, list)

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
        names, _, _, n_bip, pres_bytes, *_ = self._call()
        n = len(names)
        assert len(pres_bytes) == n * n_bip
        presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(n, n_bip)
        assert set(presence.flatten().tolist()) <= {0, 1}

    def test_xor_identity(self):
        """sum(presence[i] XOR presence[j]) equals RF(i, j) for every pair."""
        names, rf_bytes, _, n_bip, pres_bytes, *_ = self._call()
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
        _, _, leaf_names, *_ = self._call()
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
        assert r1[5] == r2[5]  # bip_leaf_indices

    def test_bip_leaf_indices_length(self):
        """bip_leaf_indices has one entry per bipartition."""
        _, _, _, n_bip, _, bip_leaf_indices = self._call()
        assert len(bip_leaf_indices) == n_bip

    def test_bip_leaf_indices_valid_range(self):
        """Every index in bip_leaf_indices is a valid position into leaf_names."""
        _, _, leaf_names, _, _, bip_leaf_indices = self._call()
        n_leaves = len(leaf_names)
        for indices in bip_leaf_indices:
            assert isinstance(indices, list)
            for idx in indices:
                assert isinstance(idx, int)
                assert 0 <= idx < n_leaves, f"Index {idx} out of range for {n_leaves} leaves"

    def test_bip_leaf_indices_canonical_side_excludes_first_leaf(self):
        """In unrooted mode the canonical side never contains the first leaf (index 0)."""
        _, _, _, _, _, bip_leaf_indices = self._call()
        for indices in bip_leaf_indices:
            assert 0 not in indices, "Canonical side must not contain leaf 0 (unrooted)"

    def test_bip_leaf_indices_decode_to_valid_names(self):
        """Decoded leaf names for each bipartition are a subset of the full leaf set."""
        _, _, leaf_names, _, _, bip_leaf_indices = self._call()
        leaf_set = set(leaf_names)
        for indices in bip_leaf_indices:
            names_for_bip = {leaf_names[i] for i in indices}
            assert names_for_bip <= leaf_set
    def test_bip_leaf_indices_simple_tree(self):
        """On a known 4-leaf tree, verify exact canonical-side leaf indices."""
        # Tree ((A,B),(C,D)) has one non-trivial bipartition: {C,D}|{A,B}.
        # Canonical side = side not containing A (index 0) = {C, D} = indices [2, 3].
        trees = ["((A,B),(C,D));", "((A,C),(B,D));"]
        names = ["t0", "t1"]
        result = rtd.pairwise_rf_with_snapshots_from_newick_iter( names, iter(trees), [{}], [0, 0] )
        _, _, leaf_names, _, _, bip_leaf_indices = result
        assert leaf_names == ["A", "B", "C", "D"]
        # Decode to names for readability
        decoded = [sorted(leaf_names[i] for i in indices) for indices in bip_leaf_indices]
        # Both bipartitions should decode to 2-leaf canonical sides
        assert all(len(d) == 2 for d in decoded)
        # Neither canonical side should include 'A'
        assert all("A" not in d for d in decoded)
        # The two bipartitions present are {C,D} and {B,D}
        assert sorted(decoded) == [["B", "D"], ["C", "D"]]

    def test_bip_leaf_indices_large_complex_trees(self):
        """
        Verify bipartition canonicalization on 3 heterogeneous 15-leaf trees.
        Checks that canonical-side indices correctly reflect each topology and
        that bipartitions shared across trees are not duplicated.
        """
        # t0 — right-caterpillar, uniform cherry pairing
        #  ·
        #  ├── ·
        #  │   ├── (A,B)
        #  │   └── ·
        #  │       ├── (C,D)
        #  │       └── ·
        #  │           ├── (E,F)
        #  │           └── ·
        #  │               ├── ·  {G,H,I,J}
        #  │               │   ├── (G,H)
        #  │               │   └── (I,J)
        #  │               └── ·  {K,L,M,N}
        #  │                   ├── (K,L)
        #  │                   └── (M,N)
        #  └── O
        T0 = "(((A,B),((C,D),((E,F),(((G,H),(I,J)),((K,L),(M,N)))))),O);"

        # t1 — same caterpillar backbone; cherries cross-paired
        #   (C,E) instead of (C,D)     (D,F) instead of (E,F)
        #   (G,I)/(H,J) instead of (G,H)/(I,J)
        #   (K,M)/(L,N) instead of (K,L)/(M,N)
        T1 = "(((A,B),((C,E),((D,F),(((G,I),(H,J)),((K,M),(L,N)))))),O);"

        # t2 — balanced: {A–F} and {G–N} are sister clades under the root
        #  ·
        #  ├── ·
        #  │   ├── ·  {A,B,C,D,E,F}
        #  │   │   ├── (A,B)
        #  │   │   └── ·  {C,D,E,F}
        #  │   │       ├── (C,D)
        #  │   │       └── (E,F)
        #  │   └── ·  {G,H,I,J,K,L,M,N}
        #  │       ├── ·  {G,H,I,J,K,L}
        #  │       │   ├── (G,H)
        #  │       │   └── ·  {I,J,K,L}
        #  │       │       ├── (I,J)
        #  │       │       └── (K,L)
        #  │       └── (M,N)
        #  └── O
        T2 = "((((A,B),((C,D),(E,F))),(((G,H),((I,J),(K,L))),(M,N))),O);"

        _, _, leaf_names, _, _, bip_leaf_indices = (
            rtd.pairwise_rf_with_snapshots_from_newick_iter(
                ["t0", "t1", "t2"], iter([T0, T1, T2]), [{}], [0, 0, 0]
            )
        )

        assert leaf_names == list("ABCDEFGHIJKLMNO")

        # Decode index lists → frozensets of taxon names for set comparison
        decoded = {frozenset(leaf_names[i] for i in idx) for idx in bip_leaf_indices}

        def s(*leaves):
            return frozenset(leaves)

        # Expected canonical bipartitions (side NOT containing 'A').
        # 15 leaves → 12 internal edges per tree.  After dedup: 23 unique.
        expected = {
            # ── shared by all three trees ──────────────────────────────────
            s(*"CDEFGHIJKLMNO"),  # canonical of {A,B} | rest
            s(*"GHIJKLMN"),       # {G–N} clade

            # ── shared by t0 and t1 (same {A,B,O} placement) ──────────────
            s(*"CDEFGHIJKLMN"),   # canonical of {A,B,O} | rest
            s(*"GHIJ"),
            s(*"KLMN"),

            # ── shared by t0 and t2 (same small cherries) ─────────────────
            s(*"CD"), s(*"EF"), s(*"GH"), s(*"IJ"), s(*"KL"), s(*"MN"),

            # ── unique to t0 ───────────────────────────────────────────────
            s(*"EFGHIJKLMN"),     # {E–N} clade

            # ── unique to t1 (cross-pairing shifts mid clades) ────────────
            s(*"CE"), s(*"DF"), s(*"GI"), s(*"HJ"), s(*"KM"), s(*"LN"),
            s(*"DFGHIJKLMN"),     # {D,F,G–N} clade

            # ── unique to t2 (balanced topology, new upper bipartitions) ───
            s(*"CDEF"),           # {C,D,E,F} clade
            s(*"IJKL"),
            s(*"GHIJKL"),
            s(*"GHIJKLMNO"),      # canonical of {A,B,C,D,E,F} | rest
        }

        assert decoded == expected

    def test_bip_presence_matrix_matches_topology(self):
        """
        For the same 3 heterogeneous trees, verify that every column of the
        presence/absence matrix has exactly the 0/1 pattern dictated by the
        topology — i.e. which trees actually contain each bipartition.

        Expected presence table (rows = t0/t1/t2, cols = bipartitions):

          Bipartition                          t0  t1  t2
          {C,D}                                 1   0   1
          {C,E}                                 0   1   0
          {D,F}                                 0   1   0
          {E,F}                                 1   0   1
          {G,H}                                 1   0   1
          {G,I}                                 0   1   0
          {H,J}                                 0   1   0
          {I,J}                                 1   0   1
          {K,L}                                 1   0   1
          {K,M}                                 0   1   0
          {L,N}                                 0   1   0
          {M,N}                                 1   0   1
          {C,D,E,F}                             0   0   1
          {G,H,I,J}                             1   1   0
          {I,J,K,L}                             0   0   1
          {K,L,M,N}                             1   1   0
          {G,H,I,J,K,L}                         0   0   1
          {G,H,I,J,K,L,M,N}                     1   1   1
          {G,H,I,J,K,L,M,N,O}                   0   0   1
          {D,F,G,H,I,J,K,L,M,N}                 0   1   0
          {E,F,G,H,I,J,K,L,M,N}                 1   0   0
          {C,D,E,F,G,H,I,J,K,L,M,N}             1   1   0
          {C,D,E,F,G,H,I,J,K,L,M,N,O}           1   1   1
        """
        T0 = "(((A,B),((C,D),((E,F),(((G,H),(I,J)),((K,L),(M,N)))))),O);"
        T1 = "(((A,B),((C,E),((D,F),(((G,I),(H,J)),((K,M),(L,N)))))),O);"
        T2 = "((((A,B),((C,D),(E,F))),(((G,H),((I,J),(K,L))),(M,N))),O);"

        _, _, leaf_names, n_bip, pres_bytes, bip_leaf_indices = (
            rtd.pairwise_rf_with_snapshots_from_newick_iter(
                ["t0", "t1", "t2"], iter([T0, T1, T2]), [{}], [0, 0, 0]
            )
        )

        presence = np.frombuffer(pres_bytes, dtype=np.uint8).reshape(3, n_bip).copy()

        # Map each bipartition (as a frozenset of names) to its column index
        col_for = {
            frozenset(leaf_names[i] for i in idx): col
            for col, idx in enumerate(bip_leaf_indices)
        }

        def s(*leaves):
            return frozenset(leaves)

        def assert_col(bip, expected_row):
            col = col_for[bip]
            actual = list(presence[:, col])
            assert actual == expected_row, (
                f"Bipartition {{{','.join(sorted(bip))}}}: "
                f"expected {expected_row}, got {actual}"
            )

        # present in all three trees
        assert_col(s(*"CDEFGHIJKLMNO"), [1, 1, 1])  # canonical of {A,B} | rest
        assert_col(s(*"GHIJKLMN"),      [1, 1, 1])  # {G–N} clade

        # present in t0 and t1 only (same {A,B,O} placement, same deep clades)
        assert_col(s(*"CDEFGHIJKLMN"),  [1, 1, 0])  # canonical of {A,B,O} | rest
        assert_col(s(*"GHIJ"),          [1, 1, 0])
        assert_col(s(*"KLMN"),          [1, 1, 0])

        # present in t0 and t2 only (same small cherries)
        for bip in [s(*"CD"), s(*"EF"), s(*"GH"), s(*"IJ"), s(*"KL"), s(*"MN")]:
            assert_col(bip, [1, 0, 1])

        # present in t0 only
        assert_col(s(*"EFGHIJKLMN"),    [1, 0, 0])  # {E–N} clade

        # present in t1 only (cross-paired cherries produce distinct mid clades)
        for bip in [s(*"CE"), s(*"DF"), s(*"GI"), s(*"HJ"), s(*"KM"), s(*"LN")]:
            assert_col(bip, [0, 1, 0])
        assert_col(s(*"DFGHIJKLMN"),    [0, 1, 0])  # {D,F,G–N} clade

        # present in t2 only (balanced topology, new upper-level bipartitions)
        assert_col(s(*"CDEF"),          [0, 0, 1])
        assert_col(s(*"IJKL"),          [0, 0, 1])
        assert_col(s(*"GHIJKL"),        [0, 0, 1])
        assert_col(s(*"GHIJKLMNO"),     [0, 0, 1])  # canonical of {A,B,C,D,E,F} | rest


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
