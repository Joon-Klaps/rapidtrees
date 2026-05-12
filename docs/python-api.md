# Python API reference

`rapidtrees` exposes six functions from its Rust core via PyO3. All accept a
Python **iterator** of newick strings, which lets the library stream through
arbitrarily large tree files without materialising all strings in memory at
once.

---

## Functions

| Function | Returns |
| --- | --- |
| `pairwise_rf_from_newick_iter` | `(names, bytes)` — RF matrix as flat `uint32` bytes, row-major |
| `pairwise_rf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — RF matrix + bipartition presence matrix + bipartition clade bitmasks |
| `pairwise_wrf_from_newick_iter` | `(names, list[float])` — Weighted RF, flat row-major |
| `pairwise_wrf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — wRF matrix + branch-length matrix + bipartition clade bitmasks |
| `pairwise_kf_from_newick_iter` | `(names, list[float])` — Kuhner-Felsenstein, flat row-major |
| `pairwise_kf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — KF matrix + branch-length matrix + bipartition clade bitmasks |

All six share the same call signature:

```python
func(
    names,           # list[str]            — one identifier per tree
    newick_iter,     # Iterator[str]        — newick strings, consumed once left to right
    translate_maps,  # list[dict[str, str]] — taxon-ID → name mappings (BEAST TRANSLATE block)
    map_indices,     # list[int]            — which translate map applies to each tree
    rooted=False,    # bool                 — compare clades (True) or bipartitions (False)
)
```

`translate_maps` and `map_indices` handle BEAST numeric taxon IDs.  If your
newick strings already use real taxon names, pass `[{}]` and `[0] * n`.

### Return types

| Function | Return type | How to decode |
| --- | --- | --- |
| `pairwise_rf_from_newick_iter` | `bytes` — flat `uint32`, row-major | `np.frombuffer(b, dtype=np.uint32).reshape(n, n)` |
| `pairwise_wrf_from_newick_iter` | `list[float]` — flat, row-major | `np.array(lst, dtype=np.float64).reshape(n, n)` |
| `pairwise_kf_from_newick_iter` | `list[float]` — flat, row-major | `np.array(lst, dtype=np.float64).reshape(n, n)` |
| `pairwise_rf_with_snapshots_from_newick_iter` | 6-tuple — see below | see below |
| `pairwise_wrf_with_snapshots_from_newick_iter` | 6-tuple — see below | see below |
| `pairwise_kf_with_snapshots_from_newick_iter` | 6-tuple — see below | see below |

### Errors raised

All functions raise `ValueError` when:
- Fewer than 2 trees are provided
- `len(names) != len(map_indices)`
- A `map_indices` value is out of bounds
- Trees have different leaf sets

---

## Examples

### Newick strings (no file I/O)

```python
import rapidtrees as rtd
import numpy as np

trees = [
    "(A:0.1,(B:0.1,C:0.1):0.1);",
    "(A:0.1,(C:0.1,B:0.1):0.1);",
    "((A:0.1,B:0.1):0.1,C:0.1);",
]
names = ["t1", "t2", "t3"]

# RF — raw uint32 bytes
tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
    names, iter(trees), [{}], [0, 0, 0]
)
n = len(tree_names)
rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)

# Weighted RF — flat list of floats
tree_names, wrf_flat = rtd.pairwise_wrf_from_newick_iter(
    names, iter(trees), [{}], [0, 0, 0]
)
wrf = np.array(wrf_flat, dtype=np.float64).reshape(n, n)

# Kuhner-Felsenstein
tree_names, kf_flat = rtd.pairwise_kf_from_newick_iter(
    names, iter(trees), [{}], [0, 0, 0]
)
kf = np.array(kf_flat, dtype=np.float64).reshape(n, n)
```

---

### BEAST `.trees` files

The Rust API does not expose file reading — parse your `.trees` file in Python
and feed the newick strings to the iterator API.  The translate maps carry the
BEAST `TRANSLATE` block so that numeric taxon IDs are resolved to real names.

```python
import re
from pathlib import Path
import rapidtrees as rtd
import numpy as np


def load_beast(path, burnin_trees=0, use_real_taxa=True):
    """Parse a BEAST .trees file. Returns (translate_map, [(name, newick)])."""
    content = Path(path).read_text()
    base = Path(path).stem

    translate = {}
    if use_real_taxa:
        in_block = False
        for line in content.splitlines():
            s = line.strip()
            if s.upper().startswith("TRANSLATE"):
                in_block = True
                continue
            if in_block:
                if s.startswith(";"):
                    break
                parts = s.rstrip(",").split()
                if len(parts) >= 2:
                    translate[parts[0]] = parts[1].strip("'")

    pairs = []
    for idx, line in enumerate(content.splitlines()):
        upper = line.upper().lstrip()
        if upper.startswith("TREE ") and " = " in line:
            header, newick = line.split(" = ", 1)
            m = re.search(r"STATE_(\d+)", header, re.IGNORECASE)
            name = f"{base}_{header.split()[-1]}" if m else f"{base}_tree_{idx}"
            newick = re.sub(r"\[&[^\]]*\]", "", newick.strip())  # strip BEAST annotations
            pairs.append((name, newick))

    return translate, pairs[burnin_trees:]


tmap, tree_pairs = load_beast("run1.trees", burnin_trees=100)
names, newicks = zip(*tree_pairs)

tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
    list(names), iter(newicks), [tmap], [0] * len(names)
)
n = len(tree_names)
rf = np.frombuffer(rf_bytes, dtype=np.uint32).reshape(n, n)
```

#### Multiple files

When combining trees from several files, give each file its own translate map
and use `map_indices` to say which map applies to each tree:

```python
all_names, all_newicks, translate_maps, map_indices = [], [], [], []

for path in ["run1.trees", "run2.trees"]:
    tmap, pairs = load_beast(path, burnin_trees=100)
    idx = len(translate_maps)
    translate_maps.append(tmap)
    for name, newick in pairs:
        all_names.append(name)
        all_newicks.append(newick)
        map_indices.append(idx)

tree_names, rf_bytes = rtd.pairwise_rf_from_newick_iter(
    all_names, iter(all_newicks), translate_maps, map_indices
)
```

---

### RF + bipartition snapshot in one pass

`pairwise_rf_with_snapshots_from_newick_iter` builds both the RF distance
matrix **and** the bipartition presence matrix in a single parse, returning a
6-tuple:

```
(tree_names, rf_bytes, leaf_names, n_bip, presence_bytes, bipartition_clade_bytes)
```

| Field | Type | Description |
| --- | --- | --- |
| `tree_names` | `list[str]` | Tree identifiers (same order as input `names`) |
| `rf_bytes` | `bytes` | Flat `uint32` RF matrix, row-major, shape `(n, n)` |
| `leaf_names` | `list[str]` | Sorted taxon names — index `i` corresponds to bit `i` in every bipartition |
| `n_bip` | `int` | Number of unique edges across all trees (internal bipartitions + pendant edges) |
| `presence_bytes` | `bytes` | Flat `uint8` presence matrix, row-major, shape `(n, n_bip)` |
| `bipartition_clade_bytes` | `bytes` | Packed bitmasks, shape `(n_bip, ceil(n_leaves/8))` — see below |

The presence matrix entry `presence[i, j]` is `1` if edge `j` appears in tree
`i`, otherwise `0`. Column order is deterministic (ascending `Bitset` order)
and stable across calls on the same tree set.

#### Edge table contents

`n_bip` counts **all** edges: pendant (leaf) edges and internal bipartitions.
Pendant edges appear as single-bit rows (one bit set), sorted before internal
bipartitions. Since all trees share the same leaf set, pendant columns are always
`1` in every row of the presence matrix.

#### Canonicalisation note

For internal bipartitions (rows with ≥ 2 bits set) the canonical side is the
half that does **not** contain the first leaf alphabetically, so bit 0 is never
set in those rows. Pendant edges are stored verbatim (no flip), so the pendant
of the first leaf has bit 0 set. The complement of an internal bipartition can
be derived as `~bip_bool[j] & True` (masked to `n_leaves` bits).

#### Bipartition clade bytes format

`bipartition_clade_bytes` is a flat `bytes` buffer of shape
`(n_bip, ceil(n_leaves / 8))`. Each row encodes the canonical leaf membership
of one bipartition as a **little-endian bitmask**: bit `i` within each row is
`1` if `leaf_names[i]` is on that bipartition's canonical side.

Decode with NumPy:

```python
import math
import numpy as np

bytes_per_bip = math.ceil(len(leaf_names) / 8)
bip_arr  = np.frombuffer(bipartition_clade_bytes, dtype=np.uint8).reshape(n_bip, bytes_per_bip)
bip_bool = np.unpackbits(bip_arr, axis=1, bitorder='little')[:, :len(leaf_names)]
# bip_bool[j, i] == 1  →  leaf_names[i] is in bipartition j
```

```python
tree_names, rf_bytes, leaf_names, n_bip, pres_bytes, bip_clade_bytes = (
    rtd.pairwise_rf_with_snapshots_from_newick_iter(
        list(names), iter(newicks), [tmap], [0] * len(names)
    )
)
n = len(tree_names)
rf       = np.frombuffer(rf_bytes,   dtype=np.uint32).reshape(n, n)
presence = np.frombuffer(pres_bytes, dtype=np.uint8 ).reshape(n, n_bip).copy()

# Verify: sum(|row_i − row_j|) == RF(i, j) for all pairs
for i in range(n):
    for j in range(n):
        assert int(np.sum(np.abs(presence[i].astype(int) - presence[j].astype(int)))) == int(rf[i, j])

# Global split frequencies across all trees (useful for Pseudo-ESS / ASDSF)
split_freq = presence.mean(axis=0)
```

#### Named presence matrix (post-hoc analysis)

Decode `bipartition_clade_bytes` to build human-readable column labels for the
presence matrix, enabling downstream analyses such as tanglegrams, identifying
unstable splits, or computing per-clade frequencies:

```python
import math
import pandas as pd

bytes_per_bip = math.ceil(len(leaf_names) / 8)
bip_arr  = np.frombuffer(bip_clade_bytes, dtype=np.uint8).reshape(n_bip, bytes_per_bip)
bip_bool = np.unpackbits(bip_arr, axis=1, bitorder='little')[:, :len(leaf_names)]

# Build a human-readable column label for each bipartition
col_labels = [
    "|".join(n for i, n in enumerate(leaf_names) if bip_bool[j, i])
    for j in range(n_bip)
]

df = pd.DataFrame(presence, index=tree_names, columns=col_labels)
# e.g.  col "C|D|E" == 1 means the split {C,D,E}|rest is present in that tree
```

---

### wRF + branch-length matrix in one pass

`pairwise_wrf_with_snapshots_from_newick_iter` builds both the pairwise wRF
distance matrix **and** a per-edge branch-length matrix in a single parse,
returning a 6-tuple:

```text
(tree_names, wrf_bytes, leaf_names, n_bip, branch_length_bytes, bipartition_clade_bytes)
```

| Field | Type | Description |
| --- | --- | --- |
| `tree_names` | `list[str]` | Tree identifiers (same order as input `names`) |
| `wrf_bytes` | `bytes` | Flat `float64` wRF matrix, row-major, shape `(n, n)` |
| `leaf_names` | `list[str]` | Sorted taxon names — index `i` corresponds to bit `i` in every bipartition |
| `n_bip` | `int` | Number of unique edges across all trees (pendant + internal) |
| `branch_length_bytes` | `bytes` | Flat `float64`, shape `(n_trees, n_bip)`, row-major |
| `bipartition_clade_bytes` | `bytes` | Packed bitmasks, shape `(n_bip, ceil(n_leaves/8))` — identical to RF snapshot |

`branch_length_bytes[i, j]` is the branch length of edge `j` in tree `i`, or
`0.0` if that edge is absent. Pendant (leaf-edge) columns are always non-zero
because every tree has every leaf.

Column order matches `bipartition_clade_bytes` (ascending `Bitset` order,
deterministic and stable across calls on the same tree set).

#### Decode and compute Fréchet ESS traces

```python
import rapidtrees as rtd
import numpy as np

tree_names, wrf_bytes, leaf_names, n_bip, bl_bytes, bip_clade_bytes = (
    rtd.pairwise_wrf_with_snapshots_from_newick_iter(
        list(names), iter(newicks), [tmap], [0] * len(names)
    )
)
n = len(tree_names)
# Actual wRF distance matrix (trees x trees)
wrf = np.frombuffer(wrf_bytes, dtype=np.float64).reshape(n, n)

# Branch length matrix (trees x bipartitions)
bl  = np.frombuffer(bl_bytes,  dtype=np.float64).reshape(n, n_bip)

# We can recompute distances relative to one reference tree using the branch-length matrix:
ref_idx = 0
wrf_trace = np.sum(np.abs(bl[ref_idx, :] - bl), axis=1)           # shape (n,)
kf_trace  = np.sqrt(np.sum((bl[ref_idx, :] - bl) ** 2, axis=1))   # shape (n,)

# Verify L1 identity: sum(|bl[i,:] - bl[j,:]|) == wrf[i,j]
for i in range(n):
    for j in range(n):
        assert abs(float(np.sum(np.abs(bl[i] - bl[j]))) - wrf[i, j]) < 1e-9
```

---

### KF + branch-length matrix in one pass

`pairwise_kf_with_snapshots_from_newick_iter` is identical to the wRF variant
except the distance matrix uses the Kuhner–Felsenstein (Branch Score) metric.
The branch-length matrix is metric-agnostic and bit-for-bit identical across
both functions for the same input.

```text
(tree_names, kf_bytes, leaf_names, n_bip, branch_length_bytes, bipartition_clade_bytes)
```

| Field | Type | Description |
| --- | --- | --- |
| `tree_names` | `list[str]` | Tree identifiers |
| `kf_bytes` | `bytes` | Flat `float64` KF matrix, row-major, shape `(n, n)` |
| `leaf_names` | `list[str]` | Sorted taxon names |
| `n_bip` | `int` | Number of unique edges |
| `branch_length_bytes` | `bytes` | Flat `float64`, shape `(n_trees, n_bip)` — same as wRF variant |
| `bipartition_clade_bytes` | `bytes` | Packed bitmasks — same as RF/wRF snapshot |

```python
tree_names, kf_bytes, leaf_names, n_bip, bl_bytes, bip_clade_bytes = (
    rtd.pairwise_kf_with_snapshots_from_newick_iter(
        list(names), iter(newicks), [tmap], [0] * len(names)
    )
)
n = len(tree_names)

# Actual KF distance matrix (trees x trees)
kf = np.frombuffer(kf_bytes, dtype=np.float64).reshape(n, n)

# Branch length matrix (trees x bipartitions)
bl  = np.frombuffer(bl_bytes, dtype=np.float64).reshape(n, n_bip)

# L2 identity: sqrt(sum((bl[i,:]-bl[j,:])**2)) == kf[i,j]
```
