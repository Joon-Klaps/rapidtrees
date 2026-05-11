# Python API reference

`rapidtrees` exposes four functions from its Rust core via PyO3. All accept a
Python **iterator** of newick strings, which lets the library stream through
arbitrarily large tree files without materialising all strings in memory at
once.

---

## Functions

| Function | Returns |
| --- | --- |
| `pairwise_rf_from_newick_iter` | `(names, bytes)` — RF matrix as flat `uint32` bytes, row-major |
| `pairwise_rf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, list[list[int]])` — RF matrix + bipartition presence matrix + bipartition leaf indices |
| `pairwise_wrf_from_newick_iter` | `(names, list[float])` — Weighted RF, flat row-major |
| `pairwise_kf_from_newick_iter` | `(names, list[float])` — Kuhner-Felsenstein, flat row-major |

All four share the same call signature:

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
(tree_names, rf_bytes, leaf_names, n_bip, presence_bytes, bipartition_leaf_indices)
```

| Field | Type | Description |
| --- | --- | --- |
| `tree_names` | `list[str]` | Tree identifiers (same order as input `names`) |
| `rf_bytes` | `bytes` | Flat `uint32` RF matrix, row-major, shape `(n, n)` |
| `leaf_names` | `list[str]` | Sorted taxon names — index `i` corresponds to bit `i` in every bipartition |
| `n_bip` | `int` | Number of unique bipartitions across all trees |
| `presence_bytes` | `bytes` | Flat `uint8` presence matrix, row-major, shape `(n, n_bip)` |
| `bipartition_leaf_indices` | `list[list[int]]` | For each bipartition column, the indices into `leaf_names` on its canonical side |

The presence matrix entry `presence[i, j]` is `1` if bipartition `j` appears
in tree `i`, otherwise `0`. Column order is deterministic (ascending `Bitset`
order) and stable across calls on the same tree set.

#### Canonicalisation note

For unrooted trees each bipartition is stored on the side that does **not**
contain the first leaf alphabetically. This means every entry in
`bipartition_leaf_indices` describes the half of the split that excludes the
first taxon — the complement can be derived via
`set(range(len(leaf_names))) - set(indices)`.

```python
tree_names, rf_bytes, leaf_names, n_bip, pres_bytes, bip_leaf_indices = (
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

`bipartition_leaf_indices` maps each column of the presence matrix to the leaf
names that make up that bipartition's canonical side.  Use it to build a
labelled pandas `DataFrame` for downstream analyses such as tanglegrams,
identifying unstable splits, or computing per-clade frequencies:

```python
import pandas as pd

# Build a human-readable column label for each bipartition
col_labels = [
    "|".join(leaf_names[i] for i in indices)
    for indices in bip_leaf_indices
]

df = pd.DataFrame(presence, index=tree_names, columns=col_labels)
# e.g.  col "C|D|E" == 1 means the split {C,D,E}|rest is present in that tree
```
