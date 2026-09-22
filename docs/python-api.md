# Python API reference

`rapidtrees` exposes iterator-based functions from its Rust core via PyO3. They
accept a Python **iterator** of newick strings, which lets callers stream source
text rather than materialising every input string solely for ingestion.
Returned distance and snapshot buffers still scale with the numbers of trees
and clades.

> **On exactness.** Splits are identified by a 128-bit fingerprint rather than
> by comparing leaf sets — that is what keeps building a tree linear in its
> taxon count instead of quadratic. Two *distinct* splits are therefore merged
> if their fingerprints collide, with probability about `e² / 2¹²⁹` for `e`
> distinct splits in the run: `1.5 × 10⁻²³` at a hundred million splits, some
> nineteen orders of magnitude below the rate at which the machine's own RAM
> flips a bit unnoticed. Results are otherwise deterministic — the same input
> gives the same split IDs and the same column order on every run and every
> machine.

---

## Functions

| Function | Returns |
| --- | --- |
| `pairwise_rf_from_newick_iter` | `(names, bytes)` — RF matrix as flat `uint32` bytes, row-major |
| `pairwise_rf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — RF matrix + bipartition presence matrix + bipartition clade bitmasks |
| `pairwise_rf_with_sparse_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, sparse)` — RF matrix + clade bitmasks + CSR presence rows |
| `pairwise_rf_with_rooted_facts_from_newick_iter` | `(names, bytes, leaf_names, n_clades, bytes, facts)` — rooted RF + compact MrHIPSTR facts |
| `pairwise_wrf_from_newick_iter` | `(names, list[float])` — Weighted RF, flat row-major |
| `pairwise_wrf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — wRF matrix + branch-length matrix + bipartition clade bitmasks |
| `pairwise_kf_from_newick_iter` | `(names, list[float])` — Kuhner-Felsenstein, flat row-major |
| `pairwise_kf_with_snapshots_from_newick_iter` | `(names, bytes, leaf_names, n_bip, bytes, bytes)` — KF matrix + branch-length matrix + bipartition clade bitmasks |

All general RF/WRF/KF endpoints share this call signature:

```python
func(
    names,            # list[str]            — one identifier per tree
    newick_iter,      # Iterator[str]        — newick strings, consumed once left to right
    translate_maps,   # list[dict[str, str]] — taxon-ID → name mappings (BEAST TRANSLATE block)
    map_indices,      # list[int]            — which translate map applies to each tree
    rooted=False,     # bool                 — compare clades (True) or bipartitions (False)
    progress=None,    # ProgressCounter | None — see "Progress reporting" below
)
```

`pairwise_rf_with_rooted_facts_from_newick_iter` is rooted by definition and
therefore omits the `rooted` argument:

```python
pairwise_rf_with_rooted_facts_from_newick_iter(
    names,
    newick_iter,
    translate_maps,
    map_indices,
    progress=None,
)
```

`translate_maps` and `map_indices` handle BEAST numeric taxon IDs.  If your
newick strings already use real taxon names, pass `[{}]` and `[0] * n`.

### Progress reporting

To display a progress bar (or any kind of UI update) while a pairwise call
runs, instantiate a `rapidtrees.ProgressCounter`, hand it to the function,
and **read it from another Python thread** while the call blocks.

```python
import threading, time, rapidtrees
from tqdm import tqdm

counter = rapidtrees.ProgressCounter()

t = threading.Thread(
    target=rapidtrees.pairwise_rf_from_newick_iter,
    args=(names, iter(newicks), translate_maps, map_indices),
    kwargs={"progress": counter},
)
t.start()

with tqdm(total=1.0, unit="frac", bar_format="{l_bar}{bar}|{n:.2f}/{total}") as bar:
    while t.is_alive():
        bar.n = counter.fraction()
        bar.refresh()
        time.sleep(0.1)
    bar.n = 1.0
    bar.refresh()

t.join()
```

`ProgressCounter` methods (all lock-free atomic reads):

| Method | Returns |
| --- | --- |
| `.value()`    | `int` — pairs completed so far |
| `.total()`    | `int` — `n*(n-1)/2` (0 before any call has started) |
| `.fraction()` | `float` in `[0.0, 1.0]`, clamped |
| `.reset()`    | clears both `value` and `total` to `0` |

Notes:

- `progress=None` (the default) holds the GIL throughout the call and adds
  zero overhead — identical to the pre-`0.6.0` behaviour.
- When a counter is supplied the GIL is released for the duration of the
  rayon loop, so a polling Python thread runs unimpeded. Rust never calls
  back into Python.
- The counter ends at exactly `value == total` (`fraction() == 1.0`) after
  the function returns.

### Return types

| Function | Return type | How to decode |
| --- | --- | --- |
| `pairwise_rf_from_newick_iter` | `bytes` — flat `uint32`, row-major | `np.frombuffer(b, dtype=np.uint32).reshape(n, n)` |
| `pairwise_wrf_from_newick_iter` | `list[float]` — flat, row-major | `np.array(lst, dtype=np.float64).reshape(n, n)` |
| `pairwise_kf_from_newick_iter` | `list[float]` — flat, row-major | `np.array(lst, dtype=np.float64).reshape(n, n)` |
| `pairwise_rf_with_snapshots_from_newick_iter` | 6-tuple — see below | see below |
| `pairwise_rf_with_sparse_snapshots_from_newick_iter` | 6-tuple — see below | see below |
| `pairwise_rf_with_rooted_facts_from_newick_iter` | 6-tuple — see below | see below |
| `pairwise_wrf_with_snapshots_from_newick_iter` | 6-tuple — see below | see below |
| `pairwise_kf_with_snapshots_from_newick_iter` | 6-tuple — see below | see below |

### Errors raised

All functions raise `ValueError` when:

- Fewer than 2 trees are provided
- `len(names) != len(map_indices)`
- A `map_indices` value is out of bounds
- Trees have different leaf sets

The rooted-facts endpoint additionally raises `ValueError` when a tree is not
strictly binary, a non-root edge lacks an explicit finite branch length, or a
cumulative root distance or calculated height is non-finite.

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
Plain Newick files need no parsing at all: read the lines, pass them straight to
the iterator API, and name the trees however you like (the CLI names them after
their line number).

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
`i`, otherwise `0`. Column order is deterministic and stable across calls on
the same tree set, so the same trees always give the same column indices.

### RF + sparse snapshot in one pass

`pairwise_rf_with_sparse_snapshots_from_newick_iter` returns the same RF bytes,
leaf names, clade count, and clade bitmasks as the dense endpoint, replacing
`presence_bytes` with a versioned CSR dictionary:

```text
(tree_names, rf_bytes, leaf_names, n_bip, bipartition_clade_bytes, sparse)
```

```python
(
    tree_names, rf_bytes, leaf_names, n_bip,
    bipartition_clade_bytes, sparse,
) = rtd.pairwise_rf_with_sparse_snapshots_from_newick_iter(
    names,
    iter(newicks),
    translate_maps,
    map_indices,
    rooted=True,
)

assert sparse["format_version"] == 1
assert sparse["encoding"] == "csr"
offsets = np.frombuffer(sparse["row_offsets"], dtype=np.uint64)
columns = np.frombuffer(sparse["column_indices"], dtype=np.uint32)

# Sorted, unique clade columns present in tree i.
i = 0
row = columns[offsets[i]:offsets[i + 1]]

# Count source-tree occurrences of every clade without a dense matrix.
counts = np.bincount(columns, minlength=n_bip)
```

The payload fields are:

| Key | Type | Description |
| --- | --- | --- |
| `format_version` | `int` | Currently `1` |
| `encoding` | `str` | Currently `"csr"` |
| `n_entries` | `int` | Number of stored tree–clade incidences (`nnz`) |
| `row_offsets` | `bytes` | Native-endian `uint64`, shape `(n_trees + 1,)` |
| `column_indices` | `bytes` | Native-endian `uint32`, shape `(n_entries,)` |

Each row is sorted and unique. Expanding it to ones at the listed columns
reconstructs the dense endpoint byte-for-byte. CSR supports both rooting modes
and variable-width rows from non-binary trees.

#### Edge table contents

`n_bip` counts **all** edges: pendant (leaf) edges and internal bipartitions.
Pendant edges appear as single-bit rows (one bit set). Since all trees share the
same leaf set, pendant columns are always `1` in every row of the presence
matrix.

> **Do not slice columns by position.** Pendants are *not* grouped before
> internal bipartitions — columns are ordered by leaf-set bit pattern, which
> interleaves them. For `(((A,B),(C,D)),(E,(F,G)));` the column sizes run
> `1 1 1 1 2 1 2 3 1 1 2 3 5`: the internal split `{C,D}` sorts ahead of the
> pendant `{E}`. Select columns by testing the row instead:
>
> ```python
> pendants = bip_bool.sum(axis=1) == 1
> internal = ~pendants
> ```

#### Canonicalisation note

For internal bipartitions (rows with ≥ 2 bits set) the canonical side is the
half that does **not** contain the first leaf alphabetically, so bit 0 is never
set in those rows. Pendant edges are stored verbatim (no flip), so the pendant
of the first leaf has bit 0 set. The complement of an internal bipartition can
be derived as `1 - bip_bool[j]`, which is already masked to `n_leaves` bits.

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

### RF + rooted facts in one pass

`pairwise_rf_with_rooted_facts_from_newick_iter` exports the compact rooted
facts needed by MrHIPSTR-style consumers without reparsing Newick in Python.
It returns rooted RF distances, the same deterministic clade catalog used by
the rooted dense-snapshot endpoint, sparse source-tree clade columns, aligned
node heights, and directly observed binary splits:

```text
(tree_names, rf_bytes, leaf_names, n_clades, clade_bytes, facts)
```

```python
(
    tree_names, rf_bytes, leaf_names, n_clades,
    clade_bytes, facts,
) = rtd.pairwise_rf_with_rooted_facts_from_newick_iter(
    names,
    iter(newicks),
    translate_maps,
    map_indices,
)

assert facts["format_version"] == 2
n_trees = len(tree_names)

clade_columns = np.frombuffer(
    facts["clade_columns"], dtype=np.uint32
).reshape(n_trees, facts["nodes_per_tree"])
node_heights = np.frombuffer(
    facts["node_heights"], dtype=np.float64
).reshape(clade_columns.shape)
root_heights = np.frombuffer(facts["root_heights"], dtype=np.float64)
split_ids = np.frombuffer(
    facts["split_ids"], dtype=np.uint32
).reshape(n_trees, facts["splits_per_tree"])
split_table = np.frombuffer(
    facts["split_table"], dtype=np.uint32
).reshape(facts["n_observed_splits"], 3)
observed_splits = split_table[split_ids]

# Source-tree occurrence count for every exported clade.
clade_counts = np.bincount(clade_columns.ravel(), minlength=n_clades)
```

The facts dictionary contains native-endian buffers:

| Key | Type | Description |
| --- | --- | --- |
| `format_version` | `int` | Currently `2` |
| `root_column` | `int` | Root sentinel; always equal to `n_clades` |
| `nodes_per_tree` | `int` | Non-root nodes per binary tree: `2 * n_leaves - 2` |
| `splits_per_tree` | `int` | Directly observed internal splits per tree: `n_leaves - 1` |
| `n_observed_splits` | `int` | Number of distinct rows in `split_table` |
| `clade_columns` | `bytes` | `uint32`, shape `(n_trees, nodes_per_tree)` |
| `node_heights` | `bytes` | `float64`, aligned with `clade_columns` |
| `root_heights` | `bytes` | `float64`, shape `(n_trees,)` |
| `split_ids` | `bytes` | `uint32`, shape `(n_trees, splits_per_tree)` |
| `split_table` | `bytes` | `uint32`, shape `(n_observed_splits, 3)` |

Each `clade_columns` row is sorted and unique and lists exactly the non-root
clades present in that source tree. Writing ones at those columns reconstructs
the corresponding rooted dense-presence row. `node_heights` is aligned
element-for-element with it, including singleton tips. `root_heights` stores
the implicit all-taxa root separately.

Rows of `split_table` are `(parent, left_child, right_child)` column triples.
Non-root values index `clade_bytes`; `root_column == n_clades` is used only as
a parent. Children are ordered by clade column, and `split_ids` maps each source
tree to only the splits directly observed in that tree.

Heights support ultrametric and non-ultrametric input and are calculated as:

```text
root_distance(root) = 0
root_distance(child) = root_distance(parent) + branch_length(child)
root_height = max(root_distance(tip))
node_height(node) = root_height - root_distance(node)
```

The endpoint requires strictly binary rooted trees and an explicit finite
branch length on every non-root edge. It has no `rooted` argument.

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

Column order matches `bipartition_clade_bytes`, and is deterministic and
stable across calls on the same tree set.

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
