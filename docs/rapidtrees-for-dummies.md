# 🌲 rapidtrees for Dummies

> No PhD in computer science required. This document explains, from the ground up, how `rapidtrees` works internally — from raw tree text to a distance number — using a concrete 7-taxon example throughout.

---

## Table of Contents

1. [What problem are we solving?](#1-what-problem-are-we-solving)
2. [Step 0 — The input: a Newick string](#step-0--the-input-a-newick-string)
3. [Step 1 — Bitsets: representing a group of taxa as a number](#step-1--bitsets-representing-a-group-of-taxa-as-a-number)
4. [Step 2 — Building bitsets via DFS](#step-2--building-bitsets-via-dfs)
5. [Step 3 — Canonicalization: making the same split look the same](#step-3--canonicalization-making-the-same-split-look-the-same)
6. [Step 4 — TreeSnapshot: the internal intermediate form](#step-4--treesnapshot-the-internal-intermediate-form)
7. [Step 5 — Computing RF distance via two-pointer merge](#step-5--computing-rf-distance-via-two-pointer-merge)
8. [Step 6 — The scaling problem at large taxa counts](#step-6--the-scaling-problem-at-large-taxa-counts)
9. [Step 7 — InternSnap: replacing bitsets with integer IDs](#step-7--internsnap-replacing-bitsets-with-integer-ids)
10. [Step 8 — Single-pair distances via the sorted merge](#step-8--single-pair-distances-via-the-sorted-merge)
11. [Step 9 — Bulk pairwise: the dense backend (RF, WRF, KF)](#step-9--bulk-pairwise-the-dense-backend-rf-wrf-kf)
12. [The full pipeline at a glance](#the-full-pipeline-at-a-glance)
13. [Cheat sheet: all representations](#cheat-sheet-all-representations)

---

## 1. What problem are we solving?

A phylogenetic tree is a hypothesis about how a set of species (called **taxa** or **leaves**) are related. When you run a Bayesian inference like BEAST, you get *thousands* of sample trees. You want to know how similar or different they are from each other.

The most common measure is the **Robinson–Foulds (RF) distance**. It counts how many "branch groupings" differ between two trees. Two identical trees have RF = 0. Two completely different trees have the maximum RF.

`rapidtrees` computes this for *all pairs* of trees in your dataset — potentially millions of comparisons — as fast as possible. The internal representations described in this document are what make that fast.

---

## Step 0 — The input: a Newick string

Trees are written in **Newick format**: nested parentheses where leaves are species names and numbers after `:` are branch lengths.

We'll use two 7-taxon example trees throughout this whole document:

```
Tree 1:  (((A,B),(C,D)),(E,(F,G)));
Tree 2:  (((A,B),(C,(D,E))),(F,G));
```

Visually:

```
Tree 1                              Tree 2

         ┌── A                               ┌── A
      ┌──┤                                ┌──┤
      │  └── B                            │  └── B
   ┌──┤                                ┌──┤
   │  │  ┌── C                         │  │  ┌── C
   │  └──┤                             │  └──┤
   │     └── D                         │     │  ┌── D
───┤                               ────┤     └──┤
   │     ┌── E                         │        └── E
   └─────┤                             │
         │  ┌── F                      │  ┌── F
         └──┤                          └──┤
            └── G                         └── G
```

These trees agree on some groupings and disagree on others. Our goal: count exactly how many groupings differ (the RF distance).

---

## Step 1 — Bitsets: representing a group of taxa as a number

A **bipartition** (also called a "split") is what an internal branch creates: it divides all taxa into two groups (everything on the left vs everything on the right of that branch).

For example, the branch connecting `{F, G}` to the rest of Tree 1 creates the bipartition:
```
{F, G}  |  {A, B, C, D, E}
```

Instead of storing this as a list of names, `rapidtrees` stores it as a compact **bitset** — a single integer where each bit represents one taxon.

Taxa are sorted alphabetically and assigned a fixed bit position **once**, used by all trees:

```
A = bit 0   (value 1)
B = bit 1   (value 2)
C = bit 2   (value 4)
D = bit 3   (value 8)
E = bit 4   (value 16)
F = bit 5   (value 32)
G = bit 6   (value 64)
```

So the group `{F, G}` becomes:

```
bit: 6 5 4 3 2 1 0
     1 1 0 0 0 0 0  =  96  (in decimal)  =  0b1100000
```

And the group `{C, D}` becomes:

```
bit: 6 5 4 3 2 1 0
     0 0 0 1 1 0 0  =  12  (in decimal)  =  0b0001100
```

> **Why bits?** Because computers can combine and compare groups of taxa with single CPU instructions. `{A,B} ∪ {C,D}` = bitwise OR. Two groups are equal iff their integers are equal.

For 7 taxa, one 64-bit integer is enough. For 2 000 taxa, you need 32 integers packed in a row (a `Vec<u64>` of length 32). This is the `Bitset` struct.

---

## Step 2 — Building bitsets via DFS

`rapidtrees` uses a **depth-first search (DFS)** — it walks the tree from the leaves upward, computing each internal node's bitset by OR-ing together its children's bitsets.

**Tree 1: `(((A,B),(C,D)),(E,(F,G)));`**

```
Start at the leaves (each leaf = one bit set):

  A  →  0b0000001  =   1
  B  →  0b0000010  =   2
  C  →  0b0000100  =   4
  D  →  0b0001000  =   8
  E  →  0b0010000  =  16
  F  →  0b0100000  =  32
  G  →  0b1000000  =  64

Work upward (OR children together):

  Node(A, B)           →    1 OR  2  =  0b0000011  =   3    ← "the A+B clade"
  Node(C, D)           →    4 OR  8  =  0b0001100  =  12    ← "the C+D clade"
  Node((A,B), (C,D))   →    3 OR 12  =  0b0001111  =  15    ← "the A+B+C+D clade"
  Node(F, G)           →   32 OR 64  =  0b1100000  =  96    ← "the F+G clade"
  Node(E, (F,G))       →   16 OR 96  =  0b1110000  = 112    ← "the E+F+G clade"
  Root                 →   15 OR 112 =  0b1111111  = 127    ← all taxa (discarded)
```

The **root node** always equals "all taxa" — it doesn't represent a bipartition and is discarded.
**Single-leaf nodes** are also discarded (trivial splits — every tree has them).

Raw internal node bitsets from Tree 1 (before canonicalization): **`{3, 12, 15, 96, 112}`**

---

## Step 3 — Canonicalization: making the same split look the same

Here's the problem. The bipartition `{A,B} | {C,D,E,F,G}` can be represented two ways:
- As the `{A,B}` side → bitset `3`
- As the `{C,D,E,F,G}` side → bitset `124`

Depending on where in the tree a node sits, the DFS might compute either one. If two trees both have this split but one computes `3` and the other computes `124`, they would look *different* in a comparison even though they're identical.

**The fix:** always store the side that does **not** contain taxon A (bit 0).

- If bit 0 is set in the raw bitset → flip to complement (the other side)
- If bit 0 is not set → keep as-is

```
Bitset   3  =  0b0000011  → bit 0 IS set → flip  →  127 XOR   3 = 0b1111100 = 124
Bitset  12  =  0b0001100  → bit 0 not set → keep →  12
Bitset  15  =  0b0001111  → bit 0 IS set → flip  →  127 XOR  15 = 0b1110000 = 112  ← same as Bitset 112!
Bitset  96  =  0b1100000  → bit 0 not set → keep →  96
Bitset 112  =  0b1110000  → bit 0 not set → keep →  112
```

> 💡 **Root-split deduplication:** Bitsets 15 and 112 are the two sides of the same root split — `{A,B,C,D} | {E,F,G}`. They canonicalize to the same value (both → 112). One copy is dropped. This is expected: a rooted binary tree has the same number of unique bipartitions (N − 2 = 5) as its unrooted equivalent, because the root placement does not add a new bipartition.

Canonical bipartitions for Tree 1: **`{12, 96, 112, 124}`**

What do these mean in plain English?

| Bitset | Binary       | Taxa set    | Meaning              |
|--------|--------------|-------------|----------------------|
| 12     | `0b0001100`  | {C, D}      | C and D share a clade |
| 96     | `0b1100000`  | {F, G}      | F and G share a clade |
| 112    | `0b1110000`  | {E, F, G}   | E, F, G share a clade |
| 124    | `0b1111100`  | {C,D,E,F,G} | A and B are sisters   |

**Tree 2** (`(((A,B),(C,(D,E))),(F,G));`) — same process:

```
Node(A, B)              →   3  → bit 0 set → flip  →  127 XOR  3 = 124
Node(D, E)              →  24  → bit 0 not set     →  24
Node(C, (D,E))          →  28  → bit 0 not set     →  28
Node((A,B), (C,(D,E)))  →  31  → bit 0 set → flip  →  127 XOR 31 =  96  ← same as Node(F,G)! Root split duplicate.
Node(F, G)              →  96  → bit 0 not set     →  96
```

Canonical bipartitions for Tree 2: **`{24, 28, 96, 124}`**

| Bitset | Taxa set     | Meaning              |
|--------|--------------|----------------------|
| 24     | {D, E}       | D and E share a clade |
| 28     | {C, D, E}    | C, D, E share a clade |
| 96     | {F, G}       | F and G share a clade ← **same as Tree 1** |
| 124    | {C,D,E,F,G}  | A and B are sisters   ← **same as Tree 1** |

---

## Step 4 — TreeSnapshot: the internal intermediate form

A `TreeSnapshot` is an **internal** intermediate form used during tree loading. It is not part of the public API — user code always works with `Snapshots`. A snapshot stores the canonicalized bipartitions for one tree as parallel sorted arrays:

```
TreeSnapshot (Tree 1)                     ← internal only, not returned to callers
├── parts:   [ Bitset(12), Bitset(96), Bitset(112), Bitset(124) ]   ← sorted!
├── lengths: [    0.20,       0.30,       0.40,        0.10     ]   ← parallel to parts
├── words:   1     (1 × u64 per Bitset, because 7 taxa ≤ 64)
└── rooted:  false

TreeSnapshot (Tree 2)
├── parts:   [ Bitset(24), Bitset(28), Bitset(96), Bitset(124) ]
├── lengths: [    0.15,       0.25,       0.30,       0.10     ]
├── words:   1
└── rooted:  false
```

Key points:
- `parts` is **sorted** — this enables the fast O(m+n) merge comparison below.
- `lengths[i]` is the branch length of `parts[i]` — **parallel arrays**, not a HashMap. 10–20× faster to access than hash lookups.
- `TreeSnapshot`s are immediately fed into `Snapshots::intern()` and then discarded.

---

## Step 5 — Computing RF distance via two-pointer merge

Both `parts` vectors are sorted. RF is computed using a **two-pointer merge** — exactly like merging two sorted lists and counting how many elements appear in both (the intersection).

```
Tree 1 parts: [  12,   96,  112,  124 ]
               ↑ i
Tree 2 parts: [  24,   28,   96,  124 ]
               ↑ j

Compare parts[i] vs parts[j]:

  12 vs  24  →  12 < 24  → unique to Tree 1  → advance i
  96 vs  24  →  96 > 24  → unique to Tree 2  → advance j
  96 vs  28  →  96 > 28  → unique to Tree 2  → advance j
  96 vs  96  →  EQUAL ✓  → shared split!  intersection = 1,  advance both
 112 vs 124  → 112 < 124 → unique to Tree 1  → advance i
 124 vs 124  →  EQUAL ✓  → shared split!  intersection = 2,  advance both

RF = |Tree 1 splits| + |Tree 2 splits| - 2 × intersection
   =      4         +        4        -    2 × 2
   =  4
```

The **shared splits** are `{F,G}` and `{C,D,E,F,G}` (= "A and B are sisters").
The **different splits** are:
- Tree 1 has: `{C,D}` and `{E,F,G}`
- Tree 2 has: `{D,E}` and `{C,D,E}`

Each comparison above (`12 vs 24`, etc.) is comparing **two `Vec<u64>`** values. With 7 taxa it's comparing 1 integer. With 2 000 taxa it's comparing 32 integers — and those integers may not be in CPU cache.

---

## Step 6 — The scaling problem at large taxa counts

Here's where things get slow.

With **2 000 taxa**:
- Each `Bitset` = 32 × 8 bytes = **256 bytes**
- A tree with 600 internal branches → 600 × 256 = **~150 KB per TreeSnapshot**
- 5 000 trees → **~750 MB** of snapshot data

750 MB does not fit in any CPU cache (L1 ≈ 32 KB, L2 ≈ 256 KB, L3 ≈ 32 MB). Every single `Bitset` comparison in the RF inner loop has to wait for data to arrive from main RAM — which takes 100–300× longer than a cache hit.

The solution is **interning**.

---

## Step 7 — InternSnap: replacing bitsets with integer IDs

After all `TreeSnapshot`s are built in a first pass, a second pass walks every bipartition in every tree and assigns each **unique** bipartition a sequential `u32` integer ID. This is called **interning** — like a dictionary lookup where every word gets a page number. The `TreeSnapshot`s are then discarded; all subsequent distance calculations work exclusively on `InternSnap`s held inside a `Snapshots` container.

**Building the global dictionary:**

```
Process Tree 1 parts in order:
  Bitset(12)  → new! → assign ID 0
  Bitset(96)  → new! → assign ID 1
  Bitset(112) → new! → assign ID 2
  Bitset(124) → new! → assign ID 3

Process Tree 2 parts:
  Bitset(24)  → new! → assign ID 4
  Bitset(28)  → new! → assign ID 5
  Bitset(96)  → seen before! → reuse ID 1  ←  same split, same ID
  Bitset(124) → seen before! → reuse ID 3  ←  same split, same ID
```

The global bipartitions table (lives once in `InternedSnapshots`):

```
ID 0  →  Bitset(12)   =  {C, D}
ID 1  →  Bitset(96)   =  {F, G}
ID 2  →  Bitset(112)  =  {E, F, G}
ID 3  →  Bitset(124)  =  {C, D, E, F, G}
ID 4  →  Bitset(24)   =  {D, E}
ID 5  →  Bitset(28)   =  {C, D, E}
```

Each tree is now represented as a **sorted `Vec<u32>`** — just small integers:

```
InternSnap (Tree 1)
├── split_ids: [ 0,    1,    2,    3  ]   ← sorted by ID
└── lengths:   [ 0.20, 0.30, 0.40, 0.10 ] ← still parallel to split_ids

InternSnap (Tree 2)
├── split_ids: [ 1,    3,    4,    5  ]   ← IDs 1 and 3 shared with Tree 1!
└── lengths:   [ 0.30, 0.10, 0.15, 0.25 ]
```

The `Snapshots` container holds:
```
Snapshots
├── snapshots:          Vec<InternSnap>        ← one per tree
├── bipartitions:       Vec<Bitset>            ← the global ID→Bitset table
└── leaf_names:         Vec<String>            ← ["A","B","C","D","E","F","G"]
```

**Memory comparison for 5 000 trees with 2 000 taxa:**

| Representation | Per bipartition | 5 000 trees × 600 splits | Fits in L3 cache? |
|---|---|---|---|
| `TreeSnapshot` (Bitset) | 256 bytes | ~750 MB | ❌ No |
| `InternSnap` (u32) | 4 bytes | ~12 MB | ✅ Yes |

The Bitsets still exist in the global table — but there's only **one copy of each unique bipartition** instead of one per tree. The per-tree data is now tiny integers.

---

## Step 8 — Single-pair distances via the sorted merge

> This section describes the single-pair functions `rf_distance`, `wrf_distance`, and `kf_distance`. For computing all pairwise distances at once (`pairwise_rf` / `pairwise_wrf` / `pairwise_kf`), see Step 9.

Exactly the same two-pointer merge as Step 5, but now each comparison is **one integer vs one integer**:

```
Tree 1 IDs: [ 0,  1,  2,  3 ]
             ↑ i
Tree 2 IDs: [ 1,  3,  4,  5 ]
             ↑ j

  0 vs 1  →  0 < 1  → advance i
  1 vs 1  →  EQUAL ✓   intersection = 1
  2 vs 3  →  2 < 3  → advance i
  3 vs 3  →  EQUAL ✓   intersection = 2

RF = 4 + 4 - 2×2 = 4   ← identical answer
```

The comparison `1 vs 1` is a **single CPU instruction** instead of comparing 32 × 64-bit words. With 5 000 trees × ~600 splits × millions of pairs, this compounds into a huge speedup because the entire working set now fits in L3 cache.

For WRF and KF, the same merge runs but accumulates branch-length differences (`|lenᵢ − lenⱼ|` for WRF) or squared differences (`(lenᵢ − lenⱼ)²` for KF) instead of counting shared IDs. The branch lengths live in the parallel `lengths` vector of each `InternSnap`.

---

## Step 9 — Bulk pairwise: the dense backend (RF, WRF, KF)

For a single pair the merge in Step 8 is fast. For *all* pairs — 5 000 trees gives ~12.5 million comparisons — a smarter layout pays off. `rapidtrees` reformulates each metric as an operation over contiguous memory rows, all powered by the same parallel upper-triangle sweep.

### RF — bit-packed presence matrix

Build a **presence matrix** with one row per tree. Each row is a bitset of width *U* (the number of unique splits across all trees), packed into `u64` words:

```text
Split IDs:      0   1   2   3   4   5
                │   │   │   │   │   │
Tree 1 row:  [  1   1   1   1   0   0  ]  packed → 0b001111
Tree 2 row:  [  0   1   0   1   1   1  ]  packed → 0b111010

RF(1, 2) = a₁  +  a₂  −  2 · popcount( row₁  AND  row₂ )
         =  4  +   4  −  2 · popcount( 0b001010 )
         =  4  +   4  −  2 · 2
         =  4
```

One extra win: splits present in *every* tree ("universal" splits — always at least the pendant/leaf edges) contribute equally to both `a` values and to the shared count, so they cancel out of RF exactly. They are dropped *before* packing, shrinking the rows. For MCMC-like sets where trees are nearly identical, most splits are universal and the bit-rows become very short — this is why RF is **~15× faster** on posterior samples.

### WRF — dense length rows

Build a **length matrix** with one row per tree. Column *u* of tree *t*'s row holds the branch length of split *u* (0.0 if absent):

```text
Split IDs:       0      1      2      3      4      5
                 │      │      │      │      │      │
Tree 1 row:  [ 0.20   0.30   0.40   0.10   0.00   0.00 ]
Tree 2 row:  [ 0.00   0.30   0.00   0.10   0.15   0.25 ]

WRF(1, 2) = |0.20−0.00| + |0.30−0.30| + |0.40−0.00| + |0.10−0.10| + |0.00−0.15| + |0.00−0.25|
           =   0.20     +    0.00     +    0.40     +    0.00     +    0.15     +    0.25
           = 1.00
```

There is no "shared term" shortcut for an absolute difference, so this is a per-element sweep — but walking two contiguous `f64` arrays is cache-friendly and vectorises well.

### KF — dot-product (Gram) form

Using the same length matrix, KF is the Euclidean distance between two length rows. Expanding the squared Euclidean distance:

```text
KF(i, j) = sqrt( ‖Wᵢ‖² + ‖Wⱼ‖² − 2 · ⟨Wᵢ, Wⱼ⟩ )
```

The **self-norms** `‖Wᵢ‖²` are computed once up front for all trees. Each pair then costs only one dot product `⟨Wᵢ, Wⱼ⟩` — analogous to RF's `popcount(AND)`. A clamp to zero before `sqrt` guards against floating-point rounding making the bracket slightly negative for near-identical trees.

---

## The full pipeline at a glance

```text
BEAST/Newick file
       │
       │  io::load_beast_trees / io::parse_trees
       ▼
 Vec<PhyloTree>  ─ (parsed tree objects from the `phylotree` crate) ─

       │
       │  TreeSnapshot::from_tree()  [for each tree, internal only]
       │    1. Sort taxa alphabetically → assign bit positions
       │    2. DFS bottom-up → compute node bitsets via OR
       │    3. Canonicalize (flip if bit 0 set)
       │    4. Sort parts vector
       ▼
 Vec<TreeSnapshot>  (internal; discarded after this step)

       │
       │  Snapshots::intern()
       │    1. Walk all parts across all trees
       │    2. Assign each unique bipartition a u32 ID
       │    3. Re-sort per-tree split_ids by ID
       ▼
 Snapshots  ← the only public representation
       │
       ├─  rf_distance / wrf_distance / kf_distance    (one pair — sorted merge)
       └─  pairwise_rf / pairwise_wrf / pairwise_kf   (all pairs — dense backend)
```

---

## Cheat sheet: all representations

### Persistent structures (live in `Snapshots`)

| Name | File | What it stores | Size per bipartition | When used |
|---|---|---|---|---|
| `Bitset` | `src/bitset.rs` | Raw packed bits, 1 bit per taxon | `ceil(n_taxa/64) × 8` bytes | Building blocks for everything |
| `TreeSnapshot` | `src/snapshot.rs` | `Vec<Bitset>` (sorted) + parallel `Vec<f64>` lengths | 256 bytes @ 2000 taxa | **Internal only** — built during loading, discarded after interning |
| `InternSnap` | `src/snapshot.rs` | `Vec<u32>` IDs (sorted); parallel `Vec<f64>` lengths (omitted on RF-only paths) | 4 bytes | **Internal** — held inside `Snapshots`, used by all distance functions |
| `Snapshots` | `src/snapshot.rs` | Global ID↔Bitset table + all `InternSnap`s + leaf names | One `Bitset` per *unique* split, shared across trees | **The public API** — returned by all loading functions |

### Transient compute structures (built and freed per pairwise call)

| Name | Metric | What it stores | Size |
|---|---|---|---|
| Packed presence matrix | `pairwise_rf` | One bit-row per tree over non-universal splits, packed into `u64` words | `n_trees × ceil(U_nonuniversal/64) × 8` bytes |
| Dense length matrix | `pairwise_wrf`, `pairwise_kf` | One `f64` row per tree over all unique splits (0.0 if absent) | `n_trees × U × 8` bytes |

> **tl;dr:** `TreeSnapshot` is a private stepping stone used during parsing. `Snapshots` is what everything else touches. The packed/dense matrices are built on the fly for each bulk pairwise call and freed immediately after.
