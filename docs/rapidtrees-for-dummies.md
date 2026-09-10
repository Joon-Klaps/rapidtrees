# 🌲 rapidtrees for Dummies

> No PhD in computer science required. This document explains, from the ground up, how `rapidtrees` works internally — from raw tree text to a distance number — using a concrete 7-taxon example throughout.

---

## Table of Contents

1. [What problem are we solving?](#1-what-problem-are-we-solving)
2. [Step 0 — The input: a Newick string](#step-0--the-input-a-newick-string)
3. [Step 1 — Splits: what a tree really *is*](#step-1--splits-what-a-tree-really-is)
4. [Step 2 — Numbering the leaves, once per run](#step-2--numbering-the-leaves-once-per-run)
5. [Step 3 — Fingerprints: a split as a 128-bit number](#step-3--fingerprints-a-split-as-a-128-bit-number)
6. [Step 4 — The DFS: folding one tree in a single pass](#step-4--the-dfs-folding-one-tree-in-a-single-pass)
7. [Step 5 — Canonicalization, for free](#step-5--canonicalization-for-free)
8. [Step 6 — Snapshot: one tree, ready to hand off](#step-6--snapshot-one-tree-ready-to-hand-off)
9. [Step 7 — Interning: from fingerprints to `u32` IDs](#step-7--interning-from-fingerprints-to-u32-ids)
10. [Step 8 — Computing RF on integers](#step-8--computing-rf-on-integers)
11. [Step 9 — Bulk pairwise: the dense and sparse backends](#step-9--bulk-pairwise-the-dense-and-sparse-backends)
12. [Step 10 — Naming the split: the clade table](#step-10--naming-the-split-the-clade-table)
13. [The full pipeline at a glance](#the-full-pipeline-at-a-glance)
14. [Cheat sheet: all representations](#cheat-sheet-all-representations)
15. [Where each piece lives in the source](#where-each-piece-lives-in-the-source)

---

## 1. What problem are we solving?

A phylogenetic tree is a hypothesis about how a set of species (called **taxa** or **leaves**) are related. When you run a Bayesian inference like BEAST, you get *thousands* of sample trees. You want to know how similar or different they are from each other.

The most common measure is the **Robinson–Foulds (RF) distance**. It counts how many "branch groupings" differ between two trees. Two identical trees have RF = 0. Two completely different trees have the maximum RF.

`rapidtrees` computes this for *all pairs* of trees in your dataset — potentially millions of comparisons — as fast as possible. The internal representations described in this document are what make that fast.

The headline constraint: with `T` trees you have `T(T−1)/2` pairs. At `T = 4000` that is eight million comparisons. Anything you do *per comparison* gets multiplied eight million times, so the whole design is about making a single comparison as close to free as possible — and doing all the real work once, up front.

---

## Step 0 — The input: a Newick string

Trees are written in **Newick format**: nested parentheses where leaves are species names and numbers after `:` are branch lengths.

We'll use two 7-taxon example trees throughout this whole document:

```text
Tree 1:  (((A,B),(C,D)),(E,(F,G)));
Tree 2:  (((A,B),(C,(D,E))),(F,G));
```

Visually:

```text
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

> **Spoiler, so you can check the arithmetic as you read:** the answer is **RF = 4**, and between them these two trees contain **13 distinct splits**. Every worked number below is real output from the code, not hand-arithmetic.

---

## Step 1 — Splits: what a tree really *is*

Cut any internal branch of an unrooted tree and it falls into two pieces. The set of leaves on each side is a **split** (also called a bipartition), written `{A,B} | {C,D,E,F,G}`.

Here is the crucial simplification: **a tree is completely described by its set of splits.** Not by its node IDs, not by its drawing, not by the order the parentheses were written. Two trees are identical exactly when they have the same set of splits, and the RF distance is simply how many splits they *don't* share:

```text
RF(t1, t2) = |splits(t1) △ splits(t2)|        ( △ = symmetric difference )
           = (splits only in t1) + (splits only in t2)
```

For our two trees:

| Split | In Tree 1? | In Tree 2? |
| --- | :---: | :---: |
| `{A,B}` | ✅ | ✅ |
| `{F,G}` | ✅ | ✅ |
| `{C,D}` | ✅ | ❌ |
| `{E,F,G}` | ✅ | ❌ |
| `{D,E}` | ❌ | ✅ |
| `{C,D,E}` | ❌ | ✅ |

Two splits in Tree 1 only, two in Tree 2 only → **RF = 4**. ✅

Two kinds of split get special treatment:

- **Pendant edges** (`{A}`, `{B}`, …) — the branch leading to a single leaf. Every tree over the same taxa has all of them, so they never contribute to RF. They are kept anyway, because the *weighted* metrics need their branch lengths.
- **Trivial splits** — a split whose one side is everything-but-one-leaf is just a pendant edge seen from the other direction, so it is dropped to avoid double-counting.

So the entire job reduces to: **turn each tree into its set of splits, then compare sets.** Everything that follows is about making both halves cheap.

---

## Step 2 — Numbering the leaves, once per run

Node IDs are assigned by the parser and differ between files — tree 1's "node 4" has nothing to do with tree 2's "node 4". Taxon *names* are the only stable identity.

So `rapidtrees` sorts the taxon names alphabetically once and assigns each a **bit index**:

```text
A→0   B→1   C→2   D→3   E→4   F→5   G→6
```

The important word is **once**. Every tree in a collection carries the same taxa (this is checked), so this table is a property of the *run*, not of any tree. Deriving it per tree would mean re-cloning and re-sorting the same seven names for every tree in the file — pure waste at 10 000 trees.

Two tables are built here and shared by every tree:

| Table | What it maps | Why it must be shared |
| --- | --- | --- |
| `leaf_index` | taxon name → bit index | Otherwise "bit 3" means different taxa in different trees |
| `labels` | bit index → random 128-bit label | Otherwise fingerprints (next step) aren't comparable |

---

## Step 3 — Fingerprints: a split as a 128-bit number

Here is the central trick of the whole library.

The obvious way to represent `{C,D}` is a **bitset**: one bit per leaf, `0b0001100`. That works, and it is what `rapidtrees` used to do. The problem is size. A bitset for `n` taxa is `⌈n/64⌉` machine words, so:

- combining two child sets costs `⌈n/64⌉` word-ORs, not one operation
- a tree with `n` leaves has ~`2n` nodes, so building one tree costs `Θ(n²/64)` work
- at 10 000 taxa that is 157 words per node, ~20 000 nodes, per tree

Instead, each taxon draws **one random 128-bit label** at the start of the run, and a group's **fingerprint** is the XOR of its members' labels:

```text
fingerprint({C,D}) = label[C] ⊕ label[D]
```

A fingerprint is a single `u128`. Combining two children is *one XOR*, regardless of taxon count. Building a tree drops from `Θ(n²/64)` to `Θ(n)`.

> 💡 **Why XOR and not a hash or an OR?** Because XOR has two properties nothing else does.
>
> 1. **It's associative and commutative**, so a subtree's fingerprint doesn't depend on the order you fold its children — exactly the property a *set* needs.
> 2. **It's its own inverse**, which makes complements free. If `total` is the XOR of all seven labels, then for any group `A`, `fingerprint(A′) = total ⊕ fingerprint(A)`. That is what Step 5 is built on, and an OR or a hash gives you nothing like it.

### The price: `rapidtrees` is no longer exact

Two *different* splits can in principle land on the same 128-bit value, and the code would then merge them. The chance is about `e² / 2¹²⁹` for `e` distinct splits in a run:

| Distinct splits `e` | Collision probability |
| --- | --- |
| 10 000 | `1.5 × 10⁻³¹` |
| 1 000 000 | `1.5 × 10⁻²⁷` |
| 100 000 000 | `1.5 × 10⁻²³` |

For scale, that last figure is roughly **nineteen orders of magnitude below** the rate at which the machine's own RAM flips a bit without telling you. Every run prints its own `e` and bound, so the guarantee is stated rather than assumed:

```text
Distinct splits e = 1056; collision bound e²/2¹²⁹ = 1.64e-33
```

The label table is seeded from a fixed constant, so this is not a source of run-to-run variation: the same input gives the same split IDs on every run and every machine.

---

## Step 4 — The DFS: folding one tree in a single pass

Each node accumulates three things — the struct is called `Acc`:

| Field | Meaning |
| --- | --- |
| `fp` | the subtree's 128-bit fingerprint |
| `first` | where this subtree's leaves start in `leaf_order` |
| `size` | how many leaves are in the subtree |

`size` is worth pausing on. It's the leaf count, carried along for free — which means the "is this a pendant?" and "is this trivial?" tests in Step 6 are integer comparisons rather than a full-width popcount over a bitset. Keeping it is what stops the `⌈n/64⌉` cost sneaking back in through the filter.

The traversal is **two iterative passes over one flat array** indexed by node ID (iterative, not recursive — a 3 000-leaf caterpillar tree would blow the stack otherwise):

**Pass 1 (pre-order)** — each leaf takes the next free slot in `leaf_order`:

```text
leaf_order = [0, 1, 2, 3, 4, 5, 6]      // for Tree 1: A B C D E F G
              ↑     ↑
              A     C
```

Because leaves are numbered *in traversal order*, every subtree occupies a **contiguous run** of `leaf_order`. That is why `first` + `size` is enough to name a leaf set later, with the tree itself long gone.

**Pass 2 (post-order)** — an internal node is just the XOR of its children:

```text
node({C,D}).fp    = acc[C].fp ⊕ acc[D].fp
node({C,D}).first = min(acc[C].first, acc[D].first)
node({C,D}).size  = acc[C].size + acc[D].size
```

One XOR, one `min`, one `+` per node. No allocation, no bitset, no popcount.

Each node then becomes a `Part` — one edge of the tree:

```rust
struct Part {
    key: Fingerprint,  // the canonical fingerprint (see Step 5)
    first: u32,        // where its leaves start in leaf_order
    size: u32,         // how many leaves
    length: f64,       // branch length of the edge above it
}
```

---

## Step 5 — Canonicalization, for free

A split has two sides, and different trees may hand you either one. Tree 1 gives you `{A,B}` from one node and `{C,D,E,F,G}` from another — **the same split**, seen from opposite ends. If they don't compare equal, RF is wrong.

The old bitset approach had to pick a side by convention (say, "the side without leaf 0") and physically flip the bitset when it guessed wrong — an `O(n)` operation, per node, per tree.

With XOR fingerprints it collapses to one line:

```rust
key = min(fp, fp ^ total)
```

Because `fingerprint(A′) = total ⊕ fingerprint(A)`, the two sides of a split give the two values `fp` and `fp ^ total` — in some order. Taking the `min` of the pair picks the same one no matter which side you started from. Both sides now produce an identical `key`, and **no leaf set was ever touched**.

Pendant edges are the one exception: they keep their raw fingerprint. Their complement is filtered out anyway, so nothing can collide with them — and at exactly two taxa, `min(h, h ^ total)` would collapse the two pendants onto each other.

---

## Step 6 — Snapshot: one tree, ready to hand off

A `Snapshot` is one tree's edges after filtering, canonicalizing and sorting:

```rust
struct Snapshot {
    parts: Vec<Part>,       // sorted by key, duplicates merged
    leaf_order: Vec<u32>,   // traversal-order leaf indices
    words: usize,           // ⌈n_leaves / 64⌉
}
```

Three things happen on the way in:

1. **Drop the trivial splits.** `size == 1` is a pendant (kept). `size >= n − 1` is a pendant's complement (dropped). Both decided by integer comparison on `size`.
2. **Sort by `key`.** Note what is being sorted: a 16-byte integer. The bitset version sorted `⌈n/64⌉`-word arrays — at 2 000 taxa, a 256-byte comparison instead of a 16-byte one.
3. **Merge duplicates.** In a rooted binary tree both children of the root canonicalize to the *same* split; without merging, RF would come out inflated by 2 versus phangorn's `RF.dist(rooted=FALSE)`. When two parts merge, their branch lengths are summed.

`Snapshot` is deliberately **short-lived**: built per tree, handed straight to the interner, dropped. Trees are processed in chunks sized to a memory budget, so peak memory stays near the *deduplicated* footprint instead of holding every tree's raw data at once.

---

## Step 7 — Interning: from fingerprints to `u32` IDs

Across 4 000 trees the *same* splits recur constantly — that's what it means for trees to be similar. So each distinct split is assigned a `u32` **ID**, once, in first-seen order.

```text
split {A,B}    → ID 0
split {C,D}    → ID 1
split {A,B,C,D}→ ID 2
...
```

Each tree then becomes an `InternSnap` — two parallel arrays, nothing more:

```rust
struct InternSnap {
    split_ids: Vec<u32>,   // sorted ascending
    lengths: Vec<f64>,     // lengths[i] belongs to split_ids[i]
}
```

The interner matches a candidate on **both** its fingerprint and the cardinality of its smaller side. The cardinality is equal for both sides of a bipartition, so it costs one integer compare and rules out a slice of the collision space for free.

For our two trees the whole run has **13 distinct splits** (7 pendants + 6 internal), so:

```text
tree 1: split_ids = [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12]
tree 2: split_ids = [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 12]
                              ↑  ↑  ↑           ↑
                       differ here — 4 positions total → RF = 4
```

Two payoffs:

- **Memory.** A tree is now `~n_bip × 4` bytes instead of `~n_bip × ⌈n/64⌉ × 8`. The working set of a large analysis fits in L2 cache rather than needing DRAM — and cache misses, not instructions, are what dominate eight million comparisons.
- **Speed.** Comparing two splits is a single integer compare.

---

## Step 8 — Computing RF on integers

Both `split_ids` lists are sorted, so counting the symmetric difference is a two-pointer merge — one linear pass, no hashing, no allocation:

```text
t1: [0, 1, 2, 3, 4, 5, 8,  9, 10, 11, 12]
t2: [0, 1, 2, 3, 5, 6, 7,  8,  9, 10, 12]
     ✓  ✓  ✓  ✓  ✗  ✓  ✗   ✗   ✓   ✓  ✗   ✓
                 └── 4 in t1 only, 4 in t2 only... 
```

More precisely: walk both lists together, count how many IDs are `shared`, then

```text
RF = (len(t1) − shared) + (len(t2) − shared)
```

With `shared = 9`, `len = 11` each: `RF = 2 + 2 = 4`. ✅

The weighted metrics use the same merge, but instead of counting they accumulate over the `lengths` arrays:

| Metric | What it accumulates |
| --- | --- |
| **RF** | count of unshared splits |
| **Weighted RF** | `Σ \|length₁ − length₂\|` over the union of splits |
| **Kuhner–Felsenstein** | `√( Σ (length₁ − length₂)² )` over the union |

A split absent from one tree contributes as if its length there were `0.0`.

---

## Step 9 — Bulk pairwise: the dense and sparse backends

One comparison is now cheap. Eight million of them still need care, and the right strategy depends on how *diverse* the tree set is.

**Sparse backend** — run the two-pointer merge per pair. Cost scales with how many splits each tree has. Best when trees are diverse, so the split table is large and each tree touches a small slice of it.

**Dense backend** — build a `(n_trees × n_splits)` presence matrix once, then compare rows with bit-parallel operations. A row is a bit-packed word array, and

```text
RF(i, j) = popcount(row_i XOR row_j)
```

which handles 64 splits per instruction. Best when trees are similar, so the split table is narrow and the matrix stays small.

`--backend auto` (the default) picks between them from the observed number of distinct splits, and refuses a dense matrix that would exceed a 5 GB budget. Both backends return **identical** matrices; the choice is purely about speed. The run log names the one that ran:

```text
Backend: dense (--backend auto)
```

Two further tricks live in the dense path:

- **Column ordering.** Splits are sorted by how often they occur, and each tree records the first and last column it touches — so a pair's comparison can skip whole stretches of the matrix.
- **Dropping dead columns.** A split present in *every* tree (or in none) can never contribute to any RF distance, so its column is removed before the loop starts.

---

## Step 10 — Naming the split: the clade table

If fingerprints replaced bitsets everywhere, what says *which taxa* a split actually names?

Not the fingerprint. It's a one-way summary — you cannot read `{C,D}` back out of a `u128`. And you need that answer: `presence.mean(axis=0)` gives you "split 4711 is in 87% of trees", which is useless until you can say that split 4711 *is* `{C,D}`. Clade-frequency and consensus-tree diagnostics live on exactly this.

So the run keeps one **clade table**: the canonical leaf set of every **distinct split** — not per node, not per tree — recovered from that split's `first`/`size` run in `leaf_order`, with the tree already dropped.

```text
Old bitset-DFS   built a leaf set once per node, per tree   (~2n × T)
Today            once per distinct split                     (e)
```

### Why leaf indices, not a bitset

The obvious storage is one bit per leaf. It was the original choice, and it is the wrong one at scale, because a bit-packed set costs `⌈n/64⌉` words **whether it names three taxa or three thousand** — so the table grows as `e × n`.

And most splits are tiny. Measured on random trees, the median clade holds **3 taxa at every taxon count**:

| taxa | splits | median clade | bit-packed | leaf indices | |
| ---: | ---: | ---: | ---: | ---: | --- |
| 250 | 3 197 | 3 | 100 K | 111 K | bit-packing marginally ahead |
| 500 | 6 450 | 3 | 403 K | 251 K | **1.6× smaller** |
| 1 000 | 12 956 | 3 | 1 620 K | 563 K | **2.9× smaller** |
| 4 000 | 51 940 | 3 | 25 564 K | 2 805 K | **9.1× smaller** |

Storing the indices costs `size × 4` bytes and stops tracking `n` altogether. They live in one flat `Vec<u32>` with an offsets array rather than a `Vec<Vec<u32>>`, because a per-split `Vec` header is 24 bytes — more than the leaf indices it would point at.

Bit packing still happens, but at the very edge: `build_bipartition_bytes` packs each row on its way out, because that is the wire format `np.unpackbits` expects. Nothing upstream of that ever sees a bit.

### The one subtlety: column order

Export columns are sorted so the same tree set always yields the same column indices. That order was defined by comparing packed `u64` words, and it has to stay *exactly* the same or every caller's column indices shift under them.

Word comparison runs low-to-high, and each word is compared by value — so the **lowest differing word** decides, and within it the **highest differing bit**. Those pull in opposite directions, which makes one case counter-intuitive:

```text
{63}  vs  {64}   →   {63} is GREATER
```

Leaf 63 is the top bit of word 0; leaf 64 is the bottom bit of word 1. Word 0 is compared first and `{64}` has nothing there, so it loses. `cmp_packed` reproduces this from two index runs without packing anything, and a test checks it against real packed bits over every subset of a 3-word universe.

### The export for our example

Decoded — note the canonical side always excludes leaf 0 (`A`):

```text
  col  0: {A}          in t1 ✓  t2 ✓      col  7: {C,D,E}    in t1 ✗  t2 ✓
  col  4: {C,D}        in t1 ✓  t2 ✗      col 10: {F,G}      in t1 ✓  t2 ✓
  col  6: {D,E}        in t1 ✗  t2 ✓      col 11: {E,F,G}    in t1 ✓  t2 ✗
                                          col 12: {C,D,E,F,G} in t1 ✓  t2 ✓
```

That last row is the `{A,B}` split — stored as its complement, because the canonical side is the one without leaf `A`. Taking that complement costs no set arithmetic either: leaves are numbered in traversal order, so a subtree owns a contiguous run of `leaf_order` and everything *else* is simply the rest of the array.

---

## The full pipeline at a glance

```text
  "(((A,B),(C,D)),(E,(F,G)));"
            │
            │  strip BEAST annotations, apply TRANSLATE, parse
            ▼
        PhyloTree                                    ── per tree, transient
            │
            │  [build]  two iterative passes:
            │           pre-order  → leaf_order slots
            │           post-order → fp = XOR of children
            │           key = min(fp, fp ^ total)
            ▼
    Snapshot { parts, leaf_order }                   ── per tree, dropped after interning
            │
            │  [intern]  dedupe by (fingerprint, smaller-side size)
            │            assign u32 IDs in first-seen order
            │            record ONE leaf set per new split
            ▼
    InternSnap { split_ids, lengths }                ── per tree, KEPT
            │
            ├──[distances]── two-pointer merge on u32  →  RF / WRF / KF matrix
            │                 (dense or sparse kernel)
            │
            └──[export]───── presence matrix, branch-length matrix,
                             bipartition bytes            →  Python / NumPy
```

The shape to remember: **all the expensive work happens once per tree, and everything in the `T(T−1)/2` loop is integer arithmetic.**

---

## Cheat sheet: all representations

### Persistent (live in `Snapshots` for the whole run)

| Structure | Type | Size | Purpose |
| --- | --- | --- | --- |
| `snapshots` | `Vec<InternSnap>` | `T × n_bip × 4 B` (+ `× 8 B` with lengths) | one tree = sorted `u32` IDs |
| `clades` | `CladeTable` | `Σ size × 4 B` | which taxa each split names — **export only** |
| `leaf_names` | `Vec<String>` | `n` strings | alphabetical taxon names |

### Transient (built per tree, then dropped)

| Structure | Lifetime | Purpose |
| --- | --- | --- |
| `PhyloTree` | one tree | parsed Newick |
| `Acc` | one node | `(fingerprint, first, size)` during the DFS |
| `Part` | one edge | `(key, first, size, length)` |
| `Snapshot` | one tree | all `Part`s, sorted and deduped |

### The two run-wide tables

| Table | Built | Purpose |
| --- | --- | --- |
| `leaf_index` | once per run | taxon name → bit index |
| `labels` | once per run | bit index → random 128-bit label |

### What `Retain` controls

Not every path needs everything, and both extras cost real work:

| Path | `lengths` | `bipartitions` |
| --- | :---: | :---: |
| `pairwise_rf_from_newick_iter` | ✗ | ✗ |
| `pairwise_wrf` / `pairwise_kf` | ✓ | ✗ |
| any `*_with_snapshots_*` | ✓ | ✓ |
| CLI | ✓ | ✓ |

---

## Where each piece lives in the source

| Module | Holds |
| --- | --- |
| `snapshot/fingerprint.rs` | the two run-wide tables, and the `Fingerprint` type |
| `snapshot/build.rs` | `Acc`, `Part`, `Snapshot` — one tree in, one snapshot out |
| `snapshot/intern.rs` | `Interner`, `InternSnap` — dedupe to `u32` IDs |
| `snapshot/export.rs` | the flat byte buffers Python reads |
| `snapshot/mod.rs` | `Snapshots`, the construction pipeline, `Retain` |
| `distances.rs` | RF / WRF / KF, dense and sparse kernels |
| `snapshot/clades.rs` | the export-only leaf-set table, and its packed ordering |
| `io.rs` | NEXUS/Newick parsing, BEAST annotation stripping |
| `api.rs` | PyO3 bindings — glue only, no computation |
