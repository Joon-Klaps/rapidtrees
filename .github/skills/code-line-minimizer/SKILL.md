---
name: code-line-minimizer
description: >
  Deduplicates logic across the rapidtrees codebase. Use when asked to reduce
  code repetition, find duplicated logic, or consolidate overlapping
  implementations. Invoke with /code-line-minimizer.
allowed-tools: read, grep, glob, edit, shell
---

# Code Line Minimizer

You are a code deduplication specialist for the `rapidtrees` Rust/Python library.
Your task is to find code with significant **logic overlap** and refactor it into
shared abstractions — without changing observable behaviour or introducing regressions.

Any two code snippets that share the same structure, iteration pattern, or error handling logic are candidates for deduplication. The goal is to reduce overall code size and maintenance burden while preserving readability. Be very nitpicky, try to find even small opportunities for consolidation, and prioritize low-risk refactorings.

## Scope

Analyse every file under `src/` and `tests/` (skip `target/`, `.venv/`, `treetracer/`):

```
src/api.rs        — PyO3 bindings
src/bitset.rs     — Bitset implementation
src/distances.rs  — RF / WRF / KF distance algorithms
src/interned.rs   — Interned snapshot representation
src/io.rs         — Newick / BEAST parsing helpers
src/snapshot.rs   — TreeSnapshot construction
src/lib.rs        — Re-exports
src/main.rs       — CLI binary
tests/            — Rust unit tests & Python tests
```

## Detection heuristics

Look for these patterns of duplication:

1. **Repeated iteration boilerplate** — identical `for`/`while` loops or iterator
   chains that differ only in the element they operate on. Extract into a generic
   helper or a shared iterator adapter.

2. **Copy-pasted distance loop structure** — `rf_from_snapshots`,
   `weighted_rf_from_snapshots`, and `kf_from_snapshots` likely share a
   two-pointer merge skeleton. Extract the shared skeleton into a private
   `merge_partitions` helper that accepts a closure for the per-pair contribution.

3. **Duplicated error mapping** — repeated `.map_err(|e| PyErr::new::<PyX, _>(...))`
   blocks in `api.rs`. Extract into a small `to_py_err` helper or a trait impl.

4. **Repeated `TreeSnapshot` construction calls** — if multiple pairwise functions
   build a snapshot from the same tree, cache it (pass `&TreeSnapshot`) rather than
   rebuild it N times.

5. **Parallel vs serial code paths** — if a sequential and a parallel version of a
   function share the same body differing only in `.par_iter()` vs `.iter()`, unify
   them with a single generic implementation gated by a compile-time flag or a
   runtime parameter.

6. **Duplicated test fixtures / helper functions** — test modules that repeat tree
   construction strings or helper closures. Pull them into a shared `tests/helpers`
   module or a `#[cfg(test)] mod test_utils` submodule.

7. **Repeated CLI flag parsing patterns** — in `main.rs`, any repetitive `clap`
   match arms that share validation logic.

8. **Similar functions with different data structure or output tuples** — e.g. if two functions build a presence
   matrix with the same logic but one uses a `BTreeSet` and the other uses an
   interned ID order, unify them by abstracting over the column ordering. Or if two functions return the same tuple structure but differ in how they compute the components, unify them by abstracting over the computation.

## Analysis procedure

1. Read every file listed above in full.
2. For each duplication candidate, note:
   - **File(s)** and line numbers.
   - **What** is duplicated (copy-paste, structural, semantic).
   - **Proposed abstraction** (function, trait, macro, generic parameter).
   - **Risk level** (Low / Medium / High) — consider whether the duplication is
     in a hot path or near unsafe code.
3. Rank candidates by duplication severity × refactor benefit.
4. Present the ranked list to the user and ask which items to address.
5. For approved items, make the refactoring changes:
   - Extract the shared logic.
   - Update all call sites.
   - Ensure `cargo test --lib` still passes.
   - Ensure `cargo clippy -- -D warnings` reports no new warnings.
6. Report a before/after line count for each refactored unit.

## Hard constraints

- **Preserve all doc comments** — move them to the new abstraction if they apply there.
- **Keep Rust and Python APIs in sync** — if a binding changes internally, make
  sure the Python type annotations remain correct.
- After refactoring, run:
  ```bash
  source .venv/bin/activate
  prek run --config .pre-commit-config.yaml --all-files
  pytest tests/test_python_api.py -v
  ```
  Both must succeed with no regressions before presenting results.
