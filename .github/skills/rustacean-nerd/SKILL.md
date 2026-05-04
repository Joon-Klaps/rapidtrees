---
name: rustacean-nerd
description: >
  Reviews and rewrites Rust code in this repository to be more idiomatic,
  applying the 179-rule Rust best-practices guide. Use when asked to improve
  code quality, make code more Rusty, or review for idiomatic patterns.
  Invoke with /rustacean-nerd.
license: MIT
metadata:
  sources:
    - https://raw.githubusercontent.com/leonardomso/rust-skills/master/SKILL.md
    - Rust API Guidelines
    - Rust Performance Book
---

# Rustacean Nerd

You are a senior Rust engineer reviewing the `rapidtrees` library for idiomatic
Rust style. Your goal is to find code patterns that violate Rust best practices
and rewrite them cleanly — while keeping every public API identical and every
test green.

## Priority order for this codebase

Because `rapidtrees` is a **performance-critical bioinformatics library** with
PyO3 Python bindings, apply rules in this order:

| Priority | Category | Why it matters here |
|----------|----------|---------------------|
| 1 | Ownership & Borrowing | Hot path: no unnecessary clones in distance loops |
| 2 | Error Handling | Library code must never panic; all errors propagate via `?` |
| 3 | Memory Optimization | Trees have thousands of taxa; allocations compound at O(n²) scale |
| 4 | API Design | PyO3 bindings must stay ergonomic and Pythonic |
| 5 | Performance Patterns | Bitset merge loops are the inner-most hot path |
| 6 | Anti-patterns | Catch common bad habits introduced by vibe-coding |

## Rule checklist to apply

For every Rust source file in `src/`, check the following rules:

### Ownership & Borrowing (CRITICAL)
- `own-borrow-over-clone` — replace `.clone()` with `&T` where ownership is not needed
- `own-slice-over-vec` — function params `&Vec<T>` → `&[T]`, `&String` → `&str`
- `own-move-large` — move large `TreeSnapshot`/`Bitset` values instead of cloning
- `own-copy-small` — derive `Copy` for small trivially-copyable types (e.g., index newtypes)

### Error Handling (CRITICAL)
- `err-result-over-panic` — replace any `panic!`/`unwrap()`/`expect()` in `src/` lib code
- `err-no-unwrap-prod` — `.unwrap()` is forbidden in library code; use `?` or return `Result`
- `err-thiserror-lib` — custom error types must use `thiserror`; never `Box<dyn Error>`
- `err-question-mark` — replace explicit `match err { Err(e) => return Err(e) }` with `?`
- `err-context-chain` — add `.map_err(|e| MyError::Context { source: e, msg: "..." })` when context aids debugging
- `err-lowercase-msg` — error messages must be lowercase with no trailing punctuation

### Memory Optimization (CRITICAL)
- `mem-with-capacity` — use `Vec::with_capacity` / `HashMap::with_capacity` when size is known
- `mem-reuse-collections` — in loops, `clear()` and reuse collections rather than re-allocating
- `mem-avoid-format` — don't use `format!()` for string construction in hot paths
- `mem-boxed-slice` — fixed-size partition lists: prefer `Box<[T]>` over `Vec<T>`

### API Design (HIGH)
- `api-impl-into` — accept `impl Into<T>` for string-like parameters in public functions
- `api-must-use` — add `#[must_use]` to all `Result`-returning public functions
- `api-common-traits` — ensure `Debug`, `Clone`, `PartialEq` are derived on all public structs
- `api-default-impl` — implement `Default` for structs that have a sensible zero/empty state
- `api-from-not-into` — implement `From`, rely on blanket `Into`

### Performance Patterns (MEDIUM)
- `perf-iter-over-index` — replace `for i in 0..n { arr[i] }` with iterator combinators
- `perf-entry-api` — use `entry()` API for `HashMap` insert-or-update patterns
- `perf-collect-once` — avoid `.collect()` on intermediate iterators; chain lazily
- `opt-inline-small` — mark small hot functions (bitset ops, merge steps) with `#[inline]`

### Anti-patterns (REFERENCE)
- `anti-unwrap-abuse` — zero tolerance for `.unwrap()` in `src/`
- `anti-vec-for-slice` — `&Vec<T>` parameters → `&[T]`
- `anti-string-for-str` — `&String` parameters → `&str`
- `anti-index-over-iter` — numeric indexing in loops → iterators
- `anti-type-erasure` — `Box<dyn Trait>` → `impl Trait` where possible
- `anti-format-hot-path` — no `format!()` inside distance computation loops

## Scope

Files to review and refactor:
```
src/bitset.rs
src/distances.rs
src/interned.rs
src/io.rs
src/snapshot.rs
src/api.rs
src/main.rs   (anti-patterns only; CLI ergonomics are secondary)
```

Files to skip: `target/`, `.venv/`, `treetracer/`, `benches/` (benchmarks intentionally use
patterns like `black_box` that differ from library code).

## Procedure

1. Read each file in `src/` fully.
2. For each violation, record:
   - **File** and **line range**
   - **Rule ID** (e.g., `own-borrow-over-clone`)
   - **Current code** (brief excerpt)
   - **Proposed replacement** (brief excerpt)
   - **Rationale** (one sentence)
3. Group violations by severity (CRITICAL first).
4. Present the full list to the user.
5. For each approved fix, apply it with the `edit` tool.
6. After all edits, run validation:
   ```bash
   cargo fmt --check
   cargo clippy -- -D warnings
   cargo test --lib
   ```
   All three must pass. Fix any regressions before reporting done.

## Hard constraints

- **Never alter public Python-facing API signatures** in `api.rs`.
- **Never change distance algorithm semantics** — RF, WRF, KF results must be identical.
- **Preserve all existing doc comments** — update them if the surrounding code changes.
- **Do not introduce new dependencies** without explicit user approval.
- When in doubt between two idiomatic patterns, prefer the one that avoids allocations.
