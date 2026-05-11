# CLAUDE.md — rapidtrees

## Project overview

`rapidtrees` is a high-performance Rust library for pairwise phylogenetic tree
distance calculations (Robinson-Foulds, Weighted RF, Kuhner-Felsenstein), exposed
to Python via PyO3/maturin. The Rust core handles all computation; Python bindings
provide the user-facing API for bioinformatics workflows.

The primary client is the `treetracer/` project, which uses `rapidtrees` to power
convergence diagnostics, ESS calculations, and tanglegram visualisations.

---

## Architecture

```
src/
  api.rs        — PyO3 bindings (keep minimal — logic lives in snapshot.rs / distances.rs)
  bitset.rs     — compact bitset for leaf sets; inner-most hot path
  snapshot.rs   — Snapshot / Snapshots types, canonicalisation, interning
  distances.rs  — RF, WRF, KF algorithms
  io.rs         — BEAST/NEXUS parsing, annotation stripping
  main.rs       — CLI binary
  lib.rs        — public re-exports
tests/
  test_python_api.py   — Python integration tests (includes phangorn R comparisons)
benches/
  memory_time_benchmark.rs
docs/
  python-api.md          — full Python API reference
  rapidtrees-for-dummies.md
```

---

## Available skills

Invoke with `/skill-name` in chat.

| Skill | When to use |
| --- | --- |
| `/rustacean-nerd` | Review and rewrite Rust code for idiomatic style (ownership, error handling, performance) |
| `/code-line-minimizer` | Find and consolidate duplicated logic across `src/` and `tests/` |
| `/rubber-duck-explanator` | Generate a plain-English Markdown diff summary of the current branch vs `master` — run last, after other skills |

---

## Code style

### Rust

- Idiomatic Rust: iterators, pattern matching, ownership semantics — no raw index loops in hot paths
- No `unsafe` without a documented, measurable justification
- No `.unwrap()` / `.expect()` in library code — propagate with `?`
- `thiserror` for error types; never `Box<dyn Error>` in library code
- `Vec::with_capacity` / `HashMap::with_capacity` whenever size is known
- Small hot functions (`Bitset` ops, merge steps) get `#[inline]`

### api.rs

- Keep `api.rs` thin: validation, serialisation, and return-value assembly only
- All computation belongs in `snapshot.rs` or `distances.rs`
- Every `#[pyfunction]` needs a Rust doc comment that is visible from Python `help()`
- Raise `PyValueError` — never panic from binding code

### Python tests

- Tests live in `tests/test_python_api.py` as pytest classes
- New functions need at least one type/shape test and one known-value regression test
- Cross-validate distance values against phangorn (R) where practical

---

## Hard rules

1. **Keep Rust and Python APIs in sync** — any change to a public function signature must update the PyO3 binding, the Python tests, and the docs (`docs/python-api.md`, `README.md`) in the same PR.
2. **Document all public API** — every `pub fn`, `pub struct`, and `#[pyfunction]` must have a doc comment.
3. **Always add tests** — new functions need unit tests; new distance metrics need a known-value regression test.
4. **Update `CHANGELOG.md`** for every user-visible change.

---

## Common commands

```bash
# Build & install (development)
pixi run pip install -e .

# Rust tests
pixi run test-rust           # cargo test --lib

# Python API tests
pixi run test-python         # pytest tests/test_python_api.py -v

# Format / lint
pixi run fmt                 # cargo fmt
pixi run clippy              # cargo clippy -- -D warnings

# Benchmarks
pixi run cargo bench
```

---

## CI requirements (must pass on every PR)

- `cargo fmt --check`
- `cargo clippy -- -D warnings`
- `cargo test --lib`
- `pytest tests/test_python_api.py -v`

Runs on `ubuntu-latest` and `macos-latest`, Python 3.10 / 3.11 / 3.12.
