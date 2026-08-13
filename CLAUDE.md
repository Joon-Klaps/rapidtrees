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

This is a single-package repo: `rapidtrees` and nothing else.

**The split with treetracer is by audience, not by subject.** `rapidtrees`
answers one question — how far apart are two trees — and is consumed by Python
(PyO3) and a CLI. Everything the *browser* needs on top of that answer lives in
`treetracer-core`, which is **not in this repo**: it is `core/` in the
[treetracer-web](https://github.com/Joon-Klaps/treetracer-web) repo, and it
depends on this crate as a pinned git dependency.

```
Cargo.toml              — the rapidtrees package + [patch.crates-io]
src/                    — the distance kernel. No wasm-bindgen here.
  api.rs        — PyO3 bindings (keep minimal — logic lives in snapshot.rs / distances.rs)
  bitset.rs     — compact bitset for leaf sets; inner-most hot path
  snapshot.rs   — Snapshot / Snapshots types, canonicalisation, interning
  distances.rs  — RF, WRF, KF algorithms
  io.rs         — BEAST/NEXUS parsing, annotation stripping
  par.rs        — rayon-or-sequential shim; `pub` so dependents share one policy
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

### Rules that keep the split honest

1. **`rapidtrees` must never depend on `wasm-bindgen`.** `treetracer-core` is the
   only crate in the graph that links it. Multiple bindgen crates would force
   every one of them onto an identical wasm-bindgen version.
2. **`rapidtrees` stays wasm-*compatible*.** No Fortran-LAPACK linear algebra, no
   C-library codec deps, and every parallel loop goes through `par.rs`. That is
   what lets it link into a wasm build at all. Enforced by
   `cargo check --lib --no-default-features --target wasm32-unknown-unknown`,
   which runs in pre-commit and in CI.
3. **Browser-only algorithms do not belong here.** They go in treetracer-web's
   `core/`. The test is whether the Python API or CLI can reach it: if not, it
   is not a `rapidtrees` concern. MCMC diagnostics — log densities, chain
   splitting, burnin, subsampling — are the browser's, which is why BEAST
   *header* semantics live there while NEXUS *primitives* (`parse_taxon_block`,
   `extract_name_state`, `strip_beast_annotations`, `rename_leaf_nodes`) stay
   here and are `pub` for it to build on.
4. **Anything `pub` in `io.rs` is treetracer-web's API.** Changing one of those
   four signatures breaks that repo's build with no compiler to warn you, since
   it is a git dependency rather than a workspace member. Treat them as
   published surface.

### Building the browser bundle

Not from this repo. In treetracer-web:

```bash
wasm-pack build core --release --target web --out-dir pkg
```

That crate pins `rapidtrees` by git branch, so **local changes here are invisible
to it until pushed.** To build treetracer-web against a working tree, uncomment
the `[patch."https://github.com/Joon-Klaps/rapidtrees"]` block at the bottom of
`core/Cargo.toml`.

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
pixi run develop

# Rust tests
pixi run pre-commit          # cargo test --lib & cargo clippy & cargo fmt

# Python API tests
pixi run test-python         # pytest tests/test_python_api.py -v

# Benchmarks
pixi run cargo bench
```

### Pixi auto-sync gotcha

Pixi's `[tool.pixi.pypi-dependencies]` declares `rapidtrees = { path = ".", editable = true }`, which means **every `pixi run X` first re-syncs the editable install from pixi's own cached wheel**. If you just edited Rust code and then run any pixi task, that sync silently overwrites your fresh `maturin develop` build with the cached (old) `.so`.

Symptoms: tests pass on the Rust side (`cargo test`) but the Python side reports `AttributeError`/missing functions or shows behaviour from a previous build, with the installed `.so` mysteriously reverting to an older size/timestamp.

Workarounds (in order of preference):
1. After `pixi run develop`, run python tests with the env interpreter directly to bypass the sync:
   ```bash
   .pixi/envs/default/bin/python -m pytest tests/test_python_api.py
   ```
2. Force a fresh install by deleting the installed package first:
   ```bash
   rm -rf .pixi/envs/default/lib/python3.12/site-packages/rapidtrees*
   pixi run maturin develop --uv
   ```
3. If pixi's cache itself is poisoned: `uv cache clean rapidtrees`.

---

## CI requirements (must pass on every PR)

- `cargo fmt --check`
- `cargo clippy -- -D warnings`
- `cargo test --lib`
- `pytest tests/test_python_api.py -v`

Runs on `ubuntu-latest` and `macos-latest`, Python 3.10 / 3.11 / 3.12.
