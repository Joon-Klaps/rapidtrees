# Contributing to rapidtrees

Thank you for your interest in contributing to `rapidtrees`! This guide will help you set up the development environment and understand our development workflow.

## 🏗️ Development Setup

### Prerequisites

- **Rust toolchain** (1.70+) — Install from [rustup.rs](https://rustup.rs/)
- **pixi** — Package manager for Python and R dependencies. Install from [pixi.sh](https://pixi.sh/)

### Quick Start

```bash
git clone https://github.com/Joon-Klaps/rapidtrees.git
cd rapidtrees

# Set up environment (creates isolated Python + R environment)
pixi install

# Install development version of the package
pixi run pip install -e .[dev]
```

This creates an isolated environment with:
- Python 3.10, 3.11, and 3.12 (for testing)
- R with the **phangorn** bioinformatics package
- Rust toolchain
- All Python dev dependencies (pytest, numpy, rpy2)

## 📋 Common Tasks

### Running Tests

```bash
# Run Rust unit tests
pixi run test-rust

# Run Python API tests (includes R comparisons with phangorn)
pixi run test-python

# Run Python tests with coverage report
pixi run test-python-cov
```

### Code Quality

```bash
# Format code (applies rustfmt)
pixi run fmt

# Lint with clippy (strict warnings-as-errors)
pixi run clippy

# Pre-commit checks
pixi run pre-commit
```

### Building

```bash
# Debug build
cargo build

# Release binary (optimized)
cargo build --release

# Python wheels (editable install for dev)
pip install -e .
```

## 📁 Project Structure

```
rapidtrees/
├── src/
│   ├── lib.rs                # Public API re-exports
│   ├── snapshot.rs           # Tree snapshot types & interning
│   ├── distances.rs          # RF, WRF, KF distance algorithms
│   ├── bitset.rs             # Efficient bitset for tree splits
│   ├── io.rs                 # BEAST/NEXUS file parsing
│   ├── main.rs               # CLI binary
│   └── api.rs                # PyO3 Python bindings
├── tests/
│   └── test_python_api.py    # Python API tests (with R comparisons)
├── benches/
│   └── memory_time_benchmark.rs  # Performance benchmarks
├── docs/
│   └── rapidtrees-for-dummies.md # Architecture guide
├── Cargo.toml                # Rust dependencies
├── pyproject.toml            # Python package metadata
└── .github/workflows/
    └── ci.yml                # CI/CD pipeline (GitHub Actions)
```

## 🔧 Rust Code Guidelines

- **Idiomatic Rust**: Prefer iterators, pattern matching, and ownership semantics
- **No `unsafe`**: Unless there's a documented performance justification (rare)
- **Error handling**: Use `?` operator and proper error types, not `.unwrap()`
- **Comments**: Only clarify non-obvious logic; clean code is self-documenting
- **Performance**: Hot paths (tree parsing, distance computation) should be benchmarked

### Example: Adding a new distance metric

1. **Add algorithm in `src/distances.rs`**:
   ```rust
   pub fn your_metric(a: &Snapshot, b: &Snapshot) -> f64 {
       // Implementation
   }

   #[cfg(test)]
   mod tests {
       #[test]
       fn test_your_metric() { /* ... */ }
   }
   ```

2. **Add PyO3 binding in `src/api.rs`**:
   ```rust
   #[pyfunction]
   fn pairwise_your_metric_from_newick_iter(...) { /* ... */ }
   ```

3. **Update `src/lib.rs`** if making public API changes

4. **Add Python tests** in `tests/test_python_api.py`

5. **Benchmark** with `benches/memory_time_benchmark.rs`

## 🐍 Python Bindings (PyO3)

- Every public Rust function exposed to Python needs:
  - Clear docstring (visible in `help()`)
  - Type-safe return types (prefer `PyBytes`, `PyList` over raw structs)
  - Meaningful Python exceptions (`PyValueError`, not panics)

Example:
```rust
/// Compute pairwise RF distances from newick strings.
///
/// Args:
///     names: Tree identifiers.
///     newick_iter: Python iterator of newick strings.
///
/// Returns:
///     (tree_names, distance_matrix_bytes)
///
/// Raises:
///     ValueError: If fewer than 2 trees or mismatched leaf sets.
#[pyfunction]
#[pyo3(signature = (names, newick_iter, ...))]
fn pairwise_rf_from_newick_iter(...) -> PyResult<(Vec<String>, Py<PyAny>)> {
    // Implementation
}
```

## 🧪 Testing Strategy

### Unit Tests (Rust)

Located inline in `#[cfg(test)]` modules:
- `bitset::tests` — bitset operations
- `distances::tests` — distance calculations with known values
- `snapshot::tests` — tree parsing and bipartition extraction
- `io::tests` — file I/O

Run with: `pixi run cargo test --lib`

### Integration Tests (Python)

`tests/test_python_api.py` tests the Python-facing API:
- Tree parsing from Newick strings
- Distance matrix computation (pairwise RF, WRF, KF)
- **Comparison with phangorn R package** for validation
- Edge cases (empty trees, single tree, mismatched taxa)

Run with: `pixi run pytest tests/test_python_api.py -v`

### Benchmarks

`benches/memory_time_benchmark.rs` measures:
- Memory usage for varying tree counts and taxa
- Wall-time and CPU-time for distance computation
- Estimated throughput for large datasets

Run with: `pixi run cargo bench`

## 🔄 CI/CD Pipeline

GitHub Actions (`.github/workflows/ci.yml`) runs on every push/PR:

1. **Rust Tests** (`ubuntu-latest`, `macos-latest`, Rust stable)
   - `cargo fmt --check`
   - `cargo clippy -- -D warnings`
   - `cargo test --lib`
   - `cargo build --release`

2. **Python Tests** (`ubuntu-latest`, `macos-latest`, Python 3.10/3.11/3.12)
   - Build PyO3 wheels with `maturin`
   - Run `pytest tests/test_python_api.py -v`
   - Includes phangorn R comparisons

**Environment**: Managed by pixi, ensuring Python + R are available in CI

## 📝 Pull Requests

Before submitting a PR:

1. **Run tests locally**: `pixi run test-rust && pixi run test-python`
2. **Format code**: `pixi run fmt`
3. **Lint**: `pixi run clippy`
4. **Update docs** if adding/changing public APIs
5. **Add tests** for new functionality
6. **Reference related issues** in your PR description

## 🚀 Release Process

1. Update version in `Cargo.toml` and `pyproject.toml`
2. Update `CHANGELOG.md`
3. Create a git tag: `git tag v0.4.1`
4. Push tag: `git push origin v0.4.1`
5. GitHub Actions publishes to PyPI and crates.io automatically

## ❓ Questions?

- Check `docs/rapidtrees-for-dummies.md` for architecture deep-dive
- Open an issue for design discussions
- See README for usage examples

---

Happy contributing! 🌲⚡

---

## 📋 Single Source of Truth for Versions

All project metadata is in **`pyproject.toml`**:
- Python package metadata (name, version, description, classifiers, URLs)
- Python, R, and Rust dependencies
- Build system configuration (maturin)
- Development tasks (under `[tool.pixi.tasks]`)

**For releases:** Update `version` in `pyproject.toml` and `Cargo.toml`, then update `CHANGELOG.md` with the new version and release notes. This ensures consistency across all project files and simplifies the release process.
