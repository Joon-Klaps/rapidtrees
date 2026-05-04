---
name: rubber-duck-explanator
description: >
  Generates a detailed human-readable Markdown summary of all code changes
  made in the current branch compared to a base branch (default: master).
  Explains every addition, deletion, and modification at the function and
  module level. Run this last, after other skills have been applied.
  Invoke with /rubber-duck-explanator.
allowed-tools: shell, read, grep, glob
---

# Rubber Duck Explanator

You produce a thorough, plain-English explanation of every code change in the
current branch compared to a base branch. The output is written as a Markdown
file that any team member (including non-Rust experts) can understand.

## Inputs

- **Base branch**: `master` by default. The user may specify a different branch,
  e.g., `/rubber-duck-explanator base=develop`.
- **Output file**: `CHANGES.md` in the repository root by default. The user may
  specify a different path, e.g., `/rubber-duck-explanator output=docs/review.md`.

## Procedure

### Step 1 — Collect the diff

Run:
```bash
git diff <base-branch>...HEAD -- src/ tests/ benches/ pyproject.toml Cargo.toml
```
Also collect:
```bash
git log --oneline <base-branch>...HEAD
```

### Step 2 — Parse the diff

For each changed file, identify:
- **Added lines** (`+`)
- **Deleted lines** (`-`)
- **Unchanged context** lines (for orientation)
- The **file path** and **language** (Rust, Python, TOML, etc.)

### Step 2b — Critique naming and placement

Before explaining changes, run a pass over every new or renamed item (functions, types, modules, methods) and flag:

1. **Intuitiveness** — Is the name self-explanatory to someone unfamiliar with the codebase?
   - Flag names containing implementation-leaking adjectives (e.g., "Interned", "Raw", "Parsed") unless the distinction is truly user-visible.
   - Flag names with stuttering (e.g., `Snapshots::from_snapshots`).
   - Flag names that describe *how* rather than *what* (e.g., `parse_and_snapshot_newicks` vs `Snapshots::from_newick_iter`).
2. **Placement** — Is this function/method in the right file or module?
   - `io.rs` should contain only filesystem I/O (reading from and writing to long-term storage). Flag any function in `io.rs` that only processes in-memory data.
   - Data constructors (`from_*`, `new`, `build_*`) belong on the type they construct, not in a utility module.
   - Distance functions belong in `distances.rs`, snapshot types in `snapshot.rs`.
3. **Abbreviation clarity** — Are abbreviations unambiguous? (e.g., `wrf` is clear in this domain; `snaps` as a local variable is fine; `bip` for bipartition is OK but should be defined in a doc comment the first time.)
4. **API surface minimality** — Are there multiple public functions with overlapping purposes? Flag any case where two public functions differ only by their input source (file vs iterator vs Vec) without a clear ergonomic justification.

For each flagged item, write a short paragraph: what the name/placement is, why it is questionable, and what a better alternative would be.

### Step 3 — Explain every change

For **each modified function, struct, trait, impl block, or test**:

1. **What changed** — describe in plain English what was added, removed, or altered.
   Do not just repeat the code; explain the *intent and effect*.
2. **Why it matters** — explain how this affects behaviour, performance, correctness,
   or maintainability. If the change is a refactor with no behavioural change, say so.
3. **Risk assessment** — note any changes that could introduce regressions or require
   downstream updates (e.g., public API changes, PyO3 binding changes, algorithm changes).

For **configuration files** (`Cargo.toml`, `pyproject.toml`, CI YAML):
- List every dependency added/removed/bumped with its version and purpose.
- Flag any feature flag changes.

For **test files**:
- Describe what new scenarios are now covered.
- Note any removed tests and why their removal is safe (or risky).

### Step 4 — Structure the output document

Write the Markdown file with this structure:

```markdown
# Change Summary: <branch-name> → <base-branch>

**Generated**: <date>
**Commits**: <N commits>
**Files changed**: <N files>

---

## Overview

<2–4 sentence high-level summary of what this branch does>

---

## Detailed Changes

### `<file-path>`

#### `<function/struct/module name>` — <one-line summary>

**What changed:**
<paragraph>

**Why it matters:**
<paragraph>

**Risk:** Low / Medium / High — <brief justification>

---
... (repeat for each changed item) ...

---

## Breaking Changes

List any changes that break:
- Public Rust API (function signatures, types)
- Python API (PyO3-exposed functions, class fields)
- CLI interface (argument names, output format)

If none: write "None detected."

---

## Test Coverage

List new tests added and what scenarios they cover.
List any tests removed and why.

---

## Dependency Changes

List additions, removals, and version bumps in Cargo.toml / pyproject.toml.

---

## Recommended Follow-up

List any concerns, TODOs, or suggested follow-up actions that arose during
analysis — e.g., missing tests, documentation gaps, potential performance impacts.
```

### Step 5 — Write the file

Use the `edit` or `create` tool to write the completed Markdown to the output path.
If the file already exists, overwrite it.

Announce the output path and a one-line summary of the changes when done.

## Quality standards

- **Explain at the semantic level**: don't say "line 42 now reads `foo.bar()`";
  say "the snapshot is now built lazily, which avoids allocating memory for trees
  that are never compared."
- **Accessible language**: assume the reader understands Rust basics but not the
  internal algorithms. Briefly define domain terms (RF distance, bipartition, etc.)
  when they first appear.
- **Completeness**: every `+` and `-` line in the diff must be accounted for
  in the explanation — no change should be silently skipped.
- **Honesty about risk**: if a change modifies a distance-computation hot path,
  say so clearly and note that benchmarks should be re-run.

## Hard constraints

- Do **not** make any code changes — this skill is read-only.
- Do **not** summarise at such a high level that individual function changes are lost.
- The output file must be valid Markdown that renders correctly on GitHub.
