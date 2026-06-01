# rapidtrees in the browser (WebAssembly)

`rapidtrees` compiles to WebAssembly and ships a tiny static site that computes
pairwise phylogenetic tree distances **entirely in the browser**. You upload a
BEAST `.trees` file, pick a metric (RF, Weighted RF, or Kuhner–Felsenstein), and
the resulting distance matrix is streamed straight to a file on your disk.
Nothing is uploaded to a server.

The site lives in [`www/`](../www) and is deployed to GitHub Pages by
[`.github/workflows/deploy-pages.yml`](../.github/workflows/deploy-pages.yml).

---

## Why a separate Newick parser?

The native (CLI / Python) builds parse Newick with the [`phylotree`] crate. That
crate pulls in C-backed dependencies (`needletail`, `ptree`, …) that do not
compile for `wasm32-unknown-unknown`, so the `wasm` feature swaps in a small
in-crate parser, [`src/phylotree-patch.rs`](../src/phylotree-patch.rs), that
mirrors exactly the slice of the `phylotree` API the core uses.

A native test (`cargo test`) cross-checks the two parsers: it asserts they
produce structurally identical trees (same bipartitions and branch lengths) for
the same Newick, which guarantees identical distance matrices. The two backends
are selected at compile time:

```rust
#[cfg(feature = "phylotree")] use phylotree::tree::{Tree as PhyloTree, TreeError};
#[cfg(all(feature = "wasm", not(feature = "phylotree")))]
use crate::phylotree_patch::{Tree as PhyloTree, TreeError};
```

The `wasm` feature also drops `rayon` (GitHub Pages cannot serve the COOP/COEP
headers `SharedArrayBuffer`-based threading needs), falling back to sequential
iterators.

[`phylotree`]: https://crates.io/crates/phylotree

---

## Building locally

You need the wasm target and [`wasm-pack`]:

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack            # or download a prebuilt binary
```

Then build the bundle and serve the site:

```bash
pixi run build-wasm                # → www/pkg/{rapidtrees.js, rapidtrees_bg.wasm}
pixi run serve-wasm                # → http://localhost:8000
```

Equivalently, without pixi:

```bash
wasm-pack build --target web --out-dir www/pkg --no-default-features --features wasm
python -m http.server --directory www
```

`www/pkg/` is generated output and is git-ignored.

> **Note:** `wasm-opt` is disabled in `Cargo.toml`
> (`[package.metadata.wasm-pack.profile.release] wasm-opt = false`) so the build
> never needs to download the binaryen toolchain. To produce a smaller bundle,
> install binaryen and flip that flag to `true`.

[`wasm-pack`]: https://rustwasm.github.io/wasm-pack/

---

## How the page works

`www/main.js`:

1. Loads the wasm module (`import init, { compute_distances } from './pkg/rapidtrees.js'`).
2. Reads the uploaded file as text and calls
   `compute_distances(content, method, rooted, burninTrees, burninStates, useRealTaxa)`.
3. Receives a `DistanceResult` exposing `names`, `n`, `dtype` (`"u32"` for RF,
   `"f64"` for WRF/KF), and `take_matrix()` — the matrix as a flat,
   little-endian `Uint8Array`.
4. Streams the matrix to disk **one CSV row at a time**, decoding each value
   with a `DataView`. With the gzip toggle on, rows are piped through the
   browser's native `CompressionStream('gzip')` so no second large buffer is
   built.

Saving uses the [File System Access API] (`showSaveFilePicker`) on Chromium
browsers, with a `Blob` + `<a download>` fallback for Firefox/Safari. The
fallback buffers the whole file in memory, so very large matrices may hit
browser memory limits there.

[File System Access API]: https://developer.mozilla.org/en-US/docs/Web/API/Window/showSaveFilePicker

---

## Deploying to GitHub Pages

The `deploy-pages.yml` workflow builds the bundle and publishes `www/` on every
push to `master`/`main` (and on manual dispatch). One-time setup: in the
repository **Settings → Pages**, set **Source** to **GitHub Actions**.

---

## Limitations

- Single-threaded (no `rayon`); large datasets compute sequentially.
- `wasm32` has a ~4 GB address-space ceiling, so extremely large tree sets may
  exhaust memory during computation.
- The in-crate parser handles standard Newick (nested clades, branch lengths
  incl. scientific notation, single-quoted labels, `[...]` comments). It keeps
  underscores literal rather than converting them to spaces.
