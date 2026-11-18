# JavaScript API reference

`npm install rapidtrees` gives the `rapidtrees` Rust library compiled to WebAssembly, for browsers and Node. It exposes three distance functions, named after the Python API, and a version function. The build is single-threaded (no rayon pool), so it runs in any browser tab and in Node without special headers.

Unlike the Python API, the JavaScript functions take the text of a whole tree file and read it with the same code as the native CLI: format detection, burn-in, `TRANSLATE` blocks and tree names all behave as they do there.

---

## Functions

| Function | Returns |
| --- | --- |
| `pairwise_rf(content, burnin_trees?, rooted?)` | `RfMatrix` — Robinson-Foulds, `matrix` is a `Uint32Array` |
| `pairwise_wrf(content, burnin_trees?, rooted?)` | `FloatMatrix` — weighted RF, `matrix` is a `Float64Array` |
| `pairwise_kf(content, burnin_trees?, rooted?)` | `FloatMatrix` — Kuhner-Felsenstein, `matrix` is a `Float64Array` |
| `rapidtrees_version()` | `string` — the version of the `rapidtrees` crate behind the build |

The three distance functions share the same call signature:

```js
pairwise_rf(
  content,       // string  — the text of a BEAST/NEXUS file (with or without TRANSLATE) or a plain Newick file
  burnin_trees,  // number  — drop this many trees from the start (default 0)
  rooted,        // boolean — compare clades (true) or bipartitions (false, the default)
)
```

Weighted RF and KF need branch lengths on every tree.

### Return types

`RfMatrix` and `FloatMatrix` have the same three getters:

| Getter | Type | Description |
| --- | --- | --- |
| `n_trees` | `number` | Trees kept after burn-in |
| `names` | `string[]` | One label per tree: the NEXUS tree name (`STATE_10000`), or `line<N>` for the Nth line of a Newick file |
| `matrix` | `Uint32Array` / `Float64Array` | The row-major `n_trees * n_trees` matrix: the distance between trees `i` and `j` is `matrix[i * n_trees + j]` |

Each read of `names` or `matrix` copies out of wasm memory, so read them once into a variable. Call `free()` on the result when you are done with it, to release the wasm memory it holds.

Unlike the CLI, the names carry no file-name prefix: `STATE_10000`, not `posterior_STATE_10000`.

### Errors thrown

All three distance functions throw an `Error` when:

- The input holds no trees
- Fewer than 2 trees are left after burn-in
- A tree fails to parse
- Trees have different leaf sets

---

## Examples

### Browser

The browser build is an ES module whose default export `init()` loads the `.wasm` file, so call it once before anything else. Bundlers such as Vite and webpack 5 pick the `.wasm` file up from the `new URL(...)` inside the module.

```js
import init, { pairwise_rf } from "rapidtrees";

await init();
const text = await fileInput.files[0].text();
const rf = pairwise_rf(text, 100); // drop the first 100 trees

const n = rf.n_trees;
const names = rf.names;
const matrix = rf.matrix;
rf.free();

console.log(`${names[0]} vs ${names[1]}: RF = ${matrix[0 * n + 1]}`);
```

### Node

Node (>= 18) resolves the same package name to a CommonJS build that loads its `.wasm` synchronously, so there is no `init()`. `import` and `require` both work.

```js
import { readFileSync } from "node:fs";
import { pairwise_kf } from "rapidtrees";

const kf = pairwise_kf(readFileSync("posterior.trees", "utf8"), 1000);
const matrix = kf.matrix;
console.log(kf.n_trees, matrix[1]);
kf.free();
```

```js
const { pairwise_rf, rapidtrees_version } = require("rapidtrees");
```

### Rows of the matrix

`matrix` is flat. To work with one row at a time without copying, take a `subarray`:

```js
const row = (matrix, n, i) => matrix.subarray(i * n, (i + 1) * n);

const distancesFromFirstTree = row(matrix, n, 0);
```

---

## Memory

The input is not streamed. The whole file goes in as one string, so the page or process holds the JS string, its UTF-8 copy in wasm linear memory and one chunk of trees while it is parsed, and wasm32 caps linear memory at 4 GiB. The native CLI and the Python package stream; use them for files that approach that size.

---

## Building and testing locally

The package is built from [`bindings/wasm`](../bindings/wasm), a separate crate that depends on `rapidtrees`; the `rapidtrees` crate itself never depends on `wasm-bindgen`. `node bindings/wasm/build.mjs` runs wasm-pack once for `--target web` and once for `--target nodejs`.

```bash
pixi run test-npm
```

builds the package, installs the packed tarball into an empty project and checks its RF, weighted RF and KF matrices against the native CLI. It needs wasm-pack, Node 18 or later and `rustup target add wasm32-unknown-unknown`.
