# rapidtrees for the browser and Node

Pairwise Robinson-Foulds (RF), weighted RF and Kuhner-Felsenstein (KF) distances between all trees of a posterior sample, computed by the [rapidtrees](https://github.com/Joon-Klaps/rapidtrees) Rust library compiled to WebAssembly. The build is single-threaded (no rayon pool), so it runs in any browser tab and in Node without special headers.

```bash
npm install rapidtrees
```

## API

All three functions take the text of a BEAST/NEXUS file (with or without a `TRANSLATE` block) or of a Newick file with one tree per line. They are named after the Python API.

| Function | Returns |
| --- | --- |
| `pairwise_rf(content, burnin_trees?, rooted?)` | `RfMatrix`: `matrix` is a `Uint32Array` |
| `pairwise_wrf(content, burnin_trees?, rooted?)` | `FloatMatrix`: `matrix` is a `Float64Array` |
| `pairwise_kf(content, burnin_trees?, rooted?)` | `FloatMatrix`: `matrix` is a `Float64Array` |
| `rapidtrees_version()` | the version of the `rapidtrees` crate behind the build, as a string |

`burnin_trees` drops that many trees from the start (default 0) and `rooted` compares clades instead of bipartitions (default `false`). A result has `n_trees`, `names` (the NEXUS tree names such as `STATE_10000`, or `line<N>` for Newick) and `matrix`, the row-major `n_trees * n_trees` matrix. Each read of `names` or `matrix` copies out of wasm memory, so read them once and call `free()` on the result when you are done. Functions throw an `Error` if no trees are left after burn-in, a tree fails to parse or the trees do not share one leaf set. Weighted RF and KF need branch lengths.

## Browser

The browser build is an ES module whose default export `init()` loads the `.wasm` file, so call it once before anything else. Bundlers such as Vite and webpack 5 pick the `.wasm` file up from the `new URL(...)` inside the module.

```js
import init, { pairwise_rf } from "rapidtrees";

await init();
const text = await fileInput.files[0].text();
const rf = pairwise_rf(text, 100);
console.log(rf.n_trees, rf.names, rf.matrix);
rf.free();
```

## Node

Node (>= 18) resolves the same name to a CommonJS build that loads its `.wasm` synchronously, so there is no `init()`. `import` and `require` both work.

```js
import { readFileSync } from "node:fs";
import { pairwise_kf } from "rapidtrees";

const kf = pairwise_kf(readFileSync("posterior.trees", "utf8"), 1000);
console.log(kf.n_trees, kf.matrix[1]);
kf.free();
```

## Memory

The input is not streamed. The whole file goes in as one string, so the page or process holds the JS string, its UTF-8 copy in wasm linear memory and the text of the trees being parsed, and wasm32 caps linear memory at 4 GiB. The native CLI and the Python package stream; use them for files that approach that size.

## Source

Built from [`bindings/wasm`](https://github.com/Joon-Klaps/rapidtrees/tree/master/bindings/wasm) in the rapidtrees repository with `node bindings/wasm/build.mjs`, which runs wasm-pack once for `--target web` and once for `--target nodejs`. MIT licence.
