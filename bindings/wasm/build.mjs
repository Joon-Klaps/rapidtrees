// Builds the npm package into bindings/wasm/pkg with wasm-pack: one build for browsers and bundlers (`--target web`) and one for Node (`--target nodejs`), joined by a hand-written package.json whose `exports` pick the right one. Run from anywhere with Node >= 18 and wasm-pack on the PATH:
//
//     node bindings/wasm/build.mjs

import { execFileSync } from "node:child_process";
import { copyFileSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const ROOT = join(HERE, "../..");
const PKG = join(HERE, "pkg");
const NAME = "rapidtrees";

const cargo = readFileSync(join(HERE, "Cargo.toml"), "utf8");
const version = /^version\s*=\s*"([^"]+)"/m.exec(cargo)?.[1];
if (!version) throw new Error("no version in bindings/wasm/Cargo.toml");

rmSync(PKG, { recursive: true, force: true });
for (const [target, dir] of [["web", "web"], ["nodejs", "node"]]) {
  execFileSync(
    "wasm-pack",
    ["build", HERE, "--release", "--target", target, "--out-dir", join(PKG, dir), "--out-name", NAME, "--no-pack"],
    { stdio: "inherit" },
  );
  // wasm-pack drops a `.gitignore` of `*` next to its output, which would hide the files from `npm publish`.
  rmSync(join(PKG, dir, ".gitignore"), { force: true });
}

const manifest = {
  name: NAME,
  version,
  description: "Fast all-pairs phylogenetic tree distances (Robinson-Foulds, weighted RF, Kuhner-Felsenstein) in the browser and Node, compiled from the rapidtrees Rust library to WebAssembly",
  license: "MIT",
  repository: { type: "git", url: "git+https://github.com/Joon-Klaps/rapidtrees.git", directory: "bindings/wasm" },
  homepage: "https://github.com/Joon-Klaps/rapidtrees#readme",
  bugs: "https://github.com/Joon-Klaps/rapidtrees/issues",
  keywords: ["phylogenetics", "tree-distance", "robinson-foulds", "bioinformatics", "wasm", "webassembly"],
  // Node resolves the `node` condition to the CommonJS build, which loads its .wasm synchronously. Everything else gets the ES module built for `--target web`, whose default export `init()` fetches the .wasm.
  main: `./node/${NAME}.js`,
  module: `./web/${NAME}.js`,
  types: `./web/${NAME}.d.ts`,
  exports: {
    ".": {
      types: `./web/${NAME}.d.ts`,
      node: `./node/${NAME}.js`,
      default: `./web/${NAME}.js`,
    },
    "./web/*": "./web/*",
    "./node/*": "./node/*",
    "./package.json": "./package.json",
  },
  files: ["web", "node", "README.md", "LICENSE"],
  sideEffects: false,
  engines: { node: ">=18" },
};
writeFileSync(join(PKG, "package.json"), `${JSON.stringify(manifest, null, 2)}\n`);
copyFileSync(join(ROOT, "LICENSE"), join(PKG, "LICENSE"));
copyFileSync(join(HERE, "README.md"), join(PKG, "README.md"));
console.log(`built ${NAME}@${version} in ${PKG}`);
