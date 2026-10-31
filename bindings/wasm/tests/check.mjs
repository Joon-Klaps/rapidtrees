// Runs inside the project that `run.mjs` installed the package into. Computes RF, weighted RF and KF with the package on the smallest BEAST file and its Newick twin in tests/data, and checks every cell against the native CLI at the relative tolerance tests/wasm/run.ts uses.

import { execFileSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { join } from "node:path";
import { pathToFileURL } from "node:url";
import * as rapidtrees from "rapidtrees";

const ROOT = process.env.RAPIDTREES_ROOT;
const SCRATCH = process.env.RAPIDTREES_SCRATCH;
const RTOL = 1e-9;
const FUNCTIONS = { rf: "pairwise_rf", weighted: "pairwise_wrf", kf: "pairwise_kf" };
const FILES = ["hiv1.trees", "hiv1.newick"];

let failures = 0;
function check(ok, message) {
  if (!ok) {
    failures += 1;
    console.error(`FAIL ${message}`);
  }
}

function native(file, metric) {
  const out = join(SCRATCH, `${file}.${metric}.tsv`);
  execFileSync(join(ROOT, "target/release/rapidtrees"), ["-i", join(ROOT, "tests/data", file), "-o", out, "--metric", metric, "--quiet"]);
  const [header, ...rows] = readFileSync(out, "utf8").trimEnd().split("\n");
  const prefix = `${file.split(".")[0]}_`;
  return {
    names: header.split("\t").slice(1).map((name) => name.replace(prefix, "")),
    values: Float64Array.from(rows.flatMap((row) => row.split("\t").slice(1).map(Number))),
  };
}

function sameCells(label, got, want) {
  if (got.length !== want.length) return check(false, `${label}: ${got.length} cells, native ${want.length}`);
  const bad = got.findIndex((g, i) => !(Math.abs(g - want[i]) <= RTOL * Math.max(1, Math.abs(want[i]))));
  check(bad < 0, `${label}: cell ${bad} is ${got[bad]}, native ${want[bad]}`);
}

check(rapidtrees.rapidtrees_version() === JSON.parse(readFileSync("node_modules/rapidtrees/package.json", "utf8")).version, "wasm and package.json versions differ");

for (const file of FILES) {
  const content = readFileSync(join(ROOT, "tests/data", file), "utf8");
  for (const [metric, fn] of Object.entries(FUNCTIONS)) {
    const result = rapidtrees[fn](content);
    const want = native(file, metric);
    const n = result.n_trees;
    check(result.matrix instanceof (metric === "rf" ? Uint32Array : Float64Array), `${file} ${metric}: wrong typed array`);
    check(result.matrix.length === n * n, `${file} ${metric}: ${result.matrix.length} cells for ${n} trees`);
    check(JSON.stringify(result.names) === JSON.stringify(want.names), `${file} ${metric}: tree names differ from native`);
    sameCells(`${file} ${metric}`, result.matrix, want.values);

    // Burn-in drops trees from the front, so the rest is the bottom-right block of the full matrix.
    const burnin = 3;
    const cut = rapidtrees[fn](content, burnin);
    const m = cut.n_trees;
    check(m === n - burnin, `${file} ${metric}: burn-in left ${m} of ${n} trees`);
    const block = Float64Array.from({ length: m * m }, (_, i) => want.values[(Math.floor(i / m) + burnin) * n + (i % m) + burnin]);
    sameCells(`${file} ${metric} burn-in ${burnin}`, cut.matrix, block);
    result.free();
    cut.free();
  }
  console.log(`${failures === 0 ? "ok  " : "FAIL"} ${file}: rf / weighted / kf against the native CLI`);
}

try {
  rapidtrees.pairwise_rf(readFileSync(join(ROOT, "tests/data", FILES[0]), "utf8"), 1000);
  check(false, "burn-in past the last tree did not throw");
} catch (e) {
  check(e instanceof Error && /no trees/.test(e.message), `unexpected error: ${e}`);
}

// The CommonJS entry that `require` resolves to is the same Node build.
const viaRequire = createRequire(import.meta.url)("rapidtrees");
check(viaRequire.rapidtrees_version() === rapidtrees.rapidtrees_version(), "require() entry differs from import");

// The browser build is a second wasm module: initialise it by hand, the way a bundler-less page would, and check it agrees with the Node build.
const web = await import(pathToFileURL(join(process.cwd(), "node_modules/rapidtrees/web/rapidtrees.js")).href);
web.initSync({ module: readFileSync(join(process.cwd(), "node_modules/rapidtrees/web/rapidtrees_bg.wasm")) });
const content = readFileSync(join(ROOT, "tests/data", FILES[0]), "utf8");
sameCells("web build kf", web.pairwise_kf(content).matrix, rapidtrees.pairwise_kf(content).matrix);
sameCells("web build rf", web.pairwise_rf(content).matrix, rapidtrees.pairwise_rf(content).matrix);
console.log(`${failures === 0 ? "ok  " : "FAIL"} web build agrees with the Node build`);

if (failures > 0) {
  console.error(`${failures} failure(s)`);
  process.exit(1);
}
