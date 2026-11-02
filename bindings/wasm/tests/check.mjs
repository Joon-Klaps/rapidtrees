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

// Where the package once read trees itself and disagreed with the CLI. It now calls the library's reader, and these keep it that way.
const QUARTETS = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"];
const MRBAYES = "#NEXUS\nbegin trees;\n   translate\n       1 A,\n       2 B,\n       3 C,\n       4 D;\n   tree gen.0 = [&U] ((1:1,2:1):1,(3:1,4:1):1);\n   tree gen.100 = [&U] ((1:1,3:1):1,(2:1,4:1):1);\nend;\n";
const READER_CASES = [
  ["two trees on one line", `${QUARTETS[0]}${QUARTETS[1]}\n`, ["line1", "line1"]],
  ["a tree wrapped over lines", `((A:1,\nB:1):1,\n(C:1,D:1):1);\n${QUARTETS[1]}\n`, ["line1", "line4"]],
  ["a byte-order mark", `\u{feff}${QUARTETS[0]}\n${QUARTETS[1]}\n`, ["line1", "line2"]],
  ["a comment line before the trees", `[written by some tool]\n${QUARTETS[0]}\n${QUARTETS[1]}\n`, ["line2", "line3"]],
  ["ape's TREE * STATE_n headers", `#NEXUS\nBEGIN TREES;\n\tTREE * STATE_1 = [&R] ${QUARTETS[0]}\n\tTREE * STATE_2 = [&R] ${QUARTETS[1]}\nEND;\n`, ["STATE_1", "STATE_2"]],
  ["MrBayes gen.N headers", MRBAYES, ["gen.0", "gen.100"]],
];
for (const [label, content, names] of READER_CASES) {
  const rf = rapidtrees.pairwise_rf(content);
  check(JSON.stringify(rf.names) === JSON.stringify(names), `${label}: names ${JSON.stringify(rf.names)}`);
  check(JSON.stringify([...rf.matrix]) === "[0,2,2,0]", `${label}: RF ${[...rf.matrix]}`);
  rf.free();
}
console.log(`${failures === 0 ? "ok  " : "FAIL"} ${READER_CASES.length} reader edge cases match the CLI's reader`);

// One tree is too few for a pairwise matrix, as in the CLI and the Python package.
const TOO_FEW = [
  ["burn-in past the last tree", readFileSync(join(ROOT, "tests/data", FILES[0]), "utf8"), 1000],
  ["a single tree", `${QUARTETS[0]}\n`, 0],
];
for (const [label, content, burnin] of TOO_FEW) {
  try {
    rapidtrees.pairwise_rf(content, burnin);
    check(false, `${label} did not throw`);
  } catch (e) {
    check(e instanceof Error && /at least 2 trees/.test(e.message), `${label}: unexpected error ${e}`);
  }
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
