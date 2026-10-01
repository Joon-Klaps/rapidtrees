// Runtime smoke test for the wasm32 build of rapidtrees.
//
// Builds `tests/wasm` for wasm32-unknown-unknown with no default features, loads it in Node, feeds it each BEAST file in `tests/data`, and checks its RF, weighted RF and KF matrices against the native CLI on the same file. Run from the repo root with Node >= 22.18 (type stripping is on by default):
//
//     node tests/wasm/run.ts

import { execFileSync } from "node:child_process";
import { mkdtempSync, readFileSync, readdirSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const ROOT = resolve(import.meta.dirname, "../..");
const DATA = join(ROOT, "tests/data");
const METRICS = ["rf", "weighted", "kf"] as const;
const RTOL = 1e-9;

interface Smoke {
  memory: WebAssembly.Memory;
  alloc(len: number): number;
  run(ptr: number, len: number): number;
  matrix(metric: number): number;
}

function cargo(args: string[]): void {
  execFileSync("cargo", args, { cwd: ROOT, stdio: "inherit" });
}

async function loadWasm(): Promise<Smoke> {
  cargo(["build", "--release", "--target", "wasm32-unknown-unknown", "--manifest-path", "tests/wasm/Cargo.toml"]);
  const bytes = readFileSync(join(ROOT, "tests/wasm/target/wasm32-unknown-unknown/release/rapidtrees_wasm_smoke.wasm"));
  const module = await WebAssembly.compile(bytes);
  // No imports means no hidden dependency on wasm-bindgen glue or a WASI host: the module runs anywhere a browser can.
  const imports = WebAssembly.Module.imports(module);
  if (imports.length > 0) {
    throw new Error(`wasm module expects host imports: ${imports.map((i) => `${i.module}.${i.name}`).join(", ")}`);
  }
  const instance = await WebAssembly.instantiate(module, {});
  return instance.exports as unknown as Smoke;
}

/** A fresh instance per file would hide leaks between runs; reusing one is closer to a long-lived browser tab. */
function runWasm(wasm: Smoke, content: Uint8Array): Float64Array[] {
  const ptr = wasm.alloc(content.length);
  new Uint8Array(wasm.memory.buffer, ptr, content.length).set(content);
  const n = wasm.run(ptr, content.length);
  if (n <= 0) throw new Error(`run() returned ${n}`);
  // Read the matrices only after `run`: it may grow memory, which detaches earlier views of `memory.buffer`.
  return METRICS.map((_, m) => new Float64Array(wasm.memory.buffer, wasm.matrix(m), n * n).slice());
}

function runNative(file: string, metric: string, out: string): Float64Array {
  execFileSync(join(ROOT, "target/release/rapidtrees"), ["-i", file, "-o", out, "--metric", metric, "--quiet"]);
  const rows = readFileSync(out, "utf8").trimEnd().split("\n").slice(1);
  return Float64Array.from(rows.flatMap((row) => row.split("\t").slice(1).map(Number)));
}

function compare(label: string, got: Float64Array, want: Float64Array): number {
  if (got.length !== want.length) {
    console.error(`  ${label}: ${got.length} cells from wasm, ${want.length} from native`);
    return 1;
  }
  const bad = got.findIndex((g, i) => !(Math.abs(g - want[i]) <= RTOL * Math.max(1, Math.abs(want[i]))));
  if (bad >= 0) {
    const n = Math.sqrt(want.length);
    console.error(`  ${label}: cell (${Math.floor(bad / n)}, ${bad % n}) is ${got[bad]} in wasm, ${want[bad]} native`);
    return 1;
  }
  return 0;
}

const files = readdirSync(DATA).filter((f) => f.endsWith(".trees")).sort();
if (files.length === 0) throw new Error(`no .trees files in ${DATA} (is git-lfs pulled?)`);

const wasm = await loadWasm();
cargo(["build", "--release", "--bin", "rapidtrees"]);

const scratch = mkdtempSync(join(tmpdir(), "rapidtrees-wasm-"));
let failures = 0;
try {
  for (const name of files) {
    const path = join(DATA, name);
    const wasmMatrices = runWasm(wasm, readFileSync(path));
    const n = Math.sqrt(wasmMatrices[0].length);
    const fileFailures = METRICS.reduce(
      (sum, metric, m) => sum + compare(metric, wasmMatrices[m], runNative(path, metric, join(scratch, `${name}.${metric}.tsv`))),
      0,
    );
    console.log(`${fileFailures === 0 ? "ok  " : "FAIL"} ${name}: ${n} trees, ${METRICS.join(" / ")}`);
    failures += fileFailures;
  }
} finally {
  rmSync(scratch, { recursive: true, force: true });
}

if (failures > 0) {
  console.error(`${failures} matrix mismatch(es) between wasm and native`);
  process.exit(1);
}
