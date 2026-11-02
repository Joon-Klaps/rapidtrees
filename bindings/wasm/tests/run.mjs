// Tests the BUILT npm package the way a user gets it: packs `bindings/wasm/pkg`, installs the tarball into an empty project and runs `check.mjs` there, which imports "rapidtrees" by name. Needs `node bindings/wasm/build.mjs` first and the native CLI (built here if missing). Run from anywhere with Node >= 18:
//
//     node bindings/wasm/tests/run.mjs

import { execFileSync } from "node:child_process";
import { copyFileSync, existsSync, mkdtempSync, readdirSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const ROOT = join(HERE, "../../..");
const PKG = join(HERE, "../pkg");
if (!existsSync(join(PKG, "package.json"))) throw new Error("run `node bindings/wasm/build.mjs` first");

execFileSync("cargo", ["build", "--release", "--bin", "rapidtrees"], { cwd: ROOT, stdio: "inherit" });

const scratch = mkdtempSync(join(tmpdir(), "rapidtrees-npm-"));
try {
  const tarball = execFileSync("npm", ["pack", PKG, "--pack-destination", scratch, "--silent"], { encoding: "utf8" }).trim().split("\n").pop();
  const project = join(scratch, "project");
  execFileSync("mkdir", [project]);
  writeFileSync(join(project, "package.json"), JSON.stringify({ name: "rapidtrees-npm-test", private: true, type: "module" }));
  execFileSync("npm", ["install", "--silent", "--no-audit", "--no-fund", join(scratch, tarball)], { cwd: project, stdio: "inherit" });
  copyFileSync(join(HERE, "check.mjs"), join(project, "check.mjs"));
  console.log(`installed ${tarball}: ${readdirSync(join(project, "node_modules/rapidtrees")).join(", ")}`);
  execFileSync("node", ["check.mjs"], {
    cwd: project,
    stdio: "inherit",
    env: { ...process.env, RAPIDTREES_ROOT: ROOT, RAPIDTREES_SCRATCH: scratch },
  });
} finally {
  rmSync(scratch, { recursive: true, force: true });
}
