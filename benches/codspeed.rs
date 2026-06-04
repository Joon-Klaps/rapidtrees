//! CodSpeed regression benches: snapshot construction + the dense pairwise
//! backends (RF, WRF, KF).
//!
//! This file is built with `codspeed-divan-compat`, which lets CodSpeed measure
//! every bench under two instruments from the *same* binary:
//!   - **simulation** — deterministic instruction counts (the CPU-speed gate);
//!   - **memory** — heap allocations tracked per bench (the memory gate that
//!     replaces the old standalone `memory_quick` bench).
//!
//! A PR that makes any bench measurably slower *or* heavier than `master` is
//! flagged. There are two bench families, so a regression is attributable:
//!   - **`construct_*`** measure `Snapshots::from_newicks` (Newick parsing +
//!     canonicalisation + interning) — the build is *inside* the measured region,
//!     so under the memory instrument this is where the persistent `Snapshots`
//!     footprint (the PR-1 interning/branch-length wins) is gated.
//!   - **`rf_* / wrf_* / kf_*`** measure only the pairwise call — the `Snapshots`
//!     are built *outside* the timed region, so simulation isolates the dense
//!     compute and memory captures just the transient result matrix.
//!
//! Two tree shapes are used, mirroring the target vs. adversarial regimes:
//!   - **similar** — one topology with a few random leaf swaps per tree (the
//!     MCMC-posterior regime rapidtrees targets, and the headline case).
//!   - **diverse** — every tree fully reshuffled, the adversarial large-`U` case.
//!
//! Inputs are deliberately small so the whole suite finishes well under a minute
//! even under the ~20–50× slowdown of instrumented execution.
//!
//! Build/run locally with `cargo codspeed build && cargo codspeed run` (add
//! `-m memory` for the memory instrument), or just `cargo bench --bench codspeed`
//! for a quick wall-clock sanity check.

use rapidtrees::Snapshots;

fn main() {
    divan::main();
}

// --- tree generation (self-contained) --------------

const TIPS: usize = 500;
const TREES: usize = 100;
const SWAPS: usize = 3;

fn build_balanced(labels: &[u32], out: &mut String) {
    if labels.len() == 1 {
        out.push_str("leaf_");
        out.push_str(&labels[0].to_string());
        out.push_str(":0.1");
        return;
    }
    let mid = labels.len() / 2;
    out.push('(');
    build_balanced(&labels[..mid], out);
    out.push(',');
    build_balanced(&labels[mid..], out);
    out.push_str("):0.1");
}

/// Small LCG so the bench needs no `rand` dependency.
fn next_rand(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

fn shuffled_labels(n: usize, seed: u64) -> Vec<u32> {
    let mut v: Vec<u32> = (0..n as u32).collect();
    let mut s = seed ^ 0x9E37_79B9_7F4A_7C15;
    for i in (1..n).rev() {
        let j = (next_rand(&mut s) as usize) % (i + 1);
        v.swap(i, j);
    }
    v
}

/// Identity labels with `swaps` random transpositions — shares most splits.
fn nearly_sorted_labels(n: usize, seed: u64, swaps: usize) -> Vec<u32> {
    let mut v: Vec<u32> = (0..n as u32).collect();
    let mut s = seed ^ 0xD1B5_4A32_D192_ED03;
    for _ in 0..swaps {
        let a = (next_rand(&mut s) as usize) % n;
        let b = (next_rand(&mut s) as usize) % n;
        v.swap(a, b);
    }
    v
}

fn newick_from_labels(labels: &[u32]) -> String {
    let mut s = String::new();
    build_balanced(labels, &mut s);
    s.push(';');
    s
}

fn similar_newicks() -> Vec<String> {
    (0..TREES)
        .map(|i| newick_from_labels(&nearly_sorted_labels(TIPS, i as u64, SWAPS)))
        .collect()
}

fn diverse_newicks() -> Vec<String> {
    (0..TREES)
        .map(|i| newick_from_labels(&shuffled_labels(TIPS, i as u64)))
        .collect()
}

fn snaps_from(newicks: &[String]) -> Snapshots {
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
    Snapshots::from_newicks(&refs, false).expect("parse failed")
}

// --- construction benches: measure Snapshots::from_newicks (parse + intern) -
// The newick strings are generated untimed; the measured region is the build,
// so the memory instrument gates the persistent `Snapshots` footprint here.

#[divan::bench]
fn construct_similar(bencher: divan::Bencher) {
    let newicks = similar_newicks();
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
    bencher.bench_local(|| Snapshots::from_newicks(&refs, false).expect("parse failed"));
}

#[divan::bench]
fn construct_diverse(bencher: divan::Bencher) {
    let newicks = diverse_newicks();
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
    bencher.bench_local(|| Snapshots::from_newicks(&refs, false).expect("parse failed"));
}

// --- pairwise benches: build snaps untimed, measure only the pairwise call --

#[divan::bench]
fn rf_similar(bencher: divan::Bencher) {
    let snaps = snaps_from(&similar_newicks());
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks());
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn wrf_similar(bencher: divan::Bencher) {
    let snaps = snaps_from(&similar_newicks());
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn wrf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks());
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn kf_similar(bencher: divan::Bencher) {
    let snaps = snaps_from(&similar_newicks());
    bencher.bench_local(|| snaps.pairwise_kf(None));
}

#[divan::bench]
fn kf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks());
    bencher.bench_local(|| snaps.pairwise_kf(None));
}
