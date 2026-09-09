//! CodSpeed regression benches: snapshot construction + the pairwise backends
//! (RF, WRF, KF).
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
//!     footprint is gated.
//!   - **`rf_* / wrf_* / kf_*`** measure only the pairwise call — the `Snapshots`
//!     are built *outside* the timed region, so simulation isolates the compute
//!     and memory captures just the transient result matrix.
//!
//! Two tree shapes are used, mirroring the target vs. adversarial regimes:
//!   - **similar** — one topology with a few random leaf swaps per tree (the
//!     MCMC-posterior regime rapidtrees targets, and the headline case).
//!   - **diverse** — every tree fully reshuffled, the adversarial large-`U` case.
//!
//! # Sizing
//!
//! `SMALL` is the historical shape, kept as a cheap regression guard. It is too
//! small to measure a kernel change: at that size most of the profile is rayon
//! workers parked in `wait_until_cold`, and a real improvement is
//! indistinguishable from scheduling noise. The larger cells exist so that a
//! kernel change is visible:
//!   - `SWEEP` — many trees at moderate taxa, where the O(t²) pair sweep
//!     dominates;
//!   - `WIDE` — many taxa at few trees, where the O(t·n²/64) per-tree build
//!     dominates. The crossover sits near `trees ≈ taxa²/1000`, so the two cells
//!     straddle it and a build-side change shows up on exactly one of them.
//!
//! Every large cell also has a `_st` variant pinned to a one-thread rayon pool.
//! With the scheduler out of the picture the measurement is of the algorithm,
//! and it is what a 1-thread-vs-1-thread comparison against other tools needs.
//!
//! `DIVERSE` is deliberately smaller than `SWEEP`: on a fully reshuffled set the
//! dense sweep costs `⌈U/64⌉` words per pair with `U ≈ trees × (taxa − 3)`, so
//! `SWEEP`'s tree count there would run for minutes.
//!
//! Each pairwise bench prints `U` and the resulting per-pair word count to
//! stderr as it builds. Those are implementation-independent and explain a
//! timing change without needing a stable machine.
//!
//! Build/run locally with `cargo codspeed build && cargo codspeed run` (add
//! `-m memory` for the memory instrument), or just `cargo bench --bench codspeed`
//! for a quick wall-clock sanity check.

use rapidtrees::Snapshots;

fn main() {
    divan::main();
}

// --- shapes -------------------------------------------------------------

/// `(taxa, trees)` per cell.
type Shape = (usize, usize);

/// Historical shape: a fast regression guard, too small to measure a kernel.
const SMALL: Shape = (500, 100);
/// Pair-sweep dominates.
const SWEEP: Shape = (500, 1500);
/// Per-tree build dominates.
const WIDE: Shape = (4000, 50);
/// Adversarial large-`U` cell; kept small because dense cost grows with trees.
const DIVERSE: Shape = (500, 300);

/// Leaf swaps per tree in the `similar` regime.
const SWAPS: usize = 3;

// --- tree generation (self-contained) --------------

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

fn similar_newicks((tips, trees): Shape) -> Vec<String> {
    (0..trees)
        .map(|i| newick_from_labels(&nearly_sorted_labels(tips, i as u64, SWAPS)))
        .collect()
}

fn diverse_newicks((tips, trees): Shape) -> Vec<String> {
    (0..trees)
        .map(|i| newick_from_labels(&shuffled_labels(tips, i as u64)))
        .collect()
}

fn snaps_from(newicks: &[String]) -> Snapshots {
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
    Snapshots::from_newicks(&refs, false).expect("parse failed")
}

/// Build the snapshots and report the diversity measures behind the timing.
fn snaps_reported(label: &str, newicks: &[String]) -> Snapshots {
    let snaps = snaps_from(newicks);
    let (distinct, universal, mean_row) = snaps.split_stats();
    let u = distinct - universal;
    eprintln!(
        "{label}: trees={} U={u} rf_words_per_pair={} mean_splits_per_tree={mean_row:.0}",
        snaps.len(),
        u.div_ceil(64),
    );
    snaps
}

/// A one-thread rayon pool: the algorithm without the scheduler.
fn single_thread_pool() -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("thread pool")
}

// --- construction benches: measure Snapshots::from_newicks (parse + intern) -
// The newick strings are generated untimed; the measured region is the build,
// so the memory instrument gates the persistent `Snapshots` footprint here.

#[divan::bench]
fn construct_similar(bencher: divan::Bencher) {
    let newicks = similar_newicks(SMALL);
    bencher.bench_local(|| snaps_from(&newicks));
}

#[divan::bench]
fn construct_diverse(bencher: divan::Bencher) {
    let newicks = diverse_newicks(SMALL);
    bencher.bench_local(|| snaps_from(&newicks));
}

#[divan::bench]
fn construct_similar_wide(bencher: divan::Bencher) {
    let newicks = similar_newicks(WIDE);
    bencher.bench_local(|| snaps_from(&newicks));
}

#[divan::bench]
fn construct_similar_wide_st(bencher: divan::Bencher) {
    let newicks = similar_newicks(WIDE);
    let pool = single_thread_pool();
    bencher.bench_local(|| pool.install(|| snaps_from(&newicks)));
}

// --- pairwise benches: build snaps untimed, measure only the pairwise call --

#[divan::bench]
fn rf_similar(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_similar", &similar_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_diverse", &diverse_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_similar_sweep(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_similar_sweep", &similar_newicks(SWEEP));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_similar_sweep_st(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_similar_sweep_st", &similar_newicks(SWEEP));
    let pool = single_thread_pool();
    bencher.bench_local(|| pool.install(|| snaps.pairwise_rf(None)));
}

#[divan::bench]
fn rf_diverse_sweep(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_diverse_sweep", &diverse_newicks(DIVERSE));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_diverse_sweep_st(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_diverse_sweep_st", &diverse_newicks(DIVERSE));
    let pool = single_thread_pool();
    bencher.bench_local(|| pool.install(|| snaps.pairwise_rf(None)));
}

#[divan::bench]
fn rf_similar_wide(bencher: divan::Bencher) {
    let snaps = snaps_reported("rf_similar_wide", &similar_newicks(WIDE));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn wrf_similar(bencher: divan::Bencher) {
    let snaps = snaps_reported("wrf_similar", &similar_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn wrf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_reported("wrf_diverse", &diverse_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn wrf_diverse_sweep(bencher: divan::Bencher) {
    let snaps = snaps_reported("wrf_diverse_sweep", &diverse_newicks(DIVERSE));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn kf_similar(bencher: divan::Bencher) {
    let snaps = snaps_reported("kf_similar", &similar_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_kf(None));
}

#[divan::bench]
fn kf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_reported("kf_diverse", &diverse_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_kf(None));
}
