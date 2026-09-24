//! CodSpeed regression benches: snapshot construction + the pairwise backends
//! (RF, WRF, KF).
//!
//! Two kinds of tree set. `posterior` is what rapidtrees is for: one random
//! base topology, and in every tree each internal edge collapsed with
//! probability [`DIVERSITY`] and re-resolved at random, the generator of the
//! manuscript benchmark (`simulate_trees.py --diversity`). Splits are held by
//! every share of the trees, as in a real BEAST posterior. `diverse` is the
//! adversarial case: independent topologies, where almost every split is
//! unique to one tree.
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

/// Share of internal edges each `posterior` tree collapses and re-resolves.
/// Normalised pairwise RF comes out near `1 − (1 − q)²`, so 0.4 lands on the
/// real BEAST posteriors the manuscript measured (0.49–0.69).
const DIVERSITY: f64 = 0.4;

/// Mean branch length; lengths are exponential, as in the simulator.
const BRANCH_SCALE: f64 = 0.1;

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

/// A rooted tree over leaves `0..n`, with any arity while it is being built.
enum Node {
    Leaf(u32),
    Inner(Vec<Node>),
}

/// Uniform in `[0, 1)` from the LCG's 31 bits.
fn uniform(state: &mut u64) -> f64 {
    next_rand(state) as f64 / (1u64 << 31) as f64
}

/// A random binary tree over `comps`, by joining random pairs until one is
/// left.
fn resolve(mut comps: Vec<Node>, state: &mut u64) -> Node {
    while comps.len() > 1 {
        let a = comps.swap_remove((next_rand(state) as usize) % comps.len());
        let b = comps.swap_remove((next_rand(state) as usize) % comps.len());
        comps.push(Node::Inner(vec![a, b]));
    }
    comps.pop().expect("at least one component")
}

/// What `node` hands its parent once its internal child edges have each been
/// collapsed with probability `q` (the child's own components are absorbed)
/// or kept (the child is re-resolved into one subtree). A kept edge's split
/// survives exactly; a collapsed one's is replaced by random ones.
fn components(node: &Node, q: f64, state: &mut u64) -> Vec<Node> {
    let Node::Inner(children) = node else {
        unreachable!("only internal nodes are expanded")
    };
    let mut comps = Vec::with_capacity(children.len());
    for child in children {
        match child {
            Node::Leaf(leaf) => comps.push(Node::Leaf(*leaf)),
            Node::Inner(_) => {
                let sub = components(child, q, state);
                if uniform(state) < q {
                    comps.extend(sub);
                } else {
                    comps.push(resolve(sub, state));
                }
            }
        }
    }
    comps
}

/// Newick for `node`, with exponential branch lengths of mean
/// [`BRANCH_SCALE`] on every edge.
fn write_newick(node: &Node, state: &mut u64, out: &mut String) {
    match node {
        Node::Leaf(leaf) => {
            out.push_str("leaf_");
            out.push_str(&leaf.to_string());
        }
        Node::Inner(children) => {
            out.push('(');
            for (k, child) in children.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                write_newick(child, state, out);
                let length = -(1.0 - uniform(state)).ln() * BRANCH_SCALE;
                out.push_str(&format!(":{length:.6}"));
            }
            out.push(')');
        }
    }
}

fn newick_from_labels(labels: &[u32]) -> String {
    let mut s = String::new();
    build_balanced(labels, &mut s);
    s.push(';');
    s
}

/// A posterior-like set: one random base topology, re-resolved per tree at
/// [`DIVERSITY`]. Seeded, so every run benchmarks the same trees.
fn posterior_newicks((tips, trees): Shape) -> Vec<String> {
    let mut state = 0x5EED_0000_0000_0040;
    let base = resolve((0..tips as u32).map(Node::Leaf).collect(), &mut state);
    (0..trees)
        .map(|_| {
            let tree = resolve(components(&base, DIVERSITY, &mut state), &mut state);
            let mut newick = String::new();
            write_newick(&tree, &mut state, &mut newick);
            newick.push(';');
            newick
        })
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
fn construct_posterior(bencher: divan::Bencher) {
    let newicks = posterior_newicks(SMALL);
    bencher.bench_local(|| snaps_from(&newicks));
}

#[divan::bench]
fn construct_diverse(bencher: divan::Bencher) {
    let newicks = diverse_newicks(SMALL);
    bencher.bench_local(|| snaps_from(&newicks));
}

#[divan::bench]
fn construct_posterior_wide(bencher: divan::Bencher) {
    let newicks = posterior_newicks(WIDE);
    bencher.bench_local(|| snaps_from(&newicks));
}

#[divan::bench]
fn construct_posterior_wide_st(bencher: divan::Bencher) {
    let newicks = posterior_newicks(WIDE);
    let pool = single_thread_pool();
    bencher.bench_local(|| pool.install(|| snaps_from(&newicks)));
}

// --- pairwise benches: build snaps untimed, measure only the pairwise call --

#[divan::bench]
fn rf_posterior(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_posterior_sweep(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SWEEP));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_posterior_sweep_st(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SWEEP));
    let pool = single_thread_pool();
    bencher.bench_local(|| pool.install(|| snaps.pairwise_rf(None)));
}

#[divan::bench]
fn rf_diverse_sweep(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks(DIVERSE));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn rf_diverse_sweep_st(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks(DIVERSE));
    let pool = single_thread_pool();
    bencher.bench_local(|| pool.install(|| snaps.pairwise_rf(None)));
}

#[divan::bench]
fn rf_posterior_wide(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(WIDE));
    bencher.bench_local(|| snaps.pairwise_rf(None));
}

#[divan::bench]
fn wrf_posterior(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn wrf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn wrf_diverse_sweep(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks(DIVERSE));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn wrf_posterior_sweep(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SWEEP));
    bencher.bench_local(|| snaps.pairwise_wrf(None));
}

#[divan::bench]
fn kf_posterior_sweep(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SWEEP));
    bencher.bench_local(|| snaps.pairwise_kf(None));
}

#[divan::bench]
fn kf_posterior(bencher: divan::Bencher) {
    let snaps = snaps_from(&posterior_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_kf(None));
}

#[divan::bench]
fn kf_diverse(bencher: divan::Bencher) {
    let snaps = snaps_from(&diverse_newicks(SMALL));
    bencher.bench_local(|| snaps.pairwise_kf(None));
}
