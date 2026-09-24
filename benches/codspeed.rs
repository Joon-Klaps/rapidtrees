//! CodSpeed regression benches: snapshot construction and the pairwise
//! backends (RF, WRF, KF).
//!
//! Every bench is named after what it measures, and that name is the only
//! place the cell is written down:
//!
//! ```text
//! {build|rf|wrf|kf}_{T}trees_{N}taxa_q{P}[_1thread]
//! ```
//!
//! `build` times `Snapshots::from_newicks` (parse + intern); `rf`, `wrf` and
//! `kf` build the snapshots untimed and time only the pairwise call. `T` trees
//! over `N` taxa are generated at diversity `q = P / 100` (below), and
//! `_1thread` runs the measured call on a one-thread rayon pool, the algorithm
//! without the scheduler. [`benches!`] turns each listed name into a bench, and
//! [`Cell::from_name`] reads the cell back out of it, so a name cannot drift
//! from what it measures: a malformed one panics with the expected pattern.
//!
//! The trees come from the manuscript simulator's generator
//! (`simulate_trees.py --diversity`): one random base topology, and in every
//! tree each internal edge collapsed with probability `q` and re-resolved at
//! random. `q40` is posterior-like, with a mean normalised pairwise RF of 0.53
//! against a median of 0.54 over the real BEAST posteriors. `q05` is
//! near-identical. `q100` re-resolves every edge, which gives independent
//! random topologies: the adversarial case, where almost every split is unique
//! to one tree.
//!
//! Build/run locally with `cargo codspeed build && cargo codspeed run` (add
//! `-m memory` for the memory instrument), or just `cargo bench --bench codspeed`
//! for a quick wall-clock sanity check.

use rapidtrees::Snapshots;

fn main() {
    divan::main();
}

// --- the grid -----------------------------------------------------------------

/// One `#[divan::bench]` per name, each measuring the cell its name spells out.
macro_rules! benches {
    ($($name:ident),* $(,)?) => {
        $(
            #[divan::bench]
            fn $name(bencher: divan::Bencher) {
                Cell::from_name(stringify!($name)).bench(bencher);
            }
        )*
    };
}

benches! {
    // Construction: a small cell, and a wide one where the per-tree build dominates.
    build_100trees_500taxa_q40,
    build_100trees_500taxa_q100,
    build_50trees_4000taxa_q40,
    build_50trees_4000taxa_q40_1thread,

    // RF: small guards, the pair sweep on a posterior, independent and
    // near-identical trees, and a wide cell.
    rf_100trees_500taxa_q40,
    rf_100trees_500taxa_q100,
    rf_1500trees_500taxa_q40,
    rf_1500trees_500taxa_q40_1thread,
    rf_1500trees_500taxa_q05,
    rf_300trees_500taxa_q100,
    rf_300trees_500taxa_q100_1thread,
    rf_50trees_4000taxa_q40,

    // Weighted RF and KF: small guards and the pair sweep on a posterior and on
    // independent trees.
    wrf_100trees_500taxa_q40,
    wrf_100trees_500taxa_q100,
    wrf_1500trees_500taxa_q40,
    wrf_300trees_500taxa_q100,
    kf_100trees_500taxa_q40,
    kf_100trees_500taxa_q100,
    kf_1500trees_500taxa_q40,
    kf_300trees_500taxa_q100,
}

// --- cells --------------------------------------------------------------------

/// What a bench times.
#[derive(Clone, Copy)]
enum Metric {
    /// `Snapshots::from_newicks`: parse and intern.
    Build,
    Rf,
    Wrf,
    Kf,
}

/// One bench, as its name spells it out.
struct Cell {
    metric: Metric,
    trees: usize,
    taxa: usize,
    /// Share of internal edges each tree collapses and re-resolves.
    q: f64,
    one_thread: bool,
}

impl Cell {
    /// Read a cell from `{metric}_{T}trees_{N}taxa_q{P}[_1thread]`, with `P`
    /// in per cent.
    fn from_name(name: &str) -> Self {
        let mut parts = name.split('_');
        let metric = match parts.next() {
            Some("build") => Metric::Build,
            Some("rf") => Metric::Rf,
            Some("wrf") => Metric::Wrf,
            Some("kf") => Metric::Kf,
            _ => malformed(name),
        };
        let mut number = |prefix: &str, suffix: &str| -> usize {
            parts
                .next()
                .and_then(|part| {
                    part.strip_prefix(prefix)?
                        .strip_suffix(suffix)?
                        .parse()
                        .ok()
                })
                .unwrap_or_else(|| malformed(name))
        };
        let (trees, taxa, percent) = (number("", "trees"), number("", "taxa"), number("q", ""));
        let one_thread = match parts.next() {
            None => false,
            Some("1thread") => true,
            Some(_) => malformed(name),
        };
        if parts.next().is_some() || percent > 100 {
            malformed(name);
        }
        Self {
            metric,
            trees,
            taxa,
            q: percent as f64 / 100.0,
            one_thread,
        }
    }

    /// Generate the trees untimed, then time what the metric names.
    fn bench(&self, bencher: divan::Bencher) {
        let newicks = tree_set(self.taxa, self.trees, self.q);
        match self.metric {
            Metric::Build => measure(bencher, self.one_thread, || snaps_from(&newicks)),
            Metric::Rf => {
                let snaps = snaps_from(&newicks);
                measure(bencher, self.one_thread, || snaps.pairwise_rf(None));
            }
            Metric::Wrf => {
                let snaps = snaps_from(&newicks);
                measure(bencher, self.one_thread, || snaps.pairwise_wrf(None));
            }
            Metric::Kf => {
                let snaps = snaps_from(&newicks);
                measure(bencher, self.one_thread, || snaps.pairwise_kf(None));
            }
        }
    }
}

fn malformed(name: &str) -> ! {
    panic!(
        "bench `{name}` must be named {{build|rf|wrf|kf}}_<T>trees_<N>taxa_q<P>[_1thread], P <= 100"
    )
}

/// Time `work`, on a one-thread rayon pool when `one_thread`.
fn measure<T: Send>(bencher: divan::Bencher, one_thread: bool, work: impl Fn() -> T + Sync) {
    if one_thread {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("thread pool");
        bencher.bench_local(|| pool.install(&work));
    } else {
        bencher.bench_local(&work);
    }
}

fn snaps_from(newicks: &[String]) -> Snapshots {
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
    Snapshots::from_newicks(&refs, false).expect("parse failed")
}

// --- tree generation (self-contained) -----------------------------------------

/// Mean branch length; lengths are exponential, as in the simulator.
const BRANCH_SCALE: f64 = 0.1;

/// `trees` trees over `taxa` leaves: one random base topology, re-resolved per
/// tree at diversity `q`. Seeded, so every run benchmarks the same trees.
fn tree_set(taxa: usize, trees: usize, q: f64) -> Vec<String> {
    let mut state = 0x5EED_0000_0000_0040;
    let base = resolve((0..taxa as u32).map(Node::Leaf).collect(), &mut state);
    (0..trees)
        .map(|_| {
            let tree = resolve(components(&base, q, &mut state), &mut state);
            let mut newick = String::new();
            write_newick(&tree, &mut state, &mut newick);
            newick.push(';');
            newick
        })
        .collect()
}

/// Small LCG so the bench needs no `rand` dependency.
fn next_rand(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

/// Uniform in `[0, 1)` from the LCG's 31 bits.
fn uniform(state: &mut u64) -> f64 {
    next_rand(state) as f64 / (1u64 << 31) as f64
}

/// A rooted tree over leaves `0..n`, with any arity while it is being built.
enum Node {
    Leaf(u32),
    Inner(Vec<Node>),
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
