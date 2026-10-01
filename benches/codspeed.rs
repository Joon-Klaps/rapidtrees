//! CodSpeed regression benches: loading a tree file, snapshot construction and
//! the pairwise backends (RF, WRF, KF).
//!
//! Every bench is named after what it measures, and that name is the only
//! place the cell is written down:
//!
//! ```text
//! {load|build|rf|wrf|kf}_{T}trees_{N}taxa_q{P}[_annotated][_1thread]
//! ```
//!
//! `load` times `io::load_beast_trees` on a NEXUS file written untimed
//! beforehand, as the CLI loads it for RF: reading the file, parsing and
//! interning. It is the only bench that sees what happens to the file's text,
//! so it is the one the memory instrument needs for changes to the reader.
//! `_annotated` gives every node BEAST-style `[&...]` annotations, which make
//! up most of the bytes of a real BEAST file. `build` times
//! `Snapshots::from_newicks` (parse + intern) on strings already in memory;
//! `rf`, `wrf` and `kf` build the snapshots untimed and time only the pairwise
//! call. `T` trees over `N` taxa are generated at diversity `q = P / 100`
//! (below), and `_1thread` runs the measured call on a one-thread rayon pool,
//! the algorithm without the scheduler. [`benches!`] turns each listed name
//! into a bench, and [`Cell::from_name`] reads the cell back out of it, so a
//! name cannot drift from what it measures: a malformed one panics with the
//! expected pattern.
//!
//! `load_real::hiv_4x21trees_162taxa` loads the four real BEAST posteriors in
//! `tests/data` (git LFS), annotations and TRANSLATE blocks included.
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

use rapidtrees::io::load_beast_trees;
use rapidtrees::{Retain, Snapshots};
use std::fmt::Write as _;
use std::path::Path;

fn main() {
    divan::main();
}

// --- the grid -----------------------------------------------------------------

/// One `#[divan::bench]` per name, each measuring the cell its name spells out,
/// in a module per group. The group is the bench path's middle segment
/// (`codspeed::rf::rf_100trees_500taxa_q40`), which is what CI shards on.
macro_rules! benches {
    ($($group:ident : [$($name:ident),* $(,)?]),* $(,)?) => {
        $(
            mod $group {
                use super::*;

                $(
                    #[divan::bench]
                    fn $name(bencher: divan::Bencher) {
                        Cell::from_name(stringify!($name)).bench(bencher);
                    }
                )*
            }
        )*
    };
}

benches! {
    // Loading a file: a posterior-sized set with and without BEAST annotations,
    // and a wide one, where the file's text dwarfs what is kept of it.
    load: [
        load_1000trees_500taxa_q40,
        load_1000trees_500taxa_q40_annotated,
        load_20trees_50000taxa_q40,
        load_20trees_50000taxa_q40_1thread,
    ],

    // Construction: a small cell, and wide ones where the per-tree build dominates.
    build: [
        build_100trees_500taxa_q40,
        build_100trees_500taxa_q100,
        build_50trees_4000taxa_q40,
        build_50trees_4000taxa_q40_1thread,
        build_20trees_50000taxa_q40,
    ],

    // RF: small guards, the pair sweep on a posterior, independent and
    // near-identical trees, and a wide cell.
    rf: [
        rf_100trees_500taxa_q40,
        rf_100trees_500taxa_q100,
        rf_1500trees_500taxa_q40,
        rf_1500trees_500taxa_q40_1thread,
        rf_1500trees_500taxa_q05,
        rf_300trees_500taxa_q100,
        rf_300trees_500taxa_q100_1thread,
        rf_50trees_4000taxa_q40,
        rf_3000trees_500taxa_q40,
        rf_20trees_50000taxa_q40,
    ],

    // Weighted RF and KF: small guards and the pair sweep on a posterior and on
    // independent trees.
    weighted: [
        wrf_100trees_500taxa_q40,
        wrf_100trees_500taxa_q100,
        wrf_1500trees_500taxa_q40,
        wrf_300trees_500taxa_q100,
        kf_100trees_500taxa_q40,
        kf_100trees_500taxa_q100,
        kf_1500trees_500taxa_q40,
        kf_300trees_500taxa_q100,
        wrf_20trees_50000taxa_q40,
    ],
}

/// The four real BEAST posteriors in `tests/data` (162 taxa, 21 trees each,
/// `[&rate=...]` annotations and TRANSLATE blocks), loaded as the CLI loads
/// them for RF.
///
/// Its module name starts with `load` so that it lands in the `load` CI shard.
mod load_real {
    use super::*;

    #[divan::bench]
    fn hiv_4x21trees_162taxa(bencher: divan::Bencher) {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let paths: Vec<_> = (1..=4).map(|k| dir.join(format!("hiv{k}.trees"))).collect();
        for path in &paths {
            // A checkout without git LFS has pointer files here, which hold no trees.
            assert!(
                load(path).len() == 21,
                "{path:?} holds no trees: fetch it with `git lfs pull`"
            );
        }
        bencher.bench_local(|| paths.iter().map(|path| load(path)).collect::<Vec<_>>());
    }
}

// --- cells --------------------------------------------------------------------

/// What a bench times.
#[derive(Clone, Copy)]
enum Metric {
    /// `io::load_beast_trees` on a NEXUS file: read, parse and intern.
    Load,
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
    /// `load` only: annotate every node as BEAST does.
    annotated: bool,
    one_thread: bool,
}

impl Cell {
    /// Read a cell from `{metric}_{T}trees_{N}taxa_q{P}[_annotated][_1thread]`,
    /// with `P` in per cent.
    fn from_name(name: &str) -> Self {
        let mut parts = name.split('_');
        let metric = match parts.next() {
            Some("load") => Metric::Load,
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
        let mut suffix = parts.next();
        let annotated = suffix == Some("annotated");
        if annotated {
            suffix = parts.next();
        }
        let one_thread = match suffix {
            None => false,
            Some("1thread") => true,
            Some(_) => malformed(name),
        };
        let loads = matches!(metric, Metric::Load);
        if parts.next().is_some() || percent > 100 || (annotated && !loads) {
            malformed(name);
        }
        Self {
            metric,
            trees,
            taxa,
            q: percent as f64 / 100.0,
            annotated,
            one_thread,
        }
    }

    /// Generate the trees untimed, then time what the metric names.
    fn bench(&self, bencher: divan::Bencher) {
        if let Metric::Load = self.metric {
            let dir = tempfile::tempdir().expect("temporary directory");
            let path = dir.path().join("trees.trees");
            let nexus = nexus_file(self.taxa, self.trees, self.q, self.annotated);
            std::fs::write(&path, nexus).expect("write the tree file");
            return measure(bencher, self.one_thread, || load(&path));
        }
        let newicks = tree_set(self.taxa, self.trees, self.q);
        match self.metric {
            Metric::Load => unreachable!("handled above"),
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
        "bench `{name}` must be named {{load|build|rf|wrf|kf}}_<T>trees_<N>taxa_q<P>[_annotated][_1thread], P <= 100, `_annotated` for `load` only"
    )
}

/// Load a tree file as the CLI does for RF: TRANSLATE applied, no burn-in,
/// unrooted, no branch lengths or bipartitions kept.
fn load(path: &Path) -> Snapshots {
    let (_, snaps) = load_beast_trees(path, 0, 0, true, false, Retain::for_distances(false));
    snaps
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

/// A BEAST-shaped NEXUS file: `T` trees over `N` taxa at diversity `q`, with
/// a TRANSLATE block, numeric leaf labels and `tree STATE_<k> = [&R] ...`
/// lines, and with `annotated`, `[&rate=...,height=...,posterior=...]` on every
/// node as BEAST writes them. Its own seed, so it leaves the trees of the other
/// benches as they were.
fn nexus_file(taxa: usize, trees: usize, q: f64, annotated: bool) -> String {
    let mut state = 0x5EED_0000_0000_0FAC;
    let base = resolve((0..taxa as u32).map(Node::Leaf).collect(), &mut state);
    let mut out = String::from("#NEXUS\n\nBegin trees;\n\tTranslate\n");
    for leaf in 0..taxa {
        let sep = if leaf + 1 < taxa { "," } else { "" };
        writeln!(out, "\t\t{} leaf_{leaf}{sep}", leaf + 1).expect("write to a String");
    }
    out.push_str("\t\t;\n");
    for k in 0..trees {
        let tree = resolve(components(&base, q, &mut state), &mut state);
        write!(out, "tree STATE_{} = [&R] ", k * 1000).expect("write to a String");
        write_beast(&tree, annotated, &mut state, &mut out);
        out.push_str(";\n");
    }
    out.push_str("End;\n");
    out
}

/// BEAST-style Newick for `node`: 1-based numeric leaf labels, exponential
/// branch lengths, and with `annotated`, an `[&...]` comment on every node
/// before its length.
fn write_beast(node: &Node, annotated: bool, state: &mut u64, out: &mut String) {
    match node {
        Node::Leaf(leaf) => write!(out, "{}", leaf + 1).expect("write to a String"),
        Node::Inner(children) => {
            out.push('(');
            for (k, child) in children.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                write_beast(child, annotated, state, out);
                if annotated {
                    let (rate, height, posterior) =
                        (uniform(state), uniform(state) * 50.0, uniform(state));
                    write!(
                        out,
                        "[&rate={rate:.10},height={height:.10},posterior={posterior:.10}]"
                    )
                    .expect("write to a String");
                }
                let length = -(1.0 - uniform(state)).ln() * BRANCH_SCALE;
                write!(out, ":{length:.6}").expect("write to a String");
            }
            out.push(')');
        }
    }
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
