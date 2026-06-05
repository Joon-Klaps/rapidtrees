use cpu_time::ProcessTime;
use memory_stats::memory_stats;
use rapidtrees::distances::{dispatch_kf, dispatch_rf, dispatch_wrf};
use rapidtrees::{Bitset, Snapshots};
use std::mem;
use std::time::{Duration, Instant};

// ---------------------------------------------------------------------------
// CLI argument parsing (no clap — must work without the `cli` feature)
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct BenchArgs {
    use_gpu: bool,
    metric: Metric,
    max_taxa: usize,
    max_trees: usize,
    iters: usize,
    swaps: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Metric {
    Rf,
    Wrf,
    Kf,
}

impl BenchArgs {
    fn parse() -> Self {
        let argv: Vec<String> = std::env::args().collect();
        let get = |flag: &str| {
            argv.iter()
                .position(|a| a == flag)
                .and_then(|i| argv.get(i + 1))
                .cloned()
        };

        let use_gpu = argv.contains(&"--gpu".to_string());

        let metric = match get("--metric").as_deref() {
            Some("wrf") | Some("weighted") => Metric::Wrf,
            Some("kf") => Metric::Kf,
            _ => Metric::Rf,
        };

        let max_taxa = get("--max-taxa")
            .and_then(|s| s.parse().ok())
            .unwrap_or(20_000);

        let max_trees = get("--max-trees")
            .and_then(|s| s.parse().ok())
            .unwrap_or(10_000);

        // Number of timed repetitions per cell. The median is reported and a
        // relative spread (stddev/mean) is shown so noisy cells are visible.
        let iters = get("--iters")
            .and_then(|s| s.parse().ok())
            .filter(|&i| i > 0)
            .unwrap_or(5);

        // Per-tree topology perturbations. 0 = every tree identical (degenerate:
        // all splits universal, RF kernel does no work). A few swaps give a
        // posterior-like set with real bipartition variation.
        let swaps = get("--swaps").and_then(|s| s.parse().ok()).unwrap_or(3);

        BenchArgs {
            use_gpu,
            metric,
            max_taxa,
            max_trees,
            iters,
            swaps,
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a balanced Newick tree whose leaves carry the labels in `order`:
/// position `i` of the balanced shape gets label `leaf_{order[i]}`.
///
/// All trees built from permutations of the same `order` share one leaf set
/// (`leaf_0..leaf_{n-1}`) but differ topologically, so a collection of them has
/// real bipartition variation — unlike `t` clones of a single tree, where every
/// split is universal and the RF kernel does zero work.
fn build_balanced_newick(order: &[usize]) -> String {
    fn recurse(order: &[usize]) -> String {
        if order.len() == 1 {
            return format!("leaf_{}", order[0]);
        }
        let mid = order.len() / 2;
        let left = recurse(&order[..mid]);
        let right = recurse(&order[mid..]);
        // Branch length keyed on the leading label so WRF/KF rows are non-trivial.
        let idx = order[0] % 97;
        format!(
            "({}:{:.4},{}:{:.4})",
            left,
            0.1 + idx as f64 * 0.05,
            right,
            0.1 + idx as f64 * 0.03
        )
    }
    format!("{};", recurse(order))
}

/// Deterministic small perturbation of `base`: apply `n_swaps` random leaf-position
/// swaps, seeded by `tree_idx` for reproducibility.
///
/// Each swap moves only the `O(log n)` bipartitions on the path between the two
/// leaves, so a few swaps yield a posterior-like set — mostly shared splits with a
/// handful unique per tree, the regime rapidtrees targets. `n_swaps == 0` reproduces
/// the old behaviour (every tree identical). Uses a tiny inline LCG to stay
/// dependency-free and reproducible.
fn perturbed_order(base: &[usize], tree_idx: usize, n_swaps: usize) -> Vec<usize> {
    let mut order = base.to_vec();
    let n = order.len();
    if n < 2 {
        return order;
    }
    let mut state = (tree_idx as u64)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(1);
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    for _ in 0..n_swaps {
        let a = next() % n;
        let b = next() % n;
        order.swap(a, b);
    }
    order
}

fn format_duration(d: Duration) -> String {
    let secs = d.as_secs_f64();
    if secs < 0.001 {
        format!("{:.2} µs", secs * 1_000_000.0)
    } else if secs < 1.0 {
        format!("{:.2} ms", secs * 1_000.0)
    } else if secs < 60.0 {
        format!("{:.2} s", secs)
    } else {
        format!("{:.2} min", secs / 60.0)
    }
}

fn format_size(bytes: usize) -> String {
    const KB: usize = 1024;
    const MB: usize = 1024 * 1024;
    const GB: usize = 1024 * 1024 * 1024;

    if bytes >= GB {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    } else if bytes >= MB {
        format!("{:.2} MB", bytes as f64 / MB as f64)
    } else if bytes >= KB {
        format!("{:.2} KB", bytes as f64 / KB as f64)
    } else {
        format!("{} B", bytes)
    }
}

fn format_count(count: u64) -> String {
    if count >= 1_000_000_000 {
        format!("{:.1}B", count as f64 / 1_000_000_000.0)
    } else if count >= 1_000_000 {
        format!("{:.1}M", count as f64 / 1_000_000.0)
    } else if count >= 1_000 {
        format!("{:.1}K", count as f64 / 1_000.0)
    } else {
        format!("{}", count)
    }
}

/// Estimate host RAM used by a `Snapshots` collection.
fn estimate_ram(snaps: &Snapshots) -> usize {
    let n_bip = snaps.bipartitions.len();
    let n_trees = snaps.len();
    let words = snaps.words_per_bitset;

    let struct_size = mem::size_of::<Snapshots>();
    let bip_vec_size =
        snaps.bipartitions.capacity() * (mem::size_of::<Bitset>() + words * mem::size_of::<u64>());
    let snap_size =
        n_trees * (24 + n_bip * mem::size_of::<u32>() + 24 + n_bip * mem::size_of::<f64>());
    let names_size = snaps
        .leaf_names
        .iter()
        .map(|n| mem::size_of::<String>() + n.len())
        .sum::<usize>();

    struct_size + bip_vec_size + snap_size + names_size
}

/// Estimate GPU VRAM for the pairwise kernel input + output buffers.
///
/// RF:      packed u32 rows  (n × ceil(n_bip/32) × 4 B)  + u32 kept (n × 4 B)
///          + output u32 matrix (n × n × 4 B)
/// WRF/KF:  f32 length rows  (n × n_bip × 4 B)
///          + output f32 matrix (n × n × 4 B)
fn estimate_vram(snaps: &Snapshots, metric: Metric) -> usize {
    let n = snaps.len();
    let n_bip = snaps.bipartitions.len();
    let words = snaps.words_per_bitset; // ceil(n_bip / 64) in u64, but GPU uses u32 words
    let gpu_words = words * 2; // u64 → two u32 words on GPU
    let output = n * n * 4;
    match metric {
        Metric::Rf => n * gpu_words * 4 + n * 4 + output,
        Metric::Wrf | Metric::Kf => n * n_bip * 4 + output,
    }
}

fn bar_log(value: f64, max: f64, width: usize) -> String {
    if value <= 0.0 {
        return "░".repeat(width);
    }
    let ratio = (value.ln() / max.ln()).clamp(0.0, 1.0);
    let filled = (ratio * width as f64).round() as usize;
    format!("{}{}", "█".repeat(filled), "░".repeat(width - filled))
}

fn mem_bar(bytes: usize, width: usize, max: usize) -> String {
    bar_log(bytes as f64, max as f64, width)
}

fn time_bar(d: Duration, width: usize, max: Duration) -> String {
    bar_log(d.as_secs_f64(), max.as_secs_f64(), width)
}

fn from_newicks_or_skip(newicks: &[String], rooted: bool, label: &str) -> Option<Snapshots> {
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();
    match Snapshots::from_newicks(&refs, rooted) {
        Ok(s) => Some(s),
        Err(e) => {
            eprintln!("Failed to parse {label}: {e}");
            None
        }
    }
}

fn run_metric(snaps: &Snapshots, metric: Metric, use_gpu: bool) {
    match metric {
        Metric::Rf => {
            let _ = dispatch_rf(snaps, None, use_gpu);
        }
        Metric::Wrf => {
            let _ = dispatch_wrf(snaps, None, use_gpu);
        }
        Metric::Kf => {
            let _ = dispatch_kf(snaps, None, use_gpu);
        }
    }
}

/// Median of a set of samples. Sorts `samples` in place.
fn median(samples: &mut [Duration]) -> Duration {
    samples.sort_unstable();
    let mid = samples.len() / 2;
    if samples.len() % 2 == 1 {
        samples[mid]
    } else {
        (samples[mid - 1] + samples[mid]) / 2
    }
}

/// Relative spread (coefficient of variation: stddev / mean) as a percentage.
/// A high value flags an unreliable cell (contention, thermal throttling, …).
fn rel_spread_pct(samples: &[Duration]) -> f64 {
    let n = samples.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let mean = samples.iter().map(Duration::as_secs_f64).sum::<f64>() / n;
    if mean <= 0.0 {
        return 0.0;
    }
    let var = samples
        .iter()
        .map(|d| {
            let x = d.as_secs_f64() - mean;
            x * x
        })
        .sum::<f64>()
        / n;
    var.sqrt() / mean * 100.0
}

/// Run `body`, sampling whole-process RSS in a background thread, and return
/// `(body_result, peak_rss_delta_over_baseline)`.
///
/// More reliable than a before/after RSS delta: the OS may reclaim or lazily
/// account pages, so a single after-reading can understate the true high-water
/// mark. We poll every millisecond and also take a final reading once `body`
/// returns so even sub-millisecond operations record their steady-state delta.
fn with_peak_rss<R>(body: impl FnOnce() -> R) -> (R, usize) {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

    let baseline = memory_stats().map(|m| m.physical_mem).unwrap_or(0);
    let stop = Arc::new(AtomicBool::new(false));
    let peak = Arc::new(AtomicUsize::new(baseline));

    let sampler = {
        let stop = Arc::clone(&stop);
        let peak = Arc::clone(&peak);
        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                if let Some(m) = memory_stats() {
                    peak.fetch_max(m.physical_mem, Ordering::Relaxed);
                }
                std::thread::sleep(Duration::from_millis(1));
            }
        })
    };

    let result = body();

    if let Some(m) = memory_stats() {
        peak.fetch_max(m.physical_mem, Ordering::Relaxed);
    }
    stop.store(true, Ordering::Relaxed);
    let _ = sampler.join();

    let peak_delta = peak.load(Ordering::Relaxed).saturating_sub(baseline);
    (result, peak_delta)
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args = BenchArgs::parse();

    let metric_name = match args.metric {
        Metric::Rf => "RF",
        Metric::Wrf => "Weighted RF",
        Metric::Kf => "KF",
    };
    let mode_name = if args.use_gpu { "GPU" } else { "CPU" };

    println!(
        "Mode: {mode_name} | Metric: {metric_name} | Iters: {}",
        args.iters
    );
    println!(
        "Max taxa: {} | Max trees: {} | Swaps/tree: {}",
        args.max_taxa, args.max_trees, args.swaps
    );
    println!(
        "Wall/CPU columns report the median over {} runs (±relative spread).",
        args.iters
    );
    if args.use_gpu {
        println!("Note: GPU falls back to CPU if no adapter found or n < 64.");
        // Self-report whether the GPU actually engages, so the .out file is
        // self-diagnosing instead of leaving us to infer it from timings.
        #[cfg(feature = "gpu")]
        match rapidtrees::gpu_adapter_label() {
            Some(label) => println!("GPU CONFIRMED: {label}"),
            None => println!("GPU NOT USED — no compatible adapter found; running on CPU"),
        }
        #[cfg(not(feature = "gpu"))]
        println!("GPU NOT USED — binary built without the `gpu` feature; running on CPU");
    }
    println!();

    let header = if args.use_gpu {
        "| Taxa (N) | Trees (T) | Combinations | Est. RAM | VRAM est / Host peak      | Wall Time | CPU Time |"
    } else {
        "| Taxa (N) | Trees (T) | Combinations | Est. RAM | Peak RAM (Δ)              | Wall Time | CPU Time |"
    };
    let sep = "|----------|-----------|--------------|----------|---------------------------|-----------|----------|";

    println!("{header}");
    println!("{sep}");

    let taxa_counts = [10, 100, 500, 1000, 2000, 5000, 10_000, 20_000]
        .into_iter()
        .filter(|&n| n <= args.max_taxa)
        .collect::<Vec<_>>();
    let tree_counts = [100, 1000, 5000]
        .into_iter()
        .filter(|&t| t <= args.max_trees)
        .collect::<Vec<_>>();

    const RAM_LIMIT: usize = 30 * 1024 * 1024 * 1024;
    const VRAM_LIMIT: usize = 16 * 1024 * 1024 * 1024; // 16 GB (conservative for P100/A100)
    const TIME_LIMIT: Duration = Duration::from_secs(1); // 1 s
    const BAR_WIDTH: usize = 10;
    const MAX_COMPARISONS: u64 = 200_000_000_000;

    for &n in &taxa_counts {
        // Base leaf order shared by every tree of this taxa count; each tree is a
        // perturbation of it, so the set has real (posterior-like) split variation.
        let base_order: Vec<usize> = (0..n).collect();
        let mk = |i: usize| build_balanced_newick(&perturbed_order(&base_order, i, args.swaps));

        for &t in &tree_counts {
            let total_comparisons = (t as u64) * (t as u64) / 2;
            let combs_str = format_count(total_comparisons);

            // --- size estimate from a small sample ---
            let subset_size = t.min(100);
            let subset_newicks: Vec<String> = (0..subset_size).map(&mk).collect();
            let Some(sample) = from_newicks_or_skip(&subset_newicks, false, "size-sample") else {
                continue;
            };

            let size_per_tree = estimate_ram(&sample) / subset_size.max(1);
            let total_est_ram = size_per_tree * t;
            let est_ram_str = format_size(total_est_ram);

            // GPU: check VRAM; CPU: check RAM
            if args.use_gpu {
                let est_vram =
                    estimate_vram(&sample, args.metric) / subset_size.max(1) * t + t * t * 4; // output matrix
                if est_vram > VRAM_LIMIT {
                    println!(
                        "| {:<8} | {:<9} | {:<12} | {:<8} | {:<25} | {:<9} | {:<8} |",
                        n,
                        t,
                        combs_str,
                        est_ram_str,
                        format!(">{}", format_size(VRAM_LIMIT)),
                        "Skipped",
                        "-"
                    );
                    continue;
                }
            } else if total_est_ram > RAM_LIMIT {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<8} | {:<25} | {:<9} | {:<8} |",
                    n, t, combs_str, est_ram_str, "Skipped (>30 GB)", "-", "-"
                );
                continue;
            }

            // --- build full snapshot, sampling peak RSS during the build ---
            let full_newicks: Vec<String> = (0..t).map(&mk).collect();

            let (full_snaps_opt, actual_ram) =
                with_peak_rss(|| from_newicks_or_skip(&full_newicks, false, "full-parse"));
            let Some(full_snaps) = full_snaps_opt else {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<8} | {:<25} | {:<9} | {:<8} |",
                    n, t, combs_str, est_ram_str, "Parse error", "-", "-"
                );
                continue;
            };

            if !args.use_gpu && actual_ram > RAM_LIMIT {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<8} | {:<25} | {:<9} | {:<8} |",
                    n,
                    t,
                    combs_str,
                    est_ram_str,
                    format!(">{}", format_size(RAM_LIMIT)),
                    "Skipped",
                    "-"
                );
                drop(full_snaps);
                continue;
            }

            // --- benchmark on a capped subset ---
            let max_subset = (MAX_COMPARISONS as f64).sqrt() as usize;
            let bench_size = t.min(max_subset);
            let bench_newicks: Vec<String> = (0..bench_size).map(&mk).collect();
            let bench_snaps =
                from_newicks_or_skip(&bench_newicks, false, "bench").expect("bench parse failed");

            // Warmup: discarded. The first call initialises the wgpu device on
            // GPU, and warms caches / branch predictors on CPU.
            run_metric(&bench_snaps, args.metric, args.use_gpu);

            // Time `iters` repetitions, sampling peak host RSS across all of them.
            let mut wall_samples: Vec<Duration> = Vec::with_capacity(args.iters);
            let mut cpu_samples: Vec<Duration> = Vec::with_capacity(args.iters);
            let (_, dispatch_host_rss) = with_peak_rss(|| {
                for _ in 0..args.iters {
                    let start_wall = Instant::now();
                    let start_cpu = ProcessTime::now();
                    run_metric(&bench_snaps, args.metric, args.use_gpu);
                    wall_samples.push(start_wall.elapsed());
                    cpu_samples.push(start_cpu.elapsed());
                }
            });

            let wall_duration = median(&mut wall_samples);
            let cpu_duration = median(&mut cpu_samples);
            let wall_spread = rel_spread_pct(&wall_samples);
            let cpu_spread = rel_spread_pct(&cpu_samples);

            let run_comparisons = (bench_size as f64) * (bench_size as f64);
            let ratio = (total_comparisons as f64) / run_comparisons;

            // Second column: VRAM estimate + host peak RSS (GPU) or peak RAM (CPU)
            let col2 = if args.use_gpu {
                let vram = estimate_vram(&full_snaps, args.metric);
                // Host peak reflects actual staging-buffer cost for the capped bench subset.
                format!("{} / {}", format_size(vram), format_size(dispatch_host_rss))
            } else {
                let ram_bar = mem_bar(actual_ram, BAR_WIDTH, RAM_LIMIT);
                format!("{} {}", ram_bar, format_size(actual_ram))
            };

            let est_wall = wall_duration.mul_f64(ratio);
            let wall_bar = time_bar(est_wall, BAR_WIDTH, TIME_LIMIT);
            let est_cpu = cpu_duration.mul_f64(ratio);
            let est_tag = if ratio > 1.01 { " (est)" } else { "" };

            let wall_str = format!(
                "{} {} ±{:.0}%{}",
                wall_bar,
                format_duration(est_wall),
                wall_spread,
                est_tag
            );
            let cpu_str = format!(
                "{} ±{:.0}%{}",
                format_duration(est_cpu),
                cpu_spread,
                est_tag
            );

            println!(
                "| {:<8} | {:<9} | {:<12} | {:<8} | {:<25} | {:<9} | {:<8} |",
                n, t, combs_str, est_ram_str, col2, wall_str, cpu_str
            );

            drop(full_snaps);
        }
    }
}
