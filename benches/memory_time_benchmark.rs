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

        BenchArgs {
            use_gpu,
            metric,
            max_taxa,
            max_trees,
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn generate_balanced_newick(start_index: usize, num_leaves: usize) -> String {
    if num_leaves == 1 {
        return format!("leaf_{}", start_index);
    }
    let left_count = num_leaves / 2;
    let right_count = num_leaves - left_count;

    let left = generate_balanced_newick(start_index, left_count);
    let right = generate_balanced_newick(start_index + left_count, right_count);

    // Use varied branch lengths so WRF/KF benchmarks are representative.
    let idx = start_index % 97;
    format!(
        "({}:{:.4},{}:{:.4})",
        left,
        0.1 + idx as f64 * 0.05,
        right,
        0.1 + idx as f64 * 0.03
    )
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

    println!("Mode: {mode_name} | Metric: {metric_name}");
    println!(
        "Max taxa: {} | Max trees: {}",
        args.max_taxa, args.max_trees
    );
    if args.use_gpu {
        println!("Note: GPU falls back to CPU if no adapter found or n < 64.");
    }
    println!();

    let header = if args.use_gpu {
        "| Taxa (N) | Trees (T) | Combinations | Est. RAM | VRAM (input+out) | Wall Time | CPU Time |"
    } else {
        "| Taxa (N) | Trees (T) | Combinations | Est. RAM | Actual RAM       | Wall Time | CPU Time |"
    };
    let sep = if args.use_gpu {
        "|----------|-----------|--------------|----------|------------------|-----------|----------|"
    } else {
        "|----------|-----------|--------------|----------|------------------|-----------|----------|"
    };

    println!("{header}");
    println!("{sep}");

    let taxa_counts = [10, 100, 500, 1000, 2000, 5000, 10_000, 20_000]
        .into_iter()
        .filter(|&n| n <= args.max_taxa)
        .collect::<Vec<_>>();
    let tree_counts = [100, 1000, 10_000]
        .into_iter()
        .filter(|&t| t <= args.max_trees)
        .collect::<Vec<_>>();

    const RAM_LIMIT: usize = 30 * 1024 * 1024 * 1024;
    const VRAM_LIMIT: usize = 16 * 1024 * 1024 * 1024; // 16 GB (conservative for P100/A100)
    const TIME_LIMIT: Duration = Duration::from_secs(60 * 60);
    const BAR_WIDTH: usize = 10;
    const MAX_COMPARISONS: u64 = 200_000_000_000;

    for &n in &taxa_counts {
        let newick = format!("{};", generate_balanced_newick(0, n));

        for &t in &tree_counts {
            let total_comparisons = (t as u64) * (t as u64) / 2;
            let combs_str = format_count(total_comparisons);

            // --- size estimate from a small sample ---
            let subset_size = t.min(100);
            let subset_newicks: Vec<String> = (0..subset_size).map(|_| newick.clone()).collect();
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
                        "| {:<8} | {:<9} | {:<12} | {:<8} | {:<16} | {:<9} | {:<8} |",
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
                    "| {:<8} | {:<9} | {:<12} | {:<8} | {:<16} | {:<9} | {:<8} |",
                    n, t, combs_str, est_ram_str, "Skipped (>30 GB)", "-", "-"
                );
                continue;
            }

            // --- build full snapshot and measure RSS delta ---
            let full_newicks: Vec<String> = (0..t).map(|_| newick.clone()).collect();

            let rss_before = memory_stats().map(|m| m.physical_mem).unwrap_or(0);
            let Some(full_snaps) = from_newicks_or_skip(&full_newicks, false, "full-parse") else {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<8} | {:<16} | {:<9} | {:<8} |",
                    n, t, combs_str, est_ram_str, "Parse error", "-", "-"
                );
                continue;
            };
            let rss_after = memory_stats().map(|m| m.physical_mem).unwrap_or(0);
            // RSS can temporarily dip due to GC/reclaim; saturating_sub gives 0 instead of panic.
            let actual_ram = rss_after.saturating_sub(rss_before);

            if !args.use_gpu && actual_ram > RAM_LIMIT {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<8} | {:<16} | {:<9} | {:<8} |",
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

            // Second column: VRAM estimate (GPU) or actual RAM (CPU)
            let col2 = if args.use_gpu {
                let vram = estimate_vram(&full_snaps, args.metric);
                format!(
                    "{} {}",
                    mem_bar(vram, BAR_WIDTH, VRAM_LIMIT),
                    format_size(vram)
                )
            } else {
                let ram_bar = mem_bar(actual_ram, BAR_WIDTH, RAM_LIMIT);
                format!("{} {}", ram_bar, format_size(actual_ram))
            };

            // --- benchmark on a capped subset ---
            let max_subset = (MAX_COMPARISONS as f64).sqrt() as usize;
            let bench_size = t.min(max_subset);
            let bench_newicks: Vec<String> = (0..bench_size).map(|_| newick.clone()).collect();
            let bench_snaps =
                from_newicks_or_skip(&bench_newicks, false, "bench").expect("bench parse failed");

            // Warmup (important for GPU: first call initialises the wgpu device)
            if args.use_gpu {
                run_metric(&bench_snaps, args.metric, true);
            }

            let start_wall = Instant::now();
            let start_cpu = ProcessTime::now();
            run_metric(&bench_snaps, args.metric, args.use_gpu);
            let wall_duration = start_wall.elapsed();
            let cpu_duration = start_cpu.elapsed();

            let run_comparisons = (bench_size as f64) * (bench_size as f64);
            let ratio = (total_comparisons as f64) / run_comparisons;

            let est_wall = wall_duration.mul_f64(ratio);
            let wall_bar = time_bar(est_wall, BAR_WIDTH, TIME_LIMIT);
            let est_cpu = cpu_duration.mul_f64(ratio);

            let wall_str = if ratio > 1.01 {
                format!("{} {} (est)", wall_bar, format_duration(est_wall))
            } else {
                format!("{} {}", wall_bar, format_duration(est_wall))
            };
            let cpu_str = if ratio > 1.01 {
                format!("{} (est)", format_duration(est_cpu))
            } else {
                format_duration(est_cpu)
            };

            println!(
                "| {:<8} | {:<9} | {:<12} | {:<8} | {:<16} | {:<9} | {:<8} |",
                n, t, combs_str, est_ram_str, col2, wall_str, cpu_str
            );

            drop(full_snaps);
        }
    }
}
