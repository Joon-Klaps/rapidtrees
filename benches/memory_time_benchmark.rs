use cpu_time::ProcessTime;
use memory_stats::memory_stats;
use phylotree::tree::Tree as PhyloTree;
use rapidtrees::distances::rf_from_snapshots;
use rapidtrees::snapshot::TreeSnapshot;
use rayon::prelude::*;
use std::mem;
use std::time::{Duration, Instant};

fn generate_balanced_newick(start_index: usize, num_leaves: usize) -> String {
    if num_leaves == 1 {
        return format!("leaf_{}", start_index);
    }
    let left_count = num_leaves / 2;
    let right_count = num_leaves - left_count;

    let left = generate_balanced_newick(start_index, left_count);
    let right = generate_balanced_newick(start_index + left_count, right_count);

    format!("({}:0.1,{}:0.1)", left, right)
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

fn estimate_size(snap: &TreeSnapshot) -> usize {
    let mut size = mem::size_of::<TreeSnapshot>();

    // parts: Vec<Bitset>
    size += snap.parts.capacity() * mem::size_of::<rapidtrees::bitset::Bitset>();
    for part in &snap.parts {
        // Bitset is a wrapper around Vec<u64>
        size += part.0.capacity() * mem::size_of::<u64>();
    }

    // lengths: Vec<f64>
    size += snap.lengths.capacity() * mem::size_of::<f64>();

    // root_children: Vec<Bitset>
    size += snap.root_children.capacity() * mem::size_of::<rapidtrees::bitset::Bitset>();
    for part in &snap.root_children {
        size += part.0.capacity() * mem::size_of::<u64>();
    }

    size
}

fn main() {
    println!(
        "| Taxa (N) | Trees (T) | Combinations | Est. Memory | Actual Memory | Wall Time | CPU Time |"
    );
    println!(
        "|----------|-----------|--------------|-------------|---------------|-----------|----------|"
    );

    let taxa_counts = [10, 100, 500, 1000, 2000, 5000];
    let tree_counts = [100, 1000, 10_000, 100_000];
    const MEMORY_LIMIT: usize = 30 * 1024 * 1024 * 1024; // 30 GB
    const MAX_COMPARISONS: usize = 200_000_000;

    for &n in &taxa_counts {
        let newick = format!("{};", generate_balanced_newick(0, n));
        let tree = PhyloTree::from_newick(&newick).expect("Failed to parse generated tree");
        let snap = TreeSnapshot::from_tree(&tree).expect("Failed to create snapshot");
        let size_per_tree = estimate_size(&snap);

        for &t in &tree_counts {
            let total_est_size = size_per_tree * t;
            let total_comparisons = (t as u64) * (t as u64);
            let combs_str = format_count(total_comparisons);
            let est_mem_str = format_size(total_est_size);

            if total_est_size > MEMORY_LIMIT {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<11} | {:<13} | {:<9} | {:<8} |",
                    n, t, combs_str, est_mem_str, "Skipped (>30GB)", "-", "-"
                );
                continue;
            }

            // Measure Memory
            let start_mem = memory_stats().map(|s| s.physical_mem).unwrap_or(0);

            // Allocate trees
            let trees: Vec<TreeSnapshot> = (0..t).map(|_| snap.clone()).collect();

            let end_mem = memory_stats().map(|s| s.physical_mem).unwrap_or(0);
            let mem_diff = end_mem.saturating_sub(start_mem);
            let actual_mem_str = format_size(mem_diff);

            if mem_diff > MEMORY_LIMIT {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<11} | {:<13} | {:<9} | {:<8} |",
                    n,
                    t,
                    combs_str,
                    est_mem_str,
                    format!("> {}", format_size(MEMORY_LIMIT)),
                    "Skipped",
                    "-"
                );
                continue;
            }

            // Benchmark Time
            // Determine subset size to keep comparisons under MAX_COMPARISONS
            // We want subset_size * subset_size <= MAX_COMPARISONS
            let max_subset = (MAX_COMPARISONS as f64).sqrt() as usize;
            let subset_size = if t > max_subset { max_subset } else { t };

            let slice = &trees[0..subset_size];

            // Pre-allocate matrix (like in main.rs)
            let mut mat = vec![vec![0.0f64; subset_size]; subset_size];

            let start_wall = Instant::now();
            let start_cpu = ProcessTime::now();

            // Parallel execution using Rayon (matching main.rs logic)
            mat.par_iter_mut().enumerate().for_each(|(i, row)| {
                for j in (i + 1)..subset_size {
                    let dist = rf_from_snapshots(&slice[i], &slice[j]);
                    row[j] = dist as f64;
                }
            });

            let wall_duration = start_wall.elapsed();
            let cpu_duration = start_cpu.elapsed();

            // Extrapolate
            let run_comparisons = (subset_size as f64) * (subset_size as f64);
            let ratio = (total_comparisons as f64) / run_comparisons;

            let est_wall = wall_duration.mul_f64(ratio);
            let est_cpu = cpu_duration.mul_f64(ratio);

            let wall_str = if ratio > 1.01 {
                format!("{} (est)", format_duration(est_wall))
            } else {
                format_duration(est_wall)
            };
            let cpu_str = if ratio > 1.01 {
                format!("{} (est)", format_duration(est_cpu))
            } else {
                format_duration(est_cpu)
            };

            println!(
                "| {:<8} | {:<9} | {:<12} | {:<11} | {:<13} | {:<9} | {:<8} |",
                n, t, combs_str, est_mem_str, actual_mem_str, wall_str, cpu_str
            );
        }
    }
}
