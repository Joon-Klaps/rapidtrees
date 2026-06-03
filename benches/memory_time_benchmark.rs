use cpu_time::ProcessTime;
use rapidtrees::{Bitset, Snapshots};
use std::alloc::{GlobalAlloc, Layout, System};
use std::mem;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

// ---------------------------------------------------------------------------
// Tracking allocator
// ---------------------------------------------------------------------------

static ALLOCATED: AtomicUsize = AtomicUsize::new(0);

struct TrackingAllocator;

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        unsafe {
            let ptr = System.alloc(layout);
            if !ptr.is_null() {
                ALLOCATED.fetch_add(layout.size(), Ordering::Relaxed);
            }
            ptr
        }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe {
            System.dealloc(ptr, layout);
            ALLOCATED.fetch_sub(layout.size(), Ordering::Relaxed);
        }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        unsafe {
            let ptr = System.alloc_zeroed(layout);
            if !ptr.is_null() {
                ALLOCATED.fetch_add(layout.size(), Ordering::Relaxed);
            }
            ptr
        }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        unsafe {
            let new_ptr = System.realloc(ptr, layout, new_size);
            if !new_ptr.is_null() {
                // Subtract old size, add new size
                ALLOCATED.fetch_sub(layout.size(), Ordering::Relaxed);
                ALLOCATED.fetch_add(new_size, Ordering::Relaxed);
            }
            new_ptr
        }
    }
}

#[global_allocator]
static GLOBAL: TrackingAllocator = TrackingAllocator;

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

/// Estimate heap memory used by a `Snapshots` collection.
///
/// Since the internal `Vec<InternSnap>` is not public, the per-tree split storage
/// is approximated: for each tree we assume `n_bip` split IDs (u32) and `n_bip`
/// branch lengths (f64). This is exact when all trees share the same topology
/// (as in this benchmark).
fn estimate_size(snaps: &Snapshots) -> usize {
    let n_bip = snaps.bipartitions.len();
    let n_trees = snaps.len();
    let words = snaps.words_per_bitset;

    // Struct overhead (Vec headers, usize fields)
    let struct_size = mem::size_of::<Snapshots>();

    // bipartitions: Vec<Bitset>, each Bitset has a Vec<u64> on the heap
    let bip_vec_size =
        snaps.bipartitions.capacity() * (mem::size_of::<Bitset>() + words * mem::size_of::<u64>());

    // snapshots (Vec<InternSnap>): per tree = split_ids Vec<u32> + lengths Vec<f64>
    // Vec header (24 bytes) + actual data for each
    let snap_size =
        n_trees * (24 + n_bip * mem::size_of::<u32>() + 24 + n_bip * mem::size_of::<f64>());

    // leaf_names: Vec<String>
    let names_size = snaps
        .leaf_names
        .iter()
        .map(|n| mem::size_of::<String>() + n.len())
        .sum::<usize>();

    struct_size + bip_vec_size + snap_size + names_size
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

fn main() {
    println!(
        "| Taxa (N) | Trees (T) | Combinations | Est. Memory | Actual Memory | Wall Time | CPU Time |"
    );
    println!(
        "|----------|-----------|--------------|-------------|---------------|-----------|----------|"
    );

    let taxa_counts = [10, 100, 500, 1000, 2000, 5000, 10_000, 20_000];
    let tree_counts = [100, 1000, 10_000];
    const MEMORY_LIMIT: usize = 30 * 1024 * 1024 * 1024; // 30 GB
    const TIME_LIMIT: Duration = Duration::from_secs(60 * 60); // 1 hour
    const BAR_WIDTH: usize = 10;
    const MAX_COMPARISONS: u64 = 200_000_000_000;

    for &n in &taxa_counts {
        let newick = format!("{};", generate_balanced_newick(0, n));

        for &t in &tree_counts {
            // --- size estimate from a small sample ---
            let subset_size = t.min(100);
            let subset_newicks: Vec<String> = (0..subset_size).map(|_| newick.clone()).collect();

            let Some(interned) = from_newicks_or_skip(&subset_newicks, false, "size-sample") else {
                continue;
            };

            let size_per_tree = estimate_size(&interned) / subset_size.max(1);
            let total_est_size = size_per_tree * t;
            let total_comparisons = (t as u64) * (t as u64) / 2;
            let combs_str = format_count(total_comparisons);
            let est_mem_str = format_size(total_est_size);

            if total_est_size > MEMORY_LIMIT {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<11} | {:<13} | {:<9} | {:<8} |",
                    n, t, combs_str, est_mem_str, "Skipped (>30GB)", "-", "-"
                );
                continue;
            }

            // --- measure actual allocation with tracking allocator ---
            let full_newicks: Vec<String> = (0..t).map(|_| newick.clone()).collect();

            let before = ALLOCATED.load(Ordering::Relaxed);
            let Some(full_interned) = from_newicks_or_skip(&full_newicks, false, "full-parse")
            else {
                println!(
                    "| {:<8} | {:<9} | {:<12} | {:<11} | {:<13} | {:<9} | {:<8} |",
                    n, t, combs_str, est_mem_str, "Parse error", "-", "-"
                );
                continue;
            };
            let after = ALLOCATED.load(Ordering::Relaxed);
            // after >= before guaranteed: we only just allocated, nothing freed yet
            let actual_mem = after.saturating_sub(before);
            let actual_mem_bar = mem_bar(actual_mem, BAR_WIDTH, MEMORY_LIMIT);
            let actual_mem_str = format!("{} {}", actual_mem_bar, format_size(actual_mem));

            if actual_mem > MEMORY_LIMIT {
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
                drop(full_interned);
                continue;
            }

            // --- benchmark pairwise distances on a capped subset ---
            let max_subset = (MAX_COMPARISONS as f64).sqrt() as usize;
            let bench_size = t.min(max_subset);
            let bench_newicks: Vec<String> = (0..bench_size).map(|_| newick.clone()).collect();
            let bench_snaps =
                from_newicks_or_skip(&bench_newicks, false, "bench").expect("bench parse failed");

            let start_wall = Instant::now();
            let start_cpu = ProcessTime::now();
            let _mat = bench_snaps.pairwise_rf();
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
                "| {:<8} | {:<9} | {:<12} | {:<11} | {:<13} | {:<9} | {:<8} |",
                n, t, combs_str, est_mem_str, actual_mem_str, wall_str, cpu_str
            );

            drop(full_interned);
        }
    }
}
