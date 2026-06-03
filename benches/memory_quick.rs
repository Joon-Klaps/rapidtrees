//! Quick memory check for the interning / branch-length changes.
//!
//! Unlike `memory_time_benchmark`, this does **no** O(N²) pairwise loop — it only
//! builds a `Snapshots` collection and reports its retained heap footprint, so it
//! finishes in a couple of seconds even for thousands of trees.
//!
//! It quantifies the two memory wins in this branch:
//!   1. removal of the `bipartition_index` (a second full copy of every unique split)
//!   2. dropping per-tree branch `lengths` on RF-only paths (+ releasing the
//!      bipartition table on the non-export RF path)
//!
//! Run with defaults, or pass `tips trees [same|distinct]`:
//!   cargo bench --bench memory_quick
//!   cargo bench --bench memory_quick -- 2000 1000 distinct
//!   cargo bench --bench memory_quick -- 2000 1000 same

use rapidtrees::{Bitset, Snapshots};
use std::alloc::{GlobalAlloc, Layout, System};
use std::mem;
use std::sync::atomic::{AtomicUsize, Ordering};

// ---------------------------------------------------------------------------
// Tracking allocator (net retained + peak)
// ---------------------------------------------------------------------------

static CURRENT: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

struct TrackingAllocator;

fn bump(delta: usize) {
    let now = CURRENT.fetch_add(delta, Ordering::Relaxed) + delta;
    PEAK.fetch_max(now, Ordering::Relaxed);
}

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            bump(layout.size());
        }
        ptr
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
        CURRENT.fetch_sub(layout.size(), Ordering::Relaxed);
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if !ptr.is_null() {
            bump(layout.size());
        }
        ptr
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new_ptr = unsafe { System.realloc(ptr, layout, new_size) };
        if !new_ptr.is_null() {
            CURRENT.fetch_sub(layout.size(), Ordering::Relaxed);
            bump(new_size);
        }
        new_ptr
    }
}

#[global_allocator]
static GLOBAL: TrackingAllocator = TrackingAllocator;

// ---------------------------------------------------------------------------
// Tree generation
// ---------------------------------------------------------------------------

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

/// Deterministic Fisher–Yates shuffle driven by a small LCG (no `rand` dep).
fn shuffled_labels(n: usize, seed: u64) -> Vec<u32> {
    let mut v: Vec<u32> = (0..n as u32).collect();
    let mut s = seed ^ 0x9E37_79B9_7F4A_7C15;
    for i in (1..n).rev() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let j = ((s >> 33) as usize) % (i + 1);
        v.swap(i, j);
    }
    v
}

fn make_newick(tips: usize, seed: u64, distinct: bool) -> String {
    let labels = if distinct {
        shuffled_labels(tips, seed)
    } else {
        (0..tips as u32).collect()
    };
    let mut s = String::new();
    build_balanced(&labels, &mut s);
    s.push(';');
    s
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

fn human(bytes: usize) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = 1024.0 * 1024.0;
    const GB: f64 = 1024.0 * 1024.0 * 1024.0;
    let b = bytes as f64;
    if b >= GB {
        format!("{:.2} GB", b / GB)
    } else if b >= MB {
        format!("{:.2} MB", b / MB)
    } else if b >= KB {
        format!("{:.2} KB", b / KB)
    } else {
        format!("{bytes} B")
    }
}

fn pct(part: usize, whole: usize) -> String {
    if whole == 0 {
        return "0.0%".to_string();
    }
    format!("{:.1}%", 100.0 * part as f64 / whole as f64)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let tips: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(1000);
    let trees: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(500);
    let distinct = args.get(3).map(|s| s != "same").unwrap_or(true);

    println!(
        "memory_quick: tips={tips} trees={trees} topology={}\n",
        if distinct { "distinct" } else { "identical" }
    );

    let newicks: Vec<String> = (0..trees)
        .map(|i| make_newick(tips, i as u64, distinct))
        .collect();
    let refs: Vec<&str> = newicks.iter().map(|s| s.as_str()).collect();

    // --- measure retained footprint of construction (default path, lengths kept) ---
    PEAK.store(CURRENT.load(Ordering::Relaxed), Ordering::Relaxed);
    let before = CURRENT.load(Ordering::Relaxed);
    let snaps = Snapshots::from_newicks(&refs, false).expect("parse failed");
    let after = CURRENT.load(Ordering::Relaxed);
    let measured = after.saturating_sub(before);
    let peak = PEAK.load(Ordering::Relaxed).saturating_sub(before);

    // --- analytic component sizes (public fields + known per-tree edge count) ---
    let words = snaps.words_per_bitset;
    let n_bip = snaps.bipartitions.len();
    let n_trees = snaps.len();
    let bitset_inline = mem::size_of::<Bitset>();
    let bitset_total = bitset_inline + words * 8; // struct + heap u64s

    // Unrooted binary tree: n pendant + (n-3) internal = 2n-3 edges per tree.
    let per_tree = if tips >= 3 { 2 * tips - 3 } else { tips };
    let total_occ = n_trees * per_tree;

    // Removed in this branch:
    //   bipartition_index = a cloned Bitset key + u32 value + map overhead per unique split.
    let index_cost = n_bip * (bitset_total + 4 + 32);
    //   lengths (RF-only paths) = f64 per split occurrence.
    let lengths_cost = total_occ * 8;
    //   bipartitions table (released on the non-export RF path only).
    let bip_cost = n_bip * bitset_total;
    let split_ids_cost = total_occ * 4;

    println!("Interned structure:");
    println!("  unique splits (n_bip):     {n_bip}");
    println!(
        "  words per bitset:          {words}  ({} per bitset)",
        human(bitset_total)
    );
    println!("  split occurrences (total): {total_occ}");
    println!();

    println!("Measured retained heap (this branch, default `from_newicks`):");
    println!("  retained: {}", human(measured));
    println!("  peak:     {}", human(peak));
    println!();

    // Footprints for the RF path, which is treetracer's hot path.
    let master_rf = measured + index_cost; // old: + duplicate index (lengths already counted)
    let new_rf = measured.saturating_sub(lengths_cost + bip_cost); // no index, no lengths, table freed
    let saved = master_rf.saturating_sub(new_rf);

    println!("RF path footprint (the hot path):");
    println!(
        "  master  (index + lengths + table): {:>12}",
        human(master_rf)
    );
    println!(
        "  this PR (no index, no lengths, table freed): {:>12}",
        human(new_rf)
    );
    println!(
        "  saved: {}  ({} of master)",
        human(saved),
        pct(saved, master_rf)
    );
    println!();

    println!("Breakdown of removed terms (vs master RF path):");
    println!(
        "  bipartition_index: {:>12}  ({} of master)",
        human(index_cost),
        pct(index_cost, master_rf)
    );
    println!(
        "  branch lengths:    {:>12}  ({} of master)",
        human(lengths_cost),
        pct(lengths_cost, master_rf)
    );
    println!(
        "  bipartition table: {:>12}  ({} of master)",
        human(bip_cost),
        pct(bip_cost, master_rf)
    );
    println!("  (split_ids kept:   {:>12})", human(split_ids_cost));
}
