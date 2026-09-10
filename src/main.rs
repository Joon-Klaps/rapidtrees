use clap::{Parser, ValueEnum};
use rapidtrees::io::{load_beast_trees, write_matrix_tsv};
use rapidtrees::{Backend, Kernel};
use std::io::{IsTerminal, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

#[derive(Parser, Debug)]
#[command(
    name = "rapidtrees",
    version,
    about = "Fast pairwise tree distance calculations (Robinson-Foulds, Weighted RF, Kuhner-Felsenstein) for phylogenetic trees"
)]
struct Args {
    /// Path to a tree file: BEAST/NEXUS .trees or plain Newick (format is auto-detected)
    #[arg(short = 'i', long = "input")]
    input: PathBuf,

    /// Burn-in by number of trees (drop first N trees)
    #[arg(short = 't', long = "burnin-trees", default_value_t = 0)]
    burnin_trees: usize,

    /// Burn-in by state (keep trees with STATE_ > value); NEXUS only, ignored for Newick
    #[arg(short = 's', long = "burnin-states", default_value_t = 0)]
    burnin_states: usize,

    /// Output path for TSV distance matrix
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Use TRANSLATE block to map taxon IDs to labels when available
    #[arg(long = "use-real-taxa", default_value_t = false)]
    use_real_taxa: bool,

    /// Distance metric to compute: rf | weighted | kf
    #[arg(long = "metric", value_enum, default_value_t = MetricArg::Rf)]
    metric: MetricArg,

    /// Compute rooted distances (compare clades) instead of unrooted (compare bipartitions)
    #[arg(long = "rooted", default_value_t = false)]
    rooted: bool,

    /// Per-pair kernel: auto | dense | sparse. `auto` picks from tree diversity
    #[arg(long = "backend", value_enum, default_value_t = Backend::Auto)]
    backend: Backend,

    /// Quiet mode: suppresses progress messages on stdout
    #[arg(short = 'q', long = "quiet", default_value_t = false)]
    quiet: bool,
}

#[derive(Copy, Clone, Debug, PartialEq, ValueEnum)]
enum MetricArg {
    Rf,
    Weighted,
    Kf,
}

fn main() {
    let args = Args::parse();

    let quiet = args.quiet;
    let t_total = Instant::now();
    let t = Instant::now();

    let (names, interned) = load_beast_trees(
        &args.input,
        args.burnin_trees,
        args.burnin_states,
        args.use_real_taxa,
        args.rooted,
    );

    log_if(
        quiet,
        format!(
            "Loaded {} trees in {:.3}s",
            names.len(),
            t.elapsed().as_secs_f64()
        ),
    );

    if names.len() < 2 {
        eprintln!("Need at least 2 trees to compute pairwise distances");
        std::process::exit(2);
    }

    log_collision_bound(quiet, interned.n_distinct_splits());

    let n_pairs = names.len() * (names.len() - 1) / 2;
    let metric_label = metric_label(args.metric);
    log_if(
        quiet,
        format!("Determining {metric_label} distances for {n_pairs} pairs"),
    );

    let t = Instant::now();
    let show_progress = !quiet && std::io::stderr().is_terminal();

    let output_path = args.output.as_path();

    // A macro rather than a generic fn: each metric keeps its own element type
    // (`u32` for RF, `f64` for the weighted pair) all the way to the writer, and
    // threading that through a function would need ten parameters to say what
    // three tokens say here.
    macro_rules! compute_and_write {
        ($metric:ident) => {{
            let dist = run_with_progress(n_pairs, show_progress, |counter| {
                interned.$metric(Some(counter), args.backend)
            });
            log_backend(quiet, args.backend, dist.kernel);
            log_computed(quiet, metric_label, &t);
            let t = Instant::now();
            let r = write_matrix_tsv(output_path, &names, &dist.matrix, interned.len());
            (r, t)
        }};
    }

    let (write_result, t) = match args.metric {
        MetricArg::Rf => compute_and_write!(pairwise_rf_with),
        MetricArg::Weighted => compute_and_write!(pairwise_wrf_with),
        MetricArg::Kf => compute_and_write!(pairwise_kf_with),
    };
    if let Err(e) = write_result {
        eprintln!("Failed to write output {}: {e}", output_path.display());
        std::process::exit(4);
    }
    log_if(
        quiet,
        format!(
            "Wrote matrix to {} in {:.3}s",
            output_path.display(),
            t.elapsed().as_secs_f64()
        ),
    );

    log_if(
        quiet,
        format!("Total runtime: {:.3}s", t_total.elapsed().as_secs_f64()),
    );
}

fn log_computed(quiet: bool, metric_label: &str, t: &Instant) {
    log_if(
        quiet,
        format!(
            "Computed {metric_label} distances in {:.3}s",
            t.elapsed().as_secs_f64()
        ),
    );
}

fn metric_label(metric: MetricArg) -> &'static str {
    match metric {
        MetricArg::Rf => "RF",
        MetricArg::Weighted => "Weighted RF",
        MetricArg::Kf => "KF",
    }
}

/// State the run's own correctness guarantee.
fn log_collision_bound(quiet: bool, distinct_splits: usize) {
    let e = distinct_splits as f64;
    // 2¹²⁹ overflows nothing here, but stays clearer written as a power.
    let bound = e * e / 2f64.powi(129);
    log_if(
        quiet,
        format!("Distinct splits e = {distinct_splits}; collision bound e²/2¹²⁹ = {bound:.2e}"),
    );
}

/// Name the kernel that ran, so `auto`'s choice is in the run log.
fn log_backend(quiet: bool, requested: Backend, chosen: Kernel) {
    let requested = format!("{requested:?}").to_lowercase();
    log_if(quiet, format!("Backend: {chosen} (--backend {requested})"));
}

fn log_if(quiet: bool, msg: String) {
    if !quiet {
        println!("{msg}");
    }
}

/// Run `work` while rendering a stderr progress bar driven by an atomic counter
/// that `work` increments as it completes pairs.
///
/// When `show_progress` is `false` the function reduces to `work(&counter)`
/// without spawning a monitor thread, so piped/redirected runs and `--quiet`
/// retain their original zero-overhead behaviour.
fn run_with_progress<R, F>(n_pairs: usize, show_progress: bool, work: F) -> R
where
    F: FnOnce(&AtomicUsize) -> R,
{
    if !show_progress {
        let counter = AtomicUsize::new(0);
        return work(&counter);
    }

    let counter = Arc::new(AtomicUsize::new(0));
    let terminate = Arc::new(AtomicBool::new(false));
    let monitor_counter = Arc::clone(&counter);
    let monitor_terminate = Arc::clone(&terminate);
    let start = Instant::now();

    let monitor = std::thread::spawn(move || {
        while !monitor_terminate.load(Ordering::Relaxed) {
            std::thread::sleep(Duration::from_millis(100));
            let done = monitor_counter.load(Ordering::Relaxed);
            render_pair_bar(done, n_pairs, start.elapsed());
        }
    });

    let result = work(&counter);

    terminate.store(true, Ordering::Relaxed);
    let _ = monitor.join();

    // Final 100%-bar so the user sees the completed state, then a newline so
    // subsequent log lines start on a clean row.
    render_pair_bar(n_pairs, n_pairs, start.elapsed());
    let _ = writeln!(std::io::stderr());
    result
}

/// Carriage-return-overwrite a single-line ASCII/Unicode progress bar on stderr.
fn render_pair_bar(done: usize, total: usize, elapsed: Duration) {
    const WIDTH: usize = 40;
    let frac = if total == 0 {
        0.0
    } else {
        (done as f64 / total as f64).clamp(0.0, 1.0)
    };
    let filled = (frac * WIDTH as f64) as usize;
    let mut bar = String::with_capacity(WIDTH * 3);
    for i in 0..WIDTH {
        bar.push(if i < filled { '█' } else { '░' });
    }
    let eta = if frac > 0.001 {
        elapsed.as_secs_f64() * (1.0 - frac) / frac
    } else {
        0.0
    };
    // \r resets the cursor; trailing spaces guard against shorter values
    // leaving stale digits behind on the line.
    let _ = write!(
        std::io::stderr(),
        "\r  [{bar}] {pct:5.1}% ({done}/{total} pairs, ETA {eta:5.1}s)   ",
        pct = frac * 100.0,
    );
    let _ = std::io::stderr().flush();
}
