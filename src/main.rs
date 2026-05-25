use clap::{Parser, ValueEnum};
use rapidtrees::io::{load_beast_trees, load_snapshots, write_matrix_tsv, write_snap};
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
    /// Path to BEAST .trees (NEXUS) file
    #[arg(
        short = 'i',
        long = "input",
        required_unless_present = "snap_input",
        conflicts_with = "snap_input"
    )]
    input: Option<PathBuf>,

    /// Path to compressed .snap file
    #[arg(
        long = "snap-input",
        required_unless_present = "input",
        conflicts_with = "input"
    )]
    snap_input: Option<PathBuf>,

    /// Burn-in by number of trees (drop first N trees)
    #[arg(short = 't', long = "burnin-trees", default_value_t = 0)]
    burnin_trees: usize,

    /// Burn-in by state (keep trees with STATE_ > value)
    #[arg(short = 's', long = "burnin-states", default_value_t = 0)]
    burnin_states: usize,

    /// Output path for TSV distance matrix
    #[arg(short = 'o', long = "output", required_unless_present = "export_snap")]
    output: Option<PathBuf>,

    #[arg(long = "export-snap")]
    export_snap: Option<PathBuf>,

    /// Use TRANSLATE block to map taxon IDs to labels when available
    #[arg(long = "use-real-taxa", default_value_t = false)]
    use_real_taxa: bool,

    /// Distance metric to compute: rf | weighted | kf
    #[arg(long = "metric", value_enum, default_value_t = MetricArg::Rf)]
    metric: MetricArg,

    /// Compute rooted distances (compare clades) instead of unrooted (compare bipartitions)
    #[arg(long = "rooted", default_value_t = false)]
    rooted: bool,

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

    if args.snap_input.is_some() && args.export_snap.is_some() {
        eprintln!("--export-snap cannot be used together with --snap-input");
        std::process::exit(7);
    }

    if args.snap_input.is_some() && (args.metric != MetricArg::Rf) {
        eprintln!("--snap-input only supports RF distances.");
        std::process::exit(8);
    }

    let quiet = args.quiet;
    let t_total = Instant::now();
    let t = Instant::now();

    let (names, interned) = match &args.snap_input {
        Some(path) => load_snapshots(path).unwrap_or_else(|e| {
            eprintln!("Failed to load trees: {e}");
            std::process::exit(1);
        }),
        None => load_beast_trees(
            args.input.as_ref().unwrap(),
            args.burnin_trees,
            args.burnin_states,
            args.use_real_taxa,
            args.rooted,
        ),
    };

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

    if let Some(snap_path) = args.export_snap {
        let t = Instant::now();
        if let Err(e) = write_snap(&snap_path, &names, &interned) {
            eprintln!("Failed to write snap {snap_path:?}: {e}");
            std::process::exit(5);
        }
        log_if(
            quiet,
            format!(
                "Exported snapshots to {snap_path:?} in {:.3}s",
                t.elapsed().as_secs_f64()
            ),
        );
        return;
    }

    let n_pairs = names.len() * (names.len() - 1) / 2;
    let metric_label = metric_label(args.metric);
    log_if(
        quiet,
        format!("Determining {metric_label} distances for {n_pairs} pairs"),
    );

    let t = Instant::now();
    let show_progress = !quiet && std::io::stderr().is_terminal();
    let mat: Vec<f64> = match args.metric {
        MetricArg::Rf => run_with_progress(n_pairs, show_progress, |counter| {
            interned
                .pairwise_rf_counted(counter)
                .into_iter()
                .map(|dist| dist as f64)
                .collect()
        }),
        MetricArg::Weighted => run_with_progress(n_pairs, show_progress, |counter| {
            interned.pairwise_wrf_counted(counter)
        }),
        MetricArg::Kf => run_with_progress(n_pairs, show_progress, |counter| {
            interned.pairwise_kf_counted(counter)
        }),
    };
    log_if(
        quiet,
        format!(
            "Computed {metric_label} distances in {:.3}s",
            t.elapsed().as_secs_f64()
        ),
    );

    let t = Instant::now();
    let output_path = args
        .output
        .as_deref()
        .expect("output is required when not exporting snap");
    if let Err(e) = write_matrix_tsv(output_path, &names, &mat, interned.len()) {
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

fn metric_label(metric: MetricArg) -> &'static str {
    match metric {
        MetricArg::Rf => "RF",
        MetricArg::Weighted => "Weighted RF",
        MetricArg::Kf => "KF",
    }
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
fn run_with_progress<T, F>(n_pairs: usize, show_progress: bool, work: F) -> Vec<T>
where
    T: Send,
    F: FnOnce(&AtomicUsize) -> Vec<T>,
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
