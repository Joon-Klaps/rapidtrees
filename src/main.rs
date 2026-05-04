use clap::{Parser, ValueEnum};
use rapidtrees::io::{load_beast_trees, load_snapshots, write_matrix_tsv, write_snap};
use std::path::PathBuf;
use std::time::Instant;

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

#[derive(Copy, Clone, Debug, ValueEnum)]
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
    let mat: Vec<f64> = match args.metric {
        MetricArg::Rf => interned
            .pairwise_rf()
            .into_iter()
            .map(|dist| dist as f64)
            .collect(),
        MetricArg::Weighted => interned.pairwise_wrf(),
        MetricArg::Kf => interned.pairwise_kf(),
    };
    log_if(
        quiet,
        format!(
            "Computed {metric_label} distances in {:.3}s",
            t.elapsed().as_secs_f64()
        ),
    );

    let t = Instant::now();
    if let Err(e) = write_matrix_tsv(&args.output, &names, &mat, interned.len()) {
        eprintln!("Failed to write output {:?}: {e}", args.output);
        std::process::exit(4);
    }
    log_if(
        quiet,
        format!(
            "Wrote matrix to {} in {:.3}s",
            args.output.display(),
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
