use clap::{Parser, ValueEnum};
use rapidtrees::distances::{pairwise_kf_matrix, pairwise_rf_matrix, pairwise_wrf_matrix};
use rapidtrees::io::{load_beast_trees, load_snapshots, write_matrix_tsv, write_snap};
use rapidtrees::snapshot::TreeSnapshot;
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

    /// Path to compressed .snap file (RF only)
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
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

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

enum DistanceInput {
    Snapshots(Vec<TreeSnapshot>),
}

fn main() {
    let args = Args::parse();

    if args.snap_input.is_some() && !matches!(args.metric, MetricArg::Rf) {
        eprintln!("--snap-input currently supports only --metric rf");
        std::process::exit(6);
    }
    if args.snap_input.is_some() && args.export_snap.is_some() {
        eprintln!("--export-snap cannot be used together with --snap-input");
        std::process::exit(7);
    }

    let quiet = args.quiet;

    let t = Instant::now();
    let (taxon_map, trees) = match &args.snap_input {
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

    let (names, snapshots): (Vec<String>, Vec<TreeSnapshot>) = trees.into_iter().unzip();

    log_if(
        quiet,
        format!(
            "Loaded {} trees from {} in {:.3}s",
            names.len(),
            if let Some(path) = &args.snap_input {
                format!("snap {path:?}")
            } else {
                "BEAST file".to_string()
            },
            t.elapsed().as_secs_f64()
        ),
    );

    // If --export-snap is specified, write the snapshots to a .snap file and exit
    if let Some(snap_path) = args.export_snap {
        let t = Instant::now();
        let mut taxa_names: Vec<String> = taxon_map.values().cloned().collect();
        taxa_names.sort_unstable();
        if let Err(e) = write_snap(&snap_path, &names, &taxa_names, &snapshots) {
            eprintln!("Failed to write snap {snap_path:?}: {e}");
            std::process::exit(5);
        }
        log_if(
            quiet,
            format!(
                "Exported snapshots to snap {snap_path:?} in {:.3}s",
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
    let mat = compute_matrix(args.metric, &DistanceInput::Snapshots(snapshots));
    log_if(
        quiet,
        format!(
            "Computed {metric_label} distances in {:.3}s",
            t.elapsed().as_secs_f64()
        ),
    );

    let t = Instant::now();
    if let Err(e) = write_matrix_tsv(&args.output, &names, &mat) {
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
}

fn compute_matrix(metric: MetricArg, input: &DistanceInput) -> Vec<Vec<f64>> {
    match (metric, input) {
        (MetricArg::Rf, DistanceInput::Snapshots(snaps)) => pairwise_rf_matrix(snaps)
            .into_iter()
            .map(|row| row.into_iter().map(|v| v as f64).collect())
            .collect(),

        (MetricArg::Weighted, DistanceInput::Snapshots(snaps)) => pairwise_wrf_matrix(snaps),
        (MetricArg::Kf, DistanceInput::Snapshots(snaps)) => pairwise_kf_matrix(snaps),
    }
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
