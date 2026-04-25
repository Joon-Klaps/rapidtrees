use clap::{Parser, ValueEnum};
use phylotree::tree::Tree;
use rapidtrees::distances::{pairwise_kf_matrix, pairwise_rf_matrix, pairwise_wrf_matrix};
use rapidtrees::io::{read_beast_trees, read_snap, write_matrix_tsv, write_snap};
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
    Presence { n_bip: usize, presence: Vec<u8> },
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

    let (names, distance_input) = match &args.snap_input {
        Some(path) => load_from_snap(path, quiet),
        None => load_from_beast(&args, quiet),
    };

    let n_pairs = names.len() * (names.len() - 1) / 2;
    let metric_label = metric_label(args.metric);
    log_if(
        quiet,
        format!("Determining {metric_label} distances for {n_pairs} pairs"),
    );

    let t = Instant::now();
    let mat = compute_matrix(args.metric, &names, &distance_input);
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

fn load_from_snap(path: &PathBuf, quiet: bool) -> (Vec<String>, DistanceInput) {
    let t = Instant::now();
    let (names, _taxa_names, n_bip, presence) = read_snap(path).unwrap_or_else(|e| {
        eprintln!("Failed to read snap {path:?}: {e}");
        std::process::exit(2);
    });
    log_if(
        quiet,
        format!(
            "Read snap with {} trees and {n_bip} bipartitions in {:.3}s",
            names.len(),
            t.elapsed().as_secs_f64()
        ),
    );
    (names, DistanceInput::Presence { n_bip, presence })
}

fn load_from_beast(args: &Args, quiet: bool) -> (Vec<String>, DistanceInput) {
    let path = args
        .input
        .as_ref()
        .expect("clap guarantees either --input or --snap-input");

    let t = Instant::now();
    let (taxons, named_trees) = read_beast_trees(
        path,
        args.burnin_trees,
        args.burnin_states,
        args.use_real_taxa,
    );
    if named_trees.is_empty() {
        eprintln!("No trees parsed from {path:?}.");
        std::process::exit(2);
    }
    log_if(
        quiet,
        format!(
            "Read {} trees ({} taxa) from BEAST file in {:.3}s",
            named_trees.len(),
            taxons.len(),
            t.elapsed().as_secs_f64()
        ),
    );

    let (names, trees): (Vec<String>, Vec<_>) = named_trees.into_iter().unzip();

    let t = Instant::now();
    let snaps: Vec<TreeSnapshot> = trees
        .iter()
        .map(|t| TreeSnapshot::from_tree(t, args.rooted))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|e| {
            eprintln!("Failed to build snapshots: {e}");
            std::process::exit(3);
        });
    log_if(
        quiet,
        format!("Built tree snapshots in {:.3}s", t.elapsed().as_secs_f64()),
    );

    if let Some(snap_path) = &args.export_snap {
        export_snap(snap_path, &names, &trees, &snaps, quiet);
    }

    (names, DistanceInput::Snapshots(snaps))
}

fn export_snap(
    snap_path: &PathBuf,
    names: &[String],
    trees: &[Tree],
    snaps: &[TreeSnapshot],
    quiet: bool,
) {
    let mut taxa_names: Vec<String> = trees[0]
        .get_leaves()
        .iter()
        .filter_map(|&id| trees[0].get(&id).ok()?.name.clone())
        .collect();
    taxa_names.sort_unstable();

    let t = Instant::now();
    if let Err(e) = write_snap(snap_path, names, &taxa_names, snaps) {
        eprintln!("Failed to write snap {snap_path:?}: {e}");
        std::process::exit(5);
    }
    log_if(
        quiet,
        format!(
            "Exported snap ({} trees, {} taxa) to {} in {:.3}s",
            snaps.len(),
            taxa_names.len(),
            snap_path.display(),
            t.elapsed().as_secs_f64()
        ),
    );
}

fn compute_matrix(metric: MetricArg, names: &[String], input: &DistanceInput) -> Vec<Vec<f64>> {
    match (metric, input) {
        (MetricArg::Rf, DistanceInput::Snapshots(snaps)) => pairwise_rf_matrix(snaps)
            .into_iter()
            .map(|row| row.into_iter().map(|v| v as f64).collect())
            .collect(),

        (MetricArg::Rf, DistanceInput::Presence { n_bip, presence }) => {
            let n = names.len();
            let mut mat = vec![vec![0.0_f64; n]; n];
            for i in 0..n {
                for j in i + 1..n {
                    let rf = presence[i * n_bip..(i + 1) * n_bip]
                        .iter()
                        .zip(&presence[j * n_bip..(j + 1) * n_bip])
                        .filter(|(a, b)| a != b)
                        .count() as f64;
                    mat[i][j] = rf;
                    mat[j][i] = rf;
                }
            }
            mat
        }

        (MetricArg::Weighted, DistanceInput::Snapshots(snaps)) => pairwise_wrf_matrix(snaps),
        (MetricArg::Kf, DistanceInput::Snapshots(snaps)) => pairwise_kf_matrix(snaps),
        _ => unreachable!("snap input is restricted to RF by argument validation"),
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
