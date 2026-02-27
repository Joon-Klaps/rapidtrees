use clap::{Parser, ValueEnum};
use rapidtrees::distances::{kf_from_snapshots, rf_from_snapshots, weighted_rf_from_snapshots};
use rapidtrees::io::{read_beast_trees, write_matrix_tsv};
use rapidtrees::snapshot::TreeSnapshot;
use rayon::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

/// Compute pairwise Robinson–Foulds distances from a BEAST/NEXUS tree file
/// and write a labeled distance matrix (TSV) where row/column names are tree names.
#[derive(Parser, Debug)]
#[command(
    name = "rapidtrees",
    version,
    about = "Fast pairwise tree distance calculations (Robinson-Foulds, Weighted RF, Kuhner-Felsenstein) for phylogenetic trees"
)]
struct Args {
    /// Path to BEAST .trees (NEXUS) file
    #[arg(short = 'i', long = "input")]
    input: PathBuf,

    /// Burn-in by number of trees (drop first N trees)
    #[arg(short = 't', long = "burnin-trees", default_value_t = 0)]
    burnin_trees: usize,

    /// Burn-in by state (keep trees with STATE_ > value)
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

    // Read trees with names
    let t0 = Instant::now();
    let (taxons, named_trees) = read_beast_trees(
        &args.input,
        args.burnin_trees,
        args.burnin_states,
        args.use_real_taxa,
    );
    if named_trees.is_empty() {
        eprintln!("No trees parsed from {:?}.", args.input);
        std::process::exit(2);
    }
    let read_s = t0.elapsed().as_secs_f64();
    log_if(!args.quiet, format!("Read beast trees in {read_s:.3}s"));
    log_if(
        !args.quiet,
        format!(
            "Read in {} taxons for {} trees",
            taxons.len(),
            named_trees.len()
        ),
    );
    let (names, trees): (Vec<String>, Vec<_>) = named_trees.into_iter().unzip();

    // Build bitset snapshots once
    let t1 = Instant::now();
    let snaps: Vec<TreeSnapshot> = trees
        .iter()
        .map(TreeSnapshot::from_tree)
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|e| {
            eprintln!("Failed to build snapshots: {e}");
            std::process::exit(3);
        });

    let snap_s = t1.elapsed().as_secs_f64();
    log_if(
        !args.quiet,
        format!("Created tree bit snapshots in {snap_s:.3}s"),
    );

    let t2 = Instant::now();
    let (metric_label, metric_fn): (&str, fn(&TreeSnapshot, &TreeSnapshot) -> f64) =
        match args.metric {
            // rf is the only one that returns usize, so cast to f64
            MetricArg::Rf => ("RF", |a, b| rf_from_snapshots(a, b) as f64),
            MetricArg::Weighted => ("Weighted", weighted_rf_from_snapshots),
            MetricArg::Kf => ("KF", kf_from_snapshots),
        };

    log_if(
        !args.quiet,
        format!(
            "Determining distances using {metric_label} for {} combinations",
            names.len() * (names.len() - 1) / 2
        ),
    );

    let n = names.len();

    // Pre-allocate matrix for better performance
    let mut mat = vec![vec![0.0f64; n]; n];

    // Compute distances in parallel using work-stealing for perfect load balance
    // Each row is independent, Rayon handles distribution automatically
    mat.par_iter_mut().enumerate().for_each(|(i, row)| {
        for j in (i + 1)..n {
            let dist = metric_fn(&snaps[i], &snaps[j]);
            row[j] = dist;
        }
    });

    // Fill lower triangle (symmetric matrix)
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        for j in (i + 1)..n {
            mat[j][i] = mat[i][j];
        }
    }

    let comp_s = t2.elapsed().as_secs_f64();
    log_if(
        !args.quiet,
        format!("Determined distances using {metric_label} in {comp_s:.3}s"),
    );

    let t3 = Instant::now();
    if let Err(e) = write_matrix_tsv(&args.output, &names, &mat) {
        eprintln!("Failed to write output {:?}: {e}", args.output);
        std::process::exit(4);
    }
    let write_s = t3.elapsed().as_secs_f64();
    log_if(
        !args.quiet,
        format!("Writing to {} in {write_s:.3}s", args.output.display()),
    );
}

fn log_if(show: bool, msg: String) {
    if show {
        println!("{}", msg);
    }
}
