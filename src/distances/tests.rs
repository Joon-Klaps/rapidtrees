//! Behavioural tests for the distance kernels: known values, metric
//! properties, and agreement with a per-pair reference on random trees.

use super::boundary::{
    RF_DENSE_SHARE, WEIGHTED_DENSE_SHARE, min_dense, parse_share, rf_dense_share,
    weighted_dense_share,
};
use super::layout::{Layout, NO_COLUMN, Trees, assign_columns};
use super::matrix::fill_symmetric;
use super::rf::distance_rf_split;
use super::treedist::TREEDIST_TREES;
use super::weighted::{PANEL, sweep, weighted_distances_split, weighted_layout};
use crate::snapshot::{InternSnap, Snapshots};
use std::collections::{BTreeMap, BTreeSet};

const T0: &str = "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";
const T1: &str = "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);";
const T2: &str = "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";

fn three_snapshots() -> Snapshots {
    Snapshots::from_newicks(&[T0, T1, T2], false).unwrap()
}

#[test]
#[allow(clippy::erasing_op)] // Keep O*n for clarity of the flattened matrix indexing
fn rf_known_values() {
    let snaps = three_snapshots();
    let n = snaps.snapshots.len(); // n = 3
    let mat = snaps.pairwise_rf(None);
    assert_eq!(mat[0 * n + 1], 4, "RF(T0,T1)");
    assert_eq!(mat[0 * n + 2], 2, "RF(T0,T2)");
    assert_eq!(mat[n + 2], 2, "RF(T1,T2)");
    for i in 0..n {
        assert_eq!(mat[i * n + i], 0, "diagonal [{i}]");
        for j in 0..n {
            assert_eq!(mat[i * n + j], mat[j * n + i], "RF symmetry [{i}][{j}]");
        }
    }
}

/// Assert `mat` has an exact-zero diagonal and is symmetric.
fn assert_symmetric_zero_diagonal(mat: &[f64], n: usize, metric: &str) {
    for (i, row) in mat.chunks(n).enumerate() {
        assert_eq!(row[i], 0.0, "{metric} diagonal [{i}]");
        for (j, v) in row.iter().enumerate() {
            assert!(
                (v - mat[j * n + i]).abs() < f64::EPSILON,
                "{metric} symmetry [{i}][{j}]"
            );
        }
    }
}

#[test]
fn weighted_symmetric_zero_diagonal() {
    let snaps = three_snapshots();
    let n = snaps.snapshots.len();
    assert_symmetric_zero_diagonal(&snaps.pairwise_wrf(None), n, "WRF");
    assert_symmetric_zero_diagonal(&snaps.pairwise_kf(None), n, "KF");
}

#[test]
fn kf_differs_from_wrf() {
    let snaps = three_snapshots();
    let wrf = snaps.pairwise_wrf(None);
    let kf = snaps.pairwise_kf(None);
    let n = snaps.snapshots.len();
    let any_different = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .any(|(i, j)| (kf[i * n + j] - wrf[i * n + j]).abs() > f64::EPSILON);
    assert!(any_different, "KF and WRF should produce different values");
}

/// Assert `d(i, k) ≤ d(i, j) + d(j, k)` for every triple, with a small
/// epsilon for the float metrics.
fn assert_triangle_inequality<T: Copy + Into<f64>>(mat: &[T], n: usize, metric: &str) {
    let d = |a: usize, b: usize| -> f64 { mat[a * n + b].into() };
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                assert!(
                    d(i, k) <= d(i, j) + d(j, k) + 1e-10,
                    "{metric} triangle inequality failed for indices {i}, {j}, {k}"
                );
            }
        }
    }
}

#[test]
fn metrics_satisfy_triangle_inequality() {
    let snaps = three_snapshots();
    let n = snaps.snapshots.len();
    assert_triangle_inequality(&snaps.pairwise_rf(None), n, "RF");
    assert_triangle_inequality(&snaps.pairwise_wrf(None), n, "WRF");
    assert_triangle_inequality(&snaps.pairwise_kf(None), n, "KF");
}

#[test]
fn split_counts_match_a_recount() {
    let snaps = three_snapshots();
    let mut recount = vec![0u32; snaps.n_distinct_splits()];
    for snap in &snaps.snapshots {
        for id in snap.ids() {
            recount[id as usize] += 1;
        }
    }
    assert_eq!(snaps.split_counts, recount);
}

#[test]
fn assign_columns_orders_by_descending_count() {
    // Split 1 is held by one tree, so `count >= 2` drops it.
    let (column_of, kept) = assign_columns(&[3, 1, 2, 3], 3, |count| count >= 2);
    assert_eq!(kept, 3);
    assert_eq!(column_of, vec![0, NO_COLUMN, 2, 1]);
}

#[test]
fn lengths_aligned_with_split_ids() {
    let snaps = three_snapshots();
    for snap in &snaps.snapshots {
        assert_eq!(snap.n_splits(), snap.lengths.len());
    }
}

/// Test oracle: all three metrics for one pair, straight from the
/// definitions over the union of both trees' splits (absent = length 0).
///
/// Deliberately naive — no filtering, no `selfᵢ + selfⱼ − 2·shared` rewrite —
/// so it cannot reproduce a backend bug by construction.
fn reference_distances(a: &InternSnap, b: &InternSnap) -> (usize, f64, f64) {
    let lengths_by_split =
        |s: &InternSnap| -> BTreeMap<u32, f64> { s.ids().zip(s.lengths.iter().copied()).collect() };
    let (ma, mb) = (lengths_by_split(a), lengths_by_split(b));

    let (mut rf, mut wrf, mut sum_sq) = (0usize, 0.0, 0.0);
    for id in ma
        .keys()
        .chain(mb.keys())
        .copied()
        .collect::<BTreeSet<u32>>()
    {
        let x = ma.get(&id).copied().unwrap_or(0.0);
        let y = mb.get(&id).copied().unwrap_or(0.0);
        if ma.contains_key(&id) != mb.contains_key(&id) {
            rf += 1;
        }
        wrf += (x - y).abs();
        sum_sq += (x - y) * (x - y);
    }
    (rf, wrf, sum_sq.sqrt())
}

/// Assert `got` is within a relative 1e-9 of the oracle's `want`.
fn assert_close(got: f64, want: f64, what: &str) {
    assert!(
        (got - want).abs() <= 1e-9 * want.max(1.0),
        "{what}: {got} vs reference {want}"
    );
}

/// Assert WRF and KF matrices agree with [`reference_distances`], and that
/// copies of the same Newick land on exactly 0.0, not a rounding crumb.
fn assert_weighted_match_reference<S: AsRef<str>>(
    snaps: &Snapshots,
    newicks: &[S],
    (wrf, kf): (&[f64], &[f64]),
    ctx: &str,
) {
    let n = snaps.snapshots.len();
    for i in 0..n {
        for j in 0..n {
            let (_, want_wrf, want_kf) =
                reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);
            let at = i * n + j;
            assert_close(wrf[at], want_wrf, &format!("WRF at [{i}][{j}], {ctx}"));
            assert_close(kf[at], want_kf, &format!("KF at [{i}][{j}], {ctx}"));
            if newicks[i].as_ref() == newicks[j].as_ref() {
                assert_eq!(wrf[at], 0.0, "WRF identical [{i}][{j}], {ctx}");
                assert_eq!(kf[at], 0.0, "KF identical [{i}][{j}], {ctx}");
            }
        }
    }
}

/// Assert an RF matrix agrees exactly with [`reference_distances`].
fn assert_rf_matches_reference(snaps: &Snapshots, rf: &[u32], ctx: &str) {
    let n = snaps.snapshots.len();
    for i in 0..n {
        for j in 0..n {
            let (want, _, _) = reference_distances(&snaps.snapshots[i], &snaps.snapshots[j]);
            assert_eq!(rf[i * n + j] as usize, want, "RF at [{i}][{j}], {ctx}");
        }
    }
}

/// Assert all three metrics agree with [`reference_distances`] on the
/// snapshots built from `newicks`.
fn assert_matches_reference<S: AsRef<str>>(snaps: &Snapshots, newicks: &[S], ctx: &str) {
    assert_rf_matches_reference(snaps, &snaps.pairwise_rf(None), ctx);
    let (wrf, kf) = (snaps.pairwise_wrf(None), snaps.pairwise_kf(None));
    assert_weighted_match_reference(snaps, newicks, (&wrf, &kf), ctx);
}

#[test]
fn pairwise_rf_identical_trees_all_zero() {
    // Every split is present in all 5 trees → every column is universal and
    // gets dropped, leaving zero-width bit rows. RF must still be 0 here.
    let t = TREEDIST_TREES[0];
    let snaps = Snapshots::from_newicks(&[t, t, t, t, t], false).unwrap();
    let mat = snaps.pairwise_rf(None);
    assert_eq!(mat.len(), 25);
    assert!(
        mat.iter().all(|&d| d == 0),
        "identical trees must have RF 0 everywhere"
    );
}

#[test]
fn pairwise_weighted_identical_trees_all_zero() {
    let t = TREEDIST_TREES[0];
    let snaps = Snapshots::from_newicks(&[t, t, t], false).unwrap();

    // Exactly 0.0 only holds because self and shared are summed in the same
    // column order; reordering either would leave a rounding crumb.
    assert!(
        snaps.pairwise_wrf(None).iter().all(|&d| d == 0.0),
        "identical trees must have WRF exactly 0 everywhere"
    );

    // Identical trees make `selfᵢ + selfⱼ − 2·dot` algebraically zero, but
    // rounding can nudge it slightly negative — the clamp must keep the
    // square root from producing NaN.
    assert!(
        snaps.pairwise_kf(None).iter().all(|&d| d == 0.0),
        "identical trees must have KF 0 everywhere (clamp must avoid NaN)"
    );
}

/// Trees sharing almost nothing, so most splits land in the "held by one
/// tree only" bucket, which gets no column and folds into its tree's `self`
/// term — which the treedist fixtures barely exercise. Branch lengths are
/// all distinct so a wrongly dropped split changes the answer visibly.
const DIVERSE_TREES: [&str; 4] = [
    "(((A:0.11,B:0.12):0.13,(C:0.14,D:0.15):0.16):0.17,((E:0.18,F:0.19):0.21,(G:0.22,H:0.23):0.24):0.25);",
    "(((A:0.31,E:0.32):0.33,(C:0.34,G:0.35):0.36):0.37,((B:0.38,F:0.39):0.41,(D:0.42,H:0.43):0.44):0.45);",
    "(((A:0.51,H:0.52):0.53,(B:0.54,G:0.55):0.56):0.57,((C:0.58,F:0.59):0.61,(D:0.62,E:0.63):0.64):0.65);",
    "(((A:0.71,D:0.72):0.73,(F:0.74,G:0.75):0.76):0.77,((B:0.78,H:0.79):0.81,(C:0.82,E:0.83):0.84):0.85);",
];

#[test]
fn pairwise_weighted_metrics_match_reference_on_diverse_trees() {
    let snaps = Snapshots::from_newicks(&DIVERSE_TREES, false).unwrap();
    assert_matches_reference(&snaps, &DIVERSE_TREES, "diverse fixtures");
}

/// Deterministic LCG, so a failure here is always reproducible.
fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

/// Random binary topology over `n_taxa` leaves, built by repeatedly joining
/// two randomly chosen subtrees.
fn random_newick(n_taxa: usize, state: &mut u64) -> String {
    let mut parts: Vec<String> = (0..n_taxa).map(|i| format!("T{i}")).collect();
    while parts.len() > 1 {
        let i = (lcg(state) as usize) % parts.len();
        let a = parts.swap_remove(i);
        let j = (lcg(state) as usize) % parts.len();
        let b = parts.swap_remove(j);
        let la = (lcg(state) % 100 + 1) as f64 / 100.0;
        let lb = (lcg(state) % 100 + 1) as f64 / 100.0;
        parts.push(format!("({a}:{la},{b}:{lb})"));
    }
    format!("{};", parts.pop().unwrap())
}

/// `n_trees` random trees over `n_taxa` leaves followed by `duplicates`
/// copies of the first ones, as Newicks and as the snapshots built from them.
fn random_snapshots(
    n_taxa: usize,
    n_trees: usize,
    duplicates: usize,
    seed: u64,
) -> (Vec<String>, Snapshots) {
    let mut state = seed;
    let mut newicks: Vec<String> = (0..n_trees)
        .map(|_| random_newick(n_taxa, &mut state))
        .collect();
    for k in 0..duplicates {
        newicks.push(newicks[k % n_trees].clone());
    }
    let refs: Vec<&str> = newicks.iter().map(String::as_str).collect();
    let snaps = Snapshots::from_newicks(&refs, false).unwrap();
    (newicks, snaps)
}

/// Differential test against [`reference_distances`] over a spread of taxon
/// counts and sharing levels.
///
/// `duplicates` matters because the column filter keys off how many trees
/// hold a split: distinct trees push splits into the dropped bucket,
/// duplicates pull them back into shared columns.
#[test]
fn pairwise_backends_match_reference_on_random_trees() {
    for &(n_taxa, n_trees, duplicates, seed) in &[
        (8usize, 6usize, 0usize, 1u64),
        (12, 8, 0, 2),
        (20, 10, 0, 3),
        (31, 9, 0, 4),
        (16, 6, 3, 5),   // forces splits into the shared-column bucket
        (24, 5, 5, 6),   // majority duplicates
        (10, 12, 11, 7), // all but one identical
    ] {
        let (newicks, snaps) = random_snapshots(n_taxa, n_trees, duplicates, seed);
        assert_matches_reference(&snaps, &newicks, &format!("seed={seed} taxa={n_taxa}"));
    }
}

/// WRF and KF with the dense/posting boundary at `min_dense`, mirroring
/// [`super::distance_wrf`] and [`super::distance_kf`]. The trees are read
/// borrowed and owned, and the two have to agree to the bit.
fn weighted_at(snaps: &Snapshots, min_dense: u32) -> (Vec<f64>, Vec<f64>) {
    let counts = &snaps.split_counts;
    let both = |trees: fn(&Snapshots) -> Trees| {
        (
            weighted_distances_split(trees(snaps), counts, None, min_dense, f64::min, |d| d),
            weighted_distances_split(
                trees(snaps),
                counts,
                None,
                min_dense,
                |a, b| a * b,
                f64::sqrt,
            ),
        )
    };
    let borrowed = both(|snaps| Trees::Borrowed(&snaps.snapshots));
    let owned = both(|snaps| Trees::Owned(snaps.snapshots.clone()));
    assert_eq!(
        borrowed, owned,
        "borrowed and owned trees disagree at min_dense={min_dense}"
    );
    borrowed
}

/// The boundary decides which path a split takes and must decide nothing
/// else: all dense (2), mixed, and all posting lists (`n + 1`) all have to
/// match the oracle.
#[test]
fn weighted_boundary_does_not_change_distances() {
    for &(n_taxa, n_trees, duplicates, seed) in &[
        (12usize, 8usize, 4usize, 11u64),
        (24, 10, 6, 12),
        (40, 6, 9, 13),
    ] {
        let (newicks, snaps) = random_snapshots(n_taxa, n_trees, duplicates, seed);
        let n = snaps.snapshots.len();
        for min_dense in [2, 3, n as u32 / 2, n as u32, n as u32 + 1] {
            let (wrf, kf) = weighted_at(&snaps, min_dense);
            let ctx = format!("seed={seed} min_dense={min_dense}");
            assert_weighted_match_reference(&snaps, &newicks, (&wrf, &kf), &ctx);
        }
    }
}

#[test]
fn share_override_accepts_only_finite_non_negative_numbers() {
    assert_eq!(parse_share(None, 0.25), 0.25);
    assert_eq!(parse_share(Some(" 0.1 "), 0.25), 0.1);
    assert_eq!(parse_share(Some("0"), 0.25), 0.0);
    assert_eq!(parse_share(Some("2"), 0.25), 2.0);
    for bad in ["", "quarter", "-0.1", "NaN", "inf"] {
        assert_eq!(parse_share(Some(bad), 0.25), 0.25, "{bad:?}");
    }
}

/// Each share is its default unless the environment overrides it. The test
/// sets no variable: the shares are read once per process, and tests share
/// that process, so a test that set one would leak into the others.
#[test]
fn dense_shares_follow_the_environment() {
    for (share, var, default) in [
        (
            rf_dense_share(),
            "RAPIDTREES_RF_DENSE_SHARE",
            RF_DENSE_SHARE,
        ),
        (
            weighted_dense_share(),
            "RAPIDTREES_WEIGHTED_DENSE_SHARE",
            WEIGHTED_DENSE_SHARE,
        ),
    ] {
        let expected = parse_share(std::env::var(var).ok().as_deref(), default);
        assert_eq!(share, expected, "{var}");
    }
}

/// The tiled sweep adds each pair's terms in the order [`sweep`] does, so it
/// gives the untiled formula to the bit. Seventeen trees leave partial
/// bands and tiles, and with 603 taxa the dense columns cross a panel
/// boundary and end in a tail.
#[test]
fn tiled_sweep_matches_the_untiled_formula_to_the_bit() {
    let (_, snaps) = random_snapshots(603, 13, 4, 51);
    let (counts, n) = (&snaps.split_counts, snaps.snapshots.len());
    let min_dense = min_dense(n, weighted_dense_share());
    type Overlap = fn(f64, f64) -> f64;
    type Finish = fn(f64) -> f64;
    let metrics: [(Overlap, Finish); 2] = [(f64::min, |d| d), (|a, b| a * b, f64::sqrt)];
    for (overlap, finish) in metrics {
        let trees = || Trees::Borrowed(&snaps.snapshots);
        let tiled = weighted_distances_split(trees(), counts, None, min_dense, overlap, finish);

        let layout = Layout::new(counts, n, min_dense, |count| count >= 2);
        let (dense, postings, self_total) = weighted_layout(trees(), counts, &layout, &overlap);
        assert!(
            dense.stride > PANEL && dense.stride % 8 != 0,
            "{}",
            dense.stride
        );
        let untiled = fill_symmetric(n, None, |i, row: &mut [f64]| {
            let (cols, lengths) = postings.of_tree(i);
            for (&col, &length) in cols.iter().zip(lengths) {
                let (trees, others) = postings.after(col, i);
                for (&j, &other) in trees.iter().zip(others) {
                    row[j as usize] += overlap(length, other);
                }
            }
            for (j, slot) in row.iter_mut().enumerate().skip(i + 1) {
                let shared = sweep(dense.row(i), dense.row(j), &overlap) + *slot;
                *slot = finish((self_total[i] + self_total[j] - 2.0 * shared).max(0.0));
            }
        });
        assert_eq!(tiled, untiled);
    }
}

#[test]
fn no_trees_give_an_empty_matrix() {
    let snaps = Snapshots::from_newicks(&[], false).unwrap();
    assert!(snaps.pairwise_rf(None).is_empty());
    assert!(snaps.pairwise_wrf(None).is_empty());
    assert!(snaps.pairwise_kf(None).is_empty());
}

/// Handing a collection to a kernel gives the matrix that borrowing it
/// gives, for every metric at the default boundaries: over several layout
/// blocks, and with enough trees for RF to give rare splits posting lists.
#[test]
fn consuming_kernels_match_borrowing_ones() {
    for &(n_taxa, n_trees, duplicates, seed) in
        &[(20usize, 12usize, 4usize, 41u64), (30, 300, 40, 42)]
    {
        let build = || random_snapshots(n_taxa, n_trees, duplicates, seed).1;
        let snaps = build();
        assert_eq!(
            snaps.pairwise_rf(None),
            build().into_pairwise_rf(None),
            "RF, seed={seed}"
        );
        assert_eq!(
            snaps.pairwise_wrf(None),
            build().into_pairwise_wrf(None),
            "WRF, seed={seed}"
        );
        assert_eq!(
            snaps.pairwise_kf(None),
            build().into_pairwise_kf(None),
            "KF, seed={seed}"
        );
    }
}

/// Same as [`weighted_boundary_does_not_change_distances`] for RF: all
/// posting lists (2), mixed, and all bit columns (`n + 1`) must match the
/// oracle exactly.
#[test]
fn rf_boundary_does_not_change_distances() {
    for &(n_taxa, n_trees, duplicates, seed) in &[
        (12usize, 8usize, 4usize, 31u64),
        (24, 10, 6, 32),
        (40, 6, 9, 33),
    ] {
        let (_, snaps) = random_snapshots(n_taxa, n_trees, duplicates, seed);
        let n = snaps.snapshots.len();
        for min_dense in [2, 3, n as u32 / 2, n as u32, n as u32 + 1] {
            let counts = &snaps.split_counts;
            let rf = distance_rf_split(Trees::Borrowed(&snaps.snapshots), counts, None, min_dense);
            let owned = distance_rf_split(
                Trees::Owned(snaps.snapshots.clone()),
                counts,
                None,
                min_dense,
            );
            assert_eq!(
                rf, owned,
                "borrowed and owned trees disagree at min_dense={min_dense}"
            );
            assert_rf_matches_reference(&snaps, &rf, &format!("seed={seed} min_dense={min_dense}"));
        }
    }
}

/// The same tree written with its children in a different order parses its
/// splits in a different order. Posted splits are summed in column order,
/// not parse order, so the copies must still cancel to exactly 0.0 on every
/// path.
#[test]
fn identical_trees_in_any_child_order_cancel_exactly() {
    let trees = [
        "(((A:0.31,B:0.17):0.23,(C:0.41,D:0.13):0.29):0.11,((E:0.37,F:0.19):0.07,(G:0.43,H:0.03):0.47):0.53);",
        "(((H:0.03,G:0.43):0.47,(F:0.19,E:0.37):0.07):0.53,((D:0.13,C:0.41):0.29,(B:0.17,A:0.31):0.23):0.11);",
        "(((A:0.31,E:0.12):0.33,(C:0.34,G:0.35):0.36):0.37,((B:0.38,F:0.39):0.41,(D:0.42,H:0.43):0.44):0.45);",
        "(((C:0.41,D:0.13):0.29,(B:0.17,A:0.31):0.23):0.11,((G:0.43,H:0.03):0.47,(E:0.37,F:0.19):0.07):0.53);",
    ];
    let snaps = Snapshots::from_newicks(&trees, false).unwrap();
    let n = snaps.snapshots.len();
    // 2: all dense. 4: pendants dense, the shared internal splits posted.
    // n + 1: everything posted.
    for min_dense in [2, 4, n as u32 + 1] {
        let (wrf, kf) = weighted_at(&snaps, min_dense);
        for (i, j) in [(0, 1), (0, 3), (1, 3)] {
            assert_eq!(wrf[i * n + j], 0.0, "WRF [{i}][{j}] min_dense={min_dense}");
            assert_eq!(kf[i * n + j], 0.0, "KF [{i}][{j}] min_dense={min_dense}");
        }
        assert!(wrf[2] > 0.0, "tree 2 differs, min_dense={min_dense}");
    }
}
