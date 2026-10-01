//! Where a split stops getting a dense column and gets a posting list.

use std::sync::OnceLock;

/// Share of the trees a split must be held by to keep its bit column in RF.
/// Rarer splits get a posting list instead.
pub(super) const RF_DENSE_SHARE: f64 = 0.03;
pub(super) const WEIGHTED_DENSE_SHARE: f64 = 0.25;

/// `default`, unless `value` is a finite, non-negative number.
///
pub(super) fn parse_share(value: Option<&str>, default: f64) -> f64 {
    value
        .and_then(|v| v.trim().parse::<f64>().ok())
        .filter(|share| share.is_finite() && *share >= 0.0)
        .unwrap_or(default)
}

/// [`RF_DENSE_SHARE`], or `RAPIDTREES_RF_DENSE_SHARE` when that is set. Read
/// once per process. The variable exists to re-tune the boundary on other data
/// or hardware; it changes the speed and never a distance.
pub(super) fn rf_dense_share() -> f64 {
    static SHARE: OnceLock<f64> = OnceLock::new();
    *SHARE.get_or_init(|| {
        parse_share(
            std::env::var("RAPIDTREES_RF_DENSE_SHARE").ok().as_deref(),
            RF_DENSE_SHARE,
        )
    })
}

/// [`WEIGHTED_DENSE_SHARE`], or `RAPIDTREES_WEIGHTED_DENSE_SHARE` when that is
/// set. See [`rf_dense_share`].
pub(super) fn weighted_dense_share() -> f64 {
    static SHARE: OnceLock<f64> = OnceLock::new();
    *SHARE.get_or_init(|| {
        parse_share(
            std::env::var("RAPIDTREES_WEIGHTED_DENSE_SHARE")
                .ok()
                .as_deref(),
            WEIGHTED_DENSE_SHARE,
        )
    })
}

/// The fewest trees that must hold a split for it to keep a dense column:
/// `share` of the `n_trees`, and never fewer than two, since a split held by
/// one tree is never shared.
pub(super) fn min_dense(n_trees: usize, share: f64) -> u32 {
    ((n_trees as f64 * share).ceil() as u32).max(2)
}
