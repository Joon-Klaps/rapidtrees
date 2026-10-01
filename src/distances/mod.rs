//! Pairwise tree distance backends for [`crate::snapshot::Snapshots`].
//!
//! All three metrics have the form `selfᵢ + selfⱼ − 2·shared`, where `shared` is
//! a popcount (RF), a running minimum (WRF), or a dot product (KF). That shape
//! lets each drop the splits that cannot affect `shared` before the O(n²) sweep.
//! A split held by only one tree is never shared, so every metric drops it; RF
//! also drops the splits held by every tree, which cancel. Both filters are pure
//! optimisations — disabling either changes no distance.
//!
//! Every metric splits its columns by how many trees hold them. A widely held
//! split keeps a dense column: a bit in RF's packed rows, an `f64` in the
//! weighted metrics' rows, swept word by word for every pair. A rarely held one
//! keeps a posting list of the trees that hold it instead, and each row adds
//! those shared terms straight into its output. A posting list of `k` trees
//! costs `k²/2` additions — HashRF's bucket — which is why only the rare splits
//! get one. A dense column costs every pair the same whether or not either
//! tree holds the split, so on a posterior, where most distinct shared splits
//! are held by a handful of trees, the lists take most of the width out of the
//! sweep. RF's bit columns are 64 to a word, so a split has to be much rarer to
//! be worth a list there (`RF_DENSE_SHARE`) than in the weighted metrics
//! (`WEIGHTED_DENSE_SHARE`).
//!
//! Seven trees, columns in descending holder count as `assign_columns` orders
//! them (`■` = the tree holds the split). The cutoff is `min_dense = 4` here
//! for illustration:
//!
//! ```text
//!             dense   │ posting │ dropped
//!             D1  D2  │ S5  S6  │ U1  U2
//! held by      7   5  │  3   2  │  1   1
//! ────────────────────┼─────────┼────────
//!     T0       ■   ■  │  ·   ·  │  ·   ·
//!     T1       ■   ■  │  ■   ·  │  ·   ·
//!     T2       ■   ■  │  ·   ·  │  ■   ·
//!     T3       ■   ·  │  ·   ■  │  ·   ·
//!     T4       ■   ■  │  ■   ·  │  ·   ·
//!     T5       ■   ·  │  ·   ·  │  ·   ■
//!     T6       ■   ■  │  ■   ■  │  ·   ·
//!                     ▲         ▲
//!                     │         └ held by ≥ 2: can be shared
//!                     └ held by ≥ min_dense
//! ```
//!
//! Dense columns are swept for every pair. `S5` becomes the posting list
//! `[T1, T4, T6]` and costs three pair updates; `S6` becomes `[T3, T6]` and
//! costs one. `U1` and `U2` fold into their tree's `self` term. RF lays its
//! columns out the same way with its own cutoff, and also drops `D1`, since a
//! split held by every tree cancels.
//!
//! There is no per-pair entry point. To compare two trees, build a two-tree
//! `Snapshots` and read the off-diagonal cell.
//!
//! # Layout of this module
//!
//! - `boundary`: the dense/posting cutoffs and their environment overrides.
//! - `layout`: the column map shared by both halves, and the trees a kernel reads.
//! - `postings`: posting lists for the rarely held splits.
//! - `matrix`: filling the symmetric output a row (or band of rows) at a time.
//! - `rf`: Robinson–Foulds over packed bit rows.
//! - `weighted`: weighted RF and Kuhner–Felsenstein over dense `f64` rows.

mod boundary;
mod layout;
mod matrix;
mod postings;
mod rf;
mod weighted;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod treedist;

pub(crate) use rf::{distance_rf, distance_rf_owned};
pub(crate) use weighted::{distance_kf, distance_kf_owned, distance_wrf, distance_wrf_owned};
