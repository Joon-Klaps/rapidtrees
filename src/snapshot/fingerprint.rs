//! The run-wide identity tables: what makes a split in one tree comparable to
//! a split in another.
//!
//! Both tables are derived once per run from the shared leaf set and borrowed
//! by every tree. Deriving either per tree would be redundant — all trees in a
//! collection carry the same taxa — and would make fingerprints incomparable
//! if the derivation ever disagreed.

use rustc_hash::FxHashMap;

/// A 128-bit XOR fingerprint of a leaf set.
///
/// Two distinct splits share a fingerprint with probability about `e² / 2¹²⁹`,
/// where `e` is the number of distinct splits the run sees.
pub(crate) type Fingerprint = u128;

/// One random 128-bit label per taxon, drawn once and shared by every tree.
///
/// The seed is fixed.
pub(super) fn taxon_labels(num_leaves: usize) -> Vec<Fingerprint> {
    // splitmix64, seeded from the digits of pi.
    let mut state = 0x243F_6A88_85A3_08D3u64;
    let mut next = move || {
        // the known-good splitmix64 recipe
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    };
    (0..num_leaves)
        .map(|_| Fingerprint::from(next()) << 64 | Fingerprint::from(next()))
        .collect()
}

/// Map each taxon name to its bit position — its index in the alphabetically
/// sorted leaf set.
///
/// Built once per run and borrowed by every tree. All trees in a collection
/// share a leaf set, so the alphabetical order is a property of the run, not of
/// any one tree: resolving it per tree would re-clone and re-sort the same names
/// once for every tree in the file.
pub(super) fn build_leaf_index(sorted_leaf_names: &[String]) -> FxHashMap<&str, usize> {
    sorted_leaf_names
        .iter()
        .enumerate()
        .map(|(bit, name)| (name.as_str(), bit))
        .collect()
}
