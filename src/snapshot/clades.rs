//! Which taxa each split names: one canonical leaf set per **distinct split**
//! in a run — not per node, not per tree.
//!
//! A split's fingerprint is a one-way summary, so `{C,D}` cannot be read back
//! out of a `u128`. This table is the readable answer, and it exists solely to
//! be exported: [`super::export`] packs it into bytes for Python and sorts
//! export columns by [`cmp_packed`]. Distances never touch it.
//!
//! The canonical side is the one without leaf 0, and indices within a clade
//! ascend — both invariants that [`cmp_packed`] and the byte packing rely on.

use std::cmp::Ordering;

/// Every distinct split's canonical leaf set, concatenated.
///
/// One flat `Vec<u32>` rather than a `Vec<Vec<u32>>`: a per-split `Vec` header
/// is 24 bytes, which at these clade sizes would cost more than the leaf
/// indices it points at.
#[derive(Debug)]
pub(crate) struct CladeTable {
    /// Ascending leaf indices, every clade back to back.
    leaves: Vec<u32>,
    /// Clade `i` is `leaves[starts[i]..starts[i + 1]]`, so this always carries
    /// one entry more than there are clades.
    starts: Vec<u32>,
}

/// Not derived: `starts` must always hold the leading `0`, and a derived
/// `Default` would leave it empty — making `len()` underflow on an otherwise
/// perfectly ordinary empty table.
impl Default for CladeTable {
    fn default() -> Self {
        Self::new()
    }
}

impl CladeTable {
    pub(crate) fn new() -> Self {
        Self {
            leaves: Vec::new(),
            starts: vec![0],
        }
    }

    /// Number of distinct splits in the table.
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.starts.len() - 1
    }

    /// Clade `i`'s leaf indices, ascending.
    #[inline]
    pub(crate) fn get(&self, i: usize) -> &[u32] {
        let (lo, hi) = (self.starts[i] as usize, self.starts[i + 1] as usize);
        &self.leaves[lo..hi]
    }

    /// Append one clade from the leaves it contains, in any order.
    pub(crate) fn push(&mut self, leaves: impl IntoIterator<Item = u32>) {
        let start = self.leaves.len();
        self.leaves.extend(leaves);
        self.leaves[start..].sort_unstable();
        self.starts.push(self.leaves.len() as u32);
    }
}

/// Order two ascending leaf-index runs exactly as their bit-packed forms would.
///
/// Packed comparison runs over `u64` words low-to-high and compares each word by
/// value, so the **lowest differing word** decides and, inside it, the **highest
/// differing bit**. Note those pull in opposite directions — leaf 63 outranks
/// leaf 64, because 64 lives in a later word that a lower word's difference
/// never reaches.
///
/// This defines export column order, which callers index against, so it must
/// match the packed form exactly rather than merely being deterministic.
pub(crate) fn cmp_packed(a: &[u32], b: &[u32]) -> Ordering {
    let Some(first_diff) = lowest_difference(a, b) else {
        return Ordering::Equal;
    };

    // Everything below `first_diff` matches, so every word below this one does.
    let word = first_diff >> 6;
    let (lo, hi) = (word << 6, (word + 1) << 6);
    let (wa, wb) = (slice_in(a, lo, hi), slice_in(b, lo, hi));

    // Within the word, the highest bit only one side has decides — which is
    // what comparing the two runs from the top down gives.
    wa.iter().rev().cmp(wb.iter().rev())
}

/// The smallest leaf index present in exactly one of two ascending runs.
fn lowest_difference(a: &[u32], b: &[u32]) -> Option<u32> {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Equal => {
                i += 1;
                j += 1;
            }
            Ordering::Less => return Some(a[i]),
            Ordering::Greater => return Some(b[j]),
        }
    }
    a.get(i).or(b.get(j)).copied()
}

/// The stretch of an ascending run lying in `lo..hi`.
fn slice_in(s: &[u32], lo: u32, hi: u32) -> &[u32] {
    let start = s.partition_point(|&x| x < lo);
    let end = s.partition_point(|&x| x < hi);
    &s[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(clades: &[&[u32]]) -> CladeTable {
        let mut t = CladeTable::new();
        for c in clades {
            t.push(c.iter().copied());
        }
        t
    }

    #[test]
    fn stores_and_returns_clades_sorted() {
        let t = table(&[&[2, 0], &[5, 1, 3]]);
        assert_eq!(t.len(), 2);
        assert_eq!(t.get(0), &[0, 2]);
        assert_eq!(t.get(1), &[1, 3, 5]);
    }

    #[test]
    fn empty_table_has_no_clades() {
        assert_eq!(CladeTable::new().len(), 0);
    }

    /// `Default` is hand-written because the derived one would leave `starts`
    /// empty, and `len()` subtracts one from it.
    #[test]
    fn default_is_an_empty_table_not_an_underflow() {
        assert_eq!(CladeTable::default().len(), 0);
    }

    /// Within one word, the highest bit only one side holds decides — so `{2}`
    /// outranks `{0,1}` even though it has fewer, smaller-looking members.
    #[test]
    fn highest_differing_bit_decides_within_a_word() {
        assert_eq!(cmp_packed(&[1], &[0, 1]), Ordering::Less);
        assert_eq!(cmp_packed(&[2], &[0, 1]), Ordering::Greater);
        assert_eq!(cmp_packed(&[0], &[1]), Ordering::Less);
        assert_eq!(cmp_packed(&[3], &[3]), Ordering::Equal);
    }

    /// The word boundary is where a "highest index wins" shortcut would break:
    /// leaf 63 sits in word 0 and leaf 64 in word 1, and word 0 is compared
    /// first, so `{63}` must outrank `{64}`.
    #[test]
    fn lower_word_outranks_higher_bit() {
        assert_eq!(cmp_packed(&[63], &[64]), Ordering::Greater);
        assert_eq!(cmp_packed(&[0, 64], &[1]), Ordering::Less);
        assert_eq!(cmp_packed(&[64], &[65]), Ordering::Less);
    }

    /// Cross-check against actually packing the bits, over every pair of subsets
    /// of a 3-word universe. This ordering pins export column order, so
    /// equivalence to the packed form has to be exhaustively true, not plausible.
    #[test]
    fn matches_bit_packed_ordering_exhaustively() {
        const N: u32 = 130;
        let interesting: Vec<u32> = vec![0, 1, 63, 64, 65, 127, 128];

        let pack = |s: &[u32]| -> Vec<u64> {
            let mut w = vec![0u64; N.div_ceil(64) as usize];
            for &i in s {
                w[(i >> 6) as usize] |= 1u64 << (i & 63);
            }
            w
        };

        // All subsets of the interesting indices, compared both ways.
        let subsets: Vec<Vec<u32>> = (0..1u32 << interesting.len())
            .map(|mask| {
                interesting
                    .iter()
                    .enumerate()
                    .filter(|(b, _)| mask >> b & 1 == 1)
                    .map(|(_, &i)| i)
                    .collect()
            })
            .collect();

        for a in &subsets {
            for b in &subsets {
                assert_eq!(
                    cmp_packed(a, b),
                    pack(a).cmp(&pack(b)),
                    "disagreed on {a:?} vs {b:?}"
                );
            }
        }
    }
}
