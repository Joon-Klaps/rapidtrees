//! Compact bitset representation for leaf sets in phylogenetic trees.
//!
//! # Overview
//! A bitset is an efficient way to represent which leaves belong to a tree partition.
//! Each bit position corresponds to a leaf index.
//!
//! # Example
//! For a tree with leaves [A, B, C, D] mapped to indices [0, 1, 2, 3]:
//! - Partition {A, C} → bitset `0b0101` (bits 0 and 2 set)
//! - Partition {B, C, D} → bitset `0b1110` (bits 1, 2, 3 set)

/// A compact bitset for representing which leaves belong to a partition.
///
/// Internally stores bits in `Vec<u64>` words to support arbitrarily large trees.
/// Each u64 word holds 64 leaf indices.
///
/// # Memory efficiency
/// - Traditional HashSet<usize>: ~24 bytes per element + overhead
/// - Bitset: 1 bit per possible element (8 bytes per 64 leaves)
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Bitset(pub Vec<u64>);

impl Bitset {
    /// Creates a new bitset with all bits set to 0.
    ///
    /// # Parameters
    /// - `words`: Number of u64 words needed. Calculate as `(num_leaves + 63) / 64`
    ///
    /// # Example
    /// ```
    /// # use rapidtrees::bitset::Bitset;
    /// // For a tree with 100 leaves, need 2 words (128 bits)
    /// let bs = Bitset::zeros(2);
    /// assert_eq!(bs.0.len(), 2);
    /// ```
    pub fn zeros(words: usize) -> Self {
        Bitset(vec![0u64; words])
    }

    /// Sets the bit at the given index to 1.
    ///
    /// Marks a leaf as present in this partition.
    ///
    /// # Parameters
    /// - `idx`: The leaf index to mark as present
    ///
    /// # Example
    /// ```
    /// # use rapidtrees::bitset::Bitset;
    /// let mut bs = Bitset::zeros(1);
    /// bs.set(0);  // Mark leaf 0 as present
    /// bs.set(5);  // Mark leaf 5 as present
    /// assert_eq!(bs.0[0], 0b00100001);
    /// ```
    #[inline]
    pub fn set(&mut self, idx: usize) {
        let word = idx >> 6; // Equivalent to idx / 64
        let bit = idx & 63; // Equivalent to idx % 64
        self.0[word] |= 1u64 << bit;
    }

    /// Performs bitwise OR with another bitset (union operation).
    ///
    /// Merges two leaf sets: `self` becomes `self ∪ other`
    ///
    /// # Example
    /// ```
    /// # use rapidtrees::bitset::Bitset;
    /// let mut left = Bitset::zeros(1);
    /// left.set(0);   // {0}
    ///
    /// let mut right = Bitset::zeros(1);
    /// right.set(1);  // {1}
    ///
    /// left.or_assign(&right);  // {0} ∪ {1} = {0, 1}
    /// assert_eq!(left.0[0], 0b11);
    /// ```
    #[inline]
    pub fn or_assign(&mut self, other: &Bitset) {
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            *a |= *b;
        }
    }

    /// Counts the number of set bits (population count).
    ///
    /// Returns how many leaves are in this partition.
    ///
    /// # Example
    /// ```
    /// # use rapidtrees::bitset::Bitset;
    /// let mut bs = Bitset::zeros(1);
    /// bs.set(0);
    /// bs.set(2);
    /// bs.set(5);
    /// assert_eq!(bs.count_ones(), 3);
    /// ```
    #[inline]
    pub fn count_ones(&self) -> usize {
        self.0.iter().map(|w| w.count_ones() as usize).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitset_basic() {
        let mut bs = Bitset::zeros(1);
        bs.set(0);
        bs.set(2);
        assert_eq!(bs.0[0], 0b0101);
    }

    #[test]
    fn test_bitset_or() {
        let mut bs1 = Bitset::zeros(1);
        bs1.set(0);
        bs1.set(1);

        let mut bs2 = Bitset::zeros(1);
        bs2.set(2);
        bs2.set(3);

        bs1.or_assign(&bs2);
        assert_eq!(bs1.0[0], 0b1111);
    }

    #[test]
    fn test_count_ones() {
        let mut bs = Bitset::zeros(1);
        bs.set(0);
        bs.set(2);
        bs.set(5);
        assert_eq!(bs.count_ones(), 3);
    }

    /// Visual example: How bitsets represent a small tree
    ///
    /// ```text
    ///           root
    ///          /    \
    ///        node1   D
    ///        /   \
    ///       A    node2
    ///            /   \
    ///           B     C
    /// ```
    ///
    /// Leaf mapping: A=0, B=1, C=2, D=3
    ///
    /// Partitions:
    /// - node2: {B, C} → `0b0110`
    /// - node1: {A, B, C} → `0b0111`
    #[test]
    fn test_mini_tree_example() {
        // Partition for node2: {B, C}
        let mut node2 = Bitset::zeros(1);
        node2.set(1); // B
        node2.set(2); // C
        assert_eq!(node2.0[0], 0b0110);
        assert_eq!(node2.count_ones(), 2);

        // Partition for node1: {A} ∪ {B, C}
        let mut node1 = Bitset::zeros(1);
        node1.set(0); // A
        node1.or_assign(&node2);
        assert_eq!(node1.0[0], 0b0111);
        assert_eq!(node1.count_ones(), 3);
    }

    #[test]
    fn test_large_tree() {
        // Test with more than 64 leaves (multiple words)
        let mut bs = Bitset::zeros(2);
        bs.set(0); // First word
        bs.set(63); // Last bit of first word
        bs.set(64); // First bit of second word
        bs.set(127); // Last bit of second word

        assert_eq!(bs.count_ones(), 4);
        assert_eq!(bs.0[0], 1u64 | (1u64 << 63));
        assert_eq!(bs.0[1], 1u64 | (1u64 << 63));
    }
}
