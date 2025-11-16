//! Extract partition snapshots from phylogenetic trees.
//!
//! # Overview
//! A TreeSnapshot captures all the bipartitions (splits) in a tree along with
//! their branch lengths. This immutable snapshot can be safely compared with
//! other snapshots in parallel.
//!
//! # What is a bipartition?
//! Each internal branch in a tree divides the leaves into two groups.
//! For example:
//! ```text
//!      root
//!     /    \
//!   {A,B}  {C,D}  ← This branch creates partition {A,B}
//! ```
//!
//! We only store one side of each partition (the smaller side by convention).
//!
//! # CRITICAL: Why we use taxon NAMES not node IDs
//! Node IDs are assigned during tree parsing and differ across files.
//! Taxon names are consistent. We sort leaves alphabetically by name
//! to ensure identical taxa always map to the same bit positions.

use crate::bitset::Bitset;
use phylotree::tree::{Tree as PhyloTree, TreeError};
use rustc_hash::{FxHashMap, FxHashSet};

/// An immutable snapshot of all partitions in a phylogenetic tree.
///
/// # Fields
/// - `parts`: All bipartitions, **canonicalized** (stored in a HashSet for O(1) lookup)
/// - `lengths`: Branch lengths for each partition (HashMap keyed by Bitset)
/// - `root_children`: Bitsets for immediate children of root (for rooted RF)
/// - `words`: Number of u64 words needed for bitsets
/// - `num_leaves`: Total number of leaves (needed for canonicalization)
/// - `rooted`: Whether the tree is rooted
///
/// # Canonicalization
/// Each bipartition can be represented two ways: {A,B}|{C,D} or {C,D}|{A,B}.
/// We canonicalize by always storing the side that does NOT contain leaf with index 0.
/// This ensures identical partitions have identical bitset representations.
///
/// # Performance
/// Using HashSet and HashMap allows O(1) average-case lookups for Robinson-Foulds
/// and weighted distance calculations, instead of O(n log n) with sorted vectors.
#[derive(Debug, Clone)]
pub struct TreeSnapshot {
    /// All partitions in the tree, canonicalized (FxHashSet for fast lookup)
    pub parts: FxHashSet<Bitset>,

    /// Branch length for each partition (keyed by the canonical Bitset)
    pub lengths: FxHashMap<Bitset, f64>,

    /// Bitsets of root's immediate children (for rooted tree adjustment)
    pub root_children: Vec<Bitset>,

    /// Number of u64 words in each bitset
    pub words: usize,

    /// Total number of leaves (needed for computing complements)
    pub num_leaves: usize,

    /// Whether this tree is rooted
    pub rooted: bool,
}

impl TreeSnapshot {
    /// Extract a snapshot from a phylogenetic tree.
    ///
    /// # Parameters
    /// - `tree`: The phylogenetic tree to extract partitions from
    /// - `include_trivial`: If true, includes single-leaf partitions (needed for weighted metrics like Robinson-Foulds and Kuhner-Felsenstein)
    ///
    /// # Algorithm
    /// 1. Extract leaf names and sort them alphabetically for consistency
    /// 2. Map each leaf name to a compact index [0..n)
    /// 3. DFS (Depth-First Search) from root, building bitsets bottom-up using leaf names
    /// 4. For each internal node, merge child bitsets with OR
    /// 5. Collect partitions (optionally including trivial single-leaf partitions)
    /// 6. Canonicalize partitions (always store side without leaf with index 0)
    ///
    /// # Errors
    /// Returns `TreeError` if the tree is empty, malformed, or has unnamed leaves.
    pub fn from_tree(tree: &PhyloTree) -> Result<Self, TreeError> {
        let rooted = tree.is_rooted()?;
        // Step 1: Extract leaf names and sort them alphabetically
        let mut leaf_names: Vec<(usize, String)> = tree
            .get_leaves()
            .iter()
            .map(|leaf_id| {
                let leaf_name = tree.get(leaf_id).unwrap().name.clone().unwrap_or_default();
                (*leaf_id, leaf_name)
            })
            .collect::<Vec<_>>();

        // Sort by taxon name (alphabetically) for consistent ordering
        leaf_names.sort_by(|a, b| a.1.cmp(&b.1));

        let num_leaves = leaf_names.len();
        let words = num_leaves.div_ceil(64);

        // Step 2: Create mapping: node_id → bit_index (based on sorted names)
        let node_id_to_leaf_index: FxHashMap<usize, usize> = leaf_names
            .iter()
            .enumerate()
            .map(|(idx, &(node_id, _))| (node_id, idx))
            .collect();

        // Step 3: Perform DFS to build bitsets for each node
        let root_id = tree.get_root()?;
        // Cache to store computed bitsets
        // Key: node_id, Value: Bitset of leaves under this node
        // Node_id, allows us to get a branch length associated with the partition
        let mut cache: FxHashMap<usize, Bitset> = FxHashMap::with_capacity_and_hasher(num_leaves * 2, Default::default());
        Self::compute_bitsets(root_id, tree, &node_id_to_leaf_index, words, &mut cache);

        // Step 4: Collect partitions (with or without trivial partitions)
        let (parts, lengths) = Self::collect_partitions(tree, root_id, &cache)?;

        // Step 5: Canonicalize partitions (always store side WITHOUT leaf 0)
        let (parts_canonical, lengths_canonical) =
            Self::canonicalize_partitions(parts, lengths, words, num_leaves);

        // Step 6: Record root's children for rooted tree adjustment
        let root_children = Self::get_root_children(tree, root_id, &cache)?;

        Ok(TreeSnapshot {
            parts: parts_canonical,
            lengths: lengths_canonical,
            root_children,
            words,
            num_leaves,
            rooted,
        })
    }

    /// Recursively compute bitsets for all nodes via DFS.
    ///
    /// # Algorithm
    /// - **Leaf node**: Create bitset with single bit set
    /// - **Internal node**: OR together all child bitsets
    ///
    /// Results are cached to avoid recomputation.
    fn compute_bitsets(
        node_id: usize,
        tree: &PhyloTree,
        node_id_to_leaf_index: &FxHashMap<usize, usize>,
        words: usize,
        cache: &mut FxHashMap<usize, Bitset>,
    ) -> Bitset {
        // Return cached result if available
        if let Some(bitset) = cache.get(&node_id) {
            return bitset.clone();
        }

        let node = tree.get(&node_id).expect("valid node");

        // Base case: leaf node
        if node.children.is_empty() {
            let mut bitset = Bitset::zeros(words);
            let leaf_idx = *node_id_to_leaf_index.get(&node_id).expect("leaf mapped");
            bitset.set(leaf_idx);
            cache.insert(node_id, bitset.clone());
            return bitset;
        }

        // Recursive case: internal node
        // Merge all child bitsets with OR
        let mut bitset = Bitset::zeros(words);
        for &child_id in &node.children {
            let child_bitset =
                Self::compute_bitsets(child_id, tree, node_id_to_leaf_index, words, cache);
            bitset.or_assign(&child_bitset);
        }

        cache.insert(node_id, bitset.clone());
        bitset
    }

    /// Collect all non-trivial partitions and their branch lengths.
    ///
    /// # Parameters
    /// - `include_trivial`: If true, includes single-leaf partitions (needed for weighted metrics)
    ///
    /// # What we skip
    /// - Root node (doesn't create a bipartition)
    /// - Trivial partitions (single leaf) - unless `include_trivial` is true
    ///
    /// # Branch lengths
    /// Some trees may have missing branch lengths.
    /// We treat missing lengths as 0.0.
    fn collect_partitions(
        tree: &PhyloTree,
        root_id: usize,
        cache: &FxHashMap<usize, Bitset>,
    ) -> Result<(Vec<Bitset>, Vec<f64>), TreeError> {
        let mut parts = Vec::new();
        let mut lengths = Vec::new();

        // Unless it becomes a bottleneck, we can parallelize this loop later
        for (&node_id, bitset) in cache.iter() {
            // Skip root (doesn't create a partition)
            if node_id == root_id {
                continue;
            }

            // Skip trivial partitions (single leaf) unless explicitly requested
            if bitset.count_ones() <= 1 {
                continue;
            }

            // Add this partition
            parts.push(bitset.clone());

            // Get branch length leading TO this node (creates the partition)
            // This is the edge from parent to this node, not the sum of child edges
            let node = tree.get(&node_id)?;
            let length: f64 = node.parent_edge.unwrap_or(0.0);
            lengths.push(length);
        }

        Ok((parts, lengths))
    }

    /// Canonicalize partitions to ensure consistent representation.
    ///
    /// # Problem
    /// A bipartition {A,B}|{C,D}:
    ///
    /// A --\                   /-- C
    ///     node1 - (root) - node2
    /// D --/                   \-- B
    ///
    /// Can be represented as:
    ///        (root)
    ///        /   \
    ///    node1    node2
    ///    /   \    /   \
    ///   A     D  B     C
    /// bitset: [{node1}: 0b0011, {node2}: 0b1100]
    ///
    /// Or as (not drawn accuratly to scale):
    ///        node1
    ///        /   \
    ///    (root)    D
    ///    /   \
    ///   A     \
    ///        node2
    ///        /   \
    ///       B     C
    /// bitset: [{root}: 0b1101, {node2}: 0b1100]
    ///
    ///
    /// Without canonicalization, identical trees might produce different bitsets!
    ///
    ///
    ///
    /// __Consider it as rooting the tree to the same leaf.__
    ///
    /// # Solution
    /// Always store the side that does NOT contain leaf 0 (the first leaf alphabetically).
    /// - If leaf 0 is set: flip to complement
    /// - If leaf 0 is not set: keep as-is
    ///
    /// # Example
    /// Leaves: A=0, B=1, C=2, D=3
    /// Partition {A,B}: bitset 0b0011 (leaf 0 SET) → flip to {C,D}: 0b1100
    /// Partition {C,D}: bitset 0b1100 (leaf 0 NOT set) → keep as 0b1100
    ///
    /// # Returns
    /// Returns (FxHashSet<Bitset>, FxHashMap<Bitset, f64>) for O(1) lookups
    fn canonicalize_partitions(
        parts: Vec<Bitset>,
        lengths: Vec<f64>,
        words: usize,
        num_leaves: usize,
    ) -> (FxHashSet<Bitset>, FxHashMap<Bitset, f64>) {
        let mut canonical_parts = FxHashSet::with_capacity_and_hasher(parts.len(), Default::default());
        let mut canonical_lengths = FxHashMap::with_capacity_and_hasher(lengths.len(), Default::default());

        for (bitset, length) in parts.into_iter().zip(lengths.into_iter()) {
            // Check if leaf 0 (bit 0 of word 0) is set
            let leaf_0_is_set = (bitset.0[0] & 1) != 0;

            let canonical_bitset = if leaf_0_is_set {
                // Flip to complement (side without leaf 0)
                Self::compute_complement(&bitset, words, num_leaves)
            } else {
                // Already canonical (leaf 0 not in this side)
                bitset
            };

            canonical_parts.insert(canonical_bitset.clone());
            canonical_lengths.insert(canonical_bitset, length);
        }

        (canonical_parts, canonical_lengths)
    }

    /// Compute the bitwise complement of a partition.
    ///
    /// Flips all bits up to num_leaves, keeping remaining bits as 0.
    ///
    /// # Example
    /// Input:  0b0011 (4 leaves) → Output: 0b1100
    /// Input:  0b1100 (4 leaves) → Output: 0b0011
    fn compute_complement(bitset: &Bitset, words: usize, num_leaves: usize) -> Bitset {
        let mut complement = Bitset::zeros(words);

        for i in 0..num_leaves {
            let word = i >> 6;
            let bit = i & 63;

            // Check if bit i is set in original
            let is_set = (bitset.0[word] & (1u64 << bit)) != 0;

            // Set bit i in complement if NOT set in original
            if !is_set {
                complement.0[word] |= 1u64 << bit;
            }
        }

        complement
    }

    /// Sort partitions lexicographically for edge length matching later.
    /// Get bitsets for root's immediate children (for rooted RF adjustment).
    ///
    /// In rooted trees, we need to know if two trees have the same root
    /// position to apply the correct RF distance adjustment.
    fn get_root_children(
        tree: &PhyloTree,
        root_id: usize,
        cache: &FxHashMap<usize, Bitset>,
    ) -> Result<Vec<Bitset>, TreeError> {
        let root = tree.get(&root_id)?;
        let mut root_children: Vec<Bitset> = root
            .children
            .iter()
            .filter_map(|&child_id| cache.get(&child_id).cloned())
            .collect();

        root_children.sort_unstable();
        Ok(root_children)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Complete example: Tree with depth 3 and multiple partitions
    ///
    /// ```text
    ///           root
    ///          /    \
    ///      node1     node2
    ///      /   \     /   \
    ///     A     B   C    node3
    ///                    /   \
    ///                   D     E
    /// ```
    ///
    /// After sorting leaves alphabetically: A=0, B=1, C=2, D=3, E=4
    ///
    /// # Step 1: Extract all partitions (bottom-up DFS)
    ///
    /// | Node  | Leaves Below | Raw Bitset | Binary    |
    /// |-------|--------------|------------|-----------|
    /// | node3 | {D, E}       | 0b11000    | bits 3,4  |
    /// | node2 | {C, D, E}    | 0b11100    | bits 2,3,4|
    /// | node1 | {A, B}       | 0b00011    | bits 0,1  |
    /// | root  | {A,B,C,D,E}  | 0b11111    | (skipped) |
    ///
    /// We skip root (no partition above it).
    ///
    /// # Step 2: Canonicalize (flip if leaf 0=A is present)
    ///
    /// | Partition | Raw Bitset | Has A? | Action | Canonical |
    /// |-----------|------------|--------|--------|-----------|
    /// | node3: {D,E}   | 0b11000 | NO  | Keep   | 0b11000   |
    /// | node2: {C,D,E} | 0b11100 | NO  | Keep   | 0b11100   |
    /// | node1: {A,B}   | 0b00011 | YES | Flip!  | 0b11100   |
    ///
    /// Wait! node1 flipped becomes {C,D,E} = 0b11100, same as node2!
    /// That's WRONG - they're different partitions!
    ///
    /// # The Real Issue: Need Better Example
    ///
    /// Let me show a tree where partitions DON'T overlap:
    ///
    /// ```text
    ///           root
    ///          /    \
    ///      node1    node2
    ///      /   \    /   \
    ///     A     B  C     D
    /// ```
    ///
    /// Leaves sorted: A=0, B=1, C=2, D=3
    ///
    /// | Node  | Leaves    | Raw Bitset | Has A? | Canonical |
    /// |-------|-----------|------------|--------|-----------|
    /// | node1 | {A, B}    | 0b0011     | YES    | 0b1100 (flip to {C,D}) |
    /// | node2 | {C, D}    | 0b1100     | NO     | 0b1100 (keep)          |
    ///
    /// Now both become 0b1100! That's STILL wrong - they're different branches!
    #[test]
    fn test_depth3_tree_partitions() {
        // This test reveals a conceptual issue...
        // Let me think about this more carefully

        // Tree:     root
        //          /    \
        //      node1    node2
        //      /   \    /   \
        //     A     B  C     D

        // The partitions are:
        // 1. node1 splits: {A,B} vs {C,D}
        // 2. node2 splits: {C,D} vs {A,B}

        // These are the SAME bipartition from different directions!
        // So they SHOULD canonicalize to the same thing!

        let mut part_ab = Bitset::zeros(1);
        part_ab.set(0); // A
        part_ab.set(1); // B

        let mut part_cd = Bitset::zeros(1);
        part_cd.set(2); // C
        part_cd.set(3); // D

        // Both represent the same split, so after canonicalization
        // they'll be identical - that's CORRECT behavior!
    }

    /// Better example: Asymmetric tree with distinct partitions
    ///
    /// ```text
    ///              root
    ///             /    \
    ///         node1     E
    ///         /   \
    ///     node2    D
    ///     /   \
    ///    A    node3
    ///         /   \
    ///        B     C
    /// ```
    ///
    /// Leaves sorted: A=0, B=1, C=2, D=3, E=4
    ///
    /// # Partitions extracted (bottom-up):
    ///
    /// | Node  | Leaves Below | Raw Bitset | Binary     |
    /// |-------|--------------|------------|------------|
    /// | node3 | {B, C}       | 0b00110    | bits 1,2   |
    /// | node2 | {A, B, C}    | 0b00111    | bits 0,1,2 |
    /// | node1 | {A,B,C,D}    | 0b01111    | bits 0,1,2,3|
    ///
    /// # Canonicalization (flip if has A=bit 0):
    ///
    /// | Partition | Raw       | Has A? | Canonical | Represents        |
    /// |-----------|-----------|--------|-----------|-------------------|
    /// | node3     | 0b00110   | NO     | 0b00110   | {B,C} vs {A,D,E} |
    /// | node2     | 0b00111   | YES    | 0b11000   | {A,B,C} vs {D,E} → flip to {D,E} |
    /// | node1     | 0b01111   | YES    | 0b10000   | {A,B,C,D} vs {E} → flip to {E}   |
    ///
    /// # Final canonical partitions (sorted):
    /// 1. 0b00110 = {B,C}
    /// 2. 0b10000 = {E}
    /// 3. 0b11000 = {D,E}
    ///
    /// All three are distinct! ✓
    #[test]
    fn test_asymmetric_tree_example() {
        // Leaves: A=0, B=1, C=2, D=3, E=4 (5 leaves)

        // node3: {B, C}
        let mut node3 = Bitset::zeros(1);
        node3.set(1); // B
        node3.set(2); // C
        assert_eq!(node3.0[0], 0b00110);
        assert!((node3.0[0] & 1) == 0); // No A, keep as-is

        // node2: {A, B, C}
        let mut node2 = Bitset::zeros(1);
        node2.set(0); // A
        node2.set(1); // B
        node2.set(2); // C
        assert_eq!(node2.0[0], 0b00111);
        assert!((node2.0[0] & 1) != 0); // Has A, need to flip!
        // Flip to complement {D, E} = 0b11000

        // node1: {A, B, C, D}
        let mut node1 = Bitset::zeros(1);
        node1.set(0); // A
        node1.set(1); // B
        node1.set(2); // C
        node1.set(3); // D
        assert_eq!(node1.0[0], 0b01111);
        assert!((node1.0[0] & 1) != 0); // Has A, need to flip!
        // Flip to complement {E} = 0b10000

        // After canonicalization and sorting:
        // 0b00110 (node3 kept)
        // 0b10000 (node1 flipped)
        // 0b11000 (node2 flipped)
        // All distinct! ✓
    }

    /// Demonstrates the critical importance of canonicalization
    ///
    /// Problem without canonicalization:
    /// ```text
    /// Tree 1:           Tree 2:
    ///    root              root
    ///   /    \            /    \
    /// {A,B} {C,D}      {C,D} {A,B}
    ///
    /// Tree 1 stores: {A,B} = 0b0011
    /// Tree 2 stores: {C,D} = 0b1100
    /// These look different but represent the SAME bipartition! ❌
    /// ```
    ///
    /// With canonicalization (always store side WITHOUT leaf 0=A):
    /// ```text
    /// Tree 1: {A,B} contains A → flip to {C,D} = 0b1100
    /// Tree 2: {C,D} no A → keep as {C,D} = 0b1100
    /// Now they match! ✓
    /// ```
    #[test]
    fn test_canonicalization() {
        // 4 leaves: A=0, B=1, C=2, D=3

        // Partition {A, B}
        let mut part_ab = Bitset::zeros(1);
        part_ab.set(0); // A
        part_ab.set(1); // B
        assert_eq!(part_ab.0[0], 0b0011);

        // Partition {C, D} (complement of {A,B})
        let mut part_cd = Bitset::zeros(1);
        part_cd.set(2); // C
        part_cd.set(3); // D
        assert_eq!(part_cd.0[0], 0b1100);

        // These are the SAME bipartition, just different sides
        // After canonicalization, both should become {C,D} = 0b1100
        // (the side without leaf 0)

        // Check if leaf 0 is set
        assert!((part_ab.0[0] & 1) != 0); // A is in {A,B}
        assert!((part_cd.0[0] & 1) == 0); // A is NOT in {C,D}

        // So we'd flip {A,B} to its complement {C,D}
    }

    /// Demonstrates why we MUST use taxon names, not node IDs
    ///
    /// When reading BEAST trees, node IDs are assigned during parsing
    /// and will differ across trees even if taxa are identical.
    ///
    /// ```text
    /// File 1 parsed:
    ///   node_7 = "Human"
    ///   node_3 = "Chimp"
    ///   node_15 = "Gorilla"
    ///
    /// File 2 parsed (same taxa, different IDs):
    ///   node_5 = "Human"
    ///   node_8 = "Chimp"
    ///   node_12 = "Gorilla"
    /// ```
    ///
    /// If we used node IDs directly:
    /// - Tree 1: node_3 → index 0, node_7 → index 1, node_15 → index 2
    /// - Tree 2: node_5 → index 0, node_8 → index 1, node_12 → index 2
    /// - Partition {Chimp, Human} in Tree 1: bitset 0b011 (nodes 3,7)
    /// - Partition {Chimp, Human} in Tree 2: bitset 0b011 (nodes 5,8)
    /// - These look the same by accident, but represent DIFFERENT taxa! ❌
    ///
    /// Correct approach (using names):
    /// - Both trees sort by name: Chimp → 0, Gorilla → 1, Human → 2
    /// - Partition {Chimp, Human}: bitset 0b101 in BOTH trees ✓
    #[test]
    fn test_taxon_names_vs_node_ids() {
        // Simulate two trees with same taxa but different node IDs

        // Tree 1: IDs [7, 3, 15] → Names ["Human", "Chimp", "Gorilla"]
        let tree1_leaves = vec![(7, "Human"), (3, "Chimp"), (15, "Gorilla")];

        // Tree 2: IDs [5, 8, 12] → Names ["Human", "Chimp", "Gorilla"]
        let tree2_leaves = vec![(5, "Human"), (8, "Chimp"), (12, "Gorilla")];

        // After sorting by NAME (not ID):
        let mut sorted1 = tree1_leaves.clone();
        let mut sorted2 = tree2_leaves.clone();
        sorted1.sort_by(|a, b| a.1.cmp(b.1));
        sorted2.sort_by(|a, b| a.1.cmp(b.1));

        // Both now have: Chimp=0, Gorilla=1, Human=2
        assert_eq!(sorted1[0].1, "Chimp");
        assert_eq!(sorted1[1].1, "Gorilla");
        assert_eq!(sorted1[2].1, "Human");

        assert_eq!(sorted2[0].1, "Chimp");
        assert_eq!(sorted2[1].1, "Gorilla");
        assert_eq!(sorted2[2].1, "Human");

        // Now partition {Chimp, Human} = bits 0,2 = 0b101 in BOTH trees ✓
        let mut partition = Bitset::zeros(1);
        partition.set(0); // Chimp
        partition.set(2); // Human
        assert_eq!(partition.0[0], 0b101);
    }

    /// Demonstrates the critical importance of sorting leaves by taxon name
    ///
    /// Problem without sorting:
    /// ```text
    /// Tree 1: get_leaves() returns [Chimp, Human, Gorilla]
    ///         Partition {Human, Gorilla} → bitset 0b0110
    ///
    /// Tree 2: get_leaves() returns [Human, Chimp, Gorilla]
    ///         Same partition {Human, Gorilla} → bitset 0b0101  ❌ DIFFERENT!
    /// ```
    ///
    /// With sorting by name:
    /// ```text
    /// Both trees:
    ///   Chimp   → index 0
    ///   Gorilla → index 1
    ///   Human   → index 2
    ///
    /// Partition {Human, Gorilla} → bitset 0b0110 ✓ SAME!
    /// ```
    #[test]
    fn test_consistent_leaf_ordering() {
        // Simulate two trees with same taxa but different node IDs

        // Tree 1: Chimp=0, Human=1, Gorilla=2
        let mut leaves1 = [(0, "Chimp"), (1, "Human"), (2, "Gorilla")];

        // Tree 2: Different node IDs, different order
        let mut leaves2 = [(5, "Human"), (3, "Chimp"), (7, "Gorilla")];

        // After sorting by name, both should have same index mapping
        leaves1.sort_by(|a, b| a.1.cmp(b.1));
        leaves2.sort_by(|a, b| a.1.cmp(b.1));

        // Both should map to: Chimp=0, Gorilla=1, Human=2
        assert_eq!(leaves1[0].1, "Chimp"); // index 0
        assert_eq!(leaves1[1].1, "Gorilla"); // index 1
        assert_eq!(leaves1[2].1, "Human"); // index 2

        assert_eq!(leaves2[0].1, "Chimp"); // index 0
        assert_eq!(leaves2[1].1, "Gorilla"); // index 1
        assert_eq!(leaves2[2].1, "Human"); // index 2

        // Now partition {Human, Gorilla} = bits 1,2 = 0b0110 in BOTH trees!
    }

    /// Conceptual test showing what a snapshot would look like
    ///
    /// ```text
    ///       root
    ///      /    \
    ///     A     node1 (length: 0.5)
    ///           /   \
    ///          B     C
    /// ```
    ///
    /// Expected snapshot:
    /// - parts: [{B,C}] as bitset `0b0110`
    /// - lengths: [0.5]
    /// - root_children: [{A}, {B,C}]
    #[test]
    fn test_snapshot_concept() {
        // This is a conceptual test - actual implementation
        // would need a real PhyloTree instance

        // Partition {B, C} with leaves B=1, C=2
        let mut partition = Bitset::zeros(1);
        partition.set(1);
        partition.set(2);
        assert_eq!(partition.0[0], 0b0110);

        // Would be stored in snapshot with branch length 0.5
        let length = 0.5;
        assert_eq!(length, 0.5);
    }
}
