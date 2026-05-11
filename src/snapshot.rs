//! Tree snapshot types and bulk-collection constructors.
//!
//! # Types
//! - [`Snapshot`]: one tree's bipartitions as sorted `Vec<Bitset>` + branch lengths.
//!   Can be compared directly with [`crate::distances::rf_distance`] etc.
//! - [`Snapshots`]: a collection of trees in an interned split-ID representation for
//!   fast, cache-friendly pairwise distance computation at scale.
//!
//! # What is a bipartition?
//! Each internal branch in a tree divides the leaves into two groups.
//! ```text
//!      root
//!     /    \
//!   {A,B}  {C,D}  ← This branch creates split {A,B}|{C,D}
//! ```
//! We store one canonical side per split (the side not containing leaf 0).
//!
//! # Why taxon names, not node IDs
//! Node IDs are assigned at parse time and differ across files.
//! Leaves are sorted alphabetically by name so identical taxa always map to
//! the same bit positions across all trees in a dataset.
//!
//! # Bulk interning
//!
//! `Snapshots::from_newick_iter` builds per-tree `Vec<Bitset>` via DFS, then
//! assigns each unique bipartition a monotonically-increasing `u32` ID.
//! Subsequent pairwise distance computation works on `Vec<u32>` — each
//! comparison is a single integer compare — and the working set of a large
//! analysis typically fits in L2 cache instead of requiring DRAM.

use crate::bitset::Bitset;
use phylotree::tree::{Tree as PhyloTree, TreeError};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};

use std::cell::RefCell;

thread_local! {
    /// Per-thread FxHashMap reused across `Snapshot::from_tree` calls.
    ///
    /// Building snapshots in a tight loop (sequential or rayon-parallel) was
    /// allocating a fresh FxHashMap per tree — at 5000 trees that's 5000
    /// hashmap allocations + drops, plus their internal node arrays. Each
    /// rayon worker keeps its own TLS slot, so this avoids cross-thread
    /// synchronization while still reusing capacity across calls.
    static BITSET_CACHE: RefCell<FxHashMap<usize, Bitset>> =
        RefCell::new(FxHashMap::default());
}

/// An immutable snapshot of all bipartitions in one phylogenetic tree.
///
/// Stores sorted `Vec<Bitset>` + parallel `Vec<f64>` branch lengths so any
/// two `Snapshot`s with the same leaf set can be directly compared via
/// [`crate::distances::rf_distance`], [`crate::distances::wrf_distance`], or
/// [`crate::distances::kf_distance`].
///
/// Passing snapshots with different leaf sets to a distance function will
/// panic — that is a programming error.
///
/// # Fields
/// - `parts`: all edges (internal bipartitions + pendant edges), canonicalized
///   and sorted ascending for O(m+n) merge
/// - `lengths`: branch lengths parallel to `parts`
/// - `leaf_names`: alphabetically sorted taxon names (defines the bit ordering)
/// - `words`: number of u64 words per bitset
///
/// # Canonicalization
/// Each internal bipartition can be represented as two complementary bitsets.
/// We always store the side that does NOT contain the first leaf (index 0)
/// so identical splits always have identical bitset representations.
/// Pendant edges are stored as-is (single-bit bitsets, no flip needed).
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// All edges (internal bipartitions + pendant edges), canonicalized and sorted ascending.
    pub(crate) parts: Vec<Bitset>,

    /// Branch lengths, parallel to `parts`.
    pub(crate) lengths: Vec<f64>,

    /// Alphabetically sorted taxon names — defines the bit ordering.
    pub leaf_names: Vec<String>,

    /// Number of u64 words per bitset.
    pub(crate) words: usize,

    /// Whether the tree was treated as rooted when snapshotted.
    pub rooted: bool,
}

impl Snapshot {
    /// Parse a Newick string and extract a snapshot of its bipartitions.
    ///
    /// The tree must be a plain Newick string (no BEAST annotations). For
    /// BEAST-format files use `Snapshots::from_newick_iter` which handles
    /// annotation stripping and taxon renaming automatically.
    ///
    /// # Errors
    /// Returns an error string if the string cannot be parsed or the tree is malformed.
    pub fn from_newick(newick: &str, rooted: bool) -> Result<Self, String> {
        let tree = PhyloTree::from_newick(newick).map_err(|e| e.to_string())?;
        Snapshot::from_tree(&tree, rooted).map_err(|e| e.to_string())
    }

    /// Extract a snapshot from an already-parsed tree.
    ///
    /// # Algorithm
    /// 1. Extract leaf names and sort them alphabetically for consistency
    /// 2. Map each leaf name to a compact index [0..n)
    /// 3. DFS (Depth-First Search) from root, building bitsets bottom-up
    /// 4. For each internal node, merge child bitsets with OR
    /// 5. Collect partitions excluding trivial single-leaf splits
    /// 6. Canonicalize partitions (always store side without leaf with index 0)
    ///
    /// # Errors
    /// Returns `TreeError` if the tree is empty, malformed, or has unnamed leaves.
    pub(crate) fn from_tree(tree: &PhyloTree, rooted_mode: bool) -> Result<Self, TreeError> {
        let rooted = tree.is_rooted()?;
        // Step 1: Extract leaf names and sort them alphabetically
        let mut leaf_id_names: Vec<(usize, String)> = tree
            .get_leaves()
            .iter()
            .map(|leaf_id| -> Result<_, TreeError> {
                let name = tree.get(leaf_id)?.name.clone().unwrap_or_default();
                Ok((*leaf_id, name))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Sort by taxon name (alphabetically) for consistent ordering
        leaf_id_names.sort_unstable_by(|a, b| a.1.cmp(&b.1));

        let sorted_leaf_names: Vec<String> = leaf_id_names.iter().map(|(_, n)| n.clone()).collect();
        let num_leaves = leaf_id_names.len();
        let words = num_leaves.div_ceil(64);

        // Step 2: Create mapping: node_id → bit_index (based on sorted names)
        let node_id_to_leaf_index: FxHashMap<usize, usize> = leaf_id_names
            .iter()
            .enumerate()
            .map(|(idx, &(node_id, _))| (node_id, idx))
            .collect();

        // Step 3: Perform DFS to build bitsets for each node.
        // Cache is borrowed from a thread-local pool: clear (drops accumulated
        // entries) and reserve (ensure capacity matches this tree's node count).
        let root_id = tree.get_root()?;
        BITSET_CACHE.with(|cell| -> Result<Snapshot, TreeError> {
            let mut cache = cell.borrow_mut();
            cache.clear();
            cache.reserve(num_leaves * 2);
            Self::compute_bitsets(root_id, tree, &node_id_to_leaf_index, words, &mut cache)?;

            // Step 4: Collect partitions (internal bipartitions + pendant edges)
            let (parts, lengths) = Self::collect_partitions(tree, root_id, &cache)?;

            // Step 5: Canonicalize partitions
            // rooted_mode=false: bipartitions (canonicalized, deduped, trivial-filtered)
            // rooted_mode=true:  clades (raw subtree bitsets, just sorted)
            let (parts_canonical, lengths_canonical) =
                Self::canonicalize_partitions(parts, lengths, words, num_leaves, rooted_mode);

            Ok(Snapshot {
                parts: parts_canonical,
                lengths: lengths_canonical,
                leaf_names: sorted_leaf_names,
                words,
                rooted,
            })
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
    ) -> Result<Bitset, TreeError> {
        if let Some(bitset) = cache.get(&node_id) {
            return Ok(bitset.clone());
        }

        let node = tree.get(&node_id)?;

        if node.children.is_empty() {
            let mut bitset = Bitset::zeros(words);
            let leaf_idx = *node_id_to_leaf_index
                .get(&node_id)
                .ok_or(TreeError::NodeNotFound(node_id))?;
            bitset.set(leaf_idx);
            cache.insert(node_id, bitset.clone());
            return Ok(bitset);
        }

        let mut bitset = Bitset::zeros(words);
        for &child_id in &node.children {
            let child_bitset =
                Self::compute_bitsets(child_id, tree, node_id_to_leaf_index, words, cache)?;
            bitset.or_assign(&child_bitset);
        }

        cache.insert(node_id, bitset.clone());
        Ok(bitset)
    }

    /// Collect all partitions (internal bipartitions and pendant edges) with branch lengths.
    ///
    /// # What we skip
    /// - Root node (doesn't create a bipartition)
    ///
    /// # Branch lengths
    /// Some trees may have missing branch lengths.
    /// We treat missing lengths as 0.0.
    fn collect_partitions(
        tree: &PhyloTree,
        root_id: usize,
        cache: &FxHashMap<usize, Bitset>,
    ) -> Result<(Vec<Bitset>, Vec<f64>), TreeError> {
        cache
            .iter()
            .filter(|(node_id, _)| **node_id != root_id)
            .map(|(node_id, bitset)| {
                let length = tree.get(node_id)?.parent_edge.unwrap_or(0.0);
                Ok((bitset.clone(), length))
            })
            .collect::<Result<Vec<_>, _>>()
            .map(|v| v.into_iter().unzip())
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
    /// Returns (Vec<Bitset>, Vec<f64>) - parallel vectors, sorted by Bitset
    fn canonicalize_partitions(
        parts: Vec<Bitset>,
        lengths: Vec<f64>,
        words: usize,
        num_leaves: usize,
        rooted_mode: bool,
    ) -> (Vec<Bitset>, Vec<f64>) {
        if rooted_mode {
            // Rooted clade mode: use raw subtree bitsets as clades.
            // No canonicalization, no trivial filter, no dedup.
            // Both root children are distinct clades; L-1 leaf clades are valid.
            let mut pairs: Vec<(Bitset, f64)> = parts.into_iter().zip(lengths).collect();
            pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            return pairs.into_iter().unzip();
        }

        // Unrooted bipartition mode: canonicalize, filter trivials, dedup.
        //
        // Three cases, decided by ones_before (count before any flip):
        //   ones_before == 1              → pendant edge: keep as single-bit bitset,
        //                                   no canonicalization needed.
        //   ones_before >= num_leaves - 1 → near-full complement of a pendant edge;
        //                                   filter out (the pendant itself is captured above).
        //   otherwise                     → internal bipartition: store the side that
        //                                   does NOT contain leaf 0.
        let mut pairs: Vec<(Bitset, f64)> = parts
            .into_iter()
            .zip(lengths)
            .filter_map(|(bitset, length)| {
                let ones_before = bitset.count_ones();

                if ones_before == 1 {
                    return Some((bitset, length));
                }

                if ones_before >= num_leaves - 1 {
                    return None;
                }

                let leaf_0_is_set = (bitset.0[0] & 1) != 0;
                let canonical_bitset = if leaf_0_is_set {
                    Self::compute_complement(&bitset, words, num_leaves)
                } else {
                    bitset
                };

                Some((canonical_bitset, length))
            })
            .collect();

        // CRITICAL: Sort by bitset for O(m+n) merge-based intersection
        pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Deduplicate consecutive identical bitsets (e.g. root bipartition
        // in rooted binary trees where both root children canonicalize to
        // the same split).  Merge by summing branch lengths.
        let mut deduped: Vec<(Bitset, f64)> = Vec::with_capacity(pairs.len());
        for (bitset, length) in pairs {
            if let Some(last) = deduped.last_mut()
                && last.0 == bitset
            {
                last.1 += length;
                continue;
            }
            deduped.push((bitset, length));
        }

        // Unzip into parallel vectors - indices now match!
        let (canonical_parts, canonical_lengths): (Vec<Bitset>, Vec<f64>) =
            deduped.into_iter().unzip();

        (canonical_parts, canonical_lengths)
    }

    /// Compute the bitwise complement of a partition.
    ///
    /// Flips all bits up to `num_leaves` using word-level NOT, then masks the
    /// trailing garbage bits in the final word.
    ///
    /// # Example
    /// Input:  0b0011 (4 leaves) → Output: 0b1100
    fn compute_complement(bitset: &Bitset, words: usize, num_leaves: usize) -> Bitset {
        let mut complement = Bitset(bitset.0.iter().map(|&w| !w).collect());
        let used = num_leaves % 64;
        if used != 0 {
            complement.0[words - 1] &= (1u64 << used) - 1;
        }
        complement
    }
}

/// One tree's bipartitions in interned form (private implementation detail of [`Snapshots`]).
///
/// `split_ids` is sorted ascending so the RF sorted-merge can run on integers.
/// `lengths[i]` is the branch length of `split_ids[i]` (parallel arrays).
#[derive(Debug, Clone)]
pub(crate) struct InternSnap {
    pub(crate) split_ids: Vec<u32>,
    pub(crate) lengths: Vec<f64>,
}

/// A bulk collection of tree snapshots in an interned split-ID representation.
///
/// All trees in the set share a single bipartition table: each unique split is
/// assigned a `u32` ID once. Pairwise distance computation then works on
/// `Vec<u32>` instead of `Vec<Bitset>`, which is faster and more cache-friendly.
///
/// # Construction
/// - [`Snapshots::from_newick_iter`] — main constructor for BEAST/NEXUS inputs
/// - [`Snapshots::from_newicks`] — convenience constructor for plain Newick slices
///
/// # Pairwise distances
/// - [`Snapshots::pairwise_rf`], [`Snapshots::pairwise_wrf`], [`Snapshots::pairwise_kf`]
#[derive(Debug)]
pub struct Snapshots {
    pub(crate) snapshots: Vec<InternSnap>,
    /// The bipartition at index `i` corresponds to split ID `i`.
    pub bipartitions: Vec<Bitset>,
    /// Reverse map: bipartition bitset → split ID.
    pub bipartition_index: FxHashMap<Bitset, u32>,
    pub words_per_bitset: usize,
    /// Alphabetically sorted taxon names shared by all trees in this set.
    pub leaf_names: Vec<String>,
}

impl Snapshots {
    /// Build a `Snapshots` collection from a lazy iterator of `(newick, translate_map)` pairs.
    ///
    /// Each `newick` may contain BEAST-format `[&...]` annotations — they are stripped
    /// automatically. The `translate_map` is applied to rename leaf labels (pass an empty
    /// map for plain Newick files with no BEAST translate block).
    ///
    /// All trees must share the same leaf set; an error is returned if any tree differs.
    ///
    /// # Errors
    /// Returns `Err(String)` if any newick fails to parse or leaf sets are inconsistent.
    pub fn from_newick_iter<'a>(
        entries: impl IntoIterator<Item = (&'a str, &'a HashMap<String, String>)>,
        rooted: bool,
    ) -> Result<Self, String> {
        let entries: Vec<_> = entries.into_iter().collect();
        if entries.is_empty() {
            return Ok(Self::empty());
        }

        // Parse the first tree to establish the reference leaf set.
        let (first_newick, first_translate) = entries[0];
        let first_clean = crate::io::strip_beast_annotations(first_newick);
        let mut first_tree = PhyloTree::from_newick(&first_clean)
            .map_err(|e| format!("Failed to parse newick at index 0: {e}"))?;
        crate::io::rename_leaf_nodes(&mut first_tree, first_translate);

        //sanity check: ensure all leaf names are unique within the first tree
        let first_leaves = first_tree.get_leaves();
        if first_leaves.len() != first_leaves.iter().collect::<HashSet<_>>().len() {
            return Err(
                "Trees have duplicate leaf names. All leaf names must be unique.".to_string(),
            );
        }

        let mut sorted_leaf_names: Vec<String> = first_leaves
            .iter()
            .filter_map(|&id| first_tree.get(&id).ok()?.name.clone())
            .collect();
        sorted_leaf_names.sort_unstable();
        sorted_leaf_names.dedup();
        let reference_leaves: HashSet<String> = sorted_leaf_names.iter().cloned().collect();

        let first_snap = Snapshot::from_tree(&first_tree, rooted)
            .map_err(|e| format!("Failed to snapshot tree at index 0: {e}"))?;
        drop(first_tree);

        // Parse all remaining trees in parallel, validating leaf sets.
        let rest: Vec<Snapshot> = entries[1..]
            .par_iter()
            .enumerate()
            .map(|(j, &(newick, translate))| {
                let i = j + 1;
                let clean = crate::io::strip_beast_annotations(newick);
                let mut tree = PhyloTree::from_newick(&clean)
                    .map_err(|e| format!("Failed to parse newick at index {i}: {e}"))?;
                crate::io::rename_leaf_nodes(&mut tree, translate);

                let leaves: HashSet<String> = tree
                    .get_leaves()
                    .iter()
                    .filter_map(|&id| tree.get(&id).ok()?.name.clone())
                    .collect();
                if leaves != reference_leaves {
                    return Err(format!(
                        "Tree {i} has a different leaf set than tree 0. All trees must share the same taxa."
                    ));
                }
                Snapshot::from_tree(&tree, rooted)
                    .map_err(|e| format!("Failed to snapshot tree at index {i}: {e}"))
            })
            .collect::<Result<_, _>>()?;

        let mut snaps = Vec::with_capacity(entries.len());
        snaps.push(first_snap);
        snaps.extend(rest);

        Ok(Self::intern(snaps, sorted_leaf_names))
    }

    /// Build a `Snapshots` collection from a slice of plain Newick strings.
    ///
    /// No BEAST annotation stripping or taxon renaming is performed — the strings
    /// must already be in standard Newick format.
    ///
    /// # Errors
    /// Returns `Err(String)` if any newick fails to parse or leaf sets differ.
    pub fn from_newicks(newicks: &[&str], rooted: bool) -> Result<Self, String> {
        let empty: HashMap<String, String> = HashMap::new();
        Self::from_newick_iter(newicks.iter().map(|&n| (n, &empty)), rooted)
    }

    /// Number of trees in this collection.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Returns `true` if the collection contains no trees.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// Map each split ID to its column position in ascending `Bitset` order.
    ///
    /// Returns a `Vec<usize>` where `result[split_id]` = column index.
    pub fn sorted_bip_id_to_col(&self) -> Vec<usize> {
        let n = self.bipartitions.len();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| self.bipartitions[a].cmp(&self.bipartitions[b]));
        let mut id_to_col = vec![0usize; n];
        for (col, &orig_id) in order.iter().enumerate() {
            id_to_col[orig_id] = col;
        }
        id_to_col
    }

    /// Build a flat row-major presence matrix `(n_trees × n_bip)` as `Vec<u8>`,
    /// together with a column-order index into `self.bipartitions`.
    ///
    /// Returns `(presence_bytes, col_to_bip_id)`:
    /// - `presence_bytes`: flat `uint8` buffer, row-major, shape `(n_trees, n_bip)`
    /// - `col_to_bip_id`: `col_to_bip_id[col]` is the index into `self.bipartitions`
    ///   for that column, in ascending `Bitset` order
    ///
    /// Columns are in ascending `Bitset` order (stable across calls on the same tree set).
    /// Each byte is `1` if the split is present in that tree, `0` otherwise.
    pub fn build_presence_matrix(&self) -> (Vec<u8>, Vec<usize>) {
        let id_to_col = self.sorted_bip_id_to_col();
        let n_trees = self.snapshots.len();
        let n_bip = self.bipartitions.len();
        let mut presence = vec![0u8; n_trees * n_bip];

        // Invert id_to_col: col_to_bip_id[col] = original bipartition ID.
        // Returns a cheap Vec<usize> instead of cloning the Bitsets themselves.
        let mut col_to_bip_id = vec![0usize; n_bip];
        for (id, &col) in id_to_col.iter().enumerate() {
            col_to_bip_id[col] = id;
        }

        if n_bip == 0 {
            return (presence, col_to_bip_id);
        }

        // Safely divide the mutable slice into row-sized chunks across threads
        presence
            .par_chunks_mut(n_bip)
            .zip(&self.snapshots) // Pair each chunk with its corresponding snapshot
            .for_each(|(row, snap)| {
                for &split_id in &snap.split_ids {
                    // Write directly into the final memory location
                    row[id_to_col[split_id as usize]] = 1;
                }
            });

        (presence, col_to_bip_id)
    }

    /// Export canonical bipartition bitmasks as a flat byte buffer.
    ///
    /// Shape: `(n_bip, ceil(n_leaves / 8))` bytes, row-major, in ascending `Bitset`
    /// (column) order — matching the column order of [`build_presence_matrix`].
    ///
    /// Bit `i` of row `j`, using **little-endian bit order within each byte**, is `1`
    /// if `leaf_names[i]` is on the canonical side of bipartition `j`.
    ///
    /// Decode on the Python side:
    /// ```python
    /// bytes_per_bip = math.ceil(len(leaf_names) / 8)
    /// bip_arr  = np.frombuffer(bip_bytes, dtype=np.uint8).reshape(n_bip, bytes_per_bip)
    /// bip_bool = np.unpackbits(bip_arr, axis=1, bitorder='little')[:, :len(leaf_names)]
    /// # bip_bool[j, i] == 1  →  leaf_names[i] is on the canonical side of bipartition j
    /// ```
    pub fn build_bipartition_bytes(&self, col_to_bip_id: &[usize]) -> Vec<u8> {
        let n_leaves = self.leaf_names.len();
        let bytes_per_bip = n_leaves.div_ceil(8);
        let mut out = Vec::with_capacity(col_to_bip_id.len() * bytes_per_bip);
        for &id in col_to_bip_id {
            let bip = &self.bipartitions[id];
            out.extend(
                bip.0
                    .iter()
                    .flat_map(|&word| word.to_le_bytes())
                    .take(bytes_per_bip),
            );
        }
        out
    }

    /// Compute all pairwise Robinson–Foulds distances as a symmetric n×n matrix.
    pub fn pairwise_rf(&self) -> Vec<usize> {
        crate::distances::pairwise_symmetric(self, crate::distances::rf_distance_fast)
    }

    /// Compute all pairwise Weighted Robinson–Foulds distances as a symmetric n×n matrix.
    pub fn pairwise_wrf(&self) -> Vec<f64> {
        crate::distances::pairwise_symmetric(self, crate::distances::wrf_distance_fast)
    }

    /// Compute all pairwise Kuhner–Felsenstein distances as a symmetric n×n matrix.
    pub fn pairwise_kf(&self) -> Vec<f64> {
        crate::distances::pairwise_symmetric(self, crate::distances::kf_distance_fast)
    }

    /// Build the interning table from a `Vec<Snapshot>` (internal helper).
    fn intern(snaps: Vec<Snapshot>, leaf_names: Vec<String>) -> Self {
        let words = snaps.first().map(|s| s.words).unwrap_or(0);
        let mut index: FxHashMap<Bitset, u32> = FxHashMap::default();
        let mut bipartitions: Vec<Bitset> = Vec::new();

        let interned: Vec<InternSnap> = snaps
            .into_iter()
            .map(|snap| {
                let mut paired: Vec<(u32, f64)> = snap
                    .parts
                    .into_iter()
                    .zip(snap.lengths)
                    .map(|(b, length)| {
                        let next_id = bipartitions.len() as u32;
                        let id = if let Some(&id) = index.get(&b) {
                            id
                        } else {
                            index.insert(b.clone(), next_id);
                            bipartitions.push(b);
                            next_id
                        };
                        (id, length)
                    })
                    .collect();

                paired.sort_unstable_by_key(|&(id, _)| id);
                let (split_ids, lengths): (Vec<u32>, Vec<f64>) = paired.into_iter().unzip();
                InternSnap { split_ids, lengths }
            })
            .collect();

        Self {
            snapshots: interned,
            bipartitions,
            bipartition_index: index,
            words_per_bitset: words,
            leaf_names,
        }
    }

    fn empty() -> Self {
        Self {
            snapshots: Vec::new(),
            bipartitions: Vec::new(),
            bipartition_index: FxHashMap::default(),
            words_per_bitset: 0,
            leaf_names: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A symmetric 4-leaf tree produces a single bipartition after canonicalization.
    ///
    /// ```text
    ///        root
    ///       /    \
    ///   node1    node2
    ///   /   \    /   \
    ///  A     B  C     D
    /// ```
    ///
    /// Leaves sorted: A=0, B=1, C=2, D=3.
    ///
    /// node1's subtree = {A,B} = 0b0011 — contains A (bit 0), so flip to {C,D} = 0b1100.
    /// node2's subtree = {C,D} = 0b1100 — no A, keep as-is.
    ///
    /// Both sides of the root bipartition canonicalize to the same bitset because they
    /// represent the **same split** viewed from either side. `canonicalize_partitions`
    /// deduplicates them into one entry and sums the branch lengths.
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

        // Verify that canonicalize_partitions deduplicates the root bipartition:
        // parts = [{A,B} = 0b0011, {C,D} = 0b1100]
        // After canonicalization: both become 0b1100
        // After dedup: only one partition remains
        let parts = vec![part_ab, part_cd];
        let lengths = vec![1.0, 2.0];
        let (canon_parts, canon_lengths) =
            Snapshot::canonicalize_partitions(parts, lengths, 1, 4, false);
        assert_eq!(
            canon_parts.len(),
            1,
            "Root bipartition should be deduplicated to 1 partition, got {}",
            canon_parts.len()
        );
        assert_eq!(
            canon_lengths[0], 3.0,
            "Deduplicated lengths should sum: 1.0 + 2.0 = 3.0"
        );
    }

    /// Test that from_tree deduplicates root bipartitions in rooted binary trees.
    ///
    /// In a rooted binary tree ((A,B),(C,D)), both root children produce the
    /// same canonical bipartition ({C,D}|{A,B}).  Without deduplication, this
    /// causes RF distances to be inflated by +2 compared to the unrooted RF
    /// (e.g. R's phangorn RF.dist with rooted=FALSE).
    #[test]
    fn test_root_bipartition_dedup() {
        use phylotree::tree::Tree as PhyloTree;

        // Tree: ((A:1,B:1):1,(C:1,D:1):1);
        // Root has two internal children - classic root bipartition duplication case
        let tree = PhyloTree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let snap = Snapshot::from_tree(&tree, false).unwrap();

        // For 4 leaves: 4 pendant edges + 1 internal bipartition = 5 entries.
        // With dedup, the duplicated root bipartition collapses to 1 internal entry.
        assert_eq!(
            snap.parts.len(),
            5,
            "Rooted 4-leaf binary tree should have 5 entries (4 pendant + 1 internal) after dedup, got {}",
            snap.parts.len()
        );
    }

    /// Test rooted vs unrooted mode partition counts and RF distances.
    #[test]
    fn test_rooted_vs_unrooted_partitions() {
        use phylotree::tree::Tree as PhyloTree;

        let rf_pair = |a: &Snapshot, b: &Snapshot| -> usize { crate::distances::rf_distance(a, b) };

        let tree1 = PhyloTree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let tree2 = PhyloTree::from_newick("((A:1,C:1):1,(B:1,D:1):1);").unwrap();

        // Unrooted mode: 4 pendant edges + 1 internal bipartition = 5 entries per tree.
        let snap1_u = Snapshot::from_tree(&tree1, false).unwrap();
        let snap2_u = Snapshot::from_tree(&tree2, false).unwrap();
        assert_eq!(
            snap1_u.parts.len(),
            5,
            "Unrooted: 4 pendant + 1 internal for 4-leaf tree"
        );
        assert_eq!(snap2_u.parts.len(), 5);
        assert_eq!(rf_pair(&snap1_u, &snap2_u), 2, "Unrooted RF = 2");

        // Rooted mode: 4 pendant edges + 2 internal clades = 6 entries per tree.
        let snap1_r = Snapshot::from_tree(&tree1, true).unwrap();
        let snap2_r = Snapshot::from_tree(&tree2, true).unwrap();
        assert_eq!(
            snap1_r.parts.len(),
            6,
            "Rooted: 4 pendant + 2 clades for 4-leaf tree"
        );
        assert_eq!(snap2_r.parts.len(), 6);
        assert_eq!(rf_pair(&snap1_r, &snap2_r), 4, "Rooted RF = 4");

        // Same topology: both modes give RF = 0
        let tree1b = PhyloTree::from_newick("((B:2,A:2):2,(D:2,C:2):2);").unwrap();
        let snap1b_u = Snapshot::from_tree(&tree1b, false).unwrap();
        let snap1b_r = Snapshot::from_tree(&tree1b, true).unwrap();
        assert_eq!(rf_pair(&snap1_u, &snap1b_u), 0, "Unrooted same topo = 0");
        assert_eq!(rf_pair(&snap1_r, &snap1b_r), 0, "Rooted same topo = 0");
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

    /// Verify that `build_bipartition_bytes` exports canonical bitmasks correctly.
    ///
    /// Pendant edges are included, so for 4 leaves the table has 6 entries:
    ///   cols 0-3: pendant edges {A}=0x01, {B}=0x02, {C}=0x04, {D}=0x08
    ///   col 4: {B,D} = bits 1,3 → 0x0A
    ///   col 5: {C,D} = bits 2,3 → 0x0C
    #[test]
    fn test_build_bipartition_bytes() {
        let snaps = Snapshots::from_newicks(
            &["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"],
            false,
        )
        .unwrap();

        let (_, col_to_bip_id) = snaps.build_presence_matrix();
        let bip_bytes = snaps.build_bipartition_bytes(&col_to_bip_id);

        // 6 bipartitions (4 pendant + 2 internal) × ceil(4/8) = 1 byte each → 6 bytes total
        assert_eq!(bip_bytes.len(), 6);

        // cols 0-3: pendant edges, sorted ascending
        assert_eq!(bip_bytes[0], 0x01, "col 0 should be pendant {{A}} = 0x01");
        assert_eq!(bip_bytes[1], 0x02, "col 1 should be pendant {{B}} = 0x02");
        assert_eq!(bip_bytes[2], 0x04, "col 2 should be pendant {{C}} = 0x04");
        assert_eq!(bip_bytes[3], 0x08, "col 3 should be pendant {{D}} = 0x08");
        // cols 4-5: internal bipartitions
        assert_eq!(bip_bytes[4], 0x0A, "col 4 should be {{B,D}} = 0x0A");
        assert_eq!(bip_bytes[5], 0x0C, "col 5 should be {{C,D}} = 0x0C");
    }

    /// Verify that `build_presence_matrix` returns col_to_bip_id in ascending Bitset order.
    ///
    /// Pendant edges appear before internal bipartitions in the sorted table.
    /// All trees share the same leaf set so pendant columns are all-1.
    #[test]
    #[allow(clippy::erasing_op)]
    fn test_build_presence_matrix_sorted() {
        let snaps = Snapshots::from_newicks(
            &["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"],
            false,
        )
        .unwrap();

        // Leaf order: A=0, B=1, C=2, D=3
        assert_eq!(snaps.leaf_names, vec!["A", "B", "C", "D"]);

        let (presence, col_to_bip_id) = snaps.build_presence_matrix();
        // 4 pendant edges + 2 internal bipartitions
        assert_eq!(col_to_bip_id.len(), 6);

        // cols 4 and 5 are the internal bipartitions in ascending order
        let col4_bits = snaps.bipartitions[col_to_bip_id[4]].0[0];
        let col5_bits = snaps.bipartitions[col_to_bip_id[5]].0[0];
        assert_eq!(col4_bits, 0b1010, "col 4 should be {{B,D}} = 0b1010");
        assert_eq!(col5_bits, 0b1100, "col 5 should be {{C,D}} = 0b1100");

        // Pendant columns 0-3 are all 1 (both trees share the same leaf set).
        // T1 has {C,D} (col5=1) but not {B,D} (col4=0).
        // T2 has {B,D} (col4=1) but not {C,D} (col5=0).
        let n = 6;
        assert_eq!(presence[0 * n + 4], 0, "T1,col4 ({{B,D}}): absent");
        assert_eq!(presence[0 * n + 5], 1, "T1,col5 ({{C,D}}): present");
        assert_eq!(presence[n + 4], 1, "T2,col4 ({{B,D}}): present");
        assert_eq!(presence[n + 5], 0, "T2,col5 ({{C,D}}): absent");
    }
}
