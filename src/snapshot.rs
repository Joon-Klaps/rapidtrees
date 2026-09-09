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
use crate::par::*;
use hashbrown::HashTable;
use phylotree::tree::{Tree as PhyloTree, TreeError};
use rustc_hash::{FxBuildHasher, FxHashMap};
use std::collections::{HashMap, HashSet};
use std::hash::BuildHasher;

/// A 128-bit XOR fingerprint of a leaf set.
///
/// Each taxon draws one random 128-bit label at the start of a run; a subtree's
/// fingerprint is the XOR of its leaves' labels, so a node costs `O(1)` to fold
/// instead of `⌈n/64⌉` words. XOR is what makes the complement free —
/// `label(A′) = total ^ label(A)` — so a bipartition's canonical form is
/// `min(h, h ^ total)`, computed without touching either leaf set.
///
/// Two distinct splits share a fingerprint with probability about `e² / 2¹²⁹`,
/// where `e` is the number of distinct splits the run sees. At `e = 10⁸` that
/// is `1.5 × 10⁻²³` — some nineteen orders of magnitude below the rate at which
/// the machine's own memory flips a bit unnoticed. The `verify` feature turns
/// the assumption into a checked one; see [`Interner::push`].
type Fingerprint = u128;

/// One random 128-bit label per taxon, drawn once and shared by every tree.
///
/// The seed is fixed: the same input must produce the same split IDs on every
/// run and every machine.
fn taxon_labels(num_leaves: usize) -> Vec<Fingerprint> {
    // splitmix64, seeded from the digits of pi.
    let mut state = 0x243F_6A88_85A3_08D3u64;
    let mut next = move || {
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

/// One edge of a tree, identified by a fingerprint instead of a leaf set.
///
/// `first`/`size` locate the subtree's leaves as a contiguous run of
/// [`Snapshot::leaf_order`], which is what lets the canonical [`Bitset`] be
/// rebuilt later — once per *unique* split — with the tree already dropped.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Part {
    /// Canonical fingerprint: `min(h, h ^ total)` for an internal bipartition,
    /// so both of its sides give the same value — that is what replaces
    /// canonicalising the bitset. Pendant edges and rooted clades keep the raw
    /// subtree fingerprint.
    key: Fingerprint,
    /// Where this subtree's leaves start in [`Snapshot::leaf_order`].
    first: u32,
    /// Leaves in the subtree, before canonicalisation. Load-bearing: it stands
    /// in for `count_ones()` in the pendant and near-full filters, which would
    /// otherwise put the full-width popcount back.
    size: u32,
    /// Length of the edge above the node.
    length: f64,
}

/// One tree's bipartitions, on the way into a [`Snapshots`] collection.
///
/// Short-lived: built per tree, handed straight to the interner, dropped. That
/// is what keeps construction peak memory near the deduplicated footprint.
/// Taxon names live once on [`Snapshots`], not per tree.
///
/// # Canonicalization
/// Each internal bipartition has two complementary sides. The fingerprint
/// canonicalises them for free by taking `min(h, h ^ total)`; the `Bitset` the
/// interner materialises is still the side without leaf 0, so the public
/// bipartition table is unchanged. Pendant edges are stored as-is.
#[derive(Debug, Clone)]
pub(crate) struct Snapshot {
    /// All edges (internal bipartitions + pendant edges), filtered and sorted
    /// by fingerprint.
    pub(crate) parts: Vec<Part>,

    /// Leaf bit indices in traversal order. Every subtree occupies a
    /// contiguous run, so `leaf_order[first..first + size]` is a part's leaves.
    pub(crate) leaf_order: Vec<u32>,

    /// Number of u64 words per bitset.
    pub(crate) words: usize,
}

impl Snapshot {
    /// Extract a snapshot from an already-parsed tree.
    ///
    /// # Algorithm
    /// 1. Extract leaf names and sort them alphabetically for consistency
    /// 2. Map each leaf name to a compact index [0..n)
    /// 3. One DFS from the root, folding a 128-bit fingerprint and a leaf count
    ///    per node — `O(1)` each, so the tree costs `Θ(n)` rather than `Θ(n²/64)`
    /// 4. Drop the trivial splits and canonicalize by fingerprint
    /// 5. Sort by fingerprint, merging any duplicate
    ///
    /// `labels` must hold one entry per taxon, indexed by the same alphabetical
    /// bit index, and must be the *same* labels for every tree in a run — see
    /// [`taxon_labels`].
    ///
    /// # Errors
    /// Returns `TreeError` if the tree is empty, malformed, or has unnamed leaves.
    pub(crate) fn from_tree(
        tree: &PhyloTree,
        rooted_mode: bool,
        labels: &[Fingerprint],
    ) -> Result<Self, TreeError> {
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

        let num_leaves = leaf_id_names.len();
        let words = num_leaves.div_ceil(64);

        // Step 2: Create mapping: node_id → bit_index (based on sorted names)
        let node_id_to_leaf_index: FxHashMap<usize, usize> = leaf_id_names
            .iter()
            .enumerate()
            .map(|(idx, &(node_id, _))| (node_id, idx))
            .collect();

        // Step 3: fingerprint every node
        let root_id = tree.get_root()?;
        let (leaf_order, parts) =
            Self::collect_partitions(tree, root_id, &node_id_to_leaf_index, labels, rooted_mode)?;

        // Steps 4-5
        // rooted_mode=false: bipartitions (canonicalized, deduped, trivial-filtered)
        // rooted_mode=true:  clades (raw subtree fingerprints, just sorted)
        let parts = Self::canonicalize_partitions(parts, num_leaves, rooted_mode);

        Ok(Snapshot {
            parts,
            leaf_order,
            words,
        })
    }

    /// Every node's subtree fingerprint and leaf run, paired with its edge.
    ///
    /// Two iterative passes over one flat arena indexed by node id. The
    /// accumulator is a 128-bit XOR plus a `u32` count, both `O(1)`, where a
    /// bitset accumulator cost `⌈n/64⌉` words per node and `Θ(n)` allocations
    /// per tree.
    ///
    /// Leaves are numbered in traversal order, so every subtree's leaves land
    /// in a contiguous run of `leaf_order`. That run is how a canonical
    /// `Bitset` is rebuilt later without keeping the tree alive.
    ///
    /// The root is skipped — it creates no bipartition. Missing branch lengths
    /// are treated as 0.0.
    fn collect_partitions(
        tree: &PhyloTree,
        root_id: usize,
        node_id_to_leaf_index: &FxHashMap<usize, usize>,
        labels: &[Fingerprint],
        rooted_mode: bool,
    ) -> Result<(Vec<u32>, Vec<Part>), TreeError> {
        let mut order = Vec::with_capacity(tree.size());
        let mut stack = vec![root_id];
        while let Some(id) = stack.pop() {
            order.push(id);
            stack.extend(tree.get(&id)?.children.iter().copied());
        }

        let mut acc = vec![Acc::default(); tree.size()];
        let mut leaf_order = Vec::with_capacity(node_id_to_leaf_index.len());

        // Pre-order: each leaf takes the next slot, which is what keeps a
        // subtree's leaves adjacent.
        for &id in &order {
            if tree.get(&id)?.children.is_empty() {
                let bit = *node_id_to_leaf_index
                    .get(&id)
                    .ok_or(TreeError::NodeNotFound(id))?;
                acc[id] = Acc {
                    fp: labels[bit],
                    first: leaf_order.len() as u32,
                    size: 1,
                };
                leaf_order.push(bit as u32);
            }
        }

        // Post-order: an internal node is the XOR of its children.
        for &id in order.iter().rev() {
            let node = tree.get(&id)?;
            if node.children.is_empty() {
                continue;
            }
            let folded = node.children.iter().fold(Acc::empty(), |a, &c| Acc {
                fp: a.fp ^ acc[c].fp,
                first: a.first.min(acc[c].first),
                size: a.size + acc[c].size,
            });
            acc[id] = folded;
        }

        let total = acc[root_id].fp;
        let parts = order
            .into_iter()
            .filter(|&id| id != root_id)
            .map(|id| {
                let Acc { fp, first, size } = acc[id];
                Ok(Part {
                    // A pendant keeps its raw fingerprint: its complement is
                    // filtered out below, so nothing else can claim the key,
                    // and at two taxa `min(h, h ^ total)` would collapse the
                    // two pendants onto each other.
                    key: if rooted_mode || size == 1 {
                        fp
                    } else {
                        fp.min(fp ^ total)
                    },
                    first,
                    size,
                    length: tree.get(&id)?.parent_edge.unwrap_or(0.0),
                })
            })
            .collect::<Result<_, TreeError>>()?;

        Ok((leaf_order, parts))
    }

    /// Drop the trivial parts, then sort by fingerprint and merge duplicates.
    ///
    /// Three cases, decided by the raw subtree size — which the DFS already
    /// accumulated, so no popcount is needed:
    ///   size == 1              → pendant edge, kept as-is.
    ///   size >= num_leaves - 1 → the complement of a pendant; dropped, since
    ///                            the pendant itself is kept above.
    ///   otherwise              → internal bipartition. `key` already carries
    ///                            `min(h, h ^ total)`, so both sides of a split
    ///                            arrive equal and there is nothing to flip.
    ///
    /// Sorting a 16-byte key replaces sorting `⌈n/64⌉`-word bitsets: at 2 000
    /// taxa, a 16-byte compare rather than a 256-byte one.
    fn canonicalize_partitions(
        mut parts: Vec<Part>,
        num_leaves: usize,
        rooted_mode: bool,
    ) -> Vec<Part> {
        if rooted_mode {
            // Rooted clade mode: raw subtree fingerprints as clades.
            // No canonicalization, no trivial filter, no dedup.
            // Both root children are distinct clades; L-1 leaf clades are valid.
            parts.sort_unstable_by_key(|p| p.key);
            return parts;
        }

        parts.retain(|p| p.size == 1 || (p.size as usize) < num_leaves - 1);
        parts.sort_unstable_by_key(|p| p.key);

        // Deduplicate identical fingerprints (e.g. the root bipartition in a
        // rooted binary tree, where both root children canonicalize to the same
        // split). Merge by summing branch lengths.
        parts.dedup_by(|later, kept| {
            later.key == kept.key && {
                kept.length += later.length;
                true
            }
        });

        parts
    }

    /// Rebuild one part's canonical `Bitset` from its leaf run.
    ///
    /// `O(subtree size)` and run once per *unique* split rather than once per
    /// occurrence — which is what the fingerprint build buys. The result is the
    /// side without leaf 0, exactly the representation the bitset-accumulator
    /// path produced, so [`Snapshots::bipartitions`] is byte-identical.
    fn materialise(&self, part: &Part, num_leaves: usize, rooted_mode: bool) -> Bitset {
        let mut bitset = Bitset::zeros(self.words);
        for &bit in &self.leaf_order[part.first as usize..][..part.size as usize] {
            bitset.set(bit as usize);
        }
        // A pendant edge is kept as-is, leaf 0 or not; its complement was
        // filtered, so no other occurrence can disagree.
        if !rooted_mode && part.size > 1 && bitset.0[0] & 1 != 0 {
            Self::complement_in_place(&mut bitset, num_leaves);
        }
        bitset
    }

    /// Flip a partition to its complement in place.
    ///
    /// Word-level NOT up to `num_leaves`, then mask the trailing garbage bits
    /// in the final word.
    ///
    /// # Example
    /// Input:  0b0011 (4 leaves) → Output: 0b1100
    fn complement_in_place(bitset: &mut Bitset, num_leaves: usize) {
        bitset.0.iter_mut().for_each(|w| *w = !*w);
        let used = num_leaves % 64;
        if used != 0
            && let Some(last) = bitset.0.last_mut()
        {
            *last &= (1u64 << used) - 1;
        }
    }
}

/// Per-node DFS accumulator: a fingerprint and the subtree's leaf run.
#[derive(Debug, Clone, Copy, Default)]
struct Acc {
    fp: Fingerprint,
    first: u32,
    size: u32,
}

impl Acc {
    /// Identity for the children fold: XOR's zero, and a `first` that any real
    /// child wins the `min` against.
    fn empty() -> Self {
        Acc {
            fp: 0,
            first: u32::MAX,
            size: 0,
        }
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
        Self::from_newick_iter_opts(entries, rooted, true)
    }

    /// Like [`Snapshots::from_newick_iter`], but lets internal callers opt out of
    /// storing branch lengths.
    ///
    /// When `store_lengths` is `false`, each `InternSnap.lengths` is left empty —
    /// saving `~n_bip × 8` bytes per tree on RF-only paths that never read lengths.
    /// The public constructor always passes `true`, so the Rust/Python API is
    /// unaffected.
    pub(crate) fn from_newick_iter_opts<'a>(
        entries: impl IntoIterator<Item = (&'a str, &'a HashMap<String, String>)>,
        rooted: bool,
        store_lengths: bool,
    ) -> Result<Self, String> {
        let entries: Vec<_> = entries.into_iter().collect();
        if entries.is_empty() {
            return Ok(Self::empty());
        }

        // Parse the first tree to establish the reference leaf set.
        let (first_newick, first_translate) = entries[0];
        let first_tree = parse_and_rename(first_newick, first_translate, 0)?;

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

        // One label set for the whole run: fingerprints are only comparable
        // across trees if every tree draws from the same table.
        let labels = taxon_labels(sorted_leaf_names.len());

        let first_snap = Snapshot::from_tree(&first_tree, rooted, &labels)
            .map_err(|e| format!("Failed to snapshot tree at index 0: {e}"))?;
        drop(first_tree);

        // Bound how many raw snapshots are alive at once. Holding every tree's
        // un-interned `Vec<Bitset>` simultaneously is the dominant memory cost at
        // construction — it can dwarf the deduplicated result and OOM the process.
        // Estimate one snapshot's raw bytes from the first tree (all trees share
        // the leaf set, so this is representative), parse the rest in chunks sized
        // to ~`CHUNK_TARGET_BYTES`, and fold each chunk into the interner — freeing
        // its bitsets — before parsing the next.
        const CHUNK_TARGET_BYTES: usize = 256 * 1024 * 1024;
        let per_snap_bytes = first_snap.parts.len() * size_of::<Part>()
            + first_snap.leaf_order.len() * size_of::<u32>();
        let chunk = (CHUNK_TARGET_BYTES / per_snap_bytes.max(1)).clamp(1, 4096);

        let mut interner = Interner::new(
            entries.len(),
            first_snap.words,
            sorted_leaf_names.len(),
            rooted,
            store_lengths,
        );
        interner.push(first_snap);

        let mut base = 1usize; // tree 0 is already interned
        for chunk_entries in entries[1..].chunks(chunk) {
            // Parse this chunk in parallel, validating leaf sets.
            let raw: Vec<Snapshot> = chunk_entries
                .par_iter()
                .enumerate()
                .map(|(k, &(newick, translate))| {
                    let i = base + k;
                    let tree = parse_and_rename(newick, translate, i)?;

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
                    Snapshot::from_tree(&tree, rooted, &labels)
                        .map_err(|e| format!("Failed to snapshot tree at index {i}: {e}"))
                })
                .collect::<Result<_, _>>()?;

            // Sequential fold: each raw snapshot's bitsets are freed right after
            // it is interned, so peak stays near the deduplicated footprint.
            for snap in raw {
                interner.push(snap);
            }
            base += chunk_entries.len();
        }

        Ok(interner.finish(sorted_leaf_names))
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

    /// Build a flat row-major branch-length matrix `(n_trees × n_bip)` as native-endian
    /// `float64` bytes, together with a column-order index into `self.bipartitions`.
    ///
    /// Returns `(branch_length_bytes, col_to_bip_id)`:
    /// - `branch_length_bytes`: flat `float64` buffer, row-major, shape `(n_trees, n_bip)`
    /// - `col_to_bip_id`: `col_to_bip_id[col]` is the index into `self.bipartitions`
    ///   for that column, in ascending `Bitset` order
    ///
    /// `branch_length_bytes[i, j]` is the branch length of edge `j` in tree `i`, or `0.0`
    /// if that edge is absent from tree `i`. Pendant (leaf) edges are always present in
    /// every tree and therefore always have a non-zero value.
    ///
    /// Column order matches [`build_presence_matrix`] (ascending `Bitset` order, stable
    /// across calls on the same tree set).
    ///
    /// Decode on the Python side:
    /// ```python
    /// bl = np.frombuffer(branch_length_bytes, dtype=np.float64).reshape(n_trees, n_bip)
    /// # Fréchet ESS traces without materialising the n×n distance matrix:
    /// wrf_trace = np.sum(np.abs(bl[ref_idx, :] - bl), axis=1)
    /// kf_trace  = np.sqrt(np.sum((bl[ref_idx, :] - bl) ** 2, axis=1))
    /// ```
    pub fn build_branch_length_matrix(&self) -> (Vec<u8>, Vec<usize>) {
        let id_to_col = self.sorted_bip_id_to_col();
        let n_trees = self.snapshots.len();
        let n_bip = self.bipartitions.len();

        let mut col_to_bip_id = vec![0usize; n_bip];
        for (id, &col) in id_to_col.iter().enumerate() {
            col_to_bip_id[col] = id;
        }

        if n_bip == 0 {
            return (Vec::new(), col_to_bip_id);
        }

        // Write bytes straight out rather than building an `f64` matrix and
        // copying it: one allocation instead of two. Absent edges stay 0.0,
        // whose native-endian encoding is the zero bytes we start from.
        let mut bytes = vec![0u8; n_trees * n_bip * 8];
        bytes
            .par_chunks_mut(n_bip * 8)
            .zip(&self.snapshots)
            .for_each(|(row, snap)| {
                for (&split_id, &length) in snap.split_ids.iter().zip(&snap.lengths) {
                    let col = id_to_col[split_id as usize];
                    row[col * 8..(col + 1) * 8].copy_from_slice(&length.to_ne_bytes());
                }
            });

        (bytes, col_to_bip_id)
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
    ///
    /// Pass `Some(counter)` to track progress: after each row `i` finishes, the
    /// counter is bumped by `n - i - 1`, reaching `n*(n-1)/2` when done. Pass
    /// `None` to skip the (negligible) counter work entirely.
    pub fn pairwise_rf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<u32> {
        crate::distances::distance_rf(self, progress)
    }

    /// Compute all pairwise Weighted Robinson–Foulds distances as a symmetric n×n matrix.
    ///
    /// See [`Self::pairwise_rf`] for the `progress` argument.
    pub fn pairwise_wrf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<f64> {
        crate::distances::distance_wrf(self, progress)
    }

    /// Compute all pairwise Kuhner–Felsenstein distances as a symmetric n×n matrix.
    ///
    /// See [`Self::pairwise_rf`] for the `progress` argument.
    pub fn pairwise_kf(&self, progress: Option<&std::sync::atomic::AtomicUsize>) -> Vec<f64> {
        crate::distances::distance_kf(self, progress)
    }

    fn empty() -> Self {
        Self {
            snapshots: Vec::new(),
            bipartitions: Vec::new(),
            words_per_bitset: 0,
            leaf_names: Vec::new(),
        }
    }
}

/// Strip BEAST annotations, parse the Newick, and apply taxon renaming.
///
/// `index` only feeds the parse-error message so it points at the offending tree.
fn parse_and_rename(
    newick: &str,
    translate: &HashMap<String, String>,
    index: usize,
) -> Result<PhyloTree, String> {
    let clean = crate::io::strip_beast_annotations(newick);
    let mut tree = PhyloTree::from_newick(&clean)
        .map_err(|e| format!("Failed to parse newick at index {index}: {e}"))?;
    crate::io::rename_leaf_nodes(&mut tree, translate);
    Ok(tree)
}

/// Incremental bipartition interner.
///
/// Deduplicates bipartitions into split IDs using a `HashTable<u32>` that stores
/// only the integer IDs — each unique bitset is kept exactly once, in
/// `bipartitions`. Snapshots are folded in one at a time via [`Interner::push`],
/// which consumes each raw `Snapshot` so its `Vec<Bitset>` is freed immediately.
/// This keeps construction peak memory near the deduplicated footprint instead
/// of holding every tree's raw bitsets at once.
///
/// IDs are assigned in first-seen order, so feeding snapshots in tree-index
/// order yields exactly the same `bipartitions` ordering and `split_ids` as
/// interning the whole `Vec<Snapshot>` at once.
struct Interner {
    hasher: FxBuildHasher,
    table: HashTable<u32>,
    /// `(canonical fingerprint, smaller side's cardinality)` per split ID —
    /// what a candidate is matched against.
    keys: Vec<(Fingerprint, u32)>,
    bipartitions: Vec<Bitset>,
    snapshots: Vec<InternSnap>,
    words: usize,
    num_leaves: usize,
    rooted: bool,
    store_lengths: bool,
}

impl Interner {
    fn new(
        n_trees: usize,
        words: usize,
        num_leaves: usize,
        rooted: bool,
        store_lengths: bool,
    ) -> Self {
        Self {
            hasher: FxBuildHasher,
            table: HashTable::new(),
            keys: Vec::new(),
            bipartitions: Vec::new(),
            snapshots: Vec::with_capacity(n_trees),
            words,
            num_leaves,
            rooted,
            store_lengths,
        }
    }

    /// Intern one raw snapshot. When `store_lengths` is `false`, branch lengths
    /// are dropped (RF-only paths never read them).
    ///
    /// A split is matched on its 128-bit fingerprint *and* the cardinality of
    /// its smaller side, which is equal for both sides of a bipartition and
    /// rules out a slice of the collision space for one integer compare. Only a
    /// fingerprint the table has never seen materialises a `Bitset`, so that
    /// `O(subtree)` walk runs once per unique split rather than once per
    /// occurrence.
    ///
    /// With the `verify` feature on, every *match* is checked against the
    /// stored bitset and a collision panics rather than silently merging two
    /// distinct splits. CI runs the suite both ways.
    fn push(&mut self, snap: Snapshot) {
        // Bind disjoint fields to locals so the closures below can borrow the
        // table mutably while reading `keys`/`hasher`.
        let hasher = &self.hasher;
        let table = &mut self.table;
        let keys = &mut self.keys;
        let bipartitions = &mut self.bipartitions;
        let (num_leaves, rooted) = (self.num_leaves, self.rooted);

        let mut paired: Vec<(u32, f64)> = snap
            .parts
            .iter()
            .map(|part| {
                let candidate = (part.key, part.size.min(num_leaves as u32 - part.size));
                let hash = hasher.hash_one(part.key);
                let id = match table.find(hash, |&id| keys[id as usize] == candidate) {
                    Some(&id) => {
                        #[cfg(feature = "verify")]
                        assert_eq!(
                            snap.materialise(part, num_leaves, rooted),
                            bipartitions[id as usize],
                            "128-bit fingerprint collision: two distinct splits share a key"
                        );
                        id
                    }
                    None => {
                        let new_id = keys.len() as u32;
                        keys.push(candidate);
                        bipartitions.push(snap.materialise(part, num_leaves, rooted));
                        // register the new ID in the hash table
                        table.insert_unique(hash, new_id, |&id| {
                            hasher.hash_one(keys[id as usize].0)
                        });
                        new_id
                    }
                };
                (id, part.length)
            })
            .collect();

        paired.sort_unstable_by_key(|&(id, _)| id);
        // On RF-only paths (store_lengths == false) skip materialising the
        // lengths column entirely instead of unzipping then discarding it.
        let (split_ids, lengths): (Vec<u32>, Vec<f64>) = if self.store_lengths {
            paired.into_iter().unzip()
        } else {
            (paired.into_iter().map(|(id, _)| id).collect(), Vec::new())
        };
        self.snapshots.push(InternSnap { split_ids, lengths });
    }

    fn finish(self, leaf_names: Vec<String>) -> Snapshots {
        Snapshots {
            snapshots: self.snapshots,
            bipartitions: self.bipartitions,
            words_per_bitset: self.words,
            leaf_names,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Decode a native-endian `f64` matrix emitted by `build_branch_length_matrix`.
    fn decode_f64(bytes: &[u8]) -> Vec<f64> {
        bytes
            .as_chunks::<8>()
            .0
            .iter()
            .copied()
            .map(f64::from_ne_bytes)
            .collect()
    }

    /// Build one snapshot with a fresh label table — the interner is not
    /// involved, so any consistent labels will do.
    fn snapshot_of(newick: &str, rooted: bool) -> Snapshot {
        let tree = PhyloTree::from_newick(newick).unwrap();
        let labels = taxon_labels(tree.get_leaves().len());
        Snapshot::from_tree(&tree, rooted, &labels).unwrap()
    }

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
    /// node1's subtree = {A,B}, node2's = {C,D}. They are the two sides of the
    /// **same split**, so `min(h, h ^ total)` gives both the same fingerprint —
    /// XOR canonicalises them without touching either leaf set. Dedup collapses
    /// them into one entry and sums the branch lengths.
    #[test]
    fn test_depth3_tree_partitions() {
        let snap = snapshot_of("((A:1,B:2):1,(C:1,D:1):2);", false);

        // 4 pendant edges + 1 internal bipartition.
        assert_eq!(snap.parts.len(), 5);

        let internal: Vec<&Part> = snap.parts.iter().filter(|p| p.size > 1).collect();
        assert_eq!(
            internal.len(),
            1,
            "both sides of the root split must share one fingerprint"
        );
        assert_eq!(
            internal[0].length, 3.0,
            "deduplicated lengths should sum: 1.0 + 2.0 = 3.0"
        );
    }

    /// Both sides of a bipartition must materialise to the same `Bitset` — the
    /// side without leaf 0 — whichever occurrence the interner happened to see
    /// first. That is what keeps the public bipartition table unchanged by the
    /// switch to fingerprints.
    #[test]
    fn test_materialised_side_excludes_leaf_zero() {
        let snaps = Snapshots::from_newicks(&["((A:1,B:1):1,(C:1,D:1):1);"], false).unwrap();
        let internal: Vec<&Bitset> = snaps
            .bipartitions
            .iter()
            .filter(|b| b.count_ones() > 1)
            .collect();
        assert_eq!(internal.len(), 1);
        assert!(!internal[0].get(0), "canonical side must exclude leaf 0");
        assert_eq!(internal[0].0[0], 0b1100, "the {{C,D}} side");
    }

    /// Test that from_tree deduplicates root bipartitions in rooted binary trees.
    ///
    /// In a rooted binary tree ((A,B),(C,D)), both root children produce the
    /// same canonical bipartition ({C,D}|{A,B}).  Without deduplication, this
    /// causes RF distances to be inflated by +2 compared to the unrooted RF
    /// (e.g. R's phangorn RF.dist with rooted=FALSE).
    #[test]
    fn test_root_bipartition_dedup() {
        // Tree: ((A:1,B:1):1,(C:1,D:1):1);
        // Root has two internal children - classic root bipartition duplication case
        let snap = snapshot_of("((A:1,B:1):1,(C:1,D:1):1);", false);

        // For 4 leaves: 4 pendant edges + 1 internal bipartition = 5 entries.
        // With dedup, the duplicated root bipartition collapses to 1 internal entry.
        assert_eq!(
            snap.parts.len(),
            5,
            "Rooted 4-leaf binary tree should have 5 entries (4 pendant + 1 internal) after dedup, got {}",
            snap.parts.len()
        );
    }

    /// Split IDs must not depend on when the process started: the label table
    /// is seeded from a constant, so two runs agree bit for bit.
    #[test]
    fn test_taxon_labels_are_deterministic_and_distinct() {
        let a = taxon_labels(64);
        assert_eq!(a, taxon_labels(64), "labels must not vary between runs");
        assert_eq!(
            a.iter().collect::<HashSet<_>>().len(),
            64,
            "labels must be distinct"
        );
        assert_eq!(&a[..8], &taxon_labels(8)[..], "a prefix is a prefix");
    }

    /// Interning is by fingerprint, so a run over the same trees must produce
    /// the same bipartition table and the same split IDs every time.
    #[test]
    fn test_interning_is_reproducible() {
        let trees = [
            "(((A:1,B:1):1,C:1):1,(D:1,E:1):1);",
            "(((A:1,D:1):1,E:1):1,(B:1,C:1):1);",
        ];
        let (first, second) = (
            Snapshots::from_newicks(&trees, false).unwrap(),
            Snapshots::from_newicks(&trees, false).unwrap(),
        );
        assert_eq!(first.bipartitions, second.bipartitions);
        for (a, b) in first.snapshots.iter().zip(&second.snapshots) {
            assert_eq!(a.split_ids, b.split_ids);
        }
    }

    /// Two taxa is the degenerate case where a pendant edge's complement is
    /// the *other* pendant edge, so `min(h, h ^ total)` would merge them. They
    /// must stay two distinct splits.
    #[test]
    fn test_two_taxa_pendants_stay_distinct() {
        let snaps = Snapshots::from_newicks(&["(A:1,B:2);", "(A:3,B:4);"], false).unwrap();
        assert_eq!(snaps.bipartitions.len(), 2);
        assert_eq!(snaps.snapshots[0].split_ids, vec![0, 1]);
    }

    /// A caterpillar tree nests as deep as it has leaves. The bitset pass is
    /// iterative so depth costs heap, not stack.
    #[test]
    fn test_deep_caterpillar_tree() {
        const LEAVES: usize = 3000;
        let mut newick = format!("l{}:0.1", LEAVES - 1);
        for i in (0..LEAVES - 1).rev() {
            newick = format!("(l{i}:0.1,{newick}):0.1");
        }
        newick.push(';');

        let snaps = Snapshots::from_newicks(&[&newick, &newick], false).unwrap();
        assert_eq!(snaps.leaf_names.len(), LEAVES);
        // LEAVES pendant edges + LEAVES-3 internal bipartitions.
        assert_eq!(snaps.snapshots[0].split_ids.len(), 2 * LEAVES - 3);
        assert_eq!(snaps.pairwise_rf(None)[1], 0, "a tree against itself");
    }

    /// Test rooted vs unrooted mode partition counts and RF distances.
    #[test]
    fn test_rooted_vs_unrooted_partitions() {
        // Two trees means a two-tree `Snapshots`, read at the off-diagonal cell.
        let rf_pair = |a: &str, b: &str, rooted: bool| -> u32 {
            Snapshots::from_newicks(&[a, b], rooted)
                .unwrap()
                .pairwise_rf(None)[1]
        };

        const TREE1: &str = "((A:1,B:1):1,(C:1,D:1):1);";
        const TREE2: &str = "((A:1,C:1):1,(B:1,D:1):1);";
        const TREE1B: &str = "((B:2,A:2):2,(D:2,C:2):2);";

        // Unrooted mode: 4 pendant edges + 1 internal bipartition = 5 entries per tree.
        let snap1_u = snapshot_of(TREE1, false);
        let snap2_u = snapshot_of(TREE2, false);
        assert_eq!(
            snap1_u.parts.len(),
            5,
            "Unrooted: 4 pendant + 1 internal for 4-leaf tree"
        );
        assert_eq!(snap2_u.parts.len(), 5);
        assert_eq!(rf_pair(TREE1, TREE2, false), 2, "Unrooted RF = 2");

        // Rooted mode: 4 pendant edges + 2 internal clades = 6 entries per tree.
        let snap1_r = snapshot_of(TREE1, true);
        let snap2_r = snapshot_of(TREE2, true);
        assert_eq!(
            snap1_r.parts.len(),
            6,
            "Rooted: 4 pendant + 2 clades for 4-leaf tree"
        );
        assert_eq!(snap2_r.parts.len(), 6);
        assert_eq!(rf_pair(TREE1, TREE2, true), 4, "Rooted RF = 4");

        // Same topology written differently: both modes give RF = 0.
        assert_eq!(rf_pair(TREE1, TREE1B, false), 0, "Unrooted same topo = 0");
        assert_eq!(rf_pair(TREE1, TREE1B, true), 0, "Rooted same topo = 0");
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

    /// `build_branch_length_matrix` returns a byte buffer of the right size
    /// and decodes to the correct number of f64 values.
    #[test]
    fn test_build_branch_length_matrix_size() {
        // 2 trees, 4 leaves → 4 pendant + 2 internal = 6 bipartitions
        let snaps = Snapshots::from_newicks(
            &["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"],
            false,
        )
        .unwrap();

        let (bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
        let n_bip = snaps.bipartitions.len();
        assert_eq!(n_bip, 6);
        assert_eq!(col_to_bip_id.len(), 6);
        // 2 trees × 6 bipartitions × 8 bytes per f64
        assert_eq!(bytes.len(), 2 * 6 * 8);
    }

    /// Pendant-edge columns are always non-zero in every tree row, because every
    /// tree has every leaf.  Internal bipartition columns are non-zero only in the
    /// trees that contain that split.
    #[test]
    fn test_build_branch_length_matrix_values() {
        // T1: ((A:0.1, B:0.2):0.3, (C:0.4, D:0.5):0.6)
        //   pendant {A}=0.1, {B}=0.2, {C}=0.4, {D}=0.5
        //   internal {C,D}=0.6  (canonical: side not containing A)
        //   internal {A,B} pendant of the root edge: store as {C,D} complement → {C,D}=0.6
        //   Actually the root branch has length 0.3 and its bipartition is {A,B}|{C,D}.
        //   Canonical side (not containing A) = {C,D}, stored with length 0.3.
        //   The inner edge (C,D) has length 0.6, canonical side {C,D} (not A) → same bitset!
        //   So T1 has two distinct bitsets for these: {A,B} edge (len 0.3) = canonical {C,D}
        //   and {C,D} edge (len 0.6) = canonical {C,D}... wait, they ARE the same bitset.
        //   Actually {A,B}|{C,D} and {C,D}|{A,B} are the same bipartition — deduplicated.
        //
        // Use simpler trees that produce clearly distinct bipartitions:
        // T1: ((A:1,B:2):10,(C:3,D:4):20)   → pendant A=1,B=2,C=3,D=4; internal {C,D}=20, {A,B}→{C,D}... hmm
        //
        // Actually for 4 leaves unrooted:
        //   ((A,B),(C,D)) has ONE internal bipartition {A,B}|{C,D}, canonical = {C,D} (excludes A).
        //   ((A,C),(B,D)) has ONE internal bipartition {A,C}|{B,D}, canonical = {B,D}.
        // T1 = ((A:1,B:2):5,(C:3,D:4):6) → internal edge = 5 or 6?
        //   The internal edge connects (A,B) cluster to (C,D) cluster.
        //   In phylotree, each child of the root carries half the root-to-clade branch.
        //   For unrooted: the branch between the two clades has no single length in Newick.
        //   In practice, rapidtrees stores the branch length of the node's edge to its parent.
        //   The root's children each have their own branch length to root.
        //   So the internal bipartition {C,D} is stored with the branch length of the
        //   (C,D)-subtree's edge to root = 6. Wait, no: the root has 2 children.
        //   Child 1 = (A:1,B:2) with branch length 5 → bipartition {A,B}|{C,D}, stored as {C,D}
        //   Child 2 = (C:3,D:4) with branch length 6 → same bipartition {A,B}|{C,D}
        //   Both children define the same bipartition! We dedup → one entry.
        //   Which branch length is stored? The first one encountered (child 1 or child 2).
        //
        // The current code stores ONE entry per unique bipartition. For the root's two
        // children that share the same bipartition, only one branch length survives.
        // This is a known property; the test just verifies the output is self-consistent:
        // the L1 identity |bl[i]-bl[j]|.sum() == wrf[i,j] must hold.
        //
        // Use the wRF distance as the ground truth.
        let trees = [
            "((A:1,B:2):5,(C:3,D:4):6);",
            "((A:7,C:8):9,(B:10,D:11):12);",
        ];
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();
        let (bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
        let n_bip = col_to_bip_id.len();
        assert_eq!(n_bip, snaps.bipartitions.len());

        // Decode bytes to f64 matrix: shape (2, n_bip)
        let floats = decode_f64(&bytes);
        assert_eq!(floats.len(), 2 * n_bip);

        let row0 = &floats[..n_bip];
        let row1 = &floats[n_bip..];

        // L1 identity: sum(|bl[0] - bl[1]|) == wrf(T0, T1)
        let l1: f64 = row0.iter().zip(row1).map(|(a, b)| (a - b).abs()).sum();
        let wrf_matrix = snaps.pairwise_wrf(None);
        // wrf_matrix is flat row-major (2×2); wrf[0,1] is at index 1
        let wrf_01 = wrf_matrix[1];
        assert!(
            (l1 - wrf_01).abs() < 1e-9,
            "L1 identity failed: |bl[0]-bl[1]|.sum()={l1:.6} vs wrf={wrf_01:.6}"
        );

        // Pendant columns must be non-zero in both rows (every tree has every leaf).
        let n_leaves = snaps.leaf_names.len();
        let bip_arr: Vec<&Bitset> = col_to_bip_id
            .iter()
            .map(|&id| &snaps.bipartitions[id])
            .collect();
        let pendant_cols: Vec<usize> = bip_arr
            .iter()
            .enumerate()
            .filter(|(_, bip)| bip.count_ones() == 1)
            .map(|(col, _)| col)
            .collect();
        assert_eq!(
            pendant_cols.len(),
            n_leaves,
            "expected one pendant col per leaf"
        );
        for &col in &pendant_cols {
            assert!(row0[col] > 0.0, "pendant col {col} is 0 in tree 0");
            assert!(row1[col] > 0.0, "pendant col {col} is 0 in tree 1");
        }
    }

    /// `build_branch_length_matrix` and `build_presence_matrix` return the same
    /// `col_to_bip_id` ordering.  A column that is non-zero in the branch-length
    /// matrix must be 1 in the presence matrix, and vice versa.
    #[test]
    fn test_branch_length_matrix_consistent_with_presence_matrix() {
        let trees = [
            "((A:1,B:2):5,(C:3,D:4):6);",
            "((A:7,C:8):9,(B:10,D:11):12);",
        ];
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();

        let (bl_bytes, bl_col_to_bip) = snaps.build_branch_length_matrix();
        let (presence, pres_col_to_bip) = snaps.build_presence_matrix();

        // Column ordering must be identical.
        assert_eq!(bl_col_to_bip, pres_col_to_bip);

        let n_bip = bl_col_to_bip.len();
        let bl = decode_f64(&bl_bytes);

        // For every (tree, bipartition) cell: bl > 0 ↔ presence == 1.
        for tree in 0..2 {
            for col in 0..n_bip {
                let has_bl = bl[tree * n_bip + col] > 0.0;
                let has_pres = presence[tree * n_bip + col] == 1;
                assert_eq!(
                    has_bl, has_pres,
                    "tree={tree} col={col}: bl>0={has_bl} but presence={has_pres}"
                );
            }
        }
    }

    /// Empty-bipartitions edge case: two identical 2-leaf trees produce a
    /// degenerate snapshot with only pendant edges (no internal bipartitions).
    /// `build_branch_length_matrix` must not panic.
    #[test]
    fn test_build_branch_length_matrix_two_leaves() {
        let snaps = Snapshots::from_newicks(&["(A:1,B:2);", "(A:3,B:4);"], false).unwrap();

        let (bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
        // 2 pendant edges, 0 internal
        assert_eq!(col_to_bip_id.len(), 2);
        assert_eq!(bytes.len(), 2 * 2 * 8); // 2 trees × 2 bips × 8 bytes

        let floats = decode_f64(&bytes);
        // Both pendant columns must be non-zero in both trees.
        assert!(floats.iter().all(|&v| v > 0.0));
    }

    /// Three-tree case: a bipartition absent from the middle tree must have a
    /// 0.0 branch length in that row and non-zero in the others.
    #[test]
    #[allow(clippy::erasing_op)]
    fn test_build_branch_length_matrix_absent_split_is_zero() {
        // T0 and T2 share bipartition {C,D}, T1 does not.
        let trees = [
            "((A:1,B:2):5,(C:3,D:4):6);",       // internal {C,D}
            "((A:7,C:8):9,(B:10,D:11):12);",    // internal {B,D}
            "((A:13,B:14):15,(C:16,D:17):18);", // internal {C,D} again
        ];
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();
        let (bytes, col_to_bip_id) = snaps.build_branch_length_matrix();
        let n_bip = col_to_bip_id.len();

        let floats = decode_f64(&bytes);

        // Find the column for the {C,D} bipartition (bits 2 and 3 set = 0b1100).
        let cd_col = col_to_bip_id.iter().position(|&id| {
            snaps.bipartitions[id].count_ones() == 2 && snaps.bipartitions[id].0[0] == 0b1100
        });
        if let Some(col) = cd_col {
            assert!(floats[0 * n_bip + col] > 0.0, "T0 should have {{C,D}}");
            assert_eq!(floats[n_bip + col], 0.0, "T1 should lack {{C,D}}");
            assert!(floats[2 * n_bip + col] > 0.0, "T2 should have {{C,D}}");
        }

        // L1 identity must hold for every pair.
        let wrf = snaps.pairwise_wrf(None);
        for i in 0..3 {
            for j in 0..3 {
                let l1: f64 = (0..n_bip)
                    .map(|c| (floats[i * n_bip + c] - floats[j * n_bip + c]).abs())
                    .sum();
                assert!(
                    (l1 - wrf[i * 3 + j]).abs() < 1e-9,
                    "L1 identity failed for ({i},{j}): l1={l1:.6} wrf={:.6}",
                    wrf[i * 3 + j]
                );
            }
        }
    }

    /// Build `Snapshots` from plain newicks with an explicit `store_lengths` flag.
    fn snaps_opts(newicks: &[&str], rooted: bool, store_lengths: bool) -> Snapshots {
        let empty: HashMap<String, String> = HashMap::new();
        Snapshots::from_newick_iter_opts(
            newicks.iter().map(|&n| (n, &empty)),
            rooted,
            store_lengths,
        )
        .unwrap()
    }

    /// RF distances are identical whether or not branch lengths are stored, and the
    /// no-lengths path leaves every `InternSnap.lengths` empty while keeping split IDs.
    #[test]
    fn test_rf_path_without_lengths_matches() {
        let trees = [
            "((A:1,B:2):5,(C:3,D:4):6);",
            "((A:7,C:8):9,(B:10,D:11):12);",
            "((A:13,D:1):2,(B:3,C:4):5);",
        ];

        let with_len = snaps_opts(&trees, false, true);
        let no_len = snaps_opts(&trees, false, false);

        // Identical interning: same bipartition count and per-tree split IDs.
        assert_eq!(with_len.bipartitions.len(), no_len.bipartitions.len());
        for (a, b) in with_len.snapshots.iter().zip(&no_len.snapshots) {
            assert_eq!(a.split_ids, b.split_ids, "split IDs must match");
        }

        // No-lengths path drops the lengths vector; default path keeps one per split.
        for snap in &no_len.snapshots {
            assert!(snap.lengths.is_empty(), "RF path must not store lengths");
        }
        for snap in &with_len.snapshots {
            assert_eq!(snap.lengths.len(), snap.split_ids.len());
        }

        assert_eq!(with_len.pairwise_rf(None), no_len.pairwise_rf(None));
    }

    /// `intern` assigns strictly-ascending, deduplicated split IDs, and identical
    /// topologies share the exact same interned IDs.
    #[test]
    fn test_intern_split_ids_sorted_and_deduped() {
        let trees = [
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,C:1):1,(B:1,D:1):1);",
        ];
        let snaps = snaps_opts(&trees, false, true);

        for snap in &snaps.snapshots {
            assert!(
                snap.split_ids.windows(2).all(|w| w[0] < w[1]),
                "split IDs must be strictly ascending and unique: {:?}",
                snap.split_ids
            );
            for &id in &snap.split_ids {
                assert!((id as usize) < snaps.bipartitions.len());
            }
        }

        assert_eq!(
            snaps.snapshots[0].split_ids, snaps.snapshots[1].split_ids,
            "identical topologies must intern to identical split IDs"
        );
    }

    /// Presence matrix and bipartition-clade bytes are unaffected by `store_lengths`
    /// (the RF-with-snapshots export path passes `false`).
    #[test]
    fn test_presence_export_independent_of_lengths() {
        let trees = [
            "((A:1,B:2):5,(C:3,D:4):6);",
            "((A:7,C:8):9,(B:10,D:11):12);",
            "((A:13,B:14):15,(C:16,D:17):18);",
        ];

        let with_len = snaps_opts(&trees, false, true);
        let no_len = snaps_opts(&trees, false, false);

        let (pres_a, cols_a) = with_len.build_presence_matrix();
        let (pres_b, cols_b) = no_len.build_presence_matrix();
        assert_eq!(pres_a, pres_b, "presence matrix must be identical");
        assert_eq!(cols_a, cols_b, "column ordering must be identical");

        let bip_a = with_len.build_bipartition_bytes(&cols_a);
        let bip_b = no_len.build_bipartition_bytes(&cols_b);
        assert_eq!(bip_a, bip_b, "bipartition clade bytes must be identical");
    }
}
