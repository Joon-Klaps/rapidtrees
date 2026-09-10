//! Behavioural tests for the snapshot pipeline.
//!
//! Kept together rather than split across the submodules: nearly all of them
//! drive the whole path from newick to interned split IDs, and they share the
//! `snapshot_of` / `snaps_opts` helpers.

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
    let mut names: Vec<String> = tree
        .get_leaves()
        .iter()
        .filter_map(|id| tree.get(id).ok()?.name.clone())
        .collect();
    names.sort_unstable();
    let labels = taxon_labels(names.len());
    Snapshot::from_tree(&tree, rooted, &labels, &build_leaf_index(&names)).unwrap()
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

/// Both sides of a bipartition must record the same leaf set — the side without
/// leaf 0 — whichever occurrence the interner happened to see first. That is
/// what keeps the exported bipartition table unchanged by the switch to
/// fingerprints.
#[test]
fn test_materialised_side_excludes_leaf_zero() {
    let snaps = Snapshots::from_newicks(&["((A:1,B:1):1,(C:1,D:1):1);"], false).unwrap();
    let internal: Vec<&[u32]> = (0..snaps.clades.len())
        .map(|id| snaps.clades.get(id))
        .filter(|c| c.len() > 1)
        .collect();
    assert_eq!(internal.len(), 1);
    assert!(
        !internal[0].contains(&0),
        "canonical side must exclude leaf 0"
    );
    assert_eq!(internal[0], &[2, 3], "the {{C,D}} side");
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
    assert_eq!(first.clades.len(), second.clades.len());
    for id in 0..first.clades.len() {
        assert_eq!(first.clades.get(id), second.clades.get(id));
    }
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
    assert_eq!(snaps.clades.len(), 2);
    assert_eq!(snaps.snapshots[0].split_ids, vec![0, 1]);
}

/// A caterpillar tree nests as deep as it has leaves. The DFS is iterative so
/// depth costs heap, not stack.
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

/// An asymmetric tree: every internal edge must yield a *distinct* canonical
/// split, and each stored side must exclude leaf 0.
///
/// ```text
///        ┌── A                 leaves sorted: A=0 B=1 C=2 D=3 E=4
///     ┌──┤
///     │  │  ┌── B              node3 = {B,C}      → kept as-is (no A)
///     │  └──┤                  node2 = {A,B,C}    → stored as {D,E}
///  ───┤     └── C              node1 = {A,B,C,D}  → stored as {E}
///     │  ┌── D
///     └──┤
///        └── E
/// ```
///
/// This used to be asserted by hand-building bitsets and flipping them. It now
/// runs through the real path, so it tests the code rather than a description
/// of it.
#[test]
fn asymmetric_tree_yields_distinct_canonical_splits() {
    let snaps = Snapshots::from_newicks(&["((A:1,(B:1,C:1):1):1,(D:1,E:1):1);"], false).unwrap();

    let internal: Vec<Vec<u32>> = (0..snaps.clades.len())
        .map(|id| snaps.clades.get(id).to_vec())
        .filter(|c| c.len() > 1)
        .collect();

    for clade in &internal {
        assert!(
            !clade.contains(&0),
            "canonical side must exclude leaf 0, got {clade:?}"
        );
    }

    let mut sorted = internal.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(
        sorted.len(),
        internal.len(),
        "every internal edge must give a distinct split, got {internal:?}"
    );

    assert!(internal.contains(&vec![1, 2]), "{{B,C}}: {internal:?}");
    assert!(internal.contains(&vec![3, 4]), "{{D,E}}: {internal:?}");
}

/// Canonicalization: both sides of one split must land on the same entry.
///
/// `((A,B),(C,D))` gives `{A,B}` from one root child and `{C,D}` from the other
/// — the *same* bipartition seen from opposite ends. Without canonicalization
/// they would be two table entries and RF would be inflated by 2 against
/// phangorn's `RF.dist(rooted = FALSE)`.
#[test]
fn both_sides_of_a_split_share_one_entry() {
    let snaps = Snapshots::from_newicks(&["((A:1,B:2):1,(C:1,D:1):2);"], false).unwrap();

    let internal: Vec<&[u32]> = (0..snaps.clades.len())
        .map(|id| snaps.clades.get(id))
        .filter(|c| c.len() > 1)
        .collect();

    assert_eq!(
        internal.len(),
        1,
        "the two root children are one split, got {internal:?}"
    );
    assert_eq!(
        internal[0],
        &[2, 3],
        "stored as {{C,D}} — the side without A"
    );
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
    // The same topology written with its leaves in a different order. The parser
    // hands out different node IDs for each, so anything keyed on node ID would
    // see two different trees.
    const AS_WRITTEN: &str = "((Human:1,Chimp:1):1,Gorilla:1);";
    const REORDERED: &str = "(Gorilla:1,(Chimp:1,Human:1):1);";

    let snaps = Snapshots::from_newicks(&[AS_WRITTEN, REORDERED], false).unwrap();

    assert_eq!(
        snaps.leaf_names,
        vec!["Chimp", "Gorilla", "Human"],
        "bit indices come from alphabetical names, not parse order"
    );
    assert_eq!(
        snaps.snapshots[0].split_ids, snaps.snapshots[1].split_ids,
        "same taxa and topology must intern identically whatever the node IDs"
    );
    assert_eq!(snaps.pairwise_rf(None)[1], 0);
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
    assert_eq!(
        snaps.clades.get(col_to_bip_id[4]),
        &[1, 3],
        "col 4 should be {{B,D}}"
    );
    assert_eq!(
        snaps.clades.get(col_to_bip_id[5]),
        &[2, 3],
        "col 5 should be {{C,D}}"
    );

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
    let n_bip = snaps.clades.len();
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
    assert_eq!(n_bip, snaps.clades.len());

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
    let bip_arr: Vec<&[u32]> = col_to_bip_id
        .iter()
        .map(|&id| snaps.clades.get(id))
        .collect();
    let pendant_cols: Vec<usize> = bip_arr
        .iter()
        .enumerate()
        .filter(|(_, bip)| bip.len() == 1)
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
    let cd_col = col_to_bip_id
        .iter()
        .position(|&id| snaps.clades.get(id) == [2, 3]);
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
    let retain = Retain {
        lengths: store_lengths,
        bipartitions: true,
    };
    Snapshots::from_newick_iter_opts(newicks.iter().map(|&n| (n, &empty)), rooted, retain).unwrap()
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
    assert_eq!(with_len.clades.len(), no_len.clades.len());
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
            assert!((id as usize) < snaps.clades.len());
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

/// Leaf-set validation runs against the shared name → bit table rather than a
/// per-tree `HashSet<String>`. These are the three ways a tree can disagree with
/// tree 0, and all of them must still be caught.
#[test]
fn leaf_set_mismatches_are_rejected() {
    const REF: &str = "((A:1,B:1):1,(C:1,D:1):1);";

    let unknown_taxon = Snapshots::from_newicks(&[REF, "((A:1,B:1):1,(C:1,Z:1):1);"], false);
    assert!(
        unknown_taxon.is_err(),
        "a taxon absent from tree 0 must be rejected"
    );

    let missing_taxon = Snapshots::from_newicks(&[REF, "((A:1,B:1):1,C:1);"], false);
    assert!(
        missing_taxon.is_err(),
        "a tree missing one of tree 0's taxa must be rejected"
    );

    // The subtle one: same leaf count, every name known, but `A` twice and no
    // `D`. A plain count-and-membership check would wave this through.
    let duplicated = Snapshots::from_newicks(&[REF, "((A:1,B:1):1,(C:1,A:1):1);"], false);
    assert!(
        duplicated.is_err(),
        "a name repeated within one tree must be rejected"
    );
}
