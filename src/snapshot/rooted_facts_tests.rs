//! Tests for optional rooted-tree fact collection, interning, and export.

use super::rooted_facts::{ROOT_ID, RawCladeRef, RawRootedFacts};
use super::*;
use crate::snapshot::fingerprint::Fingerprint;
use std::collections::HashSet;

/// Decode a native-endian `f32` buffer emitted by a snapshot exporter.
fn decode_f32(bytes: &[u8]) -> Vec<f32> {
    let (chunks, remainder) = bytes.as_chunks::<4>();
    assert!(remainder.is_empty(), "f32 buffer must contain whole values");
    chunks.iter().copied().map(f32::from_ne_bytes).collect()
}

/// Decode a native-endian `u32` buffer emitted by a snapshot exporter.
fn decode_u32(bytes: &[u8]) -> Vec<u32> {
    let (chunks, remainder) = bytes.as_chunks::<4>();
    assert!(remainder.is_empty(), "u32 buffer must contain whole values");
    chunks.iter().copied().map(u32::from_ne_bytes).collect()
}

/// Build one snapshot with a fresh label table for memory-estimate checks.
fn snapshot_of(newick: &str, rooted: bool) -> Snapshot {
    snapshot_result(newick, rooted, false).unwrap()
}

fn snapshot_result(
    newick: &str,
    rooted: bool,
    require_explicit_lengths: bool,
) -> Result<Snapshot, String> {
    let no_translate = HashMap::new();
    let mut names = newick::leaf_names(newick, &no_translate)?;
    names.sort_unstable();
    let labels = taxon_labels(names.len());
    let leaf_index = build_leaf_index(&names);
    let run = newick::RunTables {
        leaf_index: &leaf_index,
        labels: &labels,
        total: labels.iter().fold(0, |acc, &label| acc ^ label),
        rooted,
        require_explicit_lengths,
    };
    newick::snapshot(newick, &no_translate, 0, &run)
}

/// Collect optional rooted facts from one direct-parser snapshot.
fn rooted_facts_of(newick: &str) -> Result<(RawRootedFacts, Vec<Fingerprint>), String> {
    let no_translate = HashMap::new();
    let mut names = newick::leaf_names(newick, &no_translate)?;
    names.sort_unstable();
    let labels = taxon_labels(names.len());
    let snapshot = snapshot_result(newick, true, true)?;
    let facts = RawRootedFacts::from_snapshot(&snapshot)?;
    Ok((facts, labels))
}

/// Build rooted snapshots with the optional facts sidecar explicitly selected.
fn rooted_snapshots_opts(newicks: &[&str], retain_rooted_facts: bool) -> Result<Snapshots, String> {
    let empty: HashMap<String, String> = HashMap::new();
    Snapshots::from_newick_iter_opts(
        newicks.iter().map(|&newick| (newick, &empty)),
        true,
        Retain {
            lengths: false,
            bipartitions: true,
            rooted_facts: retain_rooted_facts,
        },
    )
}

/// The optional collector uses max root-to-tip distance as root height and
/// retains heights for every exact non-root clade, including singleton tips.
#[test]
fn rooted_facts_collect_non_ultrametric_heights() {
    let (facts, labels) =
        rooted_facts_of("((A:1,B:2):3,(C:4,D:5):6);").expect("collect rooted facts");

    assert_eq!(facts.root_height, 11.0_f32);
    assert_eq!(facts.nodes.len(), 6);

    let heights: std::collections::HashMap<Fingerprint, f64> = facts
        .nodes
        .iter()
        .map(|fact| (fact.clade.key, f64::from(fact.height)))
        .collect();
    assert_eq!(heights[&labels[0]], 7.0, "A");
    assert_eq!(heights[&labels[1]], 6.0, "B");
    assert_eq!(heights[&(labels[0] ^ labels[1])], 8.0, "{{A,B}}");
    assert_eq!(heights[&labels[2]], 1.0, "C");
    assert_eq!(heights[&labels[3]], 0.0, "D");
    assert_eq!(heights[&(labels[2] ^ labels[3])], 5.0, "{{C,D}}");
}

/// Quantizing branch lengths or cumulative distances before subtraction would
/// round 100000001 to 100000000 and incorrectly make the {B,C} height zero.
#[test]
fn rooted_facts_quantize_only_completed_heights() {
    let (facts, labels) =
        rooted_facts_of("(A:100000001,(B:1,C:1):100000000);").expect("collect rooted facts");
    let bc = labels[1] ^ labels[2];
    let height = facts
        .nodes
        .iter()
        .find_map(|fact| (fact.clade.key == bc).then_some(fact.height))
        .expect("{B,C} clade height");

    assert_eq!(height, 1.0_f32);
}

/// Every emitted split is a parent paired with its two immediate source-tree
/// children.  The root is represented by `None`, not by a synthetic clade.
#[test]
fn rooted_facts_collect_only_directly_observed_splits() {
    let (facts, labels) =
        rooted_facts_of("((A:1,B:2):3,(C:4,D:5):6);").expect("collect rooted facts");

    let leaf = |index: usize| RawCladeRef {
        key: labels[index],
        size: 1,
    };
    let clade = |a: usize, b: usize| RawCladeRef {
        key: labels[a] ^ labels[b],
        size: 2,
    };
    let pair = |mut children: [RawCladeRef; 2]| {
        children.sort_unstable();
        children
    };

    let expected = HashSet::from([
        (None, pair([clade(0, 1), clade(2, 3)])),
        (Some(clade(0, 1)), pair([leaf(0), leaf(1)])),
        (Some(clade(2, 3)), pair([leaf(2), leaf(3)])),
    ]);
    let observed: HashSet<_> = facts
        .splits
        .iter()
        .map(|split| (split.parent, split.children))
        .collect();

    assert_eq!(observed, expected);
}

#[test]
fn rooted_facts_require_explicit_non_root_branch_lengths() {
    let error = rooted_facts_of("(A:1,B);").expect_err("missing length must fail");
    assert!(
        error.contains("missing an explicit branch length"),
        "unexpected error: {error}"
    );
}

#[test]
fn rooted_facts_require_strictly_binary_internal_nodes() {
    let error = rooted_facts_of("(A:1,B:1,C:1);").expect_err("polytomy must fail");
    assert!(
        error.contains("must have exactly two children; found 3"),
        "unexpected error: {error}"
    );
}

#[test]
fn rooted_facts_reject_non_finite_cumulative_distances() {
    let error = rooted_facts_of("((A:1e308,B:1e308):1e308,C:1);")
        .expect_err("overflowed cumulative distance must fail");
    assert!(
        error.contains("non-finite cumulative root distance"),
        "unexpected error: {error}"
    );
}

#[test]
fn rooted_facts_reject_heights_outside_float32_range() {
    let error = rooted_facts_of("(A:1e39,B:1);")
        .expect_err("finite f64 height outside float32 range must fail");
    assert!(
        error.contains("outside the finite float32 range"),
        "unexpected error: {error}"
    );
}

/// Negative lengths are not silently rejected: RapidTrees accepts any finite
/// length here and validates representation and arithmetic.
#[test]
fn rooted_facts_accept_finite_negative_branch_lengths() {
    let (facts, labels) = rooted_facts_of("(A:-1,B:2);").expect("finite lengths");
    let heights: std::collections::HashMap<Fingerprint, f64> = facts
        .nodes
        .iter()
        .map(|fact| (fact.clade.key, f64::from(fact.height)))
        .collect();
    assert_eq!(facts.root_height, 2.0_f32);
    assert_eq!(heights[&labels[0]], 3.0);
    assert_eq!(heights[&labels[1]], 0.0);
}

/// The facts traversal must not consume call-stack depth on caterpillar trees.
#[test]
fn rooted_facts_handle_deep_caterpillar_iteratively() {
    const LEAVES: usize = 3000;
    let mut newick = format!("l{}:0.1", LEAVES - 1);
    for i in (0..LEAVES - 1).rev() {
        newick = format!("(l{i}:0.1,{newick}):0.1");
    }
    newick.push(';');

    let (facts, _) = rooted_facts_of(&newick).expect("collect deep rooted facts");
    assert_eq!(facts.nodes.len(), 2 * LEAVES - 2);
    assert_eq!(facts.splits.len(), LEAVES - 1);
}

/// Facts are resolved against the IDs assigned by the ordinary rooted
/// snapshot, remain aligned with input tree rows, and keep the root implicit.
#[test]
fn rooted_facts_sidecar_is_interned_and_tree_aligned() {
    let trees = ["((A:1,B:2):3,(C:4,D:5):6);", "((A:2,C:3):4,(B:5,D:6):7);"];
    let snaps = rooted_snapshots_opts(&trees, true).expect("facts-enabled snapshots");
    let store = snaps.rooted_facts.as_ref().expect("rooted facts sidecar");

    assert_eq!(store.trees.len(), trees.len());
    assert_eq!(
        store
            .trees
            .iter()
            .map(|facts| f64::from(facts.root_height))
            .collect::<Vec<_>>(),
        vec![11.0, 13.0]
    );

    for (row, (snapshot, facts)) in snaps.snapshots.iter().zip(&store.trees).enumerate() {
        assert_eq!(facts.node_heights.len(), 6, "row {row}");
        assert_eq!(facts.split_ids.len(), 3, "row {row}");
        assert!(facts.split_ids.is_sorted(), "row {row}");

        let root_splits = facts
            .split_ids
            .iter()
            .map(|&id| store.split_table[id as usize])
            .filter(|split| split[0] == ROOT_ID)
            .count();
        assert_eq!(root_splits, 1, "row {row}");
        for [parent, left, right] in facts
            .split_ids
            .iter()
            .map(|&id| store.split_table[id as usize])
        {
            if parent != ROOT_ID {
                assert!(snapshot.split_ids.binary_search(&parent).is_ok());
            }
            assert!(snapshot.split_ids.binary_search(&left).is_ok());
            assert!(snapshot.split_ids.binary_search(&right).is_ok());
            assert!(left < right, "children must be canonicalized");
        }
    }

    let unique_splits = store.split_table.iter().copied().collect::<HashSet<_>>();
    assert_eq!(unique_splits.len(), store.split_table.len());

    let first = &store.trees[0];
    let height_for = |wanted: &[u32]| {
        snaps.snapshots[0]
            .split_ids
            .iter()
            .zip(&first.node_heights)
            .find_map(|(&id, &height)| (snaps.clades.get(id as usize) == wanted).then_some(height))
            .unwrap_or_else(|| panic!("missing clade {wanted:?}"))
    };
    assert_eq!(height_for(&[0]), 7.0_f32);
    assert_eq!(height_for(&[0, 1]), 8.0_f32);
    assert_eq!(height_for(&[2, 3]), 5.0_f32);
}

/// Rooted facts use the stable clade-table columns, keep heights aligned with
/// sparse presence rows, and dictionary-encode directly observed child splits.
#[test]
fn rooted_fact_export_has_stable_columns_and_fixed_shapes() {
    let trees = ["((A:1,B:2):3,(C:4,D:5):6);", "((A:2,C:3):4,(B:5,D:6):7);"];
    let snaps = rooted_snapshots_opts(&trees, true).expect("facts-enabled snapshots");
    let (presence, col_to_bip_id) = snaps.build_presence_matrix();
    let facts = snaps
        .build_rooted_fact_buffers(&col_to_bip_id)
        .expect("export rooted facts");

    let n_trees = trees.len();
    let n_leaves = snaps.leaf_names.len();
    let n_clades = col_to_bip_id.len();
    assert_eq!(facts.root_column, n_clades as u32);
    assert_eq!(facts.nodes_per_tree, 2 * n_leaves - 2);
    assert_eq!(facts.splits_per_tree, n_leaves - 1);
    assert_eq!(
        facts.clade_columns.len(),
        n_trees * facts.nodes_per_tree * size_of::<u32>()
    );
    assert_eq!(
        facts.node_heights.len(),
        n_trees * facts.nodes_per_tree * size_of::<f32>()
    );
    assert_eq!(facts.root_heights.len(), n_trees * size_of::<f32>());
    assert_eq!(
        facts.split_ids.len(),
        n_trees * facts.splits_per_tree * size_of::<u32>()
    );
    assert_eq!(
        facts.split_table.len(),
        facts.n_observed_splits * 3 * size_of::<u32>()
    );

    let clade_columns = decode_u32(&facts.clade_columns);
    let node_heights = decode_f32(&facts.node_heights);
    let root_heights = decode_f32(&facts.root_heights);
    let split_ids = decode_u32(&facts.split_ids);
    let split_table = decode_u32(&facts.split_table)
        .chunks_exact(3)
        .map(|triple| [triple[0], triple[1], triple[2]])
        .collect::<Vec<_>>();
    assert_eq!(root_heights, vec![11.0_f32, 13.0_f32]);
    assert_eq!(split_table.len(), facts.n_observed_splits);
    assert!(
        split_table.windows(2).all(|pair| pair[0] < pair[1]),
        "split table must be sorted and deduplicated"
    );

    let clade_at = |column: u32| snaps.clades.get(col_to_bip_id[column as usize]);
    for tree_index in 0..n_trees {
        let node_start = tree_index * facts.nodes_per_tree;
        let tree_clades = &clade_columns[node_start..node_start + facts.nodes_per_tree];
        assert!(tree_clades.windows(2).all(|pair| pair[0] < pair[1]));
        let present_columns = presence[tree_index * n_clades..(tree_index + 1) * n_clades]
            .iter()
            .enumerate()
            .filter_map(|(column, &present)| (present != 0).then_some(column as u32))
            .collect::<Vec<_>>();
        assert_eq!(tree_clades, present_columns);
        for &column in tree_clades {
            assert!(column < facts.root_column);
        }

        let split_start = tree_index * facts.splits_per_tree;
        let tree_split_ids = &split_ids[split_start..split_start + facts.splits_per_tree];
        assert!(tree_split_ids.windows(2).all(|pair| pair[0] < pair[1]));
        let mut root_splits = 0;
        for &[parent, left, right] in tree_split_ids.iter().map(|&id| &split_table[id as usize]) {
            assert!(left < right, "child columns must be canonicalized");
            assert!(left < facts.root_column && right < facts.root_column);
            assert_eq!(presence[tree_index * n_clades + left as usize], 1);
            assert_eq!(presence[tree_index * n_clades + right as usize], 1);

            let mut child_union = clade_at(left)
                .iter()
                .chain(clade_at(right))
                .copied()
                .collect::<Vec<_>>();
            child_union.sort_unstable();
            child_union.dedup();

            let parent_clade = if parent == facts.root_column {
                root_splits += 1;
                (0..n_leaves as u32).collect::<Vec<_>>()
            } else {
                assert!(parent < facts.root_column);
                assert_eq!(presence[tree_index * n_clades + parent as usize], 1);
                clade_at(parent).to_vec()
            };
            assert_eq!(
                child_union, parent_clade,
                "row {tree_index} exports a split that is not a direct child partition"
            );
        }
        assert_eq!(root_splits, 1, "row {tree_index}");
    }

    let column_for = |wanted: &[u32]| {
        col_to_bip_id
            .iter()
            .position(|&id| snaps.clades.get(id) == wanted)
            .expect("exported clade column") as u32
    };
    let a = column_for(&[0]);
    let b = column_for(&[1]);
    let c = column_for(&[2]);
    let d = column_for(&[3]);
    let ab = column_for(&[0, 1]);
    let cd = column_for(&[2, 3]);
    let ordered = |left: u32, right: u32| [left.min(right), left.max(right)];
    let [root_left, root_right] = ordered(ab, cd);
    let [ab_left, ab_right] = ordered(a, b);
    let [cd_left, cd_right] = ordered(c, d);
    let expected_first_splits = HashSet::from([
        [facts.root_column, root_left, root_right],
        [ab, ab_left, ab_right],
        [cd, cd_left, cd_right],
    ]);
    let actual_first_splits = split_ids[..facts.splits_per_tree]
        .iter()
        .map(|&id| split_table[id as usize])
        .collect::<HashSet<_>>();
    assert_eq!(actual_first_splits, expected_first_splits);

    let first_nodes = &clade_columns[..facts.nodes_per_tree];
    let first_heights = &node_heights[..facts.nodes_per_tree];
    let height_at = |column: u32| {
        let offset = first_nodes
            .iter()
            .position(|&candidate| candidate == column)
            .expect("node column");
        first_heights[offset]
    };
    assert_eq!(height_at(a), 7.0_f32);
    assert_eq!(height_at(ab), 8.0_f32);
    assert_eq!(height_at(cd), 5.0_f32);
}

/// ID-to-column conversion is driven by the supplied mapping, and children
/// are canonicalized only after that conversion.
#[test]
fn rooted_fact_export_translates_before_canonicalizing_children() {
    let trees = ["((A:1,B:2):3,(C:4,D:5):6);"];
    let snaps = rooted_snapshots_opts(&trees, true).unwrap();
    let n_clades = snaps.n_distinct_splits();
    let reversed_columns = (0..n_clades).rev().collect::<Vec<_>>();
    let buffers = snaps
        .build_rooted_fact_buffers(&reversed_columns)
        .expect("export on supplied columns");
    let exported_nodes = decode_u32(&buffers.clade_columns);
    let exported_heights = decode_f32(&buffers.node_heights);
    let exported_split_ids = decode_u32(&buffers.split_ids);
    let exported_split_table = decode_u32(&buffers.split_table)
        .chunks_exact(3)
        .map(|triple| [triple[0], triple[1], triple[2]])
        .collect::<Vec<_>>();
    let store = snaps.rooted_facts.as_ref().unwrap();
    let stored = &store.trees[0];

    let column_for_id = |id: u32| n_clades as u32 - 1 - id;
    let mut expected_nodes = snaps.snapshots[0]
        .split_ids
        .iter()
        .zip(&stored.node_heights)
        .map(|(&id, &height)| (column_for_id(id), height))
        .collect::<Vec<_>>();
    expected_nodes.sort_unstable_by_key(|&(column, _)| column);
    assert_eq!(
        exported_nodes,
        expected_nodes
            .iter()
            .map(|&(column, _)| column)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        exported_heights,
        expected_nodes
            .iter()
            .map(|&(_, height)| height)
            .collect::<Vec<_>>()
    );

    let export_split = |[parent, left, right]: [u32; 3]| {
        let mut children = [column_for_id(left), column_for_id(right)];
        children.sort_unstable();
        [
            if parent == ROOT_ID {
                buffers.root_column
            } else {
                column_for_id(parent)
            },
            children[0],
            children[1],
        ]
    };
    let mut expected_table = store
        .split_table
        .iter()
        .copied()
        .map(export_split)
        .collect::<Vec<_>>();
    expected_table.sort_unstable();
    assert_eq!(exported_split_table, expected_table);

    let mut expected_tree_splits = stored
        .split_ids
        .iter()
        .map(|&id| export_split(store.split_table[id as usize]))
        .collect::<Vec<_>>();
    expected_tree_splits.sort_unstable();
    let actual_tree_splits = exported_split_ids
        .iter()
        .map(|&id| exported_split_table[id as usize])
        .collect::<Vec<_>>();
    assert_eq!(actual_tree_splits, expected_tree_splits);
    assert!(
        exported_split_table
            .iter()
            .all(|triple| triple[1] < triple[2])
    );
}

#[test]
fn rooted_fact_export_validates_prerequisites_and_dimensions() {
    let trees = ["((A:1,B:2):3,(C:4,D:5):6);"];

    let without_facts = rooted_snapshots_opts(&trees, false).unwrap();
    let (_, columns) = without_facts.build_presence_matrix();
    let error = without_facts
        .build_rooted_fact_buffers(&columns)
        .expect_err("missing sidecar must fail");
    assert!(error.contains("were not retained"));

    let with_facts = rooted_snapshots_opts(&trees, true).unwrap();
    let (_, columns) = with_facts.build_presence_matrix();
    let mut duplicate = columns.clone();
    duplicate[0] = duplicate[1];
    let error = with_facts
        .build_rooted_fact_buffers(&duplicate)
        .expect_err("duplicate mapping must fail");
    assert!(error.contains("more than once"));

    let mut malformed = rooted_snapshots_opts(&trees, true).unwrap();
    malformed.rooted_facts.as_mut().unwrap().trees[0]
        .node_heights
        .pop();
    let (_, columns) = malformed.build_presence_matrix();
    let error = malformed
        .build_rooted_fact_buffers(&columns)
        .expect_err("ragged node rows must fail");
    assert!(error.contains("expected 6 of each"));
}

#[test]
fn rooted_fact_export_handles_an_empty_collection() {
    let snaps = rooted_snapshots_opts(&[], true).unwrap();
    let (_, columns) = snaps.build_presence_matrix();
    let facts = snaps.build_rooted_fact_buffers(&columns).unwrap();

    assert_eq!(facts.root_column, 0);
    assert_eq!(facts.nodes_per_tree, 0);
    assert_eq!(facts.splits_per_tree, 0);
    assert_eq!(facts.n_observed_splits, 0);
    assert!(facts.clade_columns.is_empty());
    assert!(facts.node_heights.is_empty());
    assert!(facts.root_heights.is_empty());
    assert!(facts.split_ids.is_empty());
    assert!(facts.split_table.is_empty());
}

/// Opting into facts must not perturb interner IDs, clade materialization, or
/// the existing presence export. Existing constructors retain no sidecar.
#[test]
fn rooted_facts_are_optional_and_do_not_change_existing_snapshots() {
    let trees = ["((A:1,B:2):3,(C:4,D:5):6);", "((A:2,C:3):4,(B:5,D:6):7);"];
    let with_facts = rooted_snapshots_opts(&trees, true).unwrap();
    let without_facts = rooted_snapshots_opts(&trees, false).unwrap();
    let public_default = Snapshots::from_newicks(&trees, true).unwrap();

    assert!(with_facts.rooted_facts.is_some());
    assert!(without_facts.rooted_facts.is_none());
    assert!(public_default.rooted_facts.is_none());

    for (with, without) in with_facts.snapshots.iter().zip(&without_facts.snapshots) {
        assert_eq!(with.split_ids, without.split_ids);
    }
    assert_eq!(with_facts.clades.len(), without_facts.clades.len());
    for id in 0..with_facts.clades.len() {
        assert_eq!(with_facts.clades.get(id), without_facts.clades.get(id));
    }
    assert_eq!(
        with_facts.build_presence_matrix(),
        without_facts.build_presence_matrix()
    );
}

/// Strict binary/length validation belongs only to the optional facts path.
/// The established snapshot constructor keeps accepting the inputs it handled
/// before this feature existed.
#[test]
fn rooted_fact_validation_is_opt_in() {
    let missing_lengths = ["(A,B);"];
    assert!(rooted_snapshots_opts(&missing_lengths, false).is_ok());
    assert!(rooted_snapshots_opts(&missing_lengths, true).is_err());

    let polytomy = ["(A:1,B:1,C:1);"];
    assert!(rooted_snapshots_opts(&polytomy, false).is_ok());
    assert!(rooted_snapshots_opts(&polytomy, true).is_err());
}

#[test]
fn rooted_facts_errors_include_the_source_tree_index() {
    let trees = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,B:1,C:1):1,D:1);"];
    let error = rooted_snapshots_opts(&trees, true).expect_err("polytomy must fail");
    assert!(
        error.contains("tree at index 1") && error.contains("exactly two children"),
        "unexpected error: {error}"
    );
}

#[test]
fn rooted_facts_cannot_be_retained_in_unrooted_mode() {
    let empty: HashMap<String, String> = HashMap::new();
    let trees = ["((A:1,B:1):1,(C:1,D:1):1);"];
    let error = Snapshots::from_newick_iter_opts(
        trees.iter().map(|&newick| (newick, &empty)),
        false,
        Retain {
            lengths: false,
            bipartitions: true,
            rooted_facts: true,
        },
    )
    .expect_err("unrooted facts must fail");
    assert_eq!(error, "Rooted facts require rooted snapshot mode.");
}

#[test]
fn rooted_facts_are_included_in_raw_chunk_memory_estimates() {
    let newick = "((A:1,B:2):3,(C:4,D:5):6);";
    let snapshot = snapshot_of(newick, true);
    let (facts, _) = rooted_facts_of(newick).unwrap();
    let ordinary = estimated_raw_snapshot_bytes(&snapshot, None);
    let with_facts = estimated_raw_snapshot_bytes(&snapshot, Some(&facts));

    assert_eq!(
        with_facts - ordinary,
        facts.estimated_heap_bytes(),
        "chunk sizing must account for every retained raw fact element"
    );
    assert!(with_facts > ordinary);
}
