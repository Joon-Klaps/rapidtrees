//! Rectangular tree layout, for drawing trees and tanglegrams.
//!
//! Port of `treetracer/clade_freq/layout.py`, matching its conventions exactly
//! so the browser and desktop draw the same tree:
//!
//! * **x** — cumulative branch length from the root. The root sits at 0 and
//!   each tip at its total root-to-tip distance, so x is a time axis.
//! * **y** — tip rank from an in-order traversal (0, 1, 2, …); an internal node
//!   sits at the mean y of its children.
//! * **ladderize** — children sorted so the subtree with fewer tips comes
//!   first, giving an ascending staircase. Topology is unchanged; only the
//!   drawing order moves.
//!
//! Parsing reuses phylotree rather than adding a second newick parser to this
//! crate. Traversals are iterative: posterior samples of a few hundred taxa are
//! fine recursively, but a pathological ladder-shaped tree would put recursion
//! depth in the thousands, and this runs in a wasm stack.

use phylotree::tree::{NodeId, Tree as PhyloTree};

/// Node coordinates for one tree, flattened for the JS boundary.
pub struct TreeLayout {
    /// Cumulative branch length from the root, one per node.
    pub x: Vec<f64>,
    /// Vertical position, one per node.
    pub y: Vec<f64>,
    /// Index into these arrays of each node's parent; `-1` for the root.
    pub parent: Vec<i32>,
    /// Names of the tips, in drawing order (ascending y).
    pub tip_names: Vec<String>,
    /// Index into `x`/`y` for each entry of `tip_names`.
    pub tip_index: Vec<u32>,
}

impl TreeLayout {
    pub fn n_nodes(&self) -> usize {
        self.x.len()
    }
}

/// Lay out a newick string.
pub fn layout(newick: &str, ladderize: bool) -> Result<TreeLayout, String> {
    let tree = PhyloTree::from_newick(newick.trim()).map_err(|e| e.to_string())?;
    let root = tree.get_root().map_err(|e| e.to_string())?;

    // --- 1. tip counts, bottom-up -----------------------------------------
    // An explicit two-phase stack rather than recursion: push a node, then
    // push it again marked "children done" beneath its children.
    let mut tip_count: std::collections::HashMap<NodeId, usize> = std::collections::HashMap::new();
    let mut stack: Vec<(NodeId, bool)> = vec![(root, false)];
    while let Some((id, done)) = stack.pop() {
        let node = tree.get(&id).map_err(|e| e.to_string())?;
        if node.children.is_empty() {
            tip_count.insert(id, 1);
            continue;
        }
        if done {
            let total: usize = node
                .children
                .iter()
                .map(|c| tip_count.get(c).copied().unwrap_or(0))
                .sum();
            tip_count.insert(id, total);
        } else {
            stack.push((id, true));
            for &c in &node.children {
                stack.push((c, false));
            }
        }
    }

    // --- 2. child order ----------------------------------------------------
    let child_order = |id: NodeId| -> Vec<NodeId> {
        let mut kids = tree
            .get(&id)
            .map(|n| n.children.clone())
            .unwrap_or_default();
        if ladderize {
            kids.sort_by_key(|c| tip_count.get(c).copied().unwrap_or(0));
        }
        kids
    };

    // --- 3. preorder walk assigning x and slot indices ---------------------
    let mut x: Vec<f64> = Vec::new();
    let mut parent: Vec<i32> = Vec::new();
    let mut slot_of: std::collections::HashMap<NodeId, usize> = std::collections::HashMap::new();
    let mut order: Vec<NodeId> = Vec::new();

    // (node, parent slot, parent x)
    let mut walk: Vec<(NodeId, i32, f64)> = vec![(root, -1, 0.0)];
    while let Some((id, parent_slot, parent_x)) = walk.pop() {
        let node = tree.get(&id).map_err(|e| e.to_string())?;
        // The root's own parent_edge is meaningless; treat it as 0.
        let length = if parent_slot < 0 {
            0.0
        } else {
            node.parent_edge.unwrap_or(0.0)
        };
        let node_x = parent_x + length;

        let slot = x.len();
        x.push(node_x);
        parent.push(parent_slot);
        slot_of.insert(id, slot);
        order.push(id);

        // Reversed so the first child is popped first and keeps its rank.
        for &c in child_order(id).iter().rev() {
            walk.push((c, slot as i32, node_x));
        }
    }

    // --- 4. y: tips take successive ranks, internals the mean of children ---
    let mut y = vec![0.0f64; x.len()];
    let mut tip_names: Vec<String> = Vec::new();
    let mut tip_index: Vec<u32> = Vec::new();
    let mut next_rank = 0.0f64;

    let mut post: Vec<(NodeId, bool)> = vec![(root, false)];
    while let Some((id, done)) = post.pop() {
        let node = tree.get(&id).map_err(|e| e.to_string())?;
        let slot = slot_of[&id];

        if node.children.is_empty() {
            y[slot] = next_rank;
            next_rank += 1.0;
            tip_names.push(node.name.clone().unwrap_or_default());
            tip_index.push(slot as u32);
            continue;
        }

        if done {
            let kids = child_order(id);
            let sum: f64 = kids.iter().map(|c| y[slot_of[c]]).sum();
            y[slot] = sum / kids.len() as f64;
        } else {
            post.push((id, true));
            // Reversed so children are visited in drawing order, which is what
            // makes tip ranks ascend down the page.
            for &c in child_order(id).iter().rev() {
                post.push((c, false));
            }
        }
    }

    Ok(TreeLayout {
        x,
        y,
        parent,
        tip_names,
        tip_index,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn balanced_tree_ranks_tips_in_order() {
        let out = layout("((A:1,B:1):1,(C:1,D:1):1);", false).expect("layout");
        assert_eq!(out.tip_names.len(), 4);
        // Tip ranks are 0..3 in drawing order.
        let ys: Vec<f64> = out.tip_index.iter().map(|&i| out.y[i as usize]).collect();
        assert_eq!(ys, vec![0.0, 1.0, 2.0, 3.0]);
        // Root has no parent.
        assert_eq!(out.parent[0], -1);
        assert_eq!(out.x[0], 0.0);
    }

    #[test]
    fn x_is_cumulative_branch_length() {
        let out = layout("((A:1.5,B:2.5):3.0,C:0.5);", false).expect("layout");
        // Root at 0; C is a direct child at 0.5; A is 3.0 + 1.5 = 4.5.
        let pos = |name: &str| -> f64 {
            let k = out.tip_names.iter().position(|n| n == name).unwrap();
            out.x[out.tip_index[k] as usize]
        };
        assert!((pos("C") - 0.5).abs() < 1e-12, "C at {}", pos("C"));
        assert!((pos("A") - 4.5).abs() < 1e-12, "A at {}", pos("A"));
        assert!((pos("B") - 5.5).abs() < 1e-12, "B at {}", pos("B"));
    }

    #[test]
    fn internal_nodes_sit_at_the_mean_of_their_children() {
        let out = layout("((A:1,B:1):1,(C:1,D:1):1);", false).expect("layout");
        // The two cherries sit at 0.5 and 2.5; the root between them at 1.5.
        let mut internal: Vec<f64> = (0..out.n_nodes())
            .filter(|&i| !out.tip_index.contains(&(i as u32)))
            .map(|i| out.y[i])
            .collect();
        internal.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(internal, vec![0.5, 1.5, 2.5]);
    }

    #[test]
    fn ladderize_puts_the_smaller_clade_first() {
        // One tip against a clade of three.
        let plain = layout("(Z:1,((A:1,B:1):1,C:1):1);", false).expect("plain");
        let ladder = layout("(Z:1,((A:1,B:1):1,C:1):1);", true).expect("ladder");
        assert_eq!(plain.tip_names[0], "Z");
        // Ladderized, the single tip still sorts first — it is the smaller
        // subtree — but within the big clade, C (1 tip) precedes the cherry.
        let zi = ladder.tip_names.iter().position(|n| n == "Z").unwrap();
        let ci = ladder.tip_names.iter().position(|n| n == "C").unwrap();
        let ai = ladder.tip_names.iter().position(|n| n == "A").unwrap();
        assert!(zi < ci, "Z should come before C");
        assert!(ci < ai, "C (1 tip) should precede the cherry");
    }

    #[test]
    fn rejects_malformed_newick() {
        assert!(layout("((A:1,B:1;", false).is_err());
    }
}
