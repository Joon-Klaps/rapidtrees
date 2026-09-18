//! Optional per-tree facts needed by rooted summary-tree algorithms.
//!
//! This module deliberately sits beside the existing snapshot builder rather
//! than inside it.  A facts-enabled caller can run this collector while the
//! parsed [`PhyloTree`] is still alive, then hand the ordinary snapshot to
//! the unchanged interner.  Existing distance paths do neither the traversal
//! nor the allocations defined here.
//!
//! The collector retains two facts that an RF snapshot does not need:
//! exact-clade node heights and directly observed binary child splits.  Clades
//! are represented by the same run-wide fingerprints and subtree sizes used by
//! the interner, so a later step can resolve them to the IDs already assigned
//! by the ordinary rooted snapshot path.

use super::fingerprint::Fingerprint;
use phylotree::tree::{Node, Tree as PhyloTree};
use rustc_hash::FxHashMap;

/// A rooted clade on its way to the existing global clade interner.
///
/// Rooted mode keeps the raw subtree fingerprint.  `size` accompanies it so a
/// later lookup can construct exactly the same `(fingerprint, cardinality)`
/// candidate key as the existing interner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(super) struct RawCladeRef {
    pub(super) key: Fingerprint,
    pub(super) size: u32,
}

/// One non-root clade occurrence and its height in one source tree.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct RawNodeFact {
    pub(super) clade: RawCladeRef,
    pub(super) height: f64,
}

/// One directly observed binary resolution.
///
/// `parent == None` denotes the implicit all-taxa root, which is intentionally
/// absent from RapidTrees' rooted-clade table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct RawSplitFact {
    pub(super) parent: Option<RawCladeRef>,
    pub(super) children: [RawCladeRef; 2],
}

/// Facts extracted from one parsed, strictly binary rooted tree.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct RawRootedFacts {
    pub(super) nodes: Vec<RawNodeFact>,
    pub(super) root_height: f64,
    pub(super) splits: Vec<RawSplitFact>,
}

impl RawRootedFacts {
    /// Collect node heights and observed child splits without recursion.
    ///
    /// `labels` and `leaf_index` must be the run-wide tables also supplied to
    /// the ordinary rooted snapshot builder.  Reusing them is what makes the
    /// raw clade references resolvable against that snapshot's interner.
    ///
    /// Heights follow TreeTracer's existing convention:
    ///
    /// ```text
    /// root_height = max(root-to-tip distance)
    /// node_height = root_height - root-to-node distance
    /// ```
    ///
    /// Every non-root edge must have an explicit finite length, and every
    /// internal node must have exactly two children.  Finite negative lengths
    /// are accepted: this collector validates representation and arithmetic,
    /// not biological plausibility.
    pub(super) fn from_tree(
        tree: &PhyloTree,
        labels: &[Fingerprint],
        leaf_index: &FxHashMap<&str, usize>,
    ) -> Result<Self, String> {
        let root_id = tree
            .get_root()
            .map_err(|error| format!("failed to find tree root: {error}"))?;

        // Resolve every node once.  The resulting preorder guarantees that a
        // parent's distance has been calculated before any of its children.
        let mut order: Vec<(usize, &Node)> = Vec::with_capacity(tree.size());
        let mut stack = vec![root_id];
        while let Some(id) = stack.pop() {
            let node = tree
                .get(&id)
                .map_err(|error| format!("failed to access node {id}: {error}"))?;
            order.push((id, node));
            stack.extend(node.children.iter().copied());
        }

        let mut distances = vec![0.0; tree.size()];
        let mut tip_ids = Vec::with_capacity(leaf_index.len());

        for &(id, node) in &order {
            if node.children.is_empty() {
                tip_ids.push(id);
                continue;
            }
            if node.children.len() != 2 {
                return Err(format!(
                    "internal node {id} must have exactly two children; found {}",
                    node.children.len()
                ));
            }

            for &child_id in &node.children {
                let child = tree
                    .get(&child_id)
                    .map_err(|error| format!("failed to access child node {child_id}: {error}"))?;
                let length = child.parent_edge.ok_or_else(|| {
                    format!("non-root node {child_id} is missing an explicit branch length")
                })?;
                if !length.is_finite() {
                    return Err(format!(
                        "non-root node {child_id} has a non-finite branch length"
                    ));
                }
                let distance = distances[id] + length;
                if !distance.is_finite() {
                    return Err(format!(
                        "non-finite cumulative root distance at node {child_id}"
                    ));
                }
                distances[child_id] = distance;
            }
        }

        let root_height = tip_ids
            .iter()
            .map(|&id| distances[id])
            .reduce(f64::max)
            .ok_or_else(|| "tree has no tips".to_string())?;
        if !root_height.is_finite() {
            return Err("calculated a non-finite root height".to_string());
        }

        // The ordinary snapshot arena indexes its accumulator by node ID too;
        // parsed PhyloTrees use the same dense node-ID invariant here.
        let mut acc = vec![Acc::default(); tree.size()];
        for &(id, node) in order.iter().filter(|(_, node)| node.children.is_empty()) {
            let name = node
                .name
                .as_deref()
                .ok_or_else(|| format!("leaf node {id} is unnamed"))?;
            let &bit = leaf_index
                .get(name)
                .ok_or_else(|| format!("leaf {name:?} is absent from the shared taxon index"))?;
            let &fingerprint = labels
                .get(bit)
                .ok_or_else(|| format!("taxon index {bit} has no fingerprint label"))?;
            acc[id] = Acc {
                key: fingerprint,
                size: 1,
            };
        }

        for &(id, node) in order
            .iter()
            .rev()
            .filter(|(_, node)| !node.children.is_empty())
        {
            let mut folded = Acc::default();
            for &child_id in &node.children {
                folded.key ^= acc[child_id].key;
                folded.size = folded
                    .size
                    .checked_add(acc[child_id].size)
                    .ok_or_else(|| "tree contains more than u32::MAX tips".to_string())?;
            }
            acc[id] = folded;
        }

        let clade_ref = |id: usize| RawCladeRef {
            key: acc[id].key,
            size: acc[id].size,
        };

        let mut nodes = Vec::with_capacity(order.len().saturating_sub(1));
        let mut splits = Vec::with_capacity(tip_ids.len().saturating_sub(1));
        for &(id, node) in &order {
            let height = root_height - distances[id];
            if !height.is_finite() {
                return Err(format!("calculated a non-finite height at node {id}"));
            }

            if id != root_id {
                nodes.push(RawNodeFact {
                    clade: clade_ref(id),
                    height,
                });
            }

            if !node.children.is_empty() {
                let mut children = [clade_ref(node.children[0]), clade_ref(node.children[1])];
                children.sort_unstable();
                splits.push(RawSplitFact {
                    parent: (id != root_id).then(|| clade_ref(id)),
                    children,
                });
            }
        }

        let expected_nodes = tip_ids.len().saturating_mul(2).saturating_sub(2);
        let expected_splits = tip_ids.len().saturating_sub(1);
        if nodes.len() != expected_nodes || splits.len() != expected_splits {
            return Err(format!(
                "strictly binary tree with {} tips must contain {expected_nodes} non-root nodes and {expected_splits} internal splits; found {} and {}",
                tip_ids.len(),
                nodes.len(),
                splits.len()
            ));
        }

        Ok(Self {
            nodes,
            root_height,
            splits,
        })
    }
}

/// Per-node postorder accumulator for a rooted descendant set.
#[derive(Debug, Clone, Copy, Default)]
struct Acc {
    key: Fingerprint,
    size: u32,
}
