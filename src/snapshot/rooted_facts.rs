//! Optional per-tree facts needed by rooted summary-tree algorithms.
//!
//! This module deliberately sits beside the existing snapshot builder rather
//! than inside it. A facts-enabled caller derives the sidecar from the direct
//! parser's short-lived [`Snapshot`], then hands that ordinary snapshot to the
//! unchanged interner. Existing distance paths do neither the reconstruction
//! nor the allocations defined here.
//!
//! The collector retains two facts that an RF snapshot does not need:
//! exact-clade node heights and directly observed binary child splits.  Clades
//! are represented by the same run-wide fingerprints and subtree sizes used by
//! the interner, so a later step can resolve them to the IDs already assigned
//! by the ordinary rooted snapshot path.

use super::build::{Part, Snapshot};
use super::fingerprint::Fingerprint;
use rustc_hash::FxHashMap;

/// Internal parent ID for the implicit all-taxa root.
///
/// Real clade IDs are assigned by the existing `u32` interner.  The exporter
/// added in the next stage will translate this private sentinel to the public
/// rooted-clade column sentinel instead of exposing it to Python.
pub(super) const ROOT_ID: u32 = u32::MAX;

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
    pub(super) height: f32,
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
    pub(super) root_height: f32,
    pub(super) splits: Vec<RawSplitFact>,
}

/// Rooted facts after every non-root clade has been resolved by the existing
/// global interner.
///
/// Heights are aligned with the owning [`super::InternSnap::split_ids`], so the
/// same clade IDs are not retained twice. `split_ids` indexes the run-wide
/// observed-split table on [`RootedFactsStore`].
#[derive(Debug, Clone, PartialEq)]
pub(super) struct InternedRootedFacts {
    pub(super) node_heights: Vec<f32>,
    pub(super) root_height: f32,
    pub(super) split_ids: Vec<u32>,
}

/// Optional tree-aligned sidecar owned by a completed snapshot collection.
#[derive(Debug, Default, PartialEq)]
pub(super) struct RootedFactsStore {
    pub(super) trees: Vec<InternedRootedFacts>,
    /// Unique `(parent, child_a, child_b)` triples on internal clade IDs.
    pub(super) split_table: Vec<[u32; 3]>,
}

impl RootedFactsStore {
    pub(super) fn len(&self) -> usize {
        self.trees.len()
    }
}

/// Construction-only split interner. The lookup table is discarded before
/// snapshots are returned, leaving only compact per-tree IDs and one split
/// table in the retained sidecar.
#[derive(Debug, Default)]
pub(super) struct RootedFactsBuilder {
    trees: Vec<InternedRootedFacts>,
    split_table: Vec<[u32; 3]>,
    split_lookup: FxHashMap<[u32; 3], u32>,
}

impl RootedFactsBuilder {
    pub(super) fn with_capacity(n_trees: usize) -> Self {
        Self {
            trees: Vec::with_capacity(n_trees),
            split_table: Vec::new(),
            split_lookup: FxHashMap::default(),
        }
    }

    pub(super) fn len(&self) -> usize {
        self.trees.len()
    }

    pub(super) fn intern_split(&mut self, split: [u32; 3]) -> Result<u32, String> {
        if let Some(&id) = self.split_lookup.get(&split) {
            return Ok(id);
        }
        let id = u32::try_from(self.split_table.len())
            .map_err(|_| "rooted facts contain more than u32::MAX observed splits".to_string())?;
        self.split_table.push(split);
        self.split_lookup.insert(split, id);
        Ok(id)
    }

    pub(super) fn push(&mut self, facts: InternedRootedFacts) {
        self.trees.push(facts);
    }

    pub(super) fn finish(self) -> RootedFactsStore {
        RootedFactsStore {
            trees: self.trees,
            split_table: self.split_table,
        }
    }
}

impl RawRootedFacts {
    /// Heap bytes retained while this raw fact set waits for interning.
    ///
    /// The chunk-size heuristic intentionally follows the existing snapshot
    /// estimate and counts element storage rather than small `Vec` headers.
    pub(super) fn estimated_heap_bytes(&self) -> usize {
        self.nodes.len() * size_of::<RawNodeFact>() + self.splits.len() * size_of::<RawSplitFact>()
    }

    /// Collect node heights and observed child splits without reparsing Newick.
    ///
    /// In rooted mode, snapshot parts are every non-root node in postorder and
    /// each subtree occupies a contiguous leaf interval. Those two invariants
    /// are enough to reconstruct the two children of every internal node and
    /// the implicit root while the short-lived snapshot is still available.
    ///
    /// Heights follow TreeTracer's existing convention:
    ///
    /// ```text
    /// root_height = max(root-to-tip distance)
    /// node_height = root_height - root-to-node distance
    /// ```
    ///
    /// The parser has already required an explicit finite length for every
    /// non-root edge on this opt-in path. Finite negative lengths remain
    /// accepted: this collector validates representation and arithmetic, not
    /// biological plausibility.
    pub(super) fn from_snapshot(snapshot: &Snapshot) -> Result<Self, String> {
        let tip_count = snapshot.leaf_order.len();
        if tip_count == 0 {
            return Err("tree has no tips".to_string());
        }

        let expected_nodes = tip_count.saturating_mul(2).saturating_sub(2);
        let expected_splits = tip_count.saturating_sub(1);

        if tip_count == 1 {
            return Ok(Self {
                nodes: Vec::new(),
                root_height: 0.0,
                splits: Vec::new(),
            });
        }

        let mut links = Vec::with_capacity(snapshot.parts.len());
        let mut stack: Vec<usize> = Vec::with_capacity(tip_count);
        for (index, part) in snapshot.parts.iter().enumerate() {
            validate_interval(part, tip_count, index)?;
            let children = if part.size == 1 {
                None
            } else {
                let right = stack.pop().ok_or_else(|| {
                    format!("internal snapshot node {index} is missing its right child")
                })?;
                let left = stack.pop().ok_or_else(|| {
                    format!("internal snapshot node {index} is missing its left child")
                })?;
                validate_children(part, &snapshot.parts[left], &snapshot.parts[right], index)?;
                Some([left, right])
            };
            links.push(children);
            stack.push(index);
        }

        if stack.len() != 2 {
            return Err(format!(
                "root must have exactly two children; found {}",
                stack.len()
            ));
        }
        let root_children = [stack[0], stack[1]];
        validate_root_children(snapshot, root_children)?;

        let mut distances = vec![0.0; snapshot.parts.len()];
        let mut traversal = vec![(root_children[1], 0.0), (root_children[0], 0.0)];
        let mut root_height: Option<f64> = None;
        while let Some((index, parent_distance)) = traversal.pop() {
            let part = &snapshot.parts[index];
            if !part.length.is_finite() {
                return Err(format!(
                    "non-root snapshot node {index} has a non-finite branch length"
                ));
            }
            let distance = parent_distance + part.length;
            if !distance.is_finite() {
                return Err(format!(
                    "non-finite cumulative root distance at snapshot node {index}"
                ));
            }
            distances[index] = distance;
            match links[index] {
                Some([left, right]) => {
                    traversal.push((right, distance));
                    traversal.push((left, distance));
                }
                None => {
                    root_height = Some(root_height.map_or(distance, |height| height.max(distance)));
                }
            }
        }

        let root_height_f64 = root_height.ok_or_else(|| "tree has no tips".to_string())?;
        if !root_height_f64.is_finite() {
            return Err("calculated a non-finite root height".to_string());
        }

        let clade_ref = |index: usize| raw_clade(&snapshot.parts[index]);
        let mut nodes = Vec::with_capacity(expected_nodes);
        let mut splits = Vec::with_capacity(expected_splits);
        for (index, part) in snapshot.parts.iter().enumerate() {
            let height_f64 = root_height_f64 - distances[index];
            if !height_f64.is_finite() {
                return Err(format!(
                    "calculated a non-finite height at snapshot node {index}"
                ));
            }
            let height = quantize_height(height_f64).ok_or_else(|| {
                format!(
                    "calculated height at snapshot node {index} {height_f64} is outside the finite float32 range"
                )
            })?;
            nodes.push(RawNodeFact {
                clade: raw_clade(part),
                height,
            });

            if let Some([left, right]) = links[index] {
                let mut children = [clade_ref(left), clade_ref(right)];
                children.sort_unstable();
                splits.push(RawSplitFact {
                    parent: Some(raw_clade(part)),
                    children,
                });
            }
        }

        let mut root_split_children = [clade_ref(root_children[0]), clade_ref(root_children[1])];
        root_split_children.sort_unstable();
        splits.push(RawSplitFact {
            parent: None,
            children: root_split_children,
        });

        let root_height = quantize_height(root_height_f64).ok_or_else(|| {
            format!("calculated root height {root_height_f64} is outside the finite float32 range")
        })?;
        if nodes.len() != expected_nodes || splits.len() != expected_splits {
            return Err(format!(
                "strictly binary tree with {} tips must contain {expected_nodes} non-root nodes and {expected_splits} internal splits; found {} and {}",
                tip_count,
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

fn raw_clade(part: &Part) -> RawCladeRef {
    RawCladeRef {
        key: part.key,
        size: part.size,
    }
}

fn validate_interval(part: &Part, tip_count: usize, index: usize) -> Result<(), String> {
    let start = part.first as usize;
    let end = start
        .checked_add(part.size as usize)
        .ok_or_else(|| format!("snapshot node {index} has an overflowing leaf interval"))?;
    if part.size == 0 || end > tip_count {
        return Err(format!(
            "snapshot node {index} has invalid leaf interval {start}..{end} for {tip_count} tips"
        ));
    }
    Ok(())
}

fn validate_children(parent: &Part, left: &Part, right: &Part, index: usize) -> Result<(), String> {
    let left_end = left.first.checked_add(left.size);
    let right_end = right.first.checked_add(right.size);
    let parent_end = parent.first.checked_add(parent.size);
    if left.first != parent.first
        || left_end != Some(right.first)
        || right_end != parent_end
        || left.key ^ right.key != parent.key
    {
        return Err(format!(
            "internal snapshot node {index} must have exactly two children spanning its interval"
        ));
    }
    Ok(())
}

fn validate_root_children(snapshot: &Snapshot, children: [usize; 2]) -> Result<(), String> {
    let left = &snapshot.parts[children[0]];
    let right = &snapshot.parts[children[1]];
    let tip_count = u32::try_from(snapshot.leaf_order.len())
        .map_err(|_| "tree contains more than u32::MAX tips".to_string())?;
    if left.first != 0
        || left.first.checked_add(left.size) != Some(right.first)
        || right.first.checked_add(right.size) != Some(tip_count)
    {
        return Err("root children do not cover the complete contiguous leaf order".to_string());
    }
    Ok(())
}

fn quantize_height(value: f64) -> Option<f32> {
    debug_assert!(value.is_finite());
    let quantized = value as f32;
    quantized.is_finite().then_some(quantized)
}
