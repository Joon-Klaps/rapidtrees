//! Flat byte buffers for the Python side: which split is in which tree, what
//! each edge is worth, and which taxa each split names.
//!
//! Every builder here shares one column order — ascending packed leaf set,
//! stable across calls on the same tree set — so the three buffers can be
//! indexed against each other and against `leaf_names`.

use super::Snapshots;
use super::clades::cmp_packed;
use super::rooted_facts::ROOT_ID;
use crate::par::*;

/// Fixed-width, tree-major buffers for the optional rooted-tree facts.
///
/// Every non-root clade value is a column in the existing presence matrix.
/// `root_column` is the reserved value immediately after the last real clade
/// column and appears only as the parent of a root split.
///
/// This is consumed by the dedicated Python endpoint added in the next stage.
#[allow(dead_code)]
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct RootedFactBuffers {
    pub(crate) root_column: u32,
    pub(crate) nodes_per_tree: usize,
    pub(crate) splits_per_tree: usize,
    pub(crate) node_columns: Vec<u8>,
    pub(crate) node_heights: Vec<u8>,
    pub(crate) root_heights: Vec<u8>,
    pub(crate) split_columns: Vec<u8>,
}

impl Snapshots {
    /// The one column order every builder here shares, in both directions.
    ///
    /// Returns `(id_to_col, col_to_bip_id)`: `id_to_col[split_id]` is that
    /// split's column, and `col_to_bip_id[col]` is the split ID in that column.
    /// Columns ascend by packed leaf set, which is what makes the order stable across
    /// calls on the same tree set — and therefore safe to hand to Python as a
    /// column index.
    fn column_order(&self) -> (Vec<usize>, Vec<usize>) {
        let n = self.clades.len();
        let mut col_to_bip_id: Vec<usize> = (0..n).collect();
        col_to_bip_id.sort_unstable_by(|&a, &b| cmp_packed(self.clades.get(a), self.clades.get(b)));
        let mut id_to_col = vec![0usize; n];
        for (col, &orig_id) in col_to_bip_id.iter().enumerate() {
            id_to_col[orig_id] = col;
        }
        (id_to_col, col_to_bip_id)
    }

    /// Build a flat row-major presence matrix `(n_trees × n_bip)` as `Vec<u8>`,
    /// together with a column-order index into the clade table.
    ///
    /// Returns `(presence_bytes, col_to_bip_id)`:
    /// - `presence_bytes`: flat `uint8` buffer, row-major, shape `(n_trees, n_bip)`
    /// - `col_to_bip_id`: `col_to_bip_id[col]` is the split ID
    ///   for that column, in ascending packed-leaf-set order
    ///
    /// Columns are in ascending packed-leaf-set order, stable across calls on the
    /// same tree set.
    /// Each byte is `1` if the split is present in that tree, `0` otherwise.
    pub fn build_presence_matrix(&self) -> (Vec<u8>, Vec<usize>) {
        let (id_to_col, col_to_bip_id) = self.column_order();
        let n_trees = self.snapshots.len();
        let n_bip = self.clades.len();
        let mut presence = vec![0u8; n_trees * n_bip];

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

    /// Export the optional rooted-tree sidecar on existing presence columns.
    ///
    /// `col_to_bip_id` must be the complete column mapping returned by
    /// [`Snapshots::build_presence_matrix`]. It is inverted once here rather
    /// than deriving column order again, so every exported clade value indexes
    /// the accompanying presence matrix and clade table directly.
    ///
    /// The buffers are native-endian and tree-major:
    ///
    /// - `node_columns`: `u32[T, 2L - 2]`
    /// - `node_heights`: `f64[T, 2L - 2]`
    /// - `root_heights`: `f64[T]`
    /// - `split_columns`: `u32[T, L - 1, 3]`
    ///
    /// Split triples are `(parent, child_a, child_b)`. The children are sorted
    /// after ID-to-column conversion, and the implicit root is represented by
    /// `root_column == C`, where `C` is the number of exported clades.
    #[allow(dead_code)]
    pub(crate) fn build_rooted_fact_buffers(
        &self,
        col_to_bip_id: &[usize],
    ) -> Result<RootedFactBuffers, String> {
        let store = self
            .rooted_facts
            .as_ref()
            .ok_or_else(|| "rooted facts were not retained for these snapshots".to_string())?;
        let n_trees = self.snapshots.len();
        if store.len() != n_trees {
            return Err(format!(
                "rooted facts contain {} tree rows but snapshots contain {n_trees}",
                store.len()
            ));
        }

        let n_clades = self.n_distinct_splits();
        if col_to_bip_id.len() != n_clades {
            return Err(format!(
                "rooted-fact column mapping contains {} entries but snapshots contain {n_clades} distinct clades",
                col_to_bip_id.len()
            ));
        }
        let root_column = u32::try_from(n_clades).map_err(|_| {
            "rooted facts cannot export more than u32::MAX clade columns".to_string()
        })?;

        // u32::MAX is safe as an unfilled marker: a valid non-root column is
        // always strictly less than root_column, which itself is at most that
        // value.
        let mut id_to_col = vec![u32::MAX; n_clades];
        for (column, &id) in col_to_bip_id.iter().enumerate() {
            let slot = id_to_col.get_mut(id).ok_or_else(|| {
                format!(
                    "rooted-fact column mapping references clade ID {id}, but valid IDs are 0..{n_clades}"
                )
            })?;
            if *slot != u32::MAX {
                return Err(format!(
                    "rooted-fact column mapping contains clade ID {id} more than once"
                ));
            }
            *slot = u32::try_from(column).map_err(|_| {
                "rooted-fact column index cannot be represented as uint32".to_string()
            })?;
        }

        let n_leaves = self.leaf_names.len();
        let nodes_per_tree = if n_leaves == 0 {
            0
        } else {
            n_leaves
                .checked_mul(2)
                .and_then(|nodes| nodes.checked_sub(2))
                .ok_or_else(|| "rooted-fact node dimension overflows usize".to_string())?
        };
        let splits_per_tree = n_leaves.saturating_sub(1);
        let node_values = checked_product(n_trees, nodes_per_tree, "node")?;
        let split_values = checked_product(n_trees, splits_per_tree, "split")?;
        let split_scalars = checked_product(split_values, 3, "split scalar")?;

        let mut node_columns =
            vec![0; checked_bytes(node_values, size_of::<u32>(), "node column")?];
        let mut node_heights =
            vec![0; checked_bytes(node_values, size_of::<f64>(), "node height")?];
        let mut root_heights = vec![0; checked_bytes(n_trees, size_of::<f64>(), "root height")?];
        let mut split_columns =
            vec![0; checked_bytes(split_scalars, size_of::<u32>(), "split column")?];

        for (tree_index, facts) in store.trees.iter().enumerate() {
            if facts.node_ids.len() != nodes_per_tree || facts.node_heights.len() != nodes_per_tree
            {
                return Err(format!(
                    "rooted-fact row {tree_index} has {} node IDs and {} node heights; expected {nodes_per_tree} of each",
                    facts.node_ids.len(),
                    facts.node_heights.len()
                ));
            }
            if facts.splits.len() != splits_per_tree {
                return Err(format!(
                    "rooted-fact row {tree_index} has {} splits; expected {splits_per_tree}",
                    facts.splits.len()
                ));
            }
            if !facts.root_height.is_finite() {
                return Err(format!(
                    "rooted-fact row {tree_index} has a non-finite root height"
                ));
            }

            let node_base = tree_index * nodes_per_tree;
            for (offset, (&id, &height)) in
                facts.node_ids.iter().zip(&facts.node_heights).enumerate()
            {
                if !height.is_finite() {
                    return Err(format!(
                        "rooted-fact row {tree_index}, node {offset} has a non-finite height"
                    ));
                }
                let column = exported_column(id, &id_to_col, tree_index, "node")?;
                write_u32(&mut node_columns, node_base + offset, column);
                write_f64(&mut node_heights, node_base + offset, height);
            }
            write_f64(&mut root_heights, tree_index, facts.root_height);

            let expected_root_splits = usize::from(splits_per_tree > 0);
            let root_splits = facts
                .splits
                .iter()
                .filter(|split| split[0] == ROOT_ID)
                .count();
            if root_splits != expected_root_splits {
                return Err(format!(
                    "rooted-fact row {tree_index} has {root_splits} root splits; expected {expected_root_splits}"
                ));
            }

            let split_base = tree_index * splits_per_tree;
            for (offset, &[parent, child_a, child_b]) in facts.splits.iter().enumerate() {
                let parent_column = if parent == ROOT_ID {
                    root_column
                } else {
                    exported_column(parent, &id_to_col, tree_index, "split parent")?
                };
                let mut children = [
                    exported_column(child_a, &id_to_col, tree_index, "split child")?,
                    exported_column(child_b, &id_to_col, tree_index, "split child")?,
                ];
                children.sort_unstable();
                if children[0] == children[1] {
                    return Err(format!(
                        "rooted-fact row {tree_index}, split {offset} has identical child columns {}",
                        children[0]
                    ));
                }

                let triple = (split_base + offset) * 3;
                write_u32(&mut split_columns, triple, parent_column);
                write_u32(&mut split_columns, triple + 1, children[0]);
                write_u32(&mut split_columns, triple + 2, children[1]);
            }
        }

        Ok(RootedFactBuffers {
            root_column,
            nodes_per_tree,
            splits_per_tree,
            node_columns,
            node_heights,
            root_heights,
            split_columns,
        })
    }

    /// Build a flat row-major branch-length matrix `(n_trees × n_bip)` as native-endian
    /// `float64` bytes, together with a column-order index into `self.bipartitions`.
    ///
    /// Returns `(branch_length_bytes, col_to_bip_id)`:
    /// - `branch_length_bytes`: flat `float64` buffer, row-major, shape `(n_trees, n_bip)`
    /// - `col_to_bip_id`: `col_to_bip_id[col]` is the split ID
    ///   for that column, in ascending packed-leaf-set order
    ///
    /// `branch_length_bytes[i, j]` is the branch length of edge `j` in tree `i`, or `0.0`
    /// if that edge is absent from tree `i`. Pendant (leaf) edges are always present in
    /// every tree and therefore always have a non-zero value.
    ///
    /// Column order matches [`Snapshots::build_presence_matrix`], and is stable
    /// across calls on the same tree set.
    ///
    /// Decode on the Python side:
    /// ```python
    /// bl = np.frombuffer(branch_length_bytes, dtype=np.float64).reshape(n_trees, n_bip)
    /// # Fréchet ESS traces without materialising the n×n distance matrix:
    /// wrf_trace = np.sum(np.abs(bl[ref_idx, :] - bl), axis=1)
    /// kf_trace  = np.sqrt(np.sum((bl[ref_idx, :] - bl) ** 2, axis=1))
    /// ```
    pub fn build_branch_length_matrix(&self) -> (Vec<u8>, Vec<usize>) {
        let (id_to_col, col_to_bip_id) = self.column_order();
        let n_trees = self.snapshots.len();
        let n_bip = self.clades.len();

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
    /// Shape: `(n_bip, ceil(n_leaves / 8))` bytes, row-major, in the same column
    /// order as [`Snapshots::build_presence_matrix`].
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
        let bytes_per_bip = self.leaf_names.len().div_ceil(8);
        let mut out = vec![0u8; col_to_bip_id.len() * bytes_per_bip];
        // Bit packing happens here and nowhere else — the table itself holds
        // leaf indices, and this is the only place the wire format needs bits.
        for (col, &id) in col_to_bip_id.iter().enumerate() {
            let row = &mut out[col * bytes_per_bip..][..bytes_per_bip];
            for &leaf in self.clades.get(id) {
                row[leaf as usize / 8] |= 1u8 << (leaf % 8);
            }
        }
        out
    }
}

fn checked_product(left: usize, right: usize, label: &str) -> Result<usize, String> {
    left.checked_mul(right)
        .ok_or_else(|| format!("rooted-fact {label} count overflows usize"))
}

fn checked_bytes(values: usize, width: usize, label: &str) -> Result<usize, String> {
    values
        .checked_mul(width)
        .ok_or_else(|| format!("rooted-fact {label} buffer size overflows usize"))
}

fn exported_column(
    id: u32,
    id_to_col: &[u32],
    tree_index: usize,
    role: &str,
) -> Result<u32, String> {
    if id == ROOT_ID {
        return Err(format!(
            "rooted-fact row {tree_index} uses the root sentinel as a {role}"
        ));
    }
    id_to_col
        .get(id as usize)
        .copied()
        .filter(|&column| column != u32::MAX)
        .ok_or_else(|| {
            format!(
                "rooted-fact row {tree_index} references clade ID {id}, which is absent from the column mapping"
            )
        })
}

fn write_u32(bytes: &mut [u8], index: usize, value: u32) {
    let start = index * size_of::<u32>();
    bytes[start..start + size_of::<u32>()].copy_from_slice(&value.to_ne_bytes());
}

fn write_f64(bytes: &mut [u8], index: usize, value: f64) {
    let start = index * size_of::<f64>();
    bytes[start..start + size_of::<f64>()].copy_from_slice(&value.to_ne_bytes());
}
