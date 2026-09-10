//! Flat byte buffers for the Python side: which split is in which tree, what
//! each edge is worth, and which taxa each split names.
//!
//! Every builder here shares one column order — ascending packed leaf set,
//! stable across calls on the same tree set — so the three buffers can be
//! indexed against each other and against `leaf_names`.

use super::Snapshots;
use super::clades::cmp_packed;
use crate::par::*;
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
