//! WebAssembly bindings for in-browser pairwise tree-distance computation.
//!
//! Kept thin, in the same spirit as `api.rs` (the PyO3 layer): the binding does
//! validation and return-value assembly only; all real work lives in
//! [`crate::io::parse_beast_content`], [`crate::snapshot::Snapshots`] and
//! [`crate::distances`].
//!
//! [`compute_distances_core`] is plain Rust (no `wasm_bindgen`) so it can be
//! unit-tested natively. The `#[wasm_bindgen]` wrapper is gated behind the
//! `wasm` feature.

use crate::io;
use crate::snapshot::Snapshots;

/// Output of [`compute_distances_core`]: `(tree_names, matrix_bytes_le, dtype, n)`.
type ComputedMatrix = (Vec<String>, Vec<u8>, &'static str, usize);

/// Parse a BEAST `.trees` string and compute the requested pairwise distance
/// matrix, returning it as a flat little-endian byte buffer.
///
/// `method` is one of `"rf"`, `"wrf"`, `"kf"`. The matrix is row-major and
/// square (`n × n`). RF values are emitted as `u32` (4 bytes each, dtype
/// `"u32"`); WRF/KF values as `f64` (8 bytes each, dtype `"f64"`). Bytes are
/// always little-endian so the JS side can decode with a fixed `DataView`
/// byte order.
///
/// Returns `(tree_names, matrix_bytes, dtype, n)`.
///
/// # Errors
/// Returns `Err(String)` if the file contains no trees, the trees fail to parse
/// or have inconsistent leaf sets, or `method` is unrecognised.
pub(crate) fn compute_distances_core(
    content: &str,
    method: &str,
    rooted: bool,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> Result<ComputedMatrix, String> {
    let (translate_map, tree_pairs) = io::parse_beast_content(
        content,
        "upload",
        burnin_trees,
        burnin_states,
        use_real_taxa,
    );
    if tree_pairs.is_empty() {
        return Err("No trees found in the uploaded file.".to_string());
    }

    let (names, newicks): (Vec<String>, Vec<String>) = tree_pairs.into_iter().unzip();
    let entries = newicks.iter().map(|n| (n.as_str(), &translate_map));
    let snaps = Snapshots::from_newick_iter(entries, rooted)?;
    let n = names.len();

    let (bytes, dtype) = match method {
        "rf" => {
            let matrix = snaps.pairwise_rf();
            let mut bytes = Vec::with_capacity(matrix.len() * 4);
            for v in matrix {
                bytes.extend_from_slice(&(v as u32).to_le_bytes());
            }
            (bytes, "u32")
        }
        "wrf" => (f64_matrix_bytes(snaps.pairwise_wrf()), "f64"),
        "kf" => (f64_matrix_bytes(snaps.pairwise_kf()), "f64"),
        other => {
            return Err(format!(
                "Unknown method {other:?} (expected \"rf\", \"wrf\", or \"kf\")"
            ));
        }
    };

    Ok((names, bytes, dtype, n))
}

/// Serialise an `f64` distance matrix as little-endian bytes.
fn f64_matrix_bytes(matrix: Vec<f64>) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(matrix.len() * 8);
    for v in matrix {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    bytes
}

// The `#[wasm_bindgen]` items below are collected by wasm-bindgen from this
// compiled module directly; no crate-root re-export is needed.
#[cfg(feature = "wasm")]
mod bindings {
    use super::compute_distances_core;
    use wasm_bindgen::prelude::*;

    /// Install a readable panic hook so Rust panics surface in the browser
    /// console instead of an opaque `unreachable` trap. Marked as the module
    /// start function, so it runs automatically on instantiation.
    #[wasm_bindgen(start)]
    pub fn start() {
        console_error_panic_hook::set_once();
    }

    /// Result of [`compute_distances`]: tree names plus the matrix as a flat
    /// little-endian byte buffer.
    #[wasm_bindgen]
    pub struct DistanceResult {
        names: Vec<String>,
        matrix: Vec<u8>,
        dtype: String,
        n: usize,
    }

    #[wasm_bindgen]
    impl DistanceResult {
        /// Tree names, in matrix row/column order.
        #[wasm_bindgen(getter)]
        pub fn names(&self) -> Vec<String> {
            self.names.clone()
        }

        /// Element type of the matrix bytes: `"u32"` (RF) or `"f64"` (WRF/KF).
        #[wasm_bindgen(getter)]
        pub fn dtype(&self) -> String {
            self.dtype.clone()
        }

        /// Matrix dimension (number of trees); the matrix is `n × n`.
        #[wasm_bindgen(getter)]
        pub fn n(&self) -> usize {
            self.n
        }

        /// Move the matrix bytes out to JS as a `Uint8Array`, leaving this
        /// result empty. Moving (rather than cloning) avoids a second multi-GB
        /// copy inside wasm linear memory. Call exactly once.
        #[wasm_bindgen]
        pub fn take_matrix(&mut self) -> Vec<u8> {
            std::mem::take(&mut self.matrix)
        }
    }

    /// Compute the pairwise distance matrix for an uploaded BEAST `.trees` file.
    ///
    /// `method` is `"rf"`, `"wrf"`, or `"kf"`. See
    /// [`compute_distances_core`](super::compute_distances_core) for semantics.
    ///
    /// # Errors
    /// Surfaces parse/validation errors as a `JsError` (thrown in JS).
    #[wasm_bindgen]
    pub fn compute_distances(
        content: &str,
        method: &str,
        rooted: bool,
        burnin_trees: usize,
        burnin_states: usize,
        use_real_taxa: bool,
    ) -> Result<DistanceResult, JsError> {
        let (names, matrix, dtype, n) = compute_distances_core(
            content,
            method,
            rooted,
            burnin_trees,
            burnin_states,
            use_real_taxa,
        )
        .map_err(|e| JsError::new(&e))?;

        Ok(DistanceResult {
            names,
            matrix,
            dtype: dtype.to_string(),
            n,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::compute_distances_core;

    /// A tiny BEAST block: translate map + three trees with known RF distances.
    const BEAST: &str = "#NEXUS\n\
        Begin trees;\n\
        \tTranslate\n\
        \t\t1 A,\n\
        \t\t2 B,\n\
        \t\t3 C,\n\
        \t\t4 D\n\
        \t;\n\
        tree STATE_0 = ((1:1,2:1):1,(3:1,4:1):1);\n\
        tree STATE_1 = ((1:1,3:1):1,(2:1,4:1):1);\n\
        tree STATE_2 = ((1:1,2:1):1,(3:1,4:1):1);\n\
        End;\n";

    #[test]
    fn rf_matrix_has_expected_shape_and_values() {
        let (names, bytes, dtype, n) =
            compute_distances_core(BEAST, "rf", false, 0, 0, true).unwrap();
        assert_eq!(n, 3);
        assert_eq!(names.len(), 3);
        assert_eq!(dtype, "u32");
        assert_eq!(bytes.len(), 3 * 3 * 4);

        let mat: Vec<u32> = bytes
            .chunks_exact(4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect();
        let at = |i: usize, j: usize| mat[i * n + j];
        // Trees 0 and 2 are identical (RF 0); both differ from tree 1 by the
        // single non-trivial bipartition on each side → RF 2.
        assert_eq!(at(0, 0), 0);
        assert_eq!(at(0, 2), 0);
        assert_eq!(at(0, 1), 2);
        assert_eq!(at(1, 2), 2);
        // Symmetry.
        assert_eq!(at(1, 0), at(0, 1));
    }

    #[test]
    fn wrf_is_f64_and_zero_between_identical_trees() {
        let (_, bytes, dtype, n) = compute_distances_core(BEAST, "wrf", false, 0, 0, true).unwrap();
        assert_eq!(dtype, "f64");
        assert_eq!(bytes.len(), n * n * 8);
        let mat: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        let at = |i: usize, j: usize| mat[i * n + j];
        assert_eq!(at(0, 2), 0.0);
    }

    #[test]
    fn unknown_method_errors() {
        assert!(compute_distances_core(BEAST, "bogus", false, 0, 0, true).is_err());
    }

    #[test]
    fn empty_input_errors() {
        assert!(compute_distances_core("#NEXUS\n", "rf", false, 0, 0, true).is_err());
    }
}
