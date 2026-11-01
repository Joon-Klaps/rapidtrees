//! WebAssembly bindings for `rapidtrees`, built with wasm-pack into the npm package.
//!
//! One function per metric, named after the Python API: [`pairwise_rf`], [`pairwise_wrf`] and [`pairwise_kf`]. Each takes the text of a BEAST/NEXUS or plain-Newick file, drops the first `burnin_trees` trees, and returns the tree names with the flat row-major `n * n` matrix.
//!
//! The text is read by `rapidtrees::io::read_trees`, the reader the CLI uses, so format detection, burn-in, tree names and the two-tree minimum are the library's, and this crate only converts arguments and errors. The whole file is passed in as one string, so the browser build holds it in memory twice, as the JS string and as its UTF-8 copy in wasm linear memory, plus one chunk of trees while they are parsed.

use rapidtrees::io::{ReadOptions, read_trees};
use rapidtrees::{Retain, Snapshots};
use wasm_bindgen::prelude::*;

/// Result of [`pairwise_rf`]: the trees kept after burn-in and their `uint32` Robinson-Foulds matrix.
#[wasm_bindgen]
pub struct RfMatrix {
    names: Vec<String>,
    matrix: Vec<u32>,
}

/// Result of [`pairwise_wrf`] and [`pairwise_kf`]: the trees kept after burn-in and their `float64` matrix.
#[wasm_bindgen]
pub struct FloatMatrix {
    names: Vec<String>,
    matrix: Vec<f64>,
}

#[wasm_bindgen]
impl RfMatrix {
    /// Number of trees after burn-in. The matrix has `n_trees * n_trees` cells.
    #[wasm_bindgen(getter)]
    pub fn n_trees(&self) -> usize {
        self.names.len()
    }

    /// One label per tree: the NEXUS tree name (`STATE_10000`) or `line<N>` for the Nth line of a Newick file.
    #[wasm_bindgen(getter)]
    pub fn names(&self) -> Vec<String> {
        self.names.clone()
    }

    /// The matrix as a `Uint32Array`, row-major. Each read copies it out of wasm memory, so read it once.
    #[wasm_bindgen(getter)]
    pub fn matrix(&self) -> Vec<u32> {
        self.matrix.clone()
    }
}

#[wasm_bindgen]
impl FloatMatrix {
    /// Number of trees after burn-in. The matrix has `n_trees * n_trees` cells.
    #[wasm_bindgen(getter)]
    pub fn n_trees(&self) -> usize {
        self.names.len()
    }

    /// One label per tree: the NEXUS tree name (`STATE_10000`) or `line<N>` for the Nth line of a Newick file.
    #[wasm_bindgen(getter)]
    pub fn names(&self) -> Vec<String> {
        self.names.clone()
    }

    /// The matrix as a `Float64Array`, row-major. Each read copies it out of wasm memory, so read it once.
    #[wasm_bindgen(getter)]
    pub fn matrix(&self) -> Vec<f64> {
        self.matrix.clone()
    }
}

/// Pairwise Robinson-Foulds distances between all trees in `content`.
///
/// `content` is a BEAST/NEXUS or plain-Newick file, as text. `burnin_trees` drops that many trees from the start (default 0). `rooted` compares clades instead of bipartitions (default false). Throws if the input holds no trees, fewer than two are left after burn-in, a tree fails to parse or the leaf sets differ.
#[wasm_bindgen]
pub fn pairwise_rf(
    content: &str,
    burnin_trees: Option<usize>,
    rooted: Option<bool>,
) -> Result<RfMatrix, JsError> {
    let (names, snaps) = load(content, burnin_trees, rooted, false)?;
    Ok(RfMatrix {
        names,
        matrix: snaps.into_pairwise_rf(None),
    })
}

/// Pairwise weighted Robinson-Foulds distances between all trees in `content`. Arguments as for [`pairwise_rf`]; every tree needs branch lengths.
#[wasm_bindgen]
pub fn pairwise_wrf(
    content: &str,
    burnin_trees: Option<usize>,
    rooted: Option<bool>,
) -> Result<FloatMatrix, JsError> {
    let (names, snaps) = load(content, burnin_trees, rooted, true)?;
    Ok(FloatMatrix {
        names,
        matrix: snaps.into_pairwise_wrf(None),
    })
}

/// Pairwise Kuhner-Felsenstein distances between all trees in `content`. Arguments as for [`pairwise_rf`]; every tree needs branch lengths.
#[wasm_bindgen]
pub fn pairwise_kf(
    content: &str,
    burnin_trees: Option<usize>,
    rooted: Option<bool>,
) -> Result<FloatMatrix, JsError> {
    let (names, snaps) = load(content, burnin_trees, rooted, true)?;
    Ok(FloatMatrix {
        names,
        matrix: snaps.into_pairwise_kf(None),
    })
}

/// The version of the `rapidtrees` crate this build was compiled from.
#[wasm_bindgen]
pub fn rapidtrees_version() -> String {
    env!("CARGO_PKG_VERSION").to_owned()
}

/// The trees of `content`, read with what a matrix of this metric needs.
fn load(
    content: &str,
    burnin_trees: Option<usize>,
    rooted: Option<bool>,
    weighted: bool,
) -> Result<(Vec<String>, Snapshots), JsError> {
    let mut opts = ReadOptions::default();
    opts.burnin_trees = burnin_trees.unwrap_or(0);
    opts.rooted = rooted.unwrap_or(false);
    opts.retain = Retain::for_distances(weighted);
    Ok(read_trees(content.as_bytes(), &opts)?)
}
