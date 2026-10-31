//! WebAssembly bindings for `rapidtrees`, built with wasm-pack into the npm package.
//!
//! One function per metric, named after the Python API: [`pairwise_rf`], [`pairwise_wrf`] and [`pairwise_kf`]. Each takes the text of a BEAST NEXUS file (with or without a `TRANSLATE` block) or of a Newick file with one tree per line, drops the first `burnin_trees` trees, and returns the tree names with the flat row-major `n * n` matrix.
//!
//! The whole file is passed in as one string, so the browser build needs the input in memory: the JS string, its UTF-8 copy in wasm linear memory, and the trees' text while they are parsed. Only the public API that `tests/wasm` already uses is called, and nothing here adds to `rapidtrees::io`.

use rapidtrees::io::{TreeFormat, detect_format, parse_taxon_block, strip_beast_annotations};
use rapidtrees::{Retain, Snapshots};
use std::collections::HashMap;
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
/// `content` is a BEAST/NEXUS file or a Newick file with one tree per line, as text. `burnin_trees` drops that many trees from the start (default 0). `rooted` compares clades instead of bipartitions (default false). Throws if the input holds no trees, the trees fail to parse or their leaf sets differ.
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

fn load(
    content: &str,
    burnin_trees: Option<usize>,
    rooted: Option<bool>,
    weighted: bool,
) -> Result<(Vec<String>, Snapshots), JsError> {
    let (translate, trees): (_, Box<dyn Iterator<Item = (String, &str)>>) =
        match detect_format(content) {
            TreeFormat::Nexus => (parse_taxon_block(content), Box::new(nexus_trees(content))),
            TreeFormat::Newick => (HashMap::new(), Box::new(newick_trees(content))),
        };
    let (names, newicks): (Vec<String>, Vec<String>) = trees
        .skip(burnin_trees.unwrap_or(0))
        .map(|(name, body)| (name, strip_beast_annotations(body)))
        .unzip();
    if names.is_empty() {
        return Err(JsError::new("no trees left to compare after burn-in"));
    }
    let snaps = Snapshots::from_newick_iter_opts(
        newicks.iter().map(|n| (n.as_str(), &translate)),
        rooted.unwrap_or(false),
        Retain::for_distances(weighted),
    )
    .map_err(|e| JsError::new(&e))?;
    Ok((names, snaps))
}

/// `(name, Newick)` of every `tree NAME = BODY` line in the NEXUS trees block. Splits on `" = "`, not `'='`: BEAST headers carry `[&lnP=…]`.
fn nexus_trees(content: &str) -> impl Iterator<Item = (String, &str)> {
    content
        .lines()
        .skip_while(|line| !line.trim_start().to_ascii_lowercase().starts_with("tree "))
        .take_while(|line| !line.trim_start().to_ascii_lowercase().starts_with("end;"))
        .filter_map(|line| {
            let (header, body) = line.split_once(" = ")?;
            let name = header.split_whitespace().nth(1)?;
            Some((name.trim_matches('\'').to_owned(), body.trim()))
        })
}

/// `(line<N>, Newick)` of every non-empty line of a plain Newick file, `N` counting from 1.
fn newick_trees(content: &str) -> impl Iterator<Item = (String, &str)> {
    content
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty())
        .map(|(i, line)| (format!("line{}", i + 1), line.trim()))
}
