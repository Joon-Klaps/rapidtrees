//! Browser bindings.
//!
//! Built from the workspace root with:
//!
//! ```text
//! wasm-pack build treetracer-core --target web
//! ```
//!
//! No feature flags are needed any more. This crate depends on `rapidtrees`
//! with `default-features = false`, which drops the CLI (and with it `flate2`)
//! and turns off `parallel`, so rayon is not compiled in at all — see `par.rs`
//! in `rapidtrees` for why that is currently necessary. Pushing that choice
//! into the manifest means it cannot be forgotten at the command line.
//!
//! This is the only crate in the graph that links wasm-bindgen. Keeping it that
//! way avoids the schema-version lockstep that multiple bindgen crates require,
//! and lets `rapidtrees` stay a plain Rust library.
//!
//! # Why runs are pooled
//!
//! Comparing MCMC runs is the point of this tool, and that constrains the API.
//! Every tree from every run goes into **one** distance matrix and **one** MDS
//! projection. Projecting each run separately and overlaying the results would
//! be meaningless: MDS fixes coordinates only up to rotation, reflection and
//! axis order, so two independent embeddings are not in a comparable frame.
//!
//! Hence the builder. Add each run, then `finish()` once. What distinguishes
//! runs afterwards is [`TreeSet::run_ids`] — an index per tree — which the UI
//! uses to colour or facet a projection that is itself shared.
//!
//! The API is otherwise deliberately stateful: a `TreeSet` holds the parsed
//! `Snapshots` inside wasm memory and JavaScript gets a handle plus typed-array
//! views. Returning every newick across the boundary would copy the whole
//! dataset into JS and defeat the interning this crate exists to do.

use std::cell::RefCell;
use std::collections::HashMap;

use wasm_bindgen::prelude::*;

use crate::clades;
use crate::ess;
use crate::hipstr as hipstr_mod;
use crate::layout;
use crate::mds;
use rapidtrees::Snapshots;
use rapidtrees::io::load_beast_subsampled_str;

/// Default share of each chain discarded as burnin, in percent.
///
/// Exposed as a function because `wasm_bindgen` cannot export constants; the
/// point is that the UI reads its defaults from the core rather than keeping a
/// second copy that can drift.
#[wasm_bindgen(js_name = defaultBurninPercent)]
pub fn default_burnin_percent() -> f64 {
    1.0
}

/// Default number of trees kept per run after burnin.
///
/// The pairwise RF matrix is O(n²) in the *pooled* tree count, so this is the
/// setting that decides whether an analysis takes a second or a minute. 1000
/// per run keeps four runs at 4000 trees — an 8M-pair matrix, about 64 MB.
#[wasm_bindgen(js_name = defaultTargetTrees)]
pub fn default_target_trees() -> usize {
    1000
}

/// Route Rust panics to `console.error` rather than an opaque trap.
///
/// Without this a panic in release wasm surfaces as a bare `RuntimeError` with
/// no message, which is close to undebuggable from the JS side.
#[wasm_bindgen(js_name = setPanicHook)]
pub fn set_panic_hook() {
    std::panic::set_hook(Box::new(|info| {
        web_error(&info.to_string());
    }));
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console, js_name = error)]
    fn web_error(msg: &str);
}

/// Accumulates runs before they are parsed into a single comparable set.
#[wasm_bindgen(js_name = TreeSetBuilder)]
#[derive(Default)]
pub struct TreeSetBuilder {
    /// `(tree_name, newick, run_index)`, in load order.
    trees: Vec<(String, String, usize)>,
    /// One translate map per run. BEAST numbers taxa per file, and two runs can
    /// legitimately use different numbering for the same taxa, so these cannot
    /// be merged into one map.
    translates: Vec<HashMap<String, String>>,
    run_labels: Vec<String>,
    /// Which translate map each run uses. Runs no longer map one-to-one onto
    /// files, so the index has to be tracked explicitly.
    translate_of_run: Vec<usize>,
    /// Log density per tree, read from the header annotations before they were
    /// stripped. `NaN` for files that record none.
    log_density: Vec<f64>,
    /// MCMC state number per tree, so a trace can be plotted against chain
    /// position rather than sample index.
    states: Vec<u64>,
    /// Which annotation key supplied the density — `lnP`, `posterior`, … —
    /// for labelling the axis honestly.
    density_field: Option<String>,
}

#[wasm_bindgen(js_class = TreeSetBuilder)]
impl TreeSetBuilder {
    #[wasm_bindgen(constructor)]
    pub fn new() -> TreeSetBuilder {
        TreeSetBuilder::default()
    }

    /// Add one `.trees` file.
    ///
    /// `label` names the run and prefixes its tree names.
    ///
    /// Burnin and subsampling are **per run**, deliberately. Chains differ:
    /// one may have converged in a tenth of the states another needed, and
    /// runs of very different lengths should still contribute comparable
    /// numbers of trees to a pooled projection. A single global setting cannot
    /// express either.
    ///
    /// * `burnin_percent` — share of the chain to discard from the front, as a
    ///   percentage (1.0 = 1%). Proportional rather than absolute so the same
    ///   setting is meaningful across chains of different lengths.
    /// * `target_trees` — how many trees to keep after burnin, sampled evenly
    ///   across the remaining chain. `0` keeps everything. When the chain is
    ///   shorter than the target, everything after burnin is kept.
    ///
    /// Even spacing matters: taking the first N post-burnin trees would sample
    /// one narrow window of the chain and misrepresent how it mixed.
    /// `content` is the raw file bytes, not a JS string, and that is
    /// deliberate. A `.trees` file is ASCII, but a JS string stores it as
    /// UTF-16 — so `file.text()` on a 241 MB file costs 482 MB in the JS heap
    /// before wasm-bindgen transcodes it back to UTF-8 on the way in. Passing
    /// an `ArrayBuffer` skips both the inflation and the transcode.
    #[wasm_bindgen(js_name = addRun)]
    pub fn add_run(
        &mut self,
        content: &[u8],
        label: &str,
        burnin_percent: f64,
        target_trees: usize,
    ) -> Result<RunSummary, JsError> {
        let content = std::str::from_utf8(content)
            .map_err(|e| JsError::new(&format!("'{label}' is not valid UTF-8 text: {e}")))?;

        // Selection happens inside the loader, on borrowed slices, so the
        // trees being discarded are never allocated in the first place.
        let loaded = load_beast_subsampled_str(content, label, burnin_percent, target_trees, true);

        if loaded.newicks.is_empty() {
            return Err(JsError::new(&format!(
                "No trees found in '{label}'. Expected a NEXUS file containing \
                 `tree ... = ...;` statements."
            )));
        }

        // A file holding several chains becomes several runs. Eight chains
        // reported as one would hide exactly the between-chain disagreement
        // this tool exists to surface.
        let translate_slot = self.translates.len();
        let multi = loaded.chain_labels.len() > 1;
        let mut run_of_chain: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();

        let kept = loaded.newicks.len();
        for (i, ((name, newick), (density, state))) in loaded
            .names
            .into_iter()
            .zip(loaded.newicks)
            .zip(loaded.log_density.into_iter().zip(loaded.states))
            .enumerate()
        {
            let chain = loaded.chains.get(i).cloned().unwrap_or_default();
            let run = *run_of_chain.entry(chain.clone()).or_insert_with(|| {
                let label = if multi && !chain.is_empty() {
                    format!("{label}:{chain}")
                } else {
                    label.to_string()
                };
                self.run_labels.push(label);
                self.translate_of_run.push(translate_slot);
                self.run_labels.len() - 1
            });
            self.trees.push((name, newick, run));
            self.log_density.push(density);
            self.states.push(state);
        }
        if self.density_field.is_none() {
            self.density_field = loaded.density_field;
        }

        self.translates.push(loaded.translate);
        Ok(RunSummary {
            total: loaded.total,
            burnin: loaded.burnin,
            kept,
        })
    }

    #[wasm_bindgen(getter, js_name = nRuns)]
    pub fn n_runs(&self) -> usize {
        self.run_labels.len()
    }

    #[wasm_bindgen(getter, js_name = nTrees)]
    pub fn n_trees(&self) -> usize {
        self.trees.len()
    }

    /// Parse everything added so far into one comparable `TreeSet`.
    ///
    /// All runs must share a taxon set. That is a real constraint rather than
    /// an implementation limit: an RF distance between trees over different
    /// taxa is not defined, so runs with mismatched taxa cannot be placed in a
    /// common tree space at all.
    #[wasm_bindgen(js_name = finish)]
    pub fn finish(&self, rooted: bool) -> Result<TreeSet, JsError> {
        if self.trees.is_empty() {
            return Err(JsError::new("No runs added."));
        }

        let entries = self.trees.iter().map(|(_, newick, run)| {
            (
                newick.as_str(),
                &self.translates[self.translate_of_run[*run]],
            )
        });

        let snaps = Snapshots::from_newick_iter(entries, rooted).map_err(|e| {
            JsError::new(&format!(
                "{e}\n\nAll runs must describe the same taxa to be compared in a \
                 shared tree space."
            ))
        })?;

        Ok(TreeSet {
            log_density: self.log_density.clone(),
            states: self.states.clone(),
            density_field: self.density_field.clone(),
            names: self.trees.iter().map(|(n, _, _)| n.clone()).collect(),
            // Retained so `mccNewick` can hand the winning topology to a tree
            // viewer without going back to the file. Costs roughly the size of
            // the kept trees — a few tens of MB at the default 1000 per run,
            // against re-reading and re-parsing a 250 MB file.
            newicks: self.trees.iter().map(|(_, w, _)| w.clone()).collect(),
            run_ids: self.trees.iter().map(|(_, _, r)| *r as u32).collect(),
            run_labels: self.run_labels.clone(),
            translates: self.translates.clone(),
            translate_of_run: self.translate_of_run.clone(),
            snaps,
            rf_cache: RefCell::new(None),
            presence_cache: RefCell::new(None),
        })
    }
}

/// What `addRun` actually did with one file — so the UI can report
/// "1000 of 4001 trees (40 discarded as burnin)" rather than guessing.
#[wasm_bindgen]
pub struct RunSummary {
    total: usize,
    burnin: usize,
    kept: usize,
}

#[wasm_bindgen]
impl RunSummary {
    /// Trees in the file before any filtering.
    #[wasm_bindgen(getter)]
    pub fn total(&self) -> usize {
        self.total
    }
    /// Trees discarded from the front as burnin.
    #[wasm_bindgen(getter)]
    pub fn burnin(&self) -> usize {
        self.burnin
    }
    /// Trees actually contributed to the pooled set.
    #[wasm_bindgen(getter)]
    pub fn kept(&self) -> usize {
        self.kept
    }
}

/// A parsed, pooled set of trees held in wasm memory.
#[wasm_bindgen]
pub struct TreeSet {
    names: Vec<String>,
    newicks: Vec<String>,
    log_density: Vec<f64>,
    states: Vec<u64>,
    density_field: Option<String>,
    run_ids: Vec<u32>,
    run_labels: Vec<String>,
    /// Per-file translate maps, kept so a newick can be handed out with real
    /// taxon names rather than the BEAST numbering.
    translates: Vec<HashMap<String, String>>,
    /// Which translate map each run uses.
    translate_of_run: Vec<usize>,
    snaps: Snapshots,
    /// The RF matrix is O(n²) to build and wanted by both the matrix accessor
    /// and PCoA, so it is computed at most once per `TreeSet`.
    rf_cache: RefCell<Option<Vec<u32>>>,
    /// `(presence, n_bipartitions)`. Shared by clade frequencies and the MCC
    /// search, so likewise built at most once.
    presence_cache: RefCell<Option<(Vec<u8>, usize)>>,
}

#[wasm_bindgen]
impl TreeSet {
    #[wasm_bindgen(getter, js_name = nTrees)]
    pub fn n_trees(&self) -> usize {
        self.snaps.len()
    }

    #[wasm_bindgen(getter, js_name = nTaxa)]
    pub fn n_taxa(&self) -> usize {
        self.snaps.leaf_names.len()
    }

    /// Tree names, in matrix row order.
    #[wasm_bindgen(js_name = treeNames)]
    pub fn tree_names(&self) -> Vec<String> {
        self.names.clone()
    }

    /// Taxon names, alphabetically sorted (matching bipartition bit order).
    #[wasm_bindgen(js_name = taxa)]
    pub fn taxa(&self) -> Vec<String> {
        self.snaps.leaf_names.clone()
    }

    /// Log density per tree, in matrix row order.
    ///
    /// Read from each tree's `[&lnP=…]` header before annotations are
    /// stripped. `NaN` where the file records none — the caller should treat
    /// an all-NaN result as "this file has no likelihoods" rather than plot a
    /// flat line.
    #[wasm_bindgen(js_name = logDensity)]
    pub fn log_density(&self) -> Vec<f64> {
        self.log_density.clone()
    }

    /// MCMC state number per tree.
    #[wasm_bindgen(js_name = states)]
    pub fn states(&self) -> Vec<f64> {
        // f64 rather than u64: state numbers reach into the billions on long
        // chains, and BigUint64Array is awkward on the JS side for values that
        // are only ever plotted.
        self.states.iter().map(|&s| s as f64).collect()
    }

    /// Which annotation key the density came from, for axis labelling.
    #[wasm_bindgen(getter, js_name = densityField)]
    pub fn density_field(&self) -> Option<String> {
        self.density_field.clone()
    }

    /// Which run each tree came from — index into [`TreeSet::run_labels`].
    ///
    /// This is what makes a pooled projection comparable: one embedding, with
    /// runs distinguished by colour rather than by separate axes.
    #[wasm_bindgen(js_name = runIds)]
    pub fn run_ids(&self) -> Vec<u32> {
        self.run_ids.clone()
    }

    #[wasm_bindgen(js_name = runLabels)]
    pub fn run_labels(&self) -> Vec<String> {
        self.run_labels.clone()
    }

    /// Dense `n × n` Robinson–Foulds distance matrix, row-major.
    ///
    /// Spans every run, so block `(a, b)` of the matrix is the between-run
    /// comparison of runs `a` and `b`, and the diagonal blocks are within-run.
    #[wasm_bindgen(js_name = rfMatrix)]
    pub fn rf_matrix(&self) -> Vec<u32> {
        self.ensure_rf()
    }

    /// PCoA embedding of the pooled RF matrix.
    #[wasm_bindgen(js_name = pcoa)]
    pub fn pcoa(&self, k: usize) -> Result<PcoaResult, JsError> {
        let rf = self.ensure_rf();
        let out = mds::pcoa(&rf, self.snaps.len(), k)
            .ok_or_else(|| JsError::new("PCoA needs at least two trees and one axis."))?;
        Ok(PcoaResult::from(out))
    }

    /// RF distance from every tree to tree `reference` — one row of the
    /// matrix.
    ///
    /// This is the RF trace: plotted against MCMC state it shows how far the
    /// chain wandered from a chosen point in tree space, which is how
    /// `ess/rf_trace.py` diagnoses mixing. Reading a row rather than
    /// recomputing costs nothing, since the matrix is already cached.
    #[wasm_bindgen(js_name = rfTrace)]
    pub fn rf_trace(&self, reference: usize) -> Result<Vec<u32>, JsError> {
        let n = self.snaps.len();
        if reference >= n {
            return Err(JsError::new(&format!(
                "Reference tree {reference} is out of range (0..{n})."
            )));
        }
        let rf = self.ensure_rf();
        Ok(rf[reference * n..(reference + 1) * n].to_vec())
    }

    /// Pseudo-ESS across a spread of reference trees (Lanfear et al. 2016).
    ///
    /// Topology has no scalar trace to run ESS on, but RF distance to a
    /// reference tree is one. Different references give different answers, so
    /// this reports a distribution; the conservative reading of "the chain's
    /// ESS" is the minimum.
    ///
    /// References are evenly spaced rather than random, so the number is
    /// reproducible — see [`crate::ess::evenly_spaced_refs`].
    #[wasm_bindgen(js_name = pseudoEss)]
    pub fn pseudo_ess(&self, n_refs: usize) -> Result<PseudoEssResult, JsError> {
        let n = self.snaps.len();
        if n < 4 {
            return Err(JsError::new("Pseudo-ESS needs at least 4 trees."));
        }
        let rf = self.ensure_rf();
        let refs = ess::evenly_spaced_refs(n, n_refs.max(1));
        let out = ess::pseudo_ess(&rf, n, &refs);
        Ok(PseudoEssResult {
            ess_values: out.ess_values,
            ref_indices: out.ref_indices.iter().map(|&v| v as u32).collect(),
            min: out.min,
            median: out.median,
            max: out.max,
        })
    }

    /// Pseudo-ESS for one run on its own.
    ///
    /// A pooled figure across several chains measures how far apart the chains
    /// sit rather than how well any of them mixed, so the per-run breakdown is
    /// the number that actually diagnoses a chain.
    #[wasm_bindgen(js_name = pseudoEssForRun)]
    pub fn pseudo_ess_for_run(&self, run: u32, n_refs: usize) -> Result<PseudoEssResult, JsError> {
        let idx: Vec<usize> = self
            .run_ids
            .iter()
            .enumerate()
            .filter(|&(_, &r)| r == run)
            .map(|(i, _)| i)
            .collect();
        if idx.len() < 4 {
            return Err(JsError::new(
                "That run has fewer than 4 trees — too few for an ESS.",
            ));
        }
        let rf = self.ensure_rf();
        let out = ess::pseudo_ess_for(&rf, self.snaps.len(), &idx, n_refs.max(1));
        Ok(PseudoEssResult {
            ess_values: out.ess_values,
            ref_indices: out.ref_indices.iter().map(|&v| v as u32).collect(),
            min: out.min,
            median: out.median,
            max: out.max,
        })
    }

    /// Per-run frequency of every bipartition.
    ///
    /// Returns a flat `n_runs x n_bipartitions` row-major array. A clade at
    /// 0.98 in one run and 0.11 in another is the clearest single sign that
    /// two chains have not converged on the same posterior.
    #[wasm_bindgen(js_name = cladeFrequencies)]
    pub fn clade_frequencies(&self) -> CladeFreqResult {
        let (presence, n_bip) = self.ensure_presence();
        let n_runs = self.run_labels.len();
        let out =
            clades::clade_frequencies(&presence, self.snaps.len(), n_bip, &self.run_ids, n_runs);
        CladeFreqResult {
            freq: out.freq,
            n_runs,
            n_bipartitions: n_bip,
        }
    }

    /// Clade counts for a selection, paired with the taxa each clade contains.
    ///
    /// The taxon lists are what let a clade found in the frequency scatter be
    /// pointed at in a tree: a bipartition is only meaningful as the set of
    /// tips on one side of it.
    #[wasm_bindgen(js_name = cladeMembership)]
    pub fn clade_membership(&self) -> CladeMembership {
        let (_, n_bip) = self.ensure_presence();
        let col_to_bip = self.snaps.sorted_bip_id_to_col();
        // `sorted_bip_id_to_col` maps split id -> column; invert it so we can
        // walk columns in the order the presence matrix uses.
        let mut bip_of_col = vec![0usize; n_bip];
        for (bip_id, &col) in col_to_bip.iter().enumerate() {
            if col < n_bip {
                bip_of_col[col] = bip_id;
            }
        }

        let leaves = &self.snaps.leaf_names;
        let mut sizes = Vec::with_capacity(n_bip);
        let mut flat: Vec<u32> = Vec::new();
        let mut offsets: Vec<u32> = Vec::with_capacity(n_bip + 1);
        offsets.push(0);

        for &bip_id in &bip_of_col {
            let bip = &self.snaps.bipartitions[bip_id];
            let mut count = 0u32;
            for (leaf_idx, _) in leaves.iter().enumerate() {
                if bip.get(leaf_idx) {
                    flat.push(leaf_idx as u32);
                    count += 1;
                }
            }
            sizes.push(count);
            offsets.push(flat.len() as u32);
        }

        CladeMembership {
            sizes,
            members: flat,
            offsets,
            leaf_names: leaves.clone(),
        }
    }

    /// Index of the maximum-clade-credibility tree.
    #[wasm_bindgen(js_name = mccIndex)]
    pub fn mcc_index(&self) -> Result<MccInfo, JsError> {
        let (presence, n_bip) = self.ensure_presence();
        let out = clades::max_clade_credibility(&presence, self.snaps.len(), n_bip)
            .ok_or_else(|| JsError::new("No trees to summarise."))?;
        Ok(MccInfo {
            index: out.index,
            log_clade_credibility: out.log_clade_credibility,
            name: self.names.get(out.index).cloned().unwrap_or_default(),
        })
    }

    /// Build a HIPSTR tree over a selection (empty = the whole set).
    ///
    /// Unlike MCC this is not one of the sampled trees — it assembles the
    /// best-scoring combination of clades that were observed together, so it
    /// routinely scores higher than any individual sample. `majority_rule`
    /// switches to MrHIPSTR, which forces in every clade above 0.5 support.
    ///
    /// Ported from beast-mcmc's `HIPSTRTreeBuilder` (Baele et al. 2025).
    #[wasm_bindgen(js_name = hipstr)]
    pub fn hipstr(&self, indices: &[u32], majority_rule: bool) -> Result<HipstrResult, JsError> {
        let picks: Vec<usize> = if indices.is_empty() {
            (0..self.newicks.len()).collect()
        } else {
            indices.iter().map(|&v| v as usize).collect()
        };
        if picks.len() < 2 {
            return Err(JsError::new("HIPSTR needs at least two trees."));
        }

        let newicks: Vec<String> = picks
            .iter()
            .filter_map(|&i| self.newicks.get(i).cloned())
            .collect();
        let run_of: Vec<usize> = picks
            .iter()
            .map(|&i| *self.run_ids.get(i).unwrap_or(&0) as usize)
            .collect();

        let cs = hipstr_mod::build_clade_system(
            &newicks,
            &self.translates,
            &run_of
                .iter()
                .map(|&r| self.translate_of_run.get(r).copied().unwrap_or(0))
                .collect::<Vec<_>>(),
            &self.snaps.leaf_names,
        )
        .map_err(|e| JsError::new(&e))?;

        let out = hipstr_mod::hipstr(&cs, majority_rule).map_err(|e| JsError::new(&e))?;
        Ok(HipstrResult {
            newick: out.newick,
            log_credibility: out.log_credibility,
            n_trees: newicks.len(),
            n_clades: cs.n_clades(),
            mean_support: if out.clade_support.is_empty() {
                0.0
            } else {
                out.clade_support.iter().sum::<f64>() / out.clade_support.len() as f64
            },
        })
    }

    /// Newick string for one tree, with taxon names restored.
    ///
    /// The stored newick keeps BEAST's numeric tip labels — that is what the
    /// translate block is for, and applying it during snapshot construction
    /// only renames the parsed copy. Handing the raw string to a viewer showed
    /// tips as "35", "9", and made clade names impossible to match against a
    /// drawn tree.
    ///
    /// Renaming happens here, on the one tree being asked for, rather than
    /// eagerly for all of them: a 4000-tree set would otherwise pay a parse
    /// and re-emit per tree for output nobody looks at.
    #[wasm_bindgen(js_name = newick)]
    pub fn newick(&self, index: usize) -> Result<String, JsError> {
        let raw = self
            .newicks
            .get(index)
            .ok_or_else(|| JsError::new(&format!("Tree {index} is out of range.")))?;

        let run = *self.run_ids.get(index).unwrap_or(&0) as usize;
        let slot = self.translate_of_run.get(run).copied().unwrap_or(0);
        let Some(translate) = self.translates.get(slot) else {
            return Ok(raw.clone());
        };
        if translate.is_empty() {
            return Ok(raw.clone());
        }

        let mut tree = phylotree::tree::Tree::from_newick(raw.trim())
            .map_err(|e| JsError::new(&format!("Could not parse tree {index}: {e}")))?;
        rapidtrees::io::rename_leaf_nodes(&mut tree, translate);
        tree.to_newick()
            .map_err(|e| JsError::new(&format!("Could not write tree {index}: {e}")))
    }

    /// Maximum-clade-credibility tree of a *selection*.
    ///
    /// Frequencies are recomputed within the selection: a cluster in tree
    /// space is its own small posterior, and the tree that represents it is the
    /// one whose clades are commonest among its neighbours — not the global
    /// winner. This is what makes "lasso a cluster, summarise it, do it again
    /// for the other cluster, compare" a meaningful workflow.
    #[wasm_bindgen(js_name = mccForSelection)]
    pub fn mcc_for_selection(&self, indices: &[u32]) -> Result<MccInfo, JsError> {
        if indices.is_empty() {
            return Err(JsError::new("Select some trees first."));
        }
        let (presence, n_bip) = self.ensure_presence();
        let idx: Vec<usize> = indices.iter().map(|&v| v as usize).collect();
        let out = clades::max_clade_credibility_for(&presence, self.snaps.len(), n_bip, &idx)
            .ok_or_else(|| JsError::new("Could not summarise that selection."))?;
        Ok(MccInfo {
            index: out.index,
            log_clade_credibility: out.log_clade_credibility,
            name: self.names.get(out.index).cloned().unwrap_or_default(),
        })
    }

    /// Bipartition counts within a selection.
    ///
    /// Saved alongside a summary tree so two saved trees can be compared later
    /// without re-deriving anything: both index the same bipartition basis, so
    /// the comparison is column alignment.
    #[wasm_bindgen(js_name = cladeCountsForSelection)]
    pub fn clade_counts_for_selection(&self, indices: &[u32]) -> CladeCountsResult {
        let (presence, n_bip) = self.ensure_presence();
        let idx: Vec<usize> = indices.iter().map(|&v| v as usize).collect();
        let out = clades::clade_counts_for(&presence, n_bip, &idx, self.snaps.len());
        CladeCountsResult {
            counts: out.counts,
            n_trees: out.n_trees,
        }
    }

    fn ensure_presence(&self) -> (Vec<u8>, usize) {
        if let Some((p, n)) = self.presence_cache.borrow().as_ref() {
            return (p.clone(), *n);
        }
        let (presence, col_to_bip) = self.snaps.build_presence_matrix();
        let n_bip = col_to_bip.len();
        *self.presence_cache.borrow_mut() = Some((presence.clone(), n_bip));
        (presence, n_bip)
    }

    fn ensure_rf(&self) -> Vec<u32> {
        if let Some(cached) = self.rf_cache.borrow().as_ref() {
            return cached.clone();
        }
        let rf: Vec<u32> = self
            .snaps
            .pairwise_rf(None)
            .into_iter()
            .map(|v| v as u32)
            .collect();
        *self.rf_cache.borrow_mut() = Some(rf.clone());
        rf
    }
}

/// PCoA on a distance matrix supplied by the caller.
///
/// Lets the UI re-project an RF matrix already stored in OPFS — changing the
/// number of axes, say — without repeating the O(n²) pairwise pass.
#[wasm_bindgen(js_name = pcoaFromMatrix)]
pub fn pcoa_from_matrix(matrix: &[u32], n: usize, k: usize) -> Result<PcoaResult, JsError> {
    let out = mds::pcoa(matrix, n, k)
        .ok_or_else(|| JsError::new("PCoA needs at least two trees and one axis."))?;
    Ok(PcoaResult::from(out))
}

/// Gaussian kernel density estimate, for the marginal density drawn beside a
/// trace plot.
///
/// Bandwidth follows Scott's rule, matching `scipy.stats.gaussian_kde`'s
/// default — the curves are meant to agree with the desktop app's.
#[wasm_bindgen(js_name = gaussianKde)]
pub fn gaussian_kde(values: &[f64], n_grid: usize) -> Option<KdeResult> {
    crate::stats::gaussian_kde(values, n_grid).map(|(grid, density)| KdeResult { grid, density })
}

/// A kernel density estimate on a regular grid.
#[wasm_bindgen]
pub struct KdeResult {
    grid: Vec<f64>,
    density: Vec<f64>,
}

#[wasm_bindgen]
impl KdeResult {
    #[wasm_bindgen(getter)]
    pub fn grid(&self) -> Vec<f64> {
        self.grid.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn density(&self) -> Vec<f64> {
        self.density.clone()
    }
}

/// A PCoA embedding, handed to JS as typed arrays.
#[wasm_bindgen]
pub struct PcoaResult {
    coords: Vec<f64>,
    explained: Vec<f64>,
    eigenvalues: Vec<f64>,
    n: usize,
    k: usize,
}

impl From<mds::Pcoa> for PcoaResult {
    fn from(p: mds::Pcoa) -> Self {
        PcoaResult {
            coords: p.coords,
            explained: p.explained,
            eigenvalues: p.eigenvalues,
            n: p.n,
            k: p.k,
        }
    }
}

#[wasm_bindgen]
impl PcoaResult {
    /// Row-major `n × k` coordinates.
    #[wasm_bindgen(getter)]
    pub fn coords(&self) -> Vec<f64> {
        self.coords.clone()
    }
    /// Fraction of positive eigenvalue mass per axis — for axis labels.
    #[wasm_bindgen(getter)]
    pub fn explained(&self) -> Vec<f64> {
        self.explained.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn eigenvalues(&self) -> Vec<f64> {
        self.eigenvalues.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn n(&self) -> usize {
        self.n
    }
    #[wasm_bindgen(getter)]
    pub fn k(&self) -> usize {
        self.k
    }
}

/// Pseudo-ESS across a set of reference trees.
#[wasm_bindgen]
pub struct PseudoEssResult {
    ess_values: Vec<f64>,
    ref_indices: Vec<u32>,
    min: f64,
    median: f64,
    max: f64,
}

#[wasm_bindgen]
impl PseudoEssResult {
    /// One ESS estimate per reference tree.
    #[wasm_bindgen(getter, js_name = essValues)]
    pub fn ess_values(&self) -> Vec<f64> {
        self.ess_values.clone()
    }
    #[wasm_bindgen(getter, js_name = refIndices)]
    pub fn ref_indices(&self) -> Vec<u32> {
        self.ref_indices.clone()
    }
    /// The conservative summary — the tightest bound across reference choices.
    #[wasm_bindgen(getter)]
    pub fn min(&self) -> f64 {
        self.min
    }
    #[wasm_bindgen(getter)]
    pub fn median(&self) -> f64 {
        self.median
    }
    #[wasm_bindgen(getter)]
    pub fn max(&self) -> f64 {
        self.max
    }
}

/// Per-run bipartition frequencies, row-major `n_runs x n_bipartitions`.
#[wasm_bindgen]
pub struct CladeFreqResult {
    freq: Vec<f64>,
    n_runs: usize,
    n_bipartitions: usize,
}

#[wasm_bindgen]
impl CladeFreqResult {
    #[wasm_bindgen(getter)]
    pub fn freq(&self) -> Vec<f64> {
        self.freq.clone()
    }
    #[wasm_bindgen(getter, js_name = nRuns)]
    pub fn n_runs(&self) -> usize {
        self.n_runs
    }
    #[wasm_bindgen(getter, js_name = nBipartitions)]
    pub fn n_bipartitions(&self) -> usize {
        self.n_bipartitions
    }
}

/// The maximum-clade-credibility tree of a set.
#[wasm_bindgen]
pub struct MccInfo {
    index: usize,
    log_clade_credibility: f64,
    name: String,
}

#[wasm_bindgen]
impl MccInfo {
    #[wasm_bindgen(getter)]
    pub fn index(&self) -> usize {
        self.index
    }
    #[wasm_bindgen(getter, js_name = logCladeCredibility)]
    pub fn log_clade_credibility(&self) -> f64 {
        self.log_clade_credibility
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.name.clone()
    }
}

/// Bipartition counts over a selection of trees.
#[wasm_bindgen]
pub struct CladeCountsResult {
    counts: Vec<u32>,
    n_trees: usize,
}

#[wasm_bindgen]
impl CladeCountsResult {
    #[wasm_bindgen(getter)]
    pub fn counts(&self) -> Vec<u32> {
        self.counts.clone()
    }
    /// Trees actually counted — the divisor for turning counts into frequencies.
    #[wasm_bindgen(getter, js_name = nTrees)]
    pub fn n_trees(&self) -> usize {
        self.n_trees
    }
}

/// Rectangular layout of a tree, for drawing and tanglegrams.
///
/// Matches the desktop app's conventions: x is cumulative branch length from
/// the root, y is tip rank with internal nodes at the mean of their children.
#[wasm_bindgen(js_name = treeLayout)]
pub fn tree_layout(newick: &str, ladderize: bool) -> Result<LayoutResult, JsError> {
    let out = layout::layout(newick, ladderize).map_err(|e| JsError::new(&e))?;
    Ok(LayoutResult {
        x: out.x,
        y: out.y,
        parent: out.parent,
        tip_names: out.tip_names,
        tip_index: out.tip_index,
    })
}

#[wasm_bindgen]
pub struct LayoutResult {
    x: Vec<f64>,
    y: Vec<f64>,
    parent: Vec<i32>,
    tip_names: Vec<String>,
    tip_index: Vec<u32>,
}

#[wasm_bindgen]
impl LayoutResult {
    /// Cumulative branch length from the root, one per node.
    #[wasm_bindgen(getter)]
    pub fn x(&self) -> Vec<f64> {
        self.x.clone()
    }
    /// Vertical position, one per node.
    #[wasm_bindgen(getter)]
    pub fn y(&self) -> Vec<f64> {
        self.y.clone()
    }
    /// Parent slot of each node; -1 for the root.
    #[wasm_bindgen(getter)]
    pub fn parent(&self) -> Vec<i32> {
        self.parent.clone()
    }
    /// Tip names in drawing order.
    #[wasm_bindgen(getter, js_name = tipNames)]
    pub fn tip_names(&self) -> Vec<String> {
        self.tip_names.clone()
    }
    /// Node slot for each tip name.
    #[wasm_bindgen(getter, js_name = tipIndex)]
    pub fn tip_index(&self) -> Vec<u32> {
        self.tip_index.clone()
    }
}

/// Which taxa sit on the canonical side of each bipartition.
///
/// Flattened rather than nested: a 550k-bipartition analysis would otherwise
/// mean half a million JS arrays crossing the boundary. `members[offsets[j] ..
/// offsets[j+1]]` gives the leaf indices for clade `j`.
#[wasm_bindgen]
pub struct CladeMembership {
    sizes: Vec<u32>,
    members: Vec<u32>,
    offsets: Vec<u32>,
    leaf_names: Vec<String>,
}

#[wasm_bindgen]
impl CladeMembership {
    /// Number of taxa in each clade — a cheap proxy for how deep it sits.
    #[wasm_bindgen(getter)]
    pub fn sizes(&self) -> Vec<u32> {
        self.sizes.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn members(&self) -> Vec<u32> {
        self.members.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn offsets(&self) -> Vec<u32> {
        self.offsets.clone()
    }
    #[wasm_bindgen(getter, js_name = leafNames)]
    pub fn leaf_names(&self) -> Vec<String> {
        self.leaf_names.clone()
    }
}

/// A HIPSTR reconstruction.
#[wasm_bindgen]
pub struct HipstrResult {
    newick: String,
    log_credibility: f64,
    n_trees: usize,
    n_clades: usize,
    mean_support: f64,
}

#[wasm_bindgen]
impl HipstrResult {
    #[wasm_bindgen(getter)]
    pub fn newick(&self) -> String {
        self.newick.clone()
    }
    /// Sum of log clade credibilities — directly comparable with an MCC score
    /// over the same trees.
    #[wasm_bindgen(getter, js_name = logCredibility)]
    pub fn log_credibility(&self) -> f64 {
        self.log_credibility
    }
    #[wasm_bindgen(getter, js_name = nTrees)]
    pub fn n_trees(&self) -> usize {
        self.n_trees
    }
    /// Distinct clades observed across the input trees.
    #[wasm_bindgen(getter, js_name = nClades)]
    pub fn n_clades(&self) -> usize {
        self.n_clades
    }
    /// Mean posterior support of the clades in the reconstructed tree.
    #[wasm_bindgen(getter, js_name = meanSupport)]
    pub fn mean_support(&self) -> f64 {
        self.mean_support
    }
}
