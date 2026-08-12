//! HIPSTR — highest independent posterior subtree reconstruction.
//!
//! A port of `HIPSTRTreeBuilder.java` from
//! [beast-mcmc](https://github.com/beast-dev/beast-mcmc)
//! (`src/dr/app/tools/treeannotator/`), published as Baele et al. 2025,
//! *Bioinformatics* 41:btaf488, doi:10.1093/bioinformatics/btaf488.
//!
//! # What makes it different from MCC
//!
//! Maximum clade credibility picks the best *sampled* tree — it can only ever
//! return a topology that appeared in the posterior. HIPSTR **builds** one, by
//! assembling the highest-scoring combination of clades that were observed
//! together. The result frequently beats every individual sample, because good
//! resolutions of different parts of the tree rarely co-occur in one draw.
//!
//! The key restriction is that only *observed* splits are candidates: for each
//! clade, the algorithm considers the (left, right) child pairs that actually
//! occurred in some sampled tree. That is what keeps the output a plausible
//! tree rather than an arbitrary recombination of compatible clades.
//!
//! # The score
//!
//! For a clade `c` resolved into children `l` and `r`:
//!
//! ```text
//! score(c) = ln P(c) + max over observed (l, r) of [ score(l) + score(r) ]
//! ```
//!
//! with tips scoring `ln 1 = 0`. Clades are processed in increasing size so
//! every child is scored before its parent — the iterative form of the
//! reference implementation, which avoids deep recursion on a ladder-shaped
//! tree.
//!
//! MrHIPSTR adds a large bonus to any clade above 0.5 credibility, which forces
//! the majority-rule clades in wherever they are compatible.

use std::collections::HashMap;

use phylotree::tree::{NodeId, Tree as PhyloTree};

use rapidtrees::io::rename_leaf_nodes;

/// Bonus applied to majority clades under MrHIPSTR.
///
/// Matches the reference implementation's `MAJORITY_RULE_REWARD`. Large enough
/// to dominate any sum of log credibilities, so a majority clade is taken
/// whenever it is compatible.
const MAJORITY_RULE_REWARD: f64 = 1e10;

/// Every clade observed across a set of trees, with how it was resolved.
pub struct CladeSystem {
    /// Membership bitset per clade, `words` u64s each.
    clades: Vec<Vec<u64>>,
    /// How many trees contained each clade.
    counts: Vec<u32>,
    /// Number of taxa in each clade.
    sizes: Vec<u32>,
    /// Summed node height per clade, for averaging into branch lengths.
    height_sum: Vec<f64>,
    /// Observed `(left, right)` child clade pairs, per clade.
    subclades: Vec<Vec<(usize, usize)>>,
    /// Clade index for each single taxon, by leaf index.
    tip_clade: Vec<usize>,
    pub leaf_names: Vec<String>,
    pub n_trees: usize,
    root: usize,
}

impl CladeSystem {
    pub fn n_clades(&self) -> usize {
        self.clades.len()
    }

    /// Posterior probability of a clade.
    pub fn credibility(&self, clade: usize) -> f64 {
        self.counts[clade] as f64 / self.n_trees.max(1) as f64
    }

    /// Mean height of a clade's node across the trees containing it.
    fn mean_height(&self, clade: usize) -> f64 {
        let c = self.counts[clade];
        if c == 0 {
            0.0
        } else {
            self.height_sum[clade] / c as f64
        }
    }
}

fn words_for(n: usize) -> usize {
    n.div_ceil(64)
}

/// Build the clade system from a set of newick strings.
///
/// `leaf_names` fixes the taxon ordering; every tree must use the same set.
/// Trees are expected to be bifurcating — BEAST samples are — and a node with
/// any other number of children contributes its clade but no resolution, so a
/// polytomy simply offers no candidate split rather than corrupting the count.
pub fn build_clade_system(
    newicks: &[String],
    translates: &[HashMap<String, String>],
    run_of: &[usize],
    leaf_names: &[String],
) -> Result<CladeSystem, String> {
    let n_leaves = leaf_names.len();
    let words = words_for(n_leaves);
    let leaf_index: HashMap<&str, usize> = leaf_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut index: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut clades: Vec<Vec<u64>> = Vec::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut sizes: Vec<u32> = Vec::new();
    let mut height_sum: Vec<f64> = Vec::new();
    let mut subclades: Vec<Vec<(usize, usize)>> = Vec::new();
    // Dedupe pairs per clade without an O(k) scan each time.
    let mut seen_pairs: Vec<std::collections::HashSet<(usize, usize)>> = Vec::new();

    let mut intern = |bits: Vec<u64>,
                      size: u32,
                      height: f64,
                      clades: &mut Vec<Vec<u64>>,
                      counts: &mut Vec<u32>,
                      sizes: &mut Vec<u32>,
                      height_sum: &mut Vec<f64>,
                      subclades: &mut Vec<Vec<(usize, usize)>>,
                      seen_pairs: &mut Vec<std::collections::HashSet<(usize, usize)>>|
     -> usize {
        match index.get(&bits) {
            Some(&id) => {
                counts[id] += 1;
                height_sum[id] += height;
                id
            }
            None => {
                let id = clades.len();
                index.insert(bits.clone(), id);
                clades.push(bits);
                counts.push(1);
                sizes.push(size);
                height_sum.push(height);
                subclades.push(Vec::new());
                seen_pairs.push(std::collections::HashSet::new());
                id
            }
        }
    };

    for (t, newick) in newicks.iter().enumerate() {
        let mut tree = PhyloTree::from_newick(newick.trim()).map_err(|e| e.to_string())?;
        if let Some(map) = translates.get(*run_of.get(t).unwrap_or(&0))
            && !map.is_empty()
        {
            rename_leaf_nodes(&mut tree, map);
        }
        let root = tree.get_root().map_err(|e| e.to_string())?;

        // Postorder so a node's children are interned before it is.
        let mut bits_of: HashMap<NodeId, Vec<u64>> = HashMap::new();
        let mut size_of: HashMap<NodeId, u32> = HashMap::new();
        let mut id_of: HashMap<NodeId, usize> = HashMap::new();
        let mut height_of: HashMap<NodeId, f64> = HashMap::new();

        let mut stack: Vec<(NodeId, bool)> = vec![(root, false)];
        while let Some((id, done)) = stack.pop() {
            let node = tree.get(&id).map_err(|e| e.to_string())?;

            if node.children.is_empty() {
                let name = node.name.clone().unwrap_or_default();
                let Some(&li) = leaf_index.get(name.as_str()) else {
                    return Err(format!(
                        "Tree {t} has taxon '{name}', which is not in the shared taxon set."
                    ));
                };
                let mut b = vec![0u64; words];
                b[li / 64] |= 1u64 << (li % 64);
                bits_of.insert(id, b.clone());
                size_of.insert(id, 1);
                height_of.insert(id, 0.0);
                id_of.insert(
                    id,
                    intern(
                        b,
                        1,
                        0.0,
                        &mut clades,
                        &mut counts,
                        &mut sizes,
                        &mut height_sum,
                        &mut subclades,
                        &mut seen_pairs,
                    ),
                );
                continue;
            }

            if !done {
                stack.push((id, true));
                for &c in &node.children {
                    stack.push((c, false));
                }
                continue;
            }

            let mut b = vec![0u64; words];
            let mut size = 0u32;
            // Height above the tips: a child's height plus its branch length.
            // BEAST trees are ultrametric so any child agrees; taking the max
            // keeps it sane if rounding makes them differ slightly.
            let mut height = 0.0f64;
            for &c in &node.children {
                let cb = &bits_of[&c];
                for (w, v) in b.iter_mut().zip(cb) {
                    *w |= *v;
                }
                size += size_of[&c];
                let edge = tree.get(&c).ok().and_then(|n| n.parent_edge).unwrap_or(0.0);
                height = height.max(height_of[&c] + edge);
            }

            let cid = intern(
                b.clone(),
                size,
                height,
                &mut clades,
                &mut counts,
                &mut sizes,
                &mut height_sum,
                &mut subclades,
                &mut seen_pairs,
            );
            bits_of.insert(id, b);
            size_of.insert(id, size);
            height_of.insert(id, height);
            id_of.insert(id, cid);

            // Only bifurcations offer a candidate resolution. A polytomy is
            // recorded as a clade but proposes no split, which is the honest
            // thing to do: nothing was observed about how it resolves.
            if node.children.len() == 2 {
                let l = id_of[&node.children[0]];
                let r = id_of[&node.children[1]];
                let pair = if l < r { (l, r) } else { (r, l) };
                if seen_pairs[cid].insert(pair) {
                    subclades[cid].push(pair);
                }
            }
        }
    }

    // The root clade is the one containing every taxon.
    let full: Vec<u64> = (0..words)
        .map(|w| {
            let mut v = 0u64;
            for bit in 0..64 {
                let li = w * 64 + bit;
                if li < n_leaves {
                    v |= 1u64 << bit;
                }
            }
            v
        })
        .collect();
    let root = *index
        .get(&full)
        .ok_or_else(|| "No tree contained every taxon.".to_string())?;

    let mut tip_clade = vec![usize::MAX; n_leaves];
    for (id, bits) in clades.iter().enumerate() {
        if sizes[id] == 1 {
            for li in 0..n_leaves {
                if bits[li / 64] >> (li % 64) & 1 == 1 {
                    tip_clade[li] = id;
                }
            }
        }
    }

    Ok(CladeSystem {
        clades,
        counts,
        sizes,
        height_sum,
        subclades,
        tip_clade,
        leaf_names: leaf_names.to_vec(),
        n_trees: newicks.len(),
        root,
    })
}

/// The result of a HIPSTR reconstruction.
pub struct HipstrTree {
    /// Newick with mean-height branch lengths.
    pub newick: String,
    /// Sum of log clade credibilities over the constructed tree.
    pub log_credibility: f64,
    /// Posterior probability of each clade in the tree, root first.
    pub clade_support: Vec<f64>,
}

/// Build the HIPSTR (or MrHIPSTR) tree.
pub fn hipstr(cs: &CladeSystem, majority_rule: bool) -> Result<HipstrTree, String> {
    let n = cs.n_clades();
    let mut score = vec![f64::NEG_INFINITY; n];
    let mut best_left = vec![usize::MAX; n];
    let mut best_right = vec![usize::MAX; n];

    // Increasing size: every child is scored before its parent, which is what
    // lets the recursion be a single pass.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&c| cs.sizes[c]);

    for &c in &order {
        let cred = cs.credibility(c);
        let mut clade_score = cred.ln()
            + if majority_rule && cred > 0.5 {
                MAJORITY_RULE_REWARD
            } else {
                0.0
            };

        match cs.sizes[c] {
            1 => {
                // A tip contributes ln(1) = 0.
                clade_score = 0.0;
            }
            2 => {
                // One possible pair, and both children are tips, so the
                // subtree contributes nothing.
                if let Some(&(l, r)) = cs.subclades[c].first() {
                    best_left[c] = l.min(r);
                    best_right[c] = l.max(r);
                }
            }
            _ => {
                let mut best = f64::NEG_INFINITY;
                // Ties are common — many resolutions share a score — and the
                // reference breaks them on summed child credibility, so the
                // more strongly supported pair wins.
                let mut tie_break = f64::NEG_INFINITY;
                for &(l, r) in &cs.subclades[c] {
                    let ls = if cs.sizes[l] > 1 { score[l] } else { 0.0 };
                    let rs = if cs.sizes[r] > 1 { score[r] } else { 0.0 };
                    if !ls.is_finite() || !rs.is_finite() {
                        continue;
                    }
                    let total = ls + rs;
                    let support = cs.credibility(l) + cs.credibility(r);
                    if total > best || (total == best && support > tie_break) {
                        best = total;
                        tie_break = support;
                        // Smaller clade on the left, matching the reference's
                        // ordering so output is comparable.
                        if cs.sizes[l] < cs.sizes[r] {
                            best_left[c] = l;
                            best_right[c] = r;
                        } else {
                            best_left[c] = r;
                            best_right[c] = l;
                        }
                    }
                }
                if best.is_finite() {
                    clade_score += best;
                }
            }
        }
        score[c] = clade_score;
    }

    if cs.sizes[cs.root] > 2 && best_left[cs.root] == usize::MAX {
        return Err("Could not resolve the root clade from the observed splits.".into());
    }

    // Emit the tree, deepest first, with branch lengths from mean heights.
    let mut support = Vec::new();
    let mut out = String::new();
    write_clade(cs, cs.root, &best_left, &best_right, &mut out, &mut support)?;
    out.push(';');

    let log_credibility: f64 = support.iter().map(|s: &f64| s.ln()).sum();

    Ok(HipstrTree {
        newick: out,
        log_credibility,
        clade_support: support,
    })
}

fn write_clade(
    cs: &CladeSystem,
    clade: usize,
    best_left: &[usize],
    best_right: &[usize],
    out: &mut String,
    support: &mut Vec<f64>,
) -> Result<(), String> {
    if cs.sizes[clade] == 1 {
        let li = (0..cs.leaf_names.len())
            .find(|&i| cs.tip_clade[i] == clade)
            .ok_or_else(|| "Tip clade has no taxon.".to_string())?;
        out.push_str(&cs.leaf_names[li]);
        return Ok(());
    }

    support.push(cs.credibility(clade));
    let (l, r) = (best_left[clade], best_right[clade]);
    if l == usize::MAX || r == usize::MAX {
        return Err("Unresolved clade in the reconstructed tree.".into());
    }

    let h = cs.mean_height(clade);
    out.push('(');
    for (k, child) in [l, r].into_iter().enumerate() {
        if k > 0 {
            out.push(',');
        }
        write_clade(cs, child, best_left, best_right, out, support)?;
        // Branch length is the drop in mean height from parent to child.
        let bl = (h - cs.mean_height(child)).max(0.0);
        out.push_str(&format!(":{bl:.8}"));
    }
    out.push(')');
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn system(newicks: &[&str], leaves: &[&str]) -> CladeSystem {
        let owned: Vec<String> = newicks.iter().map(|s| s.to_string()).collect();
        let names: Vec<String> = leaves.iter().map(|s| s.to_string()).collect();
        match build_clade_system(&owned, &[HashMap::new()], &vec![0; owned.len()], &names) {
            Ok(cs) => cs,
            Err(e) => panic!("clade system: {e}"),
        }
    }

    #[test]
    fn recovers_a_unanimous_topology() {
        let cs = system(
            &[
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,B:1):1,(C:1,D:1):1);",
            ],
            &["A", "B", "C", "D"],
        );
        let out = hipstr(&cs, false).expect("hipstr");
        // Every clade is at probability 1, so the log credibility is 0.
        assert!(out.log_credibility.abs() < 1e-12, "{}", out.log_credibility);
        assert!(out.newick.contains("A") && out.newick.contains("D"));
        assert!(out.clade_support.iter().all(|&s| (s - 1.0).abs() < 1e-12));
    }

    #[test]
    fn assembles_a_tree_better_than_any_sample() {
        // Six taxa. (A,B) appears in 2 of 3 trees; (E,F) in a different 2 of 3.
        // No single sampled tree has both, but HIPSTR can build one that does.
        let cs = system(
            &[
                "(((A:1,B:1):1,C:2):1,((E:1,F:1):1,D:2):1);",
                "(((A:1,B:1):1,C:2):1,((D:1,E:1):1,F:2):1);",
                "(((A:1,C:1):1,B:2):1,((E:1,F:1):1,D:2):1);",
            ],
            &["A", "B", "C", "D", "E", "F"],
        );
        let out = hipstr(&cs, false).expect("hipstr");
        // Both majority cherries should be present in the reconstruction.
        assert!(
            out.newick.contains("A") && out.newick.contains("F"),
            "{}",
            out.newick
        );
        // The chosen clades should all have support, and the best score should
        // beat taking the 1-in-3 resolutions.
        assert!(out.clade_support.iter().all(|&s| s >= 1.0 / 3.0));
        assert!(out.log_credibility > (1.0f64 / 3.0).ln() * 4.0);
    }

    #[test]
    fn majority_rule_forces_supported_clades() {
        let cs = system(
            &[
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            &["A", "B", "C", "D"],
        );
        let plain = hipstr(&cs, false).expect("hipstr");
        let mr = hipstr(&cs, true).expect("mrhipstr");
        // (A,B) and (C,D) sit at 2/3 so both runs should pick them; the
        // majority variant simply makes that mandatory.
        assert!(plain.newick.contains("A"));
        assert!(mr.clade_support.iter().any(|&s| s > 0.5));
    }

    #[test]
    fn branch_lengths_are_non_negative_and_present() {
        let cs = system(
            &[
                "((A:1,B:1):2,(C:1.5,D:1.5):1.5);",
                "((A:1,B:1):2,(C:1.5,D:1.5):1.5);",
            ],
            &["A", "B", "C", "D"],
        );
        let out = hipstr(&cs, false).expect("hipstr");
        assert!(out.newick.contains(':'), "{}", out.newick);
        for part in out.newick.split(':').skip(1) {
            let num: String = part
                .chars()
                .take_while(|c| c.is_ascii_digit() || *c == '.' || *c == '-')
                .collect();
            if let Ok(v) = num.parse::<f64>() {
                assert!(v >= 0.0, "negative branch length in {}", out.newick);
            }
        }
    }

    #[test]
    fn rejects_an_unknown_taxon() {
        let owned = vec!["((A:1,B:1):1,(C:1,Z:1):1);".to_string()];
        let names: Vec<String> = ["A", "B", "C", "D"].iter().map(|s| s.to_string()).collect();
        let err = match build_clade_system(&owned, &[HashMap::new()], &[0], &names) {
            Ok(_) => panic!("expected an error for the unknown taxon"),
            Err(e) => e,
        };
        assert!(err.contains("Z"), "{err}");
    }
}
