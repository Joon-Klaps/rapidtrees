use crate::snapshot::TreeSnapshot;
use phylotree::tree::Tree;
use std::collections::{HashMap, HashSet};
use std::fs;
#[cfg(feature = "cli")]
use std::io::{self, Write};
use std::path::Path;

#[cfg(feature = "cli")]
use flate2::Compression;
#[cfg(feature = "cli")]
use flate2::write::GzEncoder;

#[cfg(feature = "cli")]
type SnapReadResult = (Vec<String>, Vec<String>, usize, Vec<u8>);

/// Strip BEAST annotations from Newick strings.
///
/// BEAST format includes annotations like :[&rate=0.123]2.45 where 2.45 is the actual branch length.
/// This function removes the [&...] annotations while preserving the branch lengths.
/// We shouldn't be needing, this TODO: update phylotree to handle BEAST annotations directly.
pub fn strip_beast_annotations(newick: &str) -> String {
    let mut result = String::with_capacity(newick.len());
    let mut in_annotation = false;
    let mut chars = newick.chars().peekable();

    while let Some(ch) = chars.next() {
        if ch == '[' && chars.peek() == Some(&'&') {
            // Start of BEAST annotation
            in_annotation = true;
        } else if ch == ']' && in_annotation {
            // End of BEAST annotation
            in_annotation = false;
        } else if !in_annotation {
            // Only copy characters outside annotations
            result.push(ch);
        }
    }

    result
}

pub fn read_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> (HashMap<String, String>, Vec<(String, Tree)>) {
    let content = match fs::read_to_string(path.as_ref()) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Failed to read {:?}: {e}", path.as_ref());
            return (HashMap::new(), Vec::new());
        }
    };

    let base_name = path
        .as_ref()
        .file_name()
        .and_then(|s| s.to_str())
        .map(|s| s.trim_end_matches(".trees"))
        .unwrap_or("unknown");

    let taxons = parse_taxon_block(&content);

    let trees = collect_tree_blocks(&content)
        .into_iter()
        .enumerate()
        //generate tree name & extract state number
        .map(|(idx, tree)| {
            let (name, state) = extract_name_state(tree.header);
            (idx, tree, state, format!("{base_name}_{name}"))
        })
        // Filter out burn-in trees based on count and/or state number if 0 we don't filter
        .filter(|(idx, _tree, state, _name)| {
            (burnin_trees == 0 && burnin_states == 0)
                || (burnin_trees > 0 && *idx >= burnin_trees)
                || (burnin_states > 0 && *state > burnin_states)
        })
        // read in the files
        .filter_map(|(idx, tree, _state, name)| {
            // Strip BEAST annotations from newick string (e.g., [&rate=...])
            // BEAST format: :[&rate=X.XX]length -> :length
            let newick = strip_beast_annotations(&tree.body);
            let mut phylo_tree = match phylotree::tree::Tree::from_newick(&newick) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!(
                        "Failed to parse tree {} at index {}: {}",
                        path.as_ref().display(),
                        idx,
                        e
                    );
                    return None;
                }
            };

            // Rename the leaves with the map
            if use_real_taxa {
                rename_leaf_nodes(&mut phylo_tree, &taxons);
            }

            Some((name, phylo_tree))
        })
        .collect::<Vec<_>>();

    (taxons, trees)
}

fn extract_name_state(header: &str) -> (String, usize) {
    // tree classic2_STATE_968940000 [...]
    // tree STATE_8500000 [...]
    // tree classic1_STATE_766540000 [...]

    // Find "STATE_" in the header (case-insensitive)
    let upper = header.to_ascii_uppercase();
    if let Some(state_pos) = upper.find("STATE_") {
        // Split on first space to get the tree name part (after "tree ")
        if let Some((_, rest)) = header.split_once(' ') {
            // Extract just the tree name (everything before the opening bracket or end of string)
            let tree_name = rest.split_whitespace().next().unwrap_or("").to_string();

            // Extract the state number (digits after "STATE_")
            let digits = header[state_pos + 6..]
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>();

            if let Ok(num) = digits.parse::<usize>() {
                return (tree_name, num);
            }
        }
    }

    ("".to_string(), 0)
}

struct TreeBlock<'a> {
    header: &'a str,
    body: String,
}

fn collect_tree_blocks(content: &str) -> Vec<TreeBlock<'_>> {
    content
        .lines()
        .skip_while(|line| !line.to_ascii_uppercase().starts_with("TREE "))
        .take_while(|line| !line.trim().to_ascii_uppercase().starts_with("END;"))
        .filter_map(|line| {
            let mut parts = line.splitn(2, " = ");
            let header = parts.next()?.trim();
            let body = parts.next()?.trim().to_string();
            Some(TreeBlock { header, body })
        })
        .collect()
}

fn parse_taxon_block(content: &str) -> HashMap<String, String> {
    content
        .lines()
        .skip_while(|line| !line.trim().to_ascii_uppercase().starts_with("TRANSLATE"))
        .skip(1)
        .take_while(|line| !line.trim().to_ascii_uppercase().starts_with(";"))
        // STRUCTURE:
        // 1 '1959.M.CD.59.ZR59',
        // 2 '1960.DRC60A',
        // 3 GU573545|Soromba_R245|Mouse|MLI|2009,
        .filter_map(|line| {
            let line = line.trim().trim_end_matches(',');
            let mut parts = line.split_whitespace();
            let id = parts.next()?.to_string();
            let label = parts.next()?.trim_matches('\'').to_string();
            Some((id, label))
        })
        .collect::<HashMap<_, _>>()
}

pub fn rename_leaf_nodes(
    phylo_tree: &mut Tree,
    translate: &std::collections::HashMap<String, String>,
) {
    for leaf_id in phylo_tree.get_leaves() {
        if let Ok(node) = phylo_tree.get_mut(&leaf_id) {
            node.name = node.name.as_ref().and_then(|n| translate.get(n).cloned());
        }
    }
}

fn build_validated_snapshots(
    named_trees: Vec<(String, Tree)>,
    rooted: bool,
) -> Result<(Vec<String>, Vec<TreeSnapshot>), String> {
    if named_trees.len() < 2 {
        return Err(format!(
            "Need at least 2 trees to compute pairwise distances, found {}",
            named_trees.len()
        ));
    }

    let (names, trees): (Vec<_>, Vec<_>) = named_trees.into_iter().unzip();

    let leaf_sets: Vec<HashSet<String>> = trees
        .iter()
        .map(|tree| {
            tree.get_leaves()
                .iter()
                .filter_map(|&id| tree.get(&id).ok()?.name.clone())
                .collect()
        })
        .collect();

    leaf_sets
        .iter()
        .enumerate()
        .skip(1)
        .try_for_each(|(i, leaves)| {
            if leaves != &leaf_sets[0] {
                Err(format!(
                    "Tree {} has a different leaf set than tree 0. All trees must have the same taxa.",
                    i
                ))
            } else {
                Ok(())
            }
        })?;

    let snapshots = trees
        .iter()
        .enumerate()
        .map(|(i, tree)| {
            TreeSnapshot::from_tree(tree, rooted)
                .map_err(|e| format!("Failed to create snapshot for tree {}: {}", i, e))
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok((names, snapshots))
}

/// Read trees from multiple files and build their [`TreeSnapshot`]s in one pass.
///
/// Tree names are prefixed with `file{idx}_` to distinguish trees across files.
/// Returns an error if any file yields no trees after burn-in, if fewer than 2
/// trees total are found, or if any tree has a different leaf set than the first.
///
/// Returns `(names, snapshots)`.
pub fn read_and_snapshot_trees(
    paths: &[String],
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> Result<(Vec<String>, Vec<TreeSnapshot>), String> {
    let named_trees: Vec<(String, Tree)> = paths
        .iter()
        .enumerate()
        .map(|(file_idx, path)| {
            let (_taxons, file_trees) =
                read_beast_trees(path, burnin_trees, burnin_states, use_real_taxa);
            if file_trees.is_empty() {
                return Err(format!("No trees found in '{}' after burnin removal", path));
            }
            Ok(file_trees
                .into_iter()
                .map(move |(name, tree)| (format!("file{}_{}", file_idx, name), tree)))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
        .collect();

    build_validated_snapshots(named_trees, rooted)
}

/// Parse newick strings (with optional BEAST annotations), apply translate maps,
/// and build [`TreeSnapshot`]s in one pass.
///
/// Validates that all trees share the same leaf set. Returns an error if fewer
/// than 2 trees are provided, lengths mismatch, or any index is out of bounds.
///
/// Returns `(snapshots, sorted_leaf_names)`.
pub fn snapshot_newicks(
    newicks: &[String],
    translate_maps: &[HashMap<String, String>],
    map_indices: &[usize],
    rooted: bool,
) -> Result<(Vec<TreeSnapshot>, Vec<String>), String> {
    if newicks.len() != map_indices.len() {
        return Err(format!(
            "newicks length ({}) must equal map_indices length ({})",
            newicks.len(),
            map_indices.len()
        ));
    }
    map_indices.iter().enumerate().try_for_each(|(i, &idx)| {
        if idx >= translate_maps.len() {
            Err(format!(
                "map_indices[{}] = {} is out of bounds ({} translate maps provided)",
                i,
                idx,
                translate_maps.len()
            ))
        } else {
            Ok(())
        }
    })?;

    let named_trees: Vec<(String, Tree)> = newicks
        .iter()
        .enumerate()
        .map(|(i, nwk)| {
            let clean = strip_beast_annotations(nwk);
            let mut tree = Tree::from_newick(&clean)
                .map_err(|e| format!("Failed to parse newick at index {}: {}", i, e))?;
            rename_leaf_nodes(&mut tree, &translate_maps[map_indices[i]]);
            Ok((format!("tree_{}", i), tree))
        })
        .collect::<Result<_, String>>()?;

    let (_, snapshots) = build_validated_snapshots(named_trees, rooted)?;

    // Re-parsing the first newick is cheaper than keeping all `Tree`s alive
    // just to recover sorted leaf names after validation.
    let first_clean = strip_beast_annotations(&newicks[0]);
    let mut first_tree = Tree::from_newick(&first_clean)
        .map_err(|e| format!("Failed to re-parse first newick for leaf names: {}", e))?;
    rename_leaf_nodes(&mut first_tree, &translate_maps[map_indices[0]]);
    let mut leaf_names: Vec<String> = first_tree
        .get_leaves()
        .iter()
        .filter_map(|&id| first_tree.get(&id).ok()?.name.clone())
        .collect();
    leaf_names.sort_unstable();

    Ok((snapshots, leaf_names))
}

/// Write a labeled square matrix as TSV to a file or stdout.
/// If `path` ends with `.gz`, the output is gzip-compressed.
/// If `path` equals `-`, the matrix is written to stdout (uncompressed).
#[cfg(feature = "cli")]
pub fn write_matrix_tsv<P: AsRef<Path>, T: std::fmt::Display>(
    path: P,
    names: &[String],
    mat: &[Vec<T>],
) -> io::Result<()> {
    use std::fs::File;
    use std::io::BufWriter;

    let p = path.as_ref();
    if p.as_os_str() == "-" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "writing to stdout is not supported by write_matrix_tsv",
        ));
    }

    let is_gz = p.to_string_lossy().ends_with(".gz");

    let mut out: Box<dyn Write> = if is_gz {
        let f = File::create(p)?;
        let enc = GzEncoder::new(f, Compression::default());
        Box::new(BufWriter::new(enc))
    } else {
        Box::new(BufWriter::new(File::create(p)?))
    };

    // Header row
    write!(&mut out, "\t")?;
    for (k, name) in names.iter().enumerate() {
        if k > 0 {
            write!(&mut out, "\t")?;
        }
        write!(&mut out, "{}", name)?;
    }
    writeln!(&mut out)?;

    // Rows
    for (i, row) in mat.iter().enumerate() {
        write!(&mut out, "{}", names[i])?;
        for val in row {
            write!(&mut out, "\t{}", val)?;
        }
        writeln!(&mut out)?;
    }

    out.flush()?;
    Ok(())
}

// ── Snap helpers (private) ────────────────────────────────────────────────

#[cfg(feature = "cli")]
fn snap_read_u64<R: std::io::Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(feature = "cli")]
fn snap_read_strings<R: std::io::Read>(r: &mut R, n: usize) -> io::Result<Vec<String>> {
    (0..n)
        .map(|_| {
            let mut len_buf = [0u8; 4];
            r.read_exact(&mut len_buf)?;
            let len = u32::from_le_bytes(len_buf) as usize;
            let mut bytes = vec![0u8; len];
            r.read_exact(&mut bytes)?;
            String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .collect()
}

// ── Public snap API ───────────────────────────────────────────────────────

/// Write tree snapshots to a gzip-compressed binary `.snap` file.
///
/// # File layout (inside gzip stream)
/// ```text
/// HEADER     4 B  magic "SNAP"
///            1 B  version (1)
///            8 B  n_trees  u64 LE
///            8 B  n_taxa   u64 LE
///            8 B  n_bip    u64 LE
/// TAXA NAMES for each taxon:  4 B len u32 LE + N B UTF-8
/// TREE NAMES for each tree:   4 B len u32 LE + N B UTF-8
/// PRESENCE   n_trees × n_bip bytes, row-major uint8
///            presence[i][j] = 1 iff bipartition j is in tree i
/// ```
///
/// Bipartition column order is deterministic: globally sorted in ascending
/// bitset order so the same tree set always produces the same columns.
///
/// Note: `sum(presence[i] XOR presence[j]) == RF(tree_i, tree_j)` exactly.
#[cfg(feature = "cli")]
pub fn write_snap<P: AsRef<Path>>(
    path: P,
    tree_names: &[String],
    taxa_names: &[String],
    snapshots: &[TreeSnapshot],
) -> io::Result<()> {
    use flate2::{Compression, write::GzEncoder};
    use std::fs::File;
    use std::io::BufWriter;

    // Each snapshot's parts list is already sorted → flatten, sort, dedup
    // gives the global sorted-unique bipartition list in one pass.
    let mut all_bips: Vec<_> = snapshots
        .iter()
        .flat_map(|s| s.parts.iter().cloned())
        .collect();
    all_bips.sort_unstable();
    all_bips.dedup();

    let n_trees = snapshots.len();
    let n_taxa = taxa_names.len();
    let n_bip = all_bips.len();

    // Build row-major presence matrix.
    // Two-pointer merge is O(|parts| + n_bip) per tree (both sides sorted).
    let mut presence = vec![0u8; n_trees * n_bip];
    for (ti, snap) in snapshots.iter().enumerate() {
        let row = &mut presence[ti * n_bip..(ti + 1) * n_bip];
        let (mut gi, mut li) = (0, 0);
        while gi < n_bip && li < snap.parts.len() {
            match all_bips[gi].cmp(&snap.parts[li]) {
                std::cmp::Ordering::Equal => {
                    row[gi] = 1;
                    gi += 1;
                    li += 1;
                }
                std::cmp::Ordering::Less => gi += 1,
                std::cmp::Ordering::Greater => li += 1,
            }
        }
    }

    let file = File::create(path.as_ref())?;
    let mut w = BufWriter::new(GzEncoder::new(file, Compression::default()));

    // Header
    w.write_all(b"SNAP")?;
    w.write_all(&[1u8])?;
    w.write_all(&(n_trees as u64).to_le_bytes())?;
    w.write_all(&(n_taxa as u64).to_le_bytes())?;
    w.write_all(&(n_bip as u64).to_le_bytes())?;

    // Names
    for name in taxa_names.iter().chain(tree_names.iter()) {
        let b = name.as_bytes();
        w.write_all(&(b.len() as u32).to_le_bytes())?;
        w.write_all(b)?;
    }

    // Presence matrix
    w.write_all(&presence)?;
    w.flush()?;
    Ok(())
}

/// Read a `.snap` file produced by [`write_snap`].
///
/// Returns `(tree_names, taxa_names, n_bip, presence_bytes)` where
/// `presence_bytes` is a flat row-major `u8` buffer of shape `n_trees × n_bip`.
///
/// RF distance between tree `i` and tree `j` is:
/// ```text
/// presence[i*n_bip .. (i+1)*n_bip]
///     .zip(presence[j*n_bip .. (j+1)*n_bip])
///     .filter(|(a,b)| a != b)
///     .count()
/// ```
#[cfg(feature = "cli")]
pub fn read_snap<P: AsRef<Path>>(path: P) -> io::Result<SnapReadResult> {
    use flate2::read::GzDecoder;
    use std::fs::File;
    use std::io::{BufReader, Read};

    let file = File::open(path.as_ref())?;
    let mut r = BufReader::new(GzDecoder::new(file));

    // Header
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != b"SNAP" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid snap magic: expected SNAP, got {:?}", &magic),
        ));
    }
    let mut ver = [0u8; 1];
    r.read_exact(&mut ver)?;
    if ver[0] != 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported snap version: {}", ver[0]),
        ));
    }

    let n_trees = snap_read_u64(&mut r)? as usize;
    let n_taxa = snap_read_u64(&mut r)? as usize;
    let n_bip = snap_read_u64(&mut r)? as usize;

    let taxa_names = snap_read_strings(&mut r, n_taxa)?;
    let tree_names = snap_read_strings(&mut r, n_trees)?;

    let mut presence = vec![0u8; n_trees * n_bip];
    r.read_exact(&mut presence)?;

    Ok((tree_names, taxa_names, n_bip, presence))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Four-leaf trees with a single topology difference — known RF = 2
    const NWK_A: &str = "(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);";
    const NWK_B: &str = "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);";

    fn identity_map(taxa: &[&str]) -> HashMap<String, String> {
        taxa.iter()
            .map(|&t| (t.to_string(), t.to_string()))
            .collect()
    }

    #[test]
    fn parse_newicks_basic() {
        let newicks = vec![NWK_A.to_string(), NWK_B.to_string()];
        let map = identity_map(&["A", "B", "C", "D"]);
        let (snaps, leaf_names) = snapshot_newicks(&newicks, &[map], &[0, 0], false).unwrap();

        assert_eq!(snaps.len(), 2);
        assert_eq!(leaf_names, vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn parse_newicks_requires_two_trees() {
        let newicks = vec![NWK_A.to_string()];
        let map = identity_map(&["A", "B", "C", "D"]);
        let result = snapshot_newicks(&newicks, &[map], &[0], false);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("at least 2"));
    }

    #[test]
    fn parse_newicks_leaf_mismatch_returns_error() {
        let nwk_different = "(A:1.0,(B:1.0,(C:1.0,E:1.0):1.0):1.0);"; // E instead of D
        let newicks = vec![NWK_A.to_string(), nwk_different.to_string()];
        let map = identity_map(&["A", "B", "C", "D", "E"]);
        let result = snapshot_newicks(&newicks, &[map], &[0, 0], false);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("different leaf set"));
    }

    #[test]
    fn parse_newicks_map_index_out_of_bounds() {
        let newicks = vec![NWK_A.to_string(), NWK_B.to_string()];
        let map = identity_map(&["A", "B", "C", "D"]);
        let result = snapshot_newicks(&newicks, &[map], &[0, 99], false);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("out of bounds"));
    }

    #[test]
    fn read_and_snapshot_trees_basic() {
        let paths = vec!["tests/data/hiv1.trees".to_string()];
        let (names, snaps) = read_and_snapshot_trees(&paths, 0, 0, false, false).unwrap();

        assert!(!names.is_empty());
        assert_eq!(names.len(), snaps.len());
        assert!(
            names[0].starts_with("file0_"),
            "names should be prefixed with 'file0_', got: {}",
            names[0]
        );
    }

    #[test]
    fn read_and_snapshot_trees_burnin_reduces_count() {
        let paths = vec!["tests/data/hiv1.trees".to_string()];
        let (names_full, _) = read_and_snapshot_trees(&paths, 0, 0, false, false).unwrap();
        let burnin = names_full.len() / 2;
        let (names_burned, _) = read_and_snapshot_trees(&paths, burnin, 0, false, false).unwrap();

        assert!(
            names_burned.len() < names_full.len(),
            "burnin should reduce tree count: full={}, burned={}",
            names_full.len(),
            names_burned.len()
        );
    }
}

#[cfg(all(test, feature = "cli"))]
mod snap_tests {
    use super::*;
    use crate::distances::pairwise_rf_matrix;
    use flate2::Compression;

    const HIV_TREES: &str = "tests/data/hiv1.trees";

    /// Build snapshots + sorted taxa names from the hiv1 fixture.
    fn hiv_setup() -> (Vec<String>, Vec<String>, Vec<TreeSnapshot>) {
        let paths = vec![HIV_TREES.to_string()];
        let (tree_names, snaps) = read_and_snapshot_trees(&paths, 0, 0, false, false).unwrap();

        // Extract sorted leaf names from the first parsed tree.
        let (_, named_trees) = read_beast_trees(HIV_TREES, 0, 0, false);
        let first_tree = &named_trees[0].1;
        let mut taxa: Vec<String> = first_tree
            .get_leaves()
            .iter()
            .filter_map(|&id| first_tree.get(&id).ok()?.name.clone())
            .collect();
        taxa.sort_unstable();

        (tree_names, taxa, snaps)
    }

    #[test]
    fn snap_names_survive_roundtrip() {
        let (tree_names, taxa_names, snaps) = hiv_setup();
        let tmp = "/tmp/rt_names.snap";

        write_snap(tmp, &tree_names, &taxa_names, &snaps).unwrap();
        let (rt_trees, rt_taxa, n_bip, presence) = read_snap(tmp).unwrap();

        assert_eq!(rt_trees, tree_names, "tree names changed");
        assert_eq!(rt_taxa, taxa_names, "taxa names changed");
        assert_eq!(
            presence.len(),
            tree_names.len() * n_bip,
            "presence buffer wrong size"
        );
    }

    #[test]
    fn snap_dimensions_are_consistent() {
        let (tree_names, taxa_names, snaps) = hiv_setup();
        let tmp = "/tmp/rt_dims.snap";

        write_snap(tmp, &tree_names, &taxa_names, &snaps).unwrap();
        let (rt_trees, rt_taxa, n_bip, presence) = read_snap(tmp).unwrap();

        assert_eq!(rt_trees.len(), snaps.len(), "n_trees mismatch");
        assert_eq!(rt_taxa.len(), taxa_names.len(), "n_taxa mismatch");
        assert!(n_bip > 0, "n_bip should be positive");
        assert_eq!(presence.len(), snaps.len() * n_bip);
    }

    /// Core invariant from the README:
    /// `sum(presence[i] XOR presence[j]) == RF(tree_i, tree_j)`
    #[test]
    fn snap_presence_rf_matches_direct() {
        let (tree_names, taxa_names, snaps) = hiv_setup();
        let tmp = "/tmp/rt_rf.snap";

        write_snap(tmp, &tree_names, &taxa_names, &snaps).unwrap();
        let (_, _, n_bip, presence) = read_snap(tmp).unwrap();

        let direct = pairwise_rf_matrix(&snaps);
        let n = snaps.len();

        for i in 0..n {
            for j in 0..n {
                let rf_pres: usize = presence[i * n_bip..(i + 1) * n_bip]
                    .iter()
                    .zip(&presence[j * n_bip..(j + 1) * n_bip])
                    .filter(|(a, b)| a != b)
                    .count();
                assert_eq!(
                    rf_pres, direct[i][j],
                    "RF mismatch at [{i}][{j}]: presence={rf_pres} direct={}",
                    direct[i][j]
                );
            }
        }
    }

    #[test]
    fn snap_diagonal_is_zero() {
        let (tree_names, taxa_names, snaps) = hiv_setup();
        let tmp = "/tmp/rt_diag.snap";

        write_snap(tmp, &tree_names, &taxa_names, &snaps).unwrap();
        let (_, _, n_bip, presence) = read_snap(tmp).unwrap();

        for i in 0..snaps.len() {
            let rf_self: usize = presence[i * n_bip..(i + 1) * n_bip]
                .iter()
                .zip(&presence[i * n_bip..(i + 1) * n_bip])
                .filter(|(a, b)| a != b)
                .count();
            assert_eq!(rf_self, 0, "self-RF at tree {i} should be 0");
        }
    }

    #[test]
    fn snap_bad_magic_returns_error() {
        let tmp = "/tmp/rt_bad_magic.snap";
        let f = std::fs::File::create(tmp).unwrap();
        let mut gz = flate2::write::GzEncoder::new(f, Compression::default());
        // Write wrong magic + plausible header bytes
        gz.write_all(b"NOPE\x01\x00\x00\x00\x00\x00\x00\x00\x00")
            .unwrap();
        gz.finish().unwrap();

        let err = read_snap(tmp).unwrap_err();
        assert!(
            err.to_string().contains("invalid snap magic"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn snap_bad_version_returns_error() {
        let tmp = "/tmp/rt_bad_ver.snap";
        let f = std::fs::File::create(tmp).unwrap();
        let mut gz = flate2::write::GzEncoder::new(f, Compression::default());
        gz.write_all(b"SNAP\x99").unwrap(); // version 153 — unsupported
        gz.finish().unwrap();

        let err = read_snap(tmp).unwrap_err();
        assert!(
            err.to_string().contains("unsupported snap version"),
            "unexpected error: {err}"
        );
    }
}
