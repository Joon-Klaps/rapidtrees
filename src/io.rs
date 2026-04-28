#[cfg(feature = "cli")]
use crate::bitset::Bitset;
use crate::snapshot::TreeSnapshot;
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use phylotree::tree::Tree;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Shared return type for functions that load trees: a taxon map and a list of
/// (tree name, snapshot) pairs.
pub type LoadedTrees = (HashMap<String, String>, Vec<(String, TreeSnapshot)>);

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

/// Parse newick strings, rename leaves, validate leaf-set consistency, and
/// create [`TreeSnapshot`]s. Returns `(snapshots, sorted_leaf_names)`.
///
/// All entries must share the same leaf set; the first tree establishes the
/// reference and subsequent trees are validated against it.
pub fn parse_and_snapshot_newicks<'a>(
    entries: impl Iterator<Item = (&'a str, &'a HashMap<String, String>)>,
    rooted: bool,
) -> Result<(Vec<TreeSnapshot>, Vec<String>), String> {
    let mut snapshots = Vec::new();
    let mut reference_leaves: Option<HashSet<String>> = None;
    let mut sorted_leaf_names = Vec::new();

    // iter & enum & map & collect is a bit cleaner than a for loop with manual indexing and error messages
    for (i, (newick, translate)) in entries.enumerate() {
        let clean = strip_beast_annotations(newick);
        let mut tree = Tree::from_newick(&clean)
            .map_err(|e| format!("Failed to parse newick at index {i}: {e}"))?;
        rename_leaf_nodes(&mut tree, translate);

        let leaves: HashSet<String> = tree
            .get_leaves()
            .iter()
            .filter_map(|&id| tree.get(&id).ok()?.name.clone())
            .collect();

        if let Some(ref_leaves) = &reference_leaves {
            if leaves != *ref_leaves {
                return Err(format!(
                    "Tree {i} has different leaf set than tree 0. All trees must have the same taxa."
                ));
            }
        } else {
            let mut sorted: Vec<String> = leaves.iter().cloned().collect();
            sorted.sort_unstable();
            sorted_leaf_names = sorted;
            reference_leaves = Some(leaves);
        }

        snapshots.push(
            TreeSnapshot::from_tree(&tree, rooted)
                .map_err(|e| format!("Failed to create tree snapshot at index {i}: {e}"))?,
        );
    }

    Ok((snapshots, sorted_leaf_names))
}

pub fn load_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> LoadedTrees {
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
    let empty_map = HashMap::new();
    let translate = if use_real_taxa { &taxons } else { &empty_map };

    let (names, newicks): (Vec<String>, Vec<String>) = collect_tree_blocks(&content)
        .into_iter()
        .enumerate()
        .map(|(idx, tree)| {
            let (name, state) = extract_name_state(tree.header);
            (idx, tree, state, format!("{base_name}_{name}"))
        })
        .filter(|(idx, _tree, state, _name)| {
            (burnin_trees == 0 && burnin_states == 0)
                || (burnin_trees > 0 && *idx >= burnin_trees)
                || (burnin_states > 0 && *state > burnin_states)
        })
        .map(|(_, tree, _, name)| (name, strip_beast_annotations(&tree.body)))
        .unzip();

    let entries = newicks.iter().map(|n| (n.as_str(), translate));
    match parse_and_snapshot_newicks(entries, rooted) {
        Ok((snapshots, _)) => (taxons, names.into_iter().zip(snapshots).collect()),
        Err(e) => {
            eprintln!("Failed to parse trees in {:?}: {e}", path.as_ref());
            (taxons, Vec::new())
        }
    }
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
    let mut all_bips: Vec<_> = snapshots
        .iter()
        .flat_map(|s| s.parts.iter().cloned())
        .collect();
    all_bips.sort_unstable();
    all_bips.dedup();

    let n_trees = snapshots.len();
    let n_taxa = taxa_names.len();
    let n_bip = all_bips.len();
    let words = snapshots.first().map_or(1, |s| s.words);

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

    let mut w = BufWriter::new(GzEncoder::new(File::create(path)?, Compression::default()));

    // Header
    w.write_all(b"SNAP")?;
    w.write_all(&[2u8])?; // version 2
    w.write_all(&(n_trees as u64).to_le_bytes())?;
    w.write_all(&(n_taxa as u64).to_le_bytes())?;
    w.write_all(&(n_bip as u64).to_le_bytes())?;
    w.write_all(&(words as u64).to_le_bytes())?; // new in v2

    // Names
    for name in taxa_names.iter().chain(tree_names.iter()) {
        let b = name.as_bytes();
        w.write_all(&(b.len() as u32).to_le_bytes())?;
        w.write_all(b)?;
    }

    // Bipartition bitsets — actual leaf content, n_bip × words × 8 bytes
    for bip in &all_bips {
        for word in &bip.0 {
            w.write_all(&word.to_le_bytes())?;
        }
    }

    // Presence matrix
    w.write_all(&presence)?;
    w.flush()?;
    Ok(())
}

/// Read a `.snap` file produced by [`write_snap`].
///
/// Returns `(tree_names, taxa_names, snapshots)`.
///
/// Snapshots are reconstructed from the on-disk presence matrix so callers can
/// feed them directly into distance functions like [`crate::distances::pairwise_rf_matrix`].
#[cfg(feature = "cli")]
pub fn load_snapshots<P: AsRef<Path>>(path: P) -> io::Result<LoadedTrees> {
    let file = File::open(path.as_ref())?;
    let mut r = BufReader::new(GzDecoder::new(file));

    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != b"SNAP" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid snap magic: expected SNAP, got {:?}", &magic),
        ));
    }

    let mut version = [0u8; 1];
    r.read_exact(&mut version)?;
    if version[0] != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported snap version: {}", version[0]),
        ));
    }

    let n_trees = snap_read_u64(&mut r)? as usize;
    let n_taxa = snap_read_u64(&mut r)? as usize;
    let n_bip = snap_read_u64(&mut r)? as usize;
    let words = snap_read_u64(&mut r)? as usize;

    let taxa_names = snap_read_strings(&mut r, n_taxa)?;
    let tree_names = snap_read_strings(&mut r, n_trees)?;

    // Number of bytes to read = n_trees * n_bip
    // Read the actual bipartition bitsets
    let all_bips: Vec<Bitset> = (0..n_bip)
        .map(|_| {
            (0..words)
                .map(|_| snap_read_u64(&mut r))
                .collect::<io::Result<Vec<u64>>>()
                .map(Bitset)
        })
        .collect::<io::Result<_>>()?;

    // Read presence matrix and reconstruct snapshots with real bitsets
    let mut presence = vec![0u8; n_trees * n_bip];
    r.read_exact(&mut presence)?;

    let snapshots: Vec<TreeSnapshot> = presence
        .chunks_exact(n_bip)
        .take(n_trees)
        .map(|row| {
            let parts: Vec<Bitset> = (0..n_bip)
                .filter(|&i| row[i] != 0)
                .map(|i| all_bips[i].clone())
                .collect();
            TreeSnapshot {
                lengths: vec![0.0; parts.len()],
                parts,
                root_children: Vec::new(),
                words,
                num_leaves: n_taxa,
                rooted: false,
            }
        })
        .collect();

    let taxa_map: HashMap<String, String> = taxa_names
        .into_iter()
        .enumerate()
        .map(|(i, name)| (i.to_string(), name))
        .collect();

    let trees: Vec<(String, TreeSnapshot)> = tree_names.into_iter().zip(snapshots).collect();

    Ok((taxa_map, trees))
}

// ── Snap helpers ──────────────────────────────────────────────────────────────

fn invalid_data(msg: impl ToString) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

fn snap_read_u64<R: io::Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn snap_read_strings<R: io::Read>(r: &mut R, n: usize) -> io::Result<Vec<String>> {
    (0..n)
        .map(|_| {
            let mut buf = [0u8; 4];
            r.read_exact(&mut buf)?;
            let mut bytes = vec![0u8; u32::from_le_bytes(buf) as usize];
            r.read_exact(&mut bytes)?;
            String::from_utf8(bytes).map_err(invalid_data) // ← passes directly now
        })
        .collect()
}

#[cfg(all(test, feature = "cli"))]
mod tests {
    use super::*;
    use crate::distances::pairwise_rf_matrix;
    use crate::snapshot::TreeSnapshot;
    use phylotree::tree::Tree as PhyloTree;

    fn make_snapshots() -> (Vec<String>, Vec<String>, Vec<TreeSnapshot>) {
        let trees = [
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,C:1):1,(B:1,D:1):1);",
            "((A:1,D:1):1,(B:1,C:1):1);",
        ];
        let snapshots: Vec<TreeSnapshot> = trees
            .iter()
            .map(|nwk| {
                TreeSnapshot::from_tree(&PhyloTree::from_newick(nwk).unwrap(), false).unwrap()
            })
            .collect();
        let tree_names = vec!["t1".to_string(), "t2".to_string(), "t3".to_string()];
        let mut taxa_names = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        taxa_names.sort_unstable();
        (tree_names, taxa_names, snapshots)
    }

    #[test]
    fn test_snap_roundtrip_preserves_rf_distances() {
        let (tree_names, taxa_names, snapshots) = make_snapshots();
        let rf_before = pairwise_rf_matrix(&snapshots);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_snap(tmp.path(), &tree_names, &taxa_names, &snapshots).unwrap();

        let (_taxon_map, loaded_trees) = load_snapshots(tmp.path()).unwrap();
        let (loaded_names, loaded_snaps): (Vec<_>, Vec<_>) = loaded_trees.into_iter().unzip();

        assert_eq!(
            loaded_names, tree_names,
            "tree names must survive roundtrip"
        );
        assert_eq!(loaded_snaps.len(), snapshots.len());

        let rf_after = pairwise_rf_matrix(&loaded_snaps);
        assert_eq!(
            rf_before, rf_after,
            "RF distances must be identical after snap roundtrip"
        );
    }

    #[test]
    fn test_snap_wrong_magic_returns_error() {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(b"NOPE").unwrap();
        assert!(load_snapshots(tmp.path()).is_err());
    }
}
