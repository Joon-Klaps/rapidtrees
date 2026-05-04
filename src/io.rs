use crate::bitset::Bitset;
use crate::snapshot::{InternSnap, Snapshots};
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use phylotree::tree::Tree;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Strip BEAST `[&...]` annotations from a Newick string.
pub(crate) fn strip_beast_annotations(newick: &str) -> String {
    let mut result = String::with_capacity(newick.len());
    let mut in_annotation = false;
    let mut chars = newick.chars().peekable();

    while let Some(ch) = chars.next() {
        if ch == '[' && chars.peek() == Some(&'&') {
            in_annotation = true;
        } else if ch == ']' && in_annotation {
            in_annotation = false;
        } else if !in_annotation {
            result.push(ch);
        }
    }

    result
}

/// Rename leaf nodes in a tree according to a translate map.
///
/// Used to apply BEAST translate blocks (numeric ID → taxon name). A no-op if
/// `translate` is empty.
pub(crate) fn rename_leaf_nodes(phylo_tree: &mut Tree, translate: &HashMap<String, String>) {
    if translate.is_empty() {
        return;
    }
    for leaf_id in phylo_tree.get_leaves() {
        if let Ok(node) = phylo_tree.get_mut(&leaf_id) {
            node.name = node.name.as_ref().and_then(|n| translate.get(n).cloned());
        }
    }
}

/// Load raw tree data from a BEAST `.trees` file without parsing.
///
/// Returns `(translate_map, [(tree_name, stripped_newick)])`. Burnin filtering
/// and annotation stripping are applied but no Newick parsing is done.
pub(crate) fn load_beast_raw<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> (HashMap<String, String>, Vec<(String, String)>) {
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
    let translate_map = if use_real_taxa {
        taxons
    } else {
        HashMap::new()
    };

    let tree_pairs: Vec<(String, String)> = collect_tree_blocks(&content)
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
        .collect();

    (translate_map, tree_pairs)
}

/// Load and parse all trees from a BEAST `.trees` file.
///
/// Returns `(tree_names, Snapshots)`. On any error, prints to stderr and
/// returns an empty `Snapshots`.
pub fn load_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
) -> (Vec<String>, Snapshots) {
    let (translate_map, tree_pairs) =
        load_beast_raw(&path, burnin_trees, burnin_states, use_real_taxa);
    let (names, newicks): (Vec<String>, Vec<String>) = tree_pairs.into_iter().unzip();
    let entries = newicks.iter().map(|n| (n.as_str(), &translate_map));
    match Snapshots::from_newick_iter(entries, rooted) {
        Ok(snaps) => (names, snaps),
        Err(e) => {
            eprintln!("Failed to parse trees in {:?}: {e}", path.as_ref());
            (Vec::new(), Snapshots::from_newicks(&[], rooted).unwrap())
        }
    }
}

/// Write a labeled square matrix as TSV to a file or stdout.
#[cfg(feature = "cli")]
pub fn write_matrix_tsv<P: AsRef<Path>, T: std::fmt::Display>(
    path: P,
    names: &[String],
    mat: &[T],
    n_trees: usize,
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

    write!(&mut out, "\t")?;
    for (k, name) in names.iter().enumerate() {
        if k > 0 {
            write!(&mut out, "\t")?;
        }
        write!(&mut out, "{}", name)?;
    }
    writeln!(&mut out)?;

    for (i, row) in mat.chunks(n_trees).enumerate() {
        write!(&mut out, "{}", names[i])?;
        for val in row {
            write!(&mut out, "\t{}", val)?;
        }
        writeln!(&mut out)?;
    }

    out.flush()?;
    Ok(())
}

// ── Snap file I/O ─────────────────────────────────────────────────────────────
//
// The `.snap` format stores only bipartition presence (no branch lengths).
// Branch lengths are zeroed on load — snap files are RF-only.

/// Write tree snapshots to a gzip-compressed binary `.snap` file.
///
/// # File layout (inside gzip stream)
/// ```text
/// HEADER     4 B  magic "SNAP"
///            1 B  version (2)
///            8 B  n_trees  u64 LE
///            8 B  n_taxa   u64 LE
///            8 B  n_bip    u64 LE
///            8 B  words    u64 LE
/// NAMES      for each taxon then each tree: 4 B len u32 LE + N B UTF-8
/// BIPARTS    n_bip × words × 8 bytes (sorted ascending by Bitset)
/// PRESENCE   n_trees × n_bip bytes, row-major uint8
/// ```
///
/// **Note:** branch lengths are not stored. A snap file loaded back via
/// [`load_snapshots`] will have zero branch lengths (RF-only).
#[cfg(feature = "cli")]
pub fn write_snap<P: AsRef<Path>>(
    path: P,
    tree_names: &[String],
    snaps: &Snapshots,
) -> io::Result<()> {
    let mut bip_order: Vec<usize> = (0..snaps.bipartitions.len()).collect();
    bip_order.sort_unstable_by(|&a, &b| snaps.bipartitions[a].cmp(&snaps.bipartitions[b]));

    let n_trees = snaps.snapshots.len();
    let n_taxa = snaps.leaf_names.len();
    let n_bip = snaps.bipartitions.len();
    let words = snaps.words_per_bitset.max(1);

    let presence = snaps.build_presence_matrix();

    let mut w = BufWriter::new(GzEncoder::new(File::create(path)?, Compression::default()));
    w.write_all(b"SNAP")?;
    w.write_all(&[2u8])?;
    w.write_all(&(n_trees as u64).to_le_bytes())?;
    w.write_all(&(n_taxa as u64).to_le_bytes())?;
    w.write_all(&(n_bip as u64).to_le_bytes())?;
    w.write_all(&(words as u64).to_le_bytes())?;

    for name in snaps.leaf_names.iter().chain(tree_names.iter()) {
        let b = name.as_bytes();
        w.write_all(&(b.len() as u32).to_le_bytes())?;
        w.write_all(b)?;
    }

    for &orig_id in &bip_order {
        for word in &snaps.bipartitions[orig_id].0 {
            w.write_all(&word.to_le_bytes())?;
        }
    }

    w.write_all(&presence)?;
    w.flush()?;
    Ok(())
}

/// Read a `.snap` file produced by [`write_snap`].
///
/// Returns `(tree_names, Snapshots)`.
///
/// **Note:** branch lengths are zeroed on load — snap files store only split
/// presence. Use the returned `Snapshots` for RF distances only.
#[cfg(feature = "cli")]
pub fn load_snapshots<P: AsRef<Path>>(path: P) -> io::Result<(Vec<String>, Snapshots)> {
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

    let all_bips: Vec<Bitset> = (0..n_bip)
        .map(|_| {
            (0..words)
                .map(|_| snap_read_u64(&mut r))
                .collect::<io::Result<Vec<u64>>>()
                .map(Bitset)
        })
        .collect::<io::Result<_>>()?;

    let mut presence = vec![0u8; n_trees * n_bip];
    r.read_exact(&mut presence)?;

    let bipartition_index: FxHashMap<Bitset, u32> = all_bips
        .iter()
        .enumerate()
        .map(|(i, bip)| (bip.clone(), i as u32))
        .collect();

    let interned_snaps: Vec<InternSnap> = presence
        .chunks_exact(n_bip)
        .take(n_trees)
        .map(|row| {
            let split_ids: Vec<u32> = (0..n_bip)
                .filter(|&i| row[i] != 0)
                .map(|i| i as u32)
                .collect();
            let lengths = vec![0.0f64; split_ids.len()];
            InternSnap { split_ids, lengths }
        })
        .collect();

    let snaps = Snapshots {
        snapshots: interned_snaps,
        bipartitions: all_bips,
        bipartition_index,
        words_per_bitset: words,
        leaf_names: taxa_names,
    };

    Ok((tree_names, snaps))
}

// ── Private helpers ────────────────────────────────────────────────────────────

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
            String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .collect()
}

fn extract_name_state(header: &str) -> (String, usize) {
    let upper = header.to_ascii_uppercase();
    if let Some(state_pos) = upper.find("STATE_")
        && let Some((_, rest)) = header.split_once(' ')
    {
        let tree_name = rest.split_whitespace().next().unwrap_or("").to_string();
        let digits = header[state_pos + 6..]
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect::<String>();
        if let Ok(num) = digits.parse::<usize>() {
            return (tree_name, num);
        }
    }
    (String::new(), 0)
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
        .filter_map(|line| {
            let line = line.trim().trim_end_matches(',');
            let mut parts = line.split_whitespace();
            let id = parts.next()?.to_string();
            let label = parts.next()?.trim_matches('\'').to_string();
            Some((id, label))
        })
        .collect::<HashMap<_, _>>()
}

#[cfg(all(test, feature = "cli"))]
mod tests {
    use super::*;
    use crate::snapshot::Snapshots;

    fn make_snapshots() -> (Vec<String>, Snapshots) {
        let trees = [
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,C:1):1,(B:1,D:1):1);",
            "((A:1,D:1):1,(B:1,C:1):1);",
        ];
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();
        let tree_names = vec!["t1".to_string(), "t2".to_string(), "t3".to_string()];
        (tree_names, snaps)
    }

    #[test]
    fn test_snap_roundtrip_preserves_rf_distances() {
        let (tree_names, snaps) = make_snapshots();
        let rf_before = snaps.pairwise_rf();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_snap(tmp.path(), &tree_names, &snaps).unwrap();

        let (loaded_names, loaded_snaps) = load_snapshots(tmp.path()).unwrap();

        assert_eq!(
            loaded_names, tree_names,
            "tree names must survive roundtrip"
        );
        assert_eq!(loaded_snaps.snapshots.len(), snaps.snapshots.len());

        let rf_after = loaded_snaps.pairwise_rf();
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
