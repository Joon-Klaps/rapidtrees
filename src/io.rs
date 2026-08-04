use crate::snapshot::Snapshots;
use phylotree::tree::Tree;
use std::collections::HashMap;
use std::fs;

use std::path::Path;

#[cfg(feature = "cli")]
use crate::bitset::Bitset;
#[cfg(feature = "cli")]
use crate::snapshot::InternSnap;
#[cfg(feature = "cli")]
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
#[cfg(feature = "cli")]
use std::fs::File;
#[cfg(feature = "cli")]
use std::io;
#[cfg(feature = "cli")]
use std::io::{BufReader, BufWriter, Read, Write};

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

    load_beast_raw_str(
        &content,
        base_name,
        burnin_trees,
        burnin_states,
        use_real_taxa,
    )
}

/// String-based counterpart to [`load_beast_raw`].
///
/// Exists for callers with no filesystem: the browser build receives the file
/// contents through the File API and never has a path. Keeping both entry
/// points on one body means the NEXUS handling — translate blocks, burnin,
/// BEAST annotation stripping — cannot drift between the desktop and web
/// builds.
pub(crate) fn load_beast_raw_str(
    content: &str,
    base_name: &str,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> (HashMap<String, String>, Vec<(String, String)>) {
    let taxons = parse_taxon_block(content);
    let translate_map = if use_real_taxa {
        taxons
    } else {
        HashMap::new()
    };

    let tree_pairs: Vec<(String, String)> = collect_tree_blocks(content)
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

    let (presence, _) = snaps.build_presence_matrix();

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
        words_per_bitset: words,
        leaf_names: taxa_names,
    };

    Ok((tree_names, snaps))
}

// ── Private helpers ────────────────────────────────────────────────────────────
#[cfg(feature = "cli")]
fn snap_read_u64<R: io::Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(feature = "cli")]
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

/// Zero-copy counterpart to [`collect_tree_blocks`]: borrows each tree body out
/// of `content` instead of allocating a `String` per tree.
///
/// Used by [`load_beast_subsampled_str`], where most trees are about to be
/// discarded and allocating them first is pure waste. On a 252 MB BEAST file
/// with 4001 trees, materialising every body costs ~240 MB before subsampling
/// even runs.
fn collect_tree_block_refs(content: &str) -> Vec<(&str, &str)> {
    content
        .lines()
        .skip_while(|line| !line.to_ascii_uppercase().starts_with("TREE "))
        .take_while(|line| !line.trim().to_ascii_uppercase().starts_with("END;"))
        .filter_map(|line| {
            let mut parts = line.splitn(2, " = ");
            let header = parts.next()?.trim();
            let body = parts.next()?.trim();
            Some((header, body))
        })
        .collect()
}

/// Load a BEAST file, applying proportional burnin and even subsampling before
/// any per-tree allocation happens.
///
/// Returns `(translate_map, kept_pairs, total_trees, burnin_count)`.
///
/// * `burnin_percent` — share of the chain dropped from the front (1.0 = 1%).
///   Proportional rather than absolute, so one setting is meaningful across
///   chains of different lengths. Never drops the whole chain.
/// * `target_trees` — how many to keep afterwards, spaced evenly across what
///   remains. `0` keeps everything; a target above the chain length keeps
///   everything.
///
/// Even spacing rather than "the first N": a contiguous slice would sample one
/// narrow window of the chain and misrepresent how it mixed.
/// Keys carrying a log density in a BEAST tree header, in preference order.
///
/// A `.trees` line looks like
/// `tree STATE_0 [&lnP=-72411.03,joint=-72411.03] = [&R] (…);`
/// and which key is present depends on the BEAST version and what was logged.
/// The order here matches `callbacks/diagnostics.py`.
const DENSITY_KEYS: [&str; 5] = ["lnP", "lnL", "posterior", "joint", "loglikelihood"];

/// Pull `key=value` pairs out of the `[&…]` block in a tree header.
///
/// Only the header is scanned, never the newick body — branch annotations like
/// `[&rate=0.65]` appear thousands of times per line and none of them are a
/// tree-level statistic.
pub(crate) fn parse_header_annotations(header: &str) -> Vec<(String, String)> {
    let Some(open) = header.find("[&") else {
        return Vec::new();
    };
    let rest = &header[open + 2..];
    let Some(close) = rest.find(']') else {
        return Vec::new();
    };
    rest[..close]
        .split(',')
        .filter_map(|kv| {
            let (k, v) = kv.split_once('=')?;
            Some((k.trim().to_string(), v.trim().to_string()))
        })
        .collect()
}

/// The log density of a tree, taking the first key that is present.
fn header_log_density(header: &str) -> (f64, Option<String>) {
    let anns = parse_header_annotations(header);
    for key in DENSITY_KEYS {
        if let Some((_, v)) = anns.iter().find(|(k, _)| k.eq_ignore_ascii_case(key))
            && let Ok(parsed) = v.parse::<f64>()
        {
            return (parsed, Some(key.to_string()));
        }
    }
    (f64::NAN, None)
}

/// Everything one `.trees` file contributes to an analysis.
pub(crate) struct LoadedTrees {
    pub translate: HashMap<String, String>,
    pub names: Vec<String>,
    pub newicks: Vec<String>,
    /// Log density per kept tree; `NaN` where the file records none.
    pub log_density: Vec<f64>,
    /// Which annotation key supplied it, for axis labelling.
    pub density_field: Option<String>,
    /// MCMC state number per kept tree, for plotting against chain position.
    pub states: Vec<u64>,
    pub total: usize,
    pub burnin: usize,
}

pub(crate) fn load_beast_subsampled_str(
    content: &str,
    base_name: &str,
    burnin_percent: f64,
    target_trees: usize,
    use_real_taxa: bool,
) -> LoadedTrees {
    let translate_map = if use_real_taxa {
        parse_taxon_block(content)
    } else {
        HashMap::new()
    };

    let blocks = collect_tree_block_refs(content);
    let total = blocks.len();
    if total == 0 {
        return LoadedTrees {
            translate: translate_map,
            names: Vec::new(),
            newicks: Vec::new(),
            log_density: Vec::new(),
            density_field: None,
            states: Vec::new(),
            total: 0,
            burnin: 0,
        };
    }

    let pct = burnin_percent.clamp(0.0, 100.0);
    let burnin = (((total as f64) * pct / 100.0).floor() as usize).min(total - 1);
    let available = total - burnin;
    let target = if target_trees == 0 {
        available
    } else {
        target_trees.min(available)
    };

    // Only now, on the survivors, do we pay for stripping and allocation.
    let mut names = Vec::with_capacity(target);
    let mut newicks = Vec::with_capacity(target);
    let mut log_density = Vec::with_capacity(target);
    let mut states = Vec::with_capacity(target);
    let mut density_field = None;

    for i in 0..target {
        let idx = burnin + (i * available) / target;
        let (header, body) = blocks[idx];
        let (name, state) = extract_name_state(header);
        // Read the header's annotations before they are stripped: `lnP` and
        // friends live there, and discarding them was what left the browser
        // build with no log-density trace at all.
        let (density, field) = header_log_density(header);
        if density_field.is_none() {
            density_field = field;
        }
        names.push(format!("{base_name}_{name}"));
        newicks.push(strip_beast_annotations(body));
        log_density.push(density);
        states.push(state as u64);
    }

    LoadedTrees {
        translate: translate_map,
        names,
        newicks,
        log_density,
        density_field,
        states,
        total,
        burnin,
    }
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

#[cfg(test)]
mod load_tests {
    use super::*;

    fn hiv2_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hiv2.trees")
    }

    // ── strip_beast_annotations ───────────────────────────────────────────────

    #[test]
    fn test_strip_annotations_removes_bracketed_content() {
        let input = "(A[&rate=0.5]:1.0,B[&rate=0.3]:2.0):0.0;";
        assert_eq!(strip_beast_annotations(input), "(A:1.0,B:2.0):0.0;");
    }

    #[test]
    fn test_strip_annotations_passthrough_plain_newick() {
        let input = "((A:1,B:1):1,(C:1,D:1):1);";
        assert_eq!(strip_beast_annotations(input), input);
    }

    #[test]
    fn test_strip_annotations_removes_leading_r_flag() {
        // [&R] marks a rooted tree in BEAST output
        let input = "[&R] ((A:1,B:1):1,(C:1,D:1):1);";
        assert_eq!(
            strip_beast_annotations(input),
            " ((A:1,B:1):1,(C:1,D:1):1);"
        );
    }

    // ── extract_name_state ────────────────────────────────────────────────────

    #[test]
    fn test_extract_name_state_standard() {
        let (name, state) = extract_name_state("tree STATE_10000 [&lnP=-123.4]");
        assert_eq!(name, "STATE_10000");
        assert_eq!(state, 10000);
    }

    #[test]
    fn test_extract_name_state_zero() {
        let (name, state) = extract_name_state("tree STATE_0 [&lnP=-123.4]");
        assert_eq!(name, "STATE_0");
        assert_eq!(state, 0);
    }

    #[test]
    fn test_extract_name_state_no_state_keyword_returns_empty() {
        let (name, state) = extract_name_state("tree my_tree");
        assert_eq!(name, "");
        assert_eq!(state, 0);
    }

    // ── parse_taxon_block ─────────────────────────────────────────────────────

    #[test]
    fn test_parse_taxon_block_basic() {
        let content =
            "Begin trees;\n\tTranslate\n\t\t1 'Alpha',\n\t\t2 'Beta',\n\t\t3 'Gamma'\n\t;\nEnd;\n";
        let map = parse_taxon_block(content);
        assert_eq!(map.len(), 3);
        assert_eq!(map.get("1").map(String::as_str), Some("Alpha"));
        assert_eq!(map.get("2").map(String::as_str), Some("Beta"));
        assert_eq!(map.get("3").map(String::as_str), Some("Gamma"));
    }

    #[test]
    fn test_parse_taxon_block_empty_when_no_translate() {
        let content = "#NEXUS\nBegin trees;\ntree t1 = (A:1,B:1);\nEnd;\n";
        assert!(parse_taxon_block(content).is_empty());
    }

    // ── collect_tree_blocks ───────────────────────────────────────────────────

    #[test]
    fn test_collect_tree_blocks_count() {
        let content = "Begin trees;\ntree t1 = (A:1,B:1);\ntree t2 = (A:2,B:2);\nEnd;\n";
        let blocks = collect_tree_blocks(content);
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn test_collect_tree_blocks_header_and_body() {
        let content = "Begin trees;\ntree STATE_0 = (A:1,B:1);\nEnd;\n";
        let blocks = collect_tree_blocks(content);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].header, "tree STATE_0");
        assert_eq!(blocks[0].body, "(A:1,B:1);");
    }

    // ── rename_leaf_nodes ─────────────────────────────────────────────────────

    #[test]
    fn test_rename_leaf_nodes_applies_translate() {
        let translate: HashMap<String, String> = [
            ("1".to_string(), "Alpha".to_string()),
            ("2".to_string(), "Beta".to_string()),
        ]
        .into();
        let mut tree = phylotree::tree::Tree::from_newick("(1:1.0,2:1.0);").unwrap();
        rename_leaf_nodes(&mut tree, &translate);
        let leaf_names: Vec<_> = tree
            .get_leaves()
            .iter()
            .filter_map(|id| tree.get(id).ok()?.name.clone())
            .collect();
        assert!(leaf_names.contains(&"Alpha".to_string()));
        assert!(leaf_names.contains(&"Beta".to_string()));
    }

    #[test]
    fn test_rename_leaf_nodes_noop_on_empty_map() {
        let mut tree = phylotree::tree::Tree::from_newick("(A:1.0,B:1.0);").unwrap();
        rename_leaf_nodes(&mut tree, &HashMap::new());
        let leaf_names: Vec<_> = tree
            .get_leaves()
            .iter()
            .filter_map(|id| tree.get(id).ok()?.name.clone())
            .collect();
        assert!(leaf_names.contains(&"A".to_string()));
        assert!(leaf_names.contains(&"B".to_string()));
    }

    // ── load_beast_raw ────────────────────────────────────────────────────────

    #[test]
    fn test_load_beast_raw_no_burnin_returns_all_trees() {
        let (_, pairs) = load_beast_raw(hiv2_path(), 0, 0, false);
        assert_eq!(pairs.len(), 21);
    }

    #[test]
    fn test_load_beast_raw_burnin_by_tree_count() {
        let (_, pairs) = load_beast_raw(hiv2_path(), 5, 0, false);
        assert_eq!(pairs.len(), 16);
    }

    #[test]
    fn test_load_beast_raw_burnin_by_state() {
        // hiv2.trees: STATE_0 .. STATE_200000 step 10000 (21 trees)
        // keep state > 50000 → STATE_60000 .. STATE_200000 = 15 trees
        let (_, pairs) = load_beast_raw(hiv2_path(), 0, 50000, false);
        assert_eq!(pairs.len(), 15);
    }

    #[test]
    fn test_load_beast_raw_use_real_taxa_populates_map() {
        let (translate, _) = load_beast_raw(hiv2_path(), 0, 0, true);
        assert!(!translate.is_empty());
        assert_eq!(
            translate.get("1").map(String::as_str),
            Some("1959.M.CD.59.ZR59")
        );
    }

    #[test]
    fn test_load_beast_raw_use_real_taxa_false_empty_map() {
        let (translate, _) = load_beast_raw(hiv2_path(), 0, 0, false);
        assert!(translate.is_empty());
    }

    #[test]
    fn test_load_beast_raw_strips_annotations() {
        let (_, pairs) = load_beast_raw(hiv2_path(), 0, 0, false);
        for (_, newick) in &pairs {
            assert!(
                !newick.contains("[&"),
                "newick must not contain BEAST annotations after stripping"
            );
        }
    }

    #[test]
    fn test_load_beast_raw_tree_names_contain_filename_stem() {
        let (_, pairs) = load_beast_raw(hiv2_path(), 0, 0, false);
        for (name, _) in &pairs {
            assert!(
                name.starts_with("hiv2_"),
                "expected name to start with 'hiv2_', got '{name}'"
            );
        }
    }

    #[test]
    fn test_load_beast_raw_nonexistent_returns_empty() {
        let (translate, pairs) = load_beast_raw("nonexistent.trees", 0, 0, false);
        assert!(pairs.is_empty());
        assert!(translate.is_empty());
    }

    // ── load_beast_trees ──────────────────────────────────────────────────────

    #[test]
    fn test_load_beast_trees_returns_correct_count() {
        let (names, snaps) = load_beast_trees(hiv2_path(), 0, 0, false, false);
        assert_eq!(names.len(), 21);
        assert_eq!(snaps.len(), 21);
        assert!(!snaps.leaf_names.is_empty());
    }

    #[test]
    fn test_load_beast_trees_burnin_reduces_count() {
        let (names, snaps) = load_beast_trees(hiv2_path(), 5, 0, false, false);
        assert_eq!(names.len(), 16);
        assert_eq!(snaps.len(), 16);
    }

    #[test]
    fn test_load_beast_trees_nonexistent_returns_empty() {
        let (names, snaps) = load_beast_trees("nonexistent.trees", 0, 0, false, false);
        assert!(names.is_empty());
        assert_eq!(snaps.len(), 0);
    }
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
        let rf_before = snaps.pairwise_rf(None);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_snap(tmp.path(), &tree_names, &snaps).unwrap();

        let (loaded_names, loaded_snaps) = load_snapshots(tmp.path()).unwrap();

        assert_eq!(
            loaded_names, tree_names,
            "tree names must survive roundtrip"
        );
        assert_eq!(loaded_snaps.snapshots.len(), snaps.snapshots.len());

        let rf_after = loaded_snaps.pairwise_rf(None);
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

#[cfg(test)]
mod annotation_tests {
    use super::*;

    #[test]
    fn reads_lnp_from_a_real_beast_header() {
        let header = "tree STATE_0 [&lnP=-72411.03428232882,joint=-72411.03428232882]";
        let (v, field) = header_log_density(header);
        assert!((v + 72411.03428232882).abs() < 1e-9, "got {v}");
        assert_eq!(field.as_deref(), Some("lnP"));
    }

    #[test]
    fn prefers_lnp_over_joint() {
        let header = "tree STATE_1 [&joint=-5.0,lnP=-1.0]";
        assert_eq!(header_log_density(header).1.as_deref(), Some("lnP"));
    }

    #[test]
    fn falls_back_through_the_key_list() {
        let (v, f) = header_log_density("tree STATE_2 [&posterior=-3.5]");
        assert_eq!(f.as_deref(), Some("posterior"));
        assert!((v + 3.5).abs() < 1e-12);
    }

    #[test]
    fn missing_annotations_give_nan() {
        assert!(header_log_density("tree STATE_3").0.is_nan());
        assert!(header_log_density("tree STATE_4 [&rate=0.5]").0.is_nan());
    }

    #[test]
    fn only_the_header_is_scanned() {
        // A branch annotation in the body must not be mistaken for a
        // tree-level statistic; callers pass the header alone.
        let anns = parse_header_annotations("tree STATE_5 [&lnP=-2.0]");
        assert_eq!(anns.len(), 1);
        assert_eq!(anns[0].0, "lnP");
    }
}
