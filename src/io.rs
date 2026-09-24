use crate::snapshot::{Retain, Snapshots};
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs;

use std::path::Path;

#[cfg(feature = "cli")]
use flate2::{Compression, write::GzEncoder};
#[cfg(feature = "cli")]
use std::io;
#[cfg(feature = "cli")]
use std::io::Write;

/// Strip BEAST `[&...]` annotations from a Newick string.
pub fn strip_beast_annotations(newick: &str) -> String {
    strip_annotations(newick).into_owned()
}

/// [`strip_beast_annotations`] without the copy when there is nothing to strip,
/// which is every tree in a plain Newick file.
pub(crate) fn strip_annotations(newick: &str) -> Cow<'_, str> {
    if !newick.contains("[&") {
        return Cow::Borrowed(newick);
    }
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

    Cow::Owned(result)
}

/// Tree-file formats [`load_beast_trees`] can read.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TreeFormat {
    /// NEXUS: optional `#NEXUS` header and `TRANSLATE` block, one
    /// `tree <name> = <newick>` line per tree.
    Nexus,
    /// Plain Newick: `;`-terminated trees and nothing else — no tree names.
    Newick,
}

/// Sniff whether `content` is a NEXUS trees file or a plain Newick file.
///
/// NEXUS is anything opening with `#NEXUS`, or holding at least one
/// `tree ... = ...` line that the NEXUS tree-line scanner recognises; everything else
/// reads as Newick. Keeping the header check means a NEXUS file whose trees
/// block is empty still reports as NEXUS, so its loader can say so instead of
/// handing the content to the Newick parser.
#[must_use]
pub fn detect_format(content: &str) -> TreeFormat {
    let nexus_header = content
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty())
        .is_some_and(|line| starts_with_ci(line, "#nexus"));

    if nexus_header || nexus_tree_lines(content).next().is_some() {
        TreeFormat::Nexus
    } else {
        TreeFormat::Newick
    }
}

/// One tree as read from a file, before burn-in filtering.
struct RawTree {
    /// Display name, already prefixed with the file stem.
    name: String,
    /// MCMC state it was sampled at; `0` when the format records none.
    state: usize,
    /// Newick string, annotations already stripped.
    newick: String,
}

/// Load raw tree data from a NEXUS or plain-Newick tree file without parsing.
///
/// The format is sniffed with [`detect_format`]. Returns
/// `(translate_map, [(tree_name, stripped_newick)])`; burnin filtering and
/// annotation stripping are applied but no Newick parsing is done.
///
/// Plain Newick files name nothing, so their trees are called
/// `<file stem>_line<n>` after the 1-based line each tree starts on. They carry
/// no `STATE_` labels either, so `burnin_states` does not apply, and they have
/// no `TRANSLATE` block, so `use_real_taxa` leaves the returned map empty.
pub(crate) fn load_beast_raw<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> (HashMap<String, String>, Vec<(String, String)>) {
    let path = path.as_ref();
    let content = match fs::read_to_string(path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Failed to read {path:?}: {e}");
            return (HashMap::new(), Vec::new());
        }
    };

    let base_name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    let format = detect_format(&content);

    // Newick trees carry no `STATE_` labels, so a state threshold cannot apply.
    let burnin_states = if format == TreeFormat::Newick && burnin_states > 0 {
        eprintln!("Ignoring burn-in by state for {path:?}: Newick trees carry no `STATE_` labels");
        0
    } else {
        burnin_states
    };

    let (translate_map, trees) = match format {
        TreeFormat::Nexus => {
            let translate = if use_real_taxa {
                parse_taxon_block(&content)
            } else {
                HashMap::new()
            };
            (translate, collect_nexus_trees(&content, base_name))
        }
        TreeFormat::Newick => (HashMap::new(), collect_newick_trees(&content, base_name)),
    };

    if trees.is_empty() && !content.trim().is_empty() {
        let missing = match format {
            TreeFormat::Nexus => "no NEXUS `tree ... = ...` lines",
            TreeFormat::Newick => "no `;`-terminated Newick trees",
        };
        eprintln!("No trees found in {path:?}: {missing}");
    }

    let tree_pairs = trees
        .into_iter()
        .enumerate()
        .filter(|(idx, tree)| keep_tree(*idx, tree.state, burnin_trees, burnin_states))
        .map(|(_, tree)| (tree.name, tree.newick))
        .collect();

    (translate_map, tree_pairs)
}

/// Burn-in predicate shared by both formats: keep tree `idx` (0-based), sampled
/// at MCMC `state`, unless a threshold excludes it.
#[inline]
fn keep_tree(idx: usize, state: usize, burnin_trees: usize, burnin_states: usize) -> bool {
    (burnin_trees == 0 && burnin_states == 0)
        || (burnin_trees > 0 && idx >= burnin_trees)
        || (burnin_states > 0 && state > burnin_states)
}

/// Load and parse all trees from a NEXUS or plain-Newick tree file.
///
/// The format is detected per file with [`detect_format`]: NEXUS trees keep
/// their `tree` header name, plain Newick trees are named `<file stem>_line<n>`
/// after the line they start on. Returns `(tree_names, Snapshots)`. On any
/// error, prints to stderr and returns an empty `Snapshots`.
///
/// `retain` says what to build besides the split IDs; pass `None` to keep
/// everything.
///
/// ```no_run
/// # use rapidtrees::{io::load_beast_trees, Retain};
/// let (names, snaps) = load_beast_trees("posterior.trees", 0, 0, false, false, None);
/// let (names, snaps) =
///     load_beast_trees("posterior.trees", 0, 0, false, false, Retain::for_distances(false));
/// ```
pub fn load_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
    retain: impl Into<Option<Retain>>,
) -> (Vec<String>, Snapshots) {
    let retain = retain.into().unwrap_or_else(Retain::everything);
    let (translate_map, tree_pairs) =
        load_beast_raw(&path, burnin_trees, burnin_states, use_real_taxa);
    let (names, newicks): (Vec<String>, Vec<String>) = tree_pairs.into_iter().unzip();
    let entries = newicks.iter().map(|n| (n.as_str(), &translate_map));
    match Snapshots::from_newick_iter_opts(entries, rooted, retain) {
        Ok(snaps) => (names, snaps),
        Err(e) => {
            eprintln!("Failed to parse trees in {:?}: {e}", path.as_ref());
            (Vec::new(), Snapshots::from_newicks(&[], rooted).unwrap())
        }
    }
}

/// Write a labeled square matrix as TSV to a file or stdout.
#[cfg(feature = "cli")]
pub fn write_matrix_tsv<P: AsRef<Path>, T: std::fmt::Display + Sync>(
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

    let rows_per_block = (BLOCK_BYTES / (n_trees.max(1) * CELL_BYTES)).max(1);
    write_rows(&mut out, names, mat, n_trees, rows_per_block)?;

    out.flush()?;
    Ok(())
}

/// Rough size of one formatted cell, tab included. It only sizes the blocks in
/// [`write_rows`], so a wrong guess costs memory or parallelism, never output.
#[cfg(feature = "cli")]
const CELL_BYTES: usize = 16;

/// Formatted text [`write_matrix_tsv`] holds at once: about 400 rows of a
/// 10 000-tree matrix.
#[cfg(feature = "cli")]
const BLOCK_BYTES: usize = 64 << 20;

/// Write the matrix body, `rows_per_block` rows at a time: each block's rows
/// are formatted in parallel, then written in order.
#[cfg(feature = "cli")]
fn write_rows<T: std::fmt::Display + Sync>(
    out: &mut dyn Write,
    names: &[String],
    mat: &[T],
    n_trees: usize,
    rows_per_block: usize,
) -> io::Result<()> {
    use crate::par::*;

    if n_trees == 0 {
        return Ok(());
    }
    for (block, rows) in mat.chunks(n_trees * rows_per_block).enumerate() {
        let first = block * rows_per_block;
        let texts = rows
            .par_chunks(n_trees)
            .enumerate()
            .map(|(k, row)| -> io::Result<Vec<u8>> {
                let mut text = Vec::with_capacity(row.len() * CELL_BYTES);
                write!(text, "{}", names[first + k])?;
                for val in row {
                    write!(text, "\t{val}")?;
                }
                text.push(b'\n');
                Ok(text)
            })
            .collect::<io::Result<Vec<_>>>()?;
        for text in &texts {
            out.write_all(text)?;
        }
    }
    Ok(())
}

/// Split a NEXUS tree header into `(tree_name, state_number)`.
pub fn extract_name_state(header: &str) -> (String, usize) {
    let upper = header.to_ascii_uppercase();
    if let Some(state_pos) = upper.find("STATE_")
        && let Some((_, rest)) = header.split_once(' ')
    {
        // NEXUS allows an optional `*` marking the default tree: `TREE * STATE_0 = ...`
        let tree_name = rest
            .split_whitespace()
            .find(|&t| t != "*")
            .unwrap_or("")
            .to_string();
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

/// Case-insensitive ASCII prefix test that allocates no uppercased copy.
#[inline]
fn starts_with_ci(s: &str, prefix: &str) -> bool {
    s.get(..prefix.len())
        .is_some_and(|p| p.eq_ignore_ascii_case(prefix))
}

/// Every `tree ... = ...` line of a NEXUS trees block, as `(header, body)`.
///
/// Lazy, so [`detect_format`] can stop at the first hit instead of collecting
/// the whole block just to ask whether one exists.
fn nexus_tree_lines(content: &str) -> impl Iterator<Item = (&str, &str)> {
    content
        .lines()
        .map(str::trim)
        .skip_while(|line| !starts_with_ci(line, "tree "))
        .take_while(|line| !starts_with_ci(line, "end;"))
        .filter_map(|line| {
            let (header, body) = line.split_once(" = ")?;
            Some((header.trim(), body.trim()))
        })
}

/// Collect a NEXUS trees block, naming each tree after its `tree` header.
fn collect_nexus_trees(content: &str, base_name: &str) -> Vec<RawTree> {
    nexus_tree_lines(content)
        .map(|(header, body)| {
            let (name, state) = extract_name_state(header);
            RawTree {
                name: format!("{base_name}_{name}"),
                state,
                newick: body.to_string(),
            }
        })
        .collect()
}

/// Collect a plain Newick file, naming each tree after the line it starts on.
///
/// Trees are delimited by `;`, so a line may hold several and one tree may wrap
/// across lines. [`strip_beast_annotations`] runs per line, which keeps the line
/// count intact and takes `[&...]` comments out of the way before the `;` split.
fn collect_newick_trees(content: &str, base_name: &str) -> Vec<RawTree> {
    let mut trees = Vec::new();
    let mut buf = String::new();
    let mut start_line = 1;

    for (idx, line) in content.lines().enumerate() {
        let stripped = strip_annotations(line.trim());
        let line = stripped.trim();
        if line.is_empty() {
            continue;
        }
        if buf.is_empty() {
            start_line = idx + 1;
        }
        buf.push_str(line);

        // Several trees may share a line; each keeps that line as its origin.
        while let Some(end) = buf.find(';') {
            let rest = buf.split_off(end + 1);
            let newick = std::mem::replace(&mut buf, rest.trim_start().to_string());
            trees.extend(newick_tree(base_name, start_line, newick));
        }
    }

    // A trailing tree with no `;` is kept, leaving the parser to complain.
    trees.extend(newick_tree(base_name, start_line, buf));
    trees
}

/// A collected chunk is a tree only if it opens with `(`, which is what keeps a
/// file of neither format from coming back as one nonsense tree.
fn newick_tree(base_name: &str, line: usize, newick: String) -> Option<RawTree> {
    newick.starts_with('(').then(|| RawTree {
        name: format!("{base_name}_line{line}"),
        state: 0,
        newick,
    })
}

/// Read a NEXUS `TRANSLATE` block into a numeric-ID → taxon-name map.
pub fn parse_taxon_block(content: &str) -> HashMap<String, String> {
    content
        .lines()
        .skip_while(|line| !starts_with_ci(line.trim(), "translate"))
        .skip(1)
        .take_while(|line| !line.trim().starts_with(';'))
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

    /// The same 21 trees as `hiv2.trees`, written one per line as plain Newick.
    fn hiv2_newick_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hiv2.newick")
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

    // ── nexus_tree_lines ──────────────────────────────────────────────────────

    #[test]
    fn test_nexus_tree_lines_count() {
        let content = "Begin trees;\ntree t1 = (A:1,B:1);\ntree t2 = (A:2,B:2);\nEnd;\n";
        assert_eq!(nexus_tree_lines(content).count(), 2);
    }

    #[test]
    fn test_nexus_tree_lines_header_and_body() {
        let content = "Begin trees;\ntree STATE_0 = (A:1,B:1);\nEnd;\n";
        let lines: Vec<_> = nexus_tree_lines(content).collect();
        assert_eq!(lines, vec![("tree STATE_0", "(A:1,B:1);")]);
    }

    #[test]
    fn test_nexus_tree_lines_indented_and_starred() {
        let content = "Begin trees;\n\tTREE * STATE_10 = (A:1,B:1);\n\tEnd;\n";
        let lines: Vec<_> = nexus_tree_lines(content).collect();
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].1, "(A:1,B:1);");
        assert_eq!(extract_name_state(lines[0].0), ("STATE_10".into(), 10));
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

    /// NEXUS bodies come back as written. The parse strips `[&...]` annotations,
    /// and it runs in parallel, so the serial reader leaves them in place.
    #[test]
    fn test_load_beast_raw_leaves_annotations_to_the_parse() {
        let (_, pairs) = load_beast_raw(hiv2_path(), 0, 0, false);
        assert!(pairs.iter().any(|(_, newick)| newick.contains("[&")));
        for (_, newick) in &pairs {
            assert!(
                !strip_annotations(newick).contains("[&"),
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

    /// Writes `content` to a temporary `.trees` file and loads it.
    fn load_raw_from_str(content: &str) -> Vec<(String, String)> {
        load_raw_from_str_burnin(content, ".trees", 0, 0)
    }

    /// Writes `content` to a temporary file named `*<suffix>` and loads it.
    fn load_raw_from_str_burnin(
        content: &str,
        suffix: &str,
        burnin_trees: usize,
        burnin_states: usize,
    ) -> Vec<(String, String)> {
        use std::io::Write;
        let mut tmp = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        tmp.write_all(content.as_bytes()).unwrap();
        load_beast_raw(tmp.path(), burnin_trees, burnin_states, false).1
    }

    #[test]
    fn test_load_beast_raw_reads_indented_starred_tree_lines() {
        // ape/R writes tab-indented `TREE * STATE_n = ...` lines; see issue #21.
        let pairs = load_raw_from_str(
            "#NEXUS\nBEGIN TREES;\n\tTREE * STATE_1 = [&R] (A:1,B:1);\n\tTREE * STATE_2 = [&R] (A:2,B:2);\nEND;\n",
        );
        assert_eq!(pairs.len(), 2);
        assert!(pairs[0].0.ends_with("_STATE_1"), "got {}", pairs[0].0);
        assert_eq!(strip_annotations(&pairs[0].1).trim(), "(A:1,B:1);");
    }

    #[test]
    fn test_load_beast_raw_without_tree_lines_returns_empty() {
        // Non-empty NEXUS with no `tree` lines: warns on stderr, yields nothing.
        let pairs = load_raw_from_str("#NEXUS\nBEGIN TAXA;\n\tDIMENSIONS NTAX=2;\nEND;\n");
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_load_beast_raw_empty_file_returns_empty() {
        assert!(load_raw_from_str("").is_empty());
    }

    // ── load_beast_trees ──────────────────────────────────────────────────────

    #[test]
    fn test_load_beast_trees_returns_correct_count() {
        let (names, snaps) = load_beast_trees(hiv2_path(), 0, 0, false, false, None);
        assert_eq!(names.len(), 21);
        assert_eq!(snaps.len(), 21);
        assert!(!snaps.leaf_names.is_empty());
    }

    #[test]
    fn test_load_beast_trees_burnin_reduces_count() {
        let (names, snaps) = load_beast_trees(hiv2_path(), 5, 0, false, false, None);
        assert_eq!(names.len(), 16);
        assert_eq!(snaps.len(), 16);
    }

    #[test]
    fn test_load_beast_trees_nonexistent_returns_empty() {
        let (names, snaps) = load_beast_trees("nonexistent.trees", 0, 0, false, false, None);
        assert!(names.is_empty());
        assert_eq!(snaps.len(), 0);
    }

    #[test]
    fn test_load_beast_trees_for_distances_matches_everything() {
        // The CLI loads with `Retain::for_distances`, which skips what a matrix
        // never reads. That must change no distance, and must actually skip it.
        let load =
            |retain: Option<Retain>| load_beast_trees(hiv2_path(), 0, 0, false, false, retain).1;
        let all = load(None);
        let rf = load(Some(Retain::for_distances(false)));
        let weighted = load(Some(Retain::for_distances(true)));

        assert_eq!(rf.pairwise_rf(None), all.pairwise_rf(None));
        assert_eq!(weighted.pairwise_wrf(None), all.pairwise_wrf(None));
        assert_eq!(weighted.pairwise_kf(None), all.pairwise_kf(None));

        assert_ne!(all.clades.len(), 0);
        assert_eq!((rf.clades.len(), weighted.clades.len()), (0, 0));
        for (r, w) in rf.snapshots.iter().zip(&weighted.snapshots) {
            assert!(r.lengths.is_empty(), "RF path must not store lengths");
            assert_eq!(w.lengths.len(), w.split_ids.len());
        }
    }

    // ── detect_format ─────────────────────────────────────────────────────────

    #[test]
    fn test_detect_format_nexus_header() {
        assert_eq!(
            detect_format("#NEXUS\nBEGIN TREES;\ntree t1 = (A:1,B:1);\nEND;\n"),
            TreeFormat::Nexus
        );
    }

    #[test]
    fn test_detect_format_nexus_without_header() {
        // Some writers omit `#NEXUS`; the trees block still gives it away.
        assert_eq!(
            detect_format("Begin trees;\n\ttree t1 = (A:1,B:1);\nEnd;\n"),
            TreeFormat::Nexus
        );
    }

    #[test]
    fn test_detect_format_plain_newick() {
        assert_eq!(
            detect_format("(A:1,B:1);\n(A:2,B:2);\n"),
            TreeFormat::Newick
        );
    }

    #[test]
    fn test_detect_format_newick_with_rooting_flag() {
        assert_eq!(detect_format("[&R] (A:1,B:1);\n"), TreeFormat::Newick);
    }

    #[test]
    fn test_detect_format_newick_after_blank_lines() {
        assert_eq!(detect_format("\n\n   \n(A:1,B:1);\n"), TreeFormat::Newick);
    }

    #[test]
    fn test_detect_format_unrecognised_reads_as_newick() {
        // Nothing NEXUS about it, so it falls through to the Newick reader,
        // which drops any block that does not open with `(`.
        assert_eq!(detect_format("not a tree file\n"), TreeFormat::Newick);
        assert!(collect_newick_trees("not a tree file\n", "f").is_empty());
    }

    #[test]
    fn test_detect_format_nexus_header_survives_empty_trees_block() {
        // No `tree` lines, but the header keeps it on the NEXUS error path.
        assert_eq!(
            detect_format("#NEXUS\nBEGIN TAXA;\n\tDIMENSIONS NTAX=2;\nEND;\n"),
            TreeFormat::Nexus
        );
    }

    #[test]
    fn test_detect_format_real_files() {
        let nexus = std::fs::read_to_string(hiv2_path()).unwrap();
        let newick = std::fs::read_to_string(hiv2_newick_path()).unwrap();
        assert_eq!(detect_format(&nexus), TreeFormat::Nexus);
        assert_eq!(detect_format(&newick), TreeFormat::Newick);
    }

    // ── collect_newick_trees ──────────────────────────────────────────────────

    /// `(line number, newick)` for each tree, which is all these tests assert on.
    fn newick_blocks(content: &str) -> Vec<(String, String)> {
        collect_newick_trees(content, "f")
            .into_iter()
            .map(|tree| (tree.name, tree.newick))
            .collect()
    }

    #[test]
    fn test_collect_newick_trees_one_tree_per_line() {
        assert_eq!(
            newick_blocks("(A:1,B:1);\n(A:2,B:2);\n"),
            vec![
                ("f_line1".to_string(), "(A:1,B:1);".to_string()),
                ("f_line2".to_string(), "(A:2,B:2);".to_string()),
            ]
        );
    }

    #[test]
    fn test_collect_newick_trees_skips_blank_lines() {
        // Line numbers track the file, not the tree index.
        let blocks = newick_blocks("\n(A:1,B:1);\n\n(A:2,B:2);\n");
        assert_eq!(
            blocks.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>(),
            ["f_line2", "f_line4"]
        );
    }

    #[test]
    fn test_collect_newick_trees_several_trees_on_one_line() {
        let blocks = newick_blocks("(A:1,B:1);(A:2,B:2);\n");
        assert_eq!(blocks.len(), 2);
        assert!(blocks.iter().all(|(name, _)| name == "f_line1"));
    }

    #[test]
    fn test_collect_newick_trees_tree_wrapped_over_lines() {
        // A tree split across lines is joined and reported at its first line.
        assert_eq!(
            newick_blocks("((A:1,\nB:1):1,\n(C:1,D:1):1);\n"),
            vec![(
                "f_line1".to_string(),
                "((A:1,B:1):1,(C:1,D:1):1);".to_string()
            )]
        );
    }

    #[test]
    fn test_collect_newick_trees_strips_annotations_before_splitting() {
        // Annotations go first, so a `;` inside one cannot end the tree early.
        assert_eq!(
            newick_blocks("(A[&note=a;b]:1,B:1);\n"),
            vec![("f_line1".to_string(), "(A:1,B:1);".to_string())]
        );
    }

    #[test]
    fn test_collect_newick_trees_drops_leading_rooting_flag() {
        assert_eq!(
            newick_blocks("[&R] (A:1,B:1);\n"),
            vec![("f_line1".to_string(), "(A:1,B:1);".to_string())]
        );
    }

    #[test]
    fn test_collect_newick_trees_keeps_unterminated_tail() {
        let blocks = newick_blocks("(A:1,B:1);\n(A:2,B:2)\n");
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[1], ("f_line2".to_string(), "(A:2,B:2)".to_string()));
    }

    #[test]
    fn test_collect_newick_trees_have_no_state() {
        // No `STATE_` labels in Newick, so every tree sits at state 0.
        let trees = collect_newick_trees("(A:1,B:1);\n(A:2,B:2);\n", "f");
        assert!(trees.iter().all(|tree| tree.state == 0));
    }

    // ── load_beast_raw: plain Newick ──────────────────────────────────────────

    #[test]
    fn test_load_newick_raw_names_trees_by_line() {
        let pairs = load_raw_from_str_burnin("(A:1,B:1);\n(A:2,B:2);\n", ".newick", 0, 0);
        assert_eq!(pairs.len(), 2);
        assert!(pairs[0].0.ends_with("_line1"), "got {}", pairs[0].0);
        assert!(pairs[1].0.ends_with("_line2"), "got {}", pairs[1].0);
    }

    #[test]
    fn test_load_newick_raw_strips_annotations() {
        let pairs = load_raw_from_str_burnin("(A[&rate=0.5]:1,B:1);\n", ".newick", 0, 0);
        assert_eq!(pairs[0].1, "(A:1,B:1);");
    }

    #[test]
    fn test_load_newick_raw_burnin_by_tree_count() {
        let pairs =
            load_raw_from_str_burnin("(A:1,B:1);\n(A:2,B:2);\n(A:3,B:3);\n", ".newick", 2, 0);
        assert_eq!(pairs.len(), 1);
        assert!(pairs[0].0.ends_with("_line3"), "got {}", pairs[0].0);
    }

    #[test]
    fn test_load_newick_raw_burnin_by_state_is_ignored() {
        // Newick carries no STATE_ labels, so burnin-by-state cannot apply.
        let pairs = load_raw_from_str_burnin("(A:1,B:1);\n(A:2,B:2);\n", ".newick", 0, 500);
        assert_eq!(pairs.len(), 2);
    }

    #[test]
    fn test_load_newick_raw_hiv2_returns_all_trees() {
        let (translate, pairs) = load_beast_raw(hiv2_newick_path(), 0, 0, true);
        assert_eq!(pairs.len(), 21);
        assert!(translate.is_empty(), "Newick files have no TRANSLATE block");
        assert!(pairs[0].0.starts_with("hiv2_line"), "got {}", pairs[0].0);
        for (_, newick) in &pairs {
            assert!(!newick.contains("[&"), "annotations must be stripped");
        }
    }

    #[test]
    fn test_load_beast_trees_newick_matches_nexus() {
        // hiv2.newick holds the same 21 trees as hiv2.trees, so both files must
        // yield identical RF distances regardless of format.
        let (nexus_names, nexus_snaps) = load_beast_trees(hiv2_path(), 0, 0, false, false, None);
        let (newick_names, newick_snaps) =
            load_beast_trees(hiv2_newick_path(), 0, 0, false, false, None);
        assert_eq!(newick_names.len(), nexus_names.len());
        assert_eq!(newick_snaps.len(), nexus_snaps.len());
        assert_eq!(
            newick_snaps.pairwise_rf(None),
            nexus_snaps.pairwise_rf(None)
        );
    }
}

#[cfg(all(test, feature = "cli"))]
mod write_tests {
    use super::{write_matrix_tsv, write_rows};
    use std::io::Read;

    /// The matrix body as the writer produced it before it went parallel: one
    /// row after another, one cell at a time.
    fn serial_body<T: std::fmt::Display>(names: &[String], mat: &[T], n: usize) -> String {
        let mut text = String::new();
        for (i, row) in mat.chunks(n).enumerate() {
            text.push_str(&names[i]);
            for val in row {
                text.push_str(&format!("\t{val}"));
            }
            text.push('\n');
        }
        text
    }

    fn tree_names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("tree_{i}")).collect()
    }

    #[test]
    fn rows_match_serial_formatting_across_block_boundaries() {
        let n = 7;
        let names = tree_names(n);
        let weighted: Vec<f64> = (0..n * n).map(|k| k as f64 / 3.0).collect();
        let rf: Vec<u32> = (0..n * n).map(|k| k as u32 * 2).collect();
        let (want_weighted, want_rf) = (
            serial_body(&names, &weighted, n),
            serial_body(&names, &rf, n),
        );

        for rows_per_block in [1, 2, 3, n, n + 5] {
            let mut got = Vec::new();
            write_rows(&mut got, &names, &weighted, n, rows_per_block).unwrap();
            assert_eq!(
                String::from_utf8(got).unwrap(),
                want_weighted,
                "f64, rows_per_block={rows_per_block}"
            );

            let mut got = Vec::new();
            write_rows(&mut got, &names, &rf, n, rows_per_block).unwrap();
            assert_eq!(
                String::from_utf8(got).unwrap(),
                want_rf,
                "u32, rows_per_block={rows_per_block}"
            );
        }
    }

    #[test]
    fn writes_header_and_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("matrix.tsv");
        let mat: Vec<u32> = vec![0, 4, 2, 4, 0, 6, 2, 6, 0];
        write_matrix_tsv(&path, &tree_names(3), &mat, 3).unwrap();
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "\ttree_0\ttree_1\ttree_2\ntree_0\t0\t4\t2\ntree_1\t4\t0\t6\ntree_2\t2\t6\t0\n"
        );
    }

    #[test]
    fn gz_output_decompresses_to_the_plain_text() {
        let dir = tempfile::tempdir().unwrap();
        let (plain, gz) = (dir.path().join("m.tsv"), dir.path().join("m.tsv.gz"));
        let names = tree_names(4);
        let mat: Vec<f64> = (0..16).map(|k| k as f64 * 0.125).collect();
        write_matrix_tsv(&plain, &names, &mat, 4).unwrap();
        write_matrix_tsv(&gz, &names, &mat, 4).unwrap();

        let mut decoded = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&gz).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(decoded, std::fs::read_to_string(&plain).unwrap());
    }
}
