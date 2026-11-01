//! Reading a whole tree file into [`Snapshots`].

use super::format::TreeFormat;
use super::reader::{TreeReader, TreeRecord, read_ahead};
use crate::snapshot::{Retain, Snapshots};
use std::collections::HashMap;
use std::io::{self, BufRead};
use std::path::Path;

/// Fewest trees [`read_trees`] and [`read_tree_file`] hand back: a pairwise
/// distance needs two. The Python bindings ask for as many.
pub const MIN_TREES: usize = 2;

/// How [`read_trees`] and [`read_tree_file`] read their input.
///
/// The default reads every tree, unrooted, with leaves named as the trees name
/// them, and keeps everything. Set the fields that differ:
///
/// ```
/// # use rapidtrees::{Retain, io::ReadOptions};
/// let mut opts = ReadOptions::default();
/// opts.burnin_trees = 100;
/// opts.retain = Retain::for_distances(false);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct ReadOptions {
    /// Drop the first `burnin_trees` trees.
    pub burnin_trees: usize,
    /// NEXUS only: drop the trees sampled at MCMC state `burnin_states` or
    /// before, read from the `STATE_<n>` in their names. With both burn-ins
    /// set, a tree is kept if it is past either one.
    pub burnin_states: usize,
    /// NEXUS only: name leaves by their `TRANSLATE` labels instead of the IDs
    /// the trees use.
    pub use_real_taxa: bool,
    /// Compare clades instead of bipartitions.
    pub rooted: bool,
    /// What to build besides the split IDs; see [`Retain`].
    pub retain: Retain,
}

/// Why [`read_trees`] or [`read_tree_file`] returned no trees.
///
/// A binding matches on it to raise errors of its own, and needs a last arm
/// for the variants a later release may add:
///
/// ```
/// # use rapidtrees::io::{ReadOptions, ReadTreesError, read_trees};
/// let one_tree = "((A:1,B:1):1,(C:1,D:1):1);\n";
/// let message = match read_trees(one_tree.as_bytes(), &ReadOptions::default()) {
///     Ok(_) => unreachable!("one tree is too few"),
///     Err(ReadTreesError::TooFewTrees { found }) => format!("only {found} tree"),
///     Err(other) => other.to_string(),
/// };
/// assert_eq!(message, "only 1 tree");
/// ```
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReadTreesError {
    /// The input could not be opened or read, or is not UTF-8.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Nothing in the input reads as a tree of the format it was sniffed as.
    #[error("no trees found: expected {}", .0.tree_shape())]
    NoTrees(TreeFormat),
    /// Fewer than [`MIN_TREES`] trees were left after burn-in.
    #[error("need at least {MIN_TREES} trees to compute pairwise distances, found {found}")]
    TooFewTrees {
        /// How many were left.
        found: usize,
    },
    /// A tree did not parse, or its taxa differ from the first tree's.
    #[error("{0}")]
    Parse(String),
}

/// Read the NEXUS or plain-Newick trees of `input` into [`Snapshots`].
///
/// `input` is any [`BufRead`]: a byte slice, a [`BufReader`](io::BufReader)
/// over a file or a socket, or a `&mut` one to keep it. The format is sniffed
/// from its first non-blank line. NEXUS trees are named after their `tree`
/// header (`STATE_0`), plain-Newick trees after the line they start on
/// (`line1`). Trees are read on the calling thread a chunk at a time, and each
/// chunk's text is dropped once parsed, so memory does not grow with the input.
/// Returns `(tree_names, snapshots)`, in input order.
///
/// # Errors
///
/// [`ReadTreesError::Io`] if `input` fails or is not UTF-8,
/// [`ReadTreesError::NoTrees`] if it holds no trees,
/// [`ReadTreesError::TooFewTrees`] if burn-in leaves fewer than [`MIN_TREES`],
/// and [`ReadTreesError::Parse`] if a tree does not parse or its taxa differ
/// from the first tree's.
///
/// ```
/// # use rapidtrees::io::{ReadOptions, read_trees};
/// let text = "((A:1,B:1):1,(C:1,D:1):1);\n((A:1,C:1):1,(B:1,D:1):1);\n";
/// let (names, snaps) = read_trees(text.as_bytes(), &ReadOptions::default())?;
/// assert_eq!(names, ["line1", "line2"]);
/// assert_eq!(snaps.pairwise_rf(None), [0, 2, 2, 0]);
/// # Ok::<(), rapidtrees::io::ReadTreesError>(())
/// ```
pub fn read_trees<R: BufRead>(
    input: R,
    opts: &ReadOptions,
) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
    let (translate, reader) = TreeReader::new(input, String::new(), "the input", None, opts)?;
    enough_trees(build(reader, &translate, opts)?)
}

/// [`read_trees`] on the file at `path`, with each tree name prefixed by the
/// file's stem (`posterior_STATE_0`). With the `parallel` feature, outside
/// wasm32, the file is read ahead on a background thread, so reading it
/// overlaps parsing it.
///
/// # Errors
///
/// As [`read_trees`], and [`ReadTreesError::Io`] if the file cannot be opened.
///
/// ```no_run
/// # use rapidtrees::io::{ReadOptions, read_tree_file};
/// let mut opts = ReadOptions::default();
/// opts.burnin_trees = 100;
/// let (names, snaps) = read_tree_file("posterior.trees", &opts)?;
/// let rf = snaps.into_pairwise_rf(None);
/// # Ok::<(), rapidtrees::io::ReadTreesError>(())
/// ```
pub fn read_tree_file<P: AsRef<Path>>(
    path: P,
    opts: &ReadOptions,
) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
    enough_trees(load_file(path.as_ref(), opts)?)
}

/// Load and parse all trees from a NEXUS or plain-Newick tree file.
///
/// [`read_tree_file`] with positional arguments, and an empty result where
/// that returns an error: the error is printed to stderr instead. Unlike
/// [`read_tree_file`], it loads a file of a single tree.
///
/// `retain` says what to build besides the split IDs; pass `None` to keep
/// everything.
#[deprecated(note = "use `read_tree_file`, which returns the error instead of printing it")]
pub fn load_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
    rooted: bool,
    retain: impl Into<Option<Retain>>,
) -> (Vec<String>, Snapshots) {
    let opts = ReadOptions {
        burnin_trees,
        burnin_states,
        use_real_taxa,
        rooted,
        retain: retain.into().unwrap_or_default(),
    };
    let path = path.as_ref();
    load_file(path, &opts).unwrap_or_else(|e| {
        eprintln!("Failed to read {path:?}: {e}");
        (Vec::new(), Snapshots::empty())
    })
}

/// The trees of the file at `path`, however few.
fn load_file(path: &Path, opts: &ReadOptions) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
    let (translate, reader) = TreeReader::open(path, opts)?;
    build(read_ahead(reader), &translate, opts)
}

/// Parse `trees` into [`Snapshots`], stopping at the first error, whether it
/// comes from reading or from parsing.
fn build(
    trees: impl Iterator<Item = Result<TreeRecord, ReadTreesError>>,
    translate: &HashMap<String, String>,
    opts: &ReadOptions,
) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
    let mut names = Vec::new();
    let mut failed = None;
    // Fused: the builder asks again after a short chunk, and `map_while`
    // would read on past the error that ended it.
    let entries = trees
        .map_while(|tree| match tree {
            Ok(mut tree) => {
                names.push(std::mem::take(&mut tree.name));
                Some((tree, translate))
            }
            Err(e) => {
                failed = Some(e);
                None
            }
        })
        .fuse();
    let built = Snapshots::from_newick_iter_opts(entries, opts.rooted, opts.retain);
    // A read error cut the trees short, so it outranks whatever the builder
    // made of the ones before it.
    if let Some(e) = failed {
        return Err(e);
    }
    let snaps = built.map_err(ReadTreesError::Parse)?;
    Ok((names, snaps))
}

/// `loaded`, if it holds at least [`MIN_TREES`] trees.
fn enough_trees(
    (names, snaps): (Vec<String>, Snapshots),
) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
    if names.len() < MIN_TREES {
        return Err(ReadTreesError::TooFewTrees { found: names.len() });
    }
    Ok((names, snaps))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{hiv2_newick_path, hiv2_path};

    /// What these tests read with: RF only, TRANSLATE resolved.
    fn rf_opts() -> ReadOptions {
        ReadOptions {
            use_real_taxa: true,
            retain: Retain::for_distances(false),
            ..ReadOptions::default()
        }
    }

    fn read(content: &[u8]) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
        read_trees(content, &rf_opts())
    }

    /// Writes `content` to a temporary `*<suffix>` file and reads it back.
    fn read_file(content: &[u8], suffix: &str) -> Result<(Vec<String>, Snapshots), ReadTreesError> {
        use std::io::Write;
        let mut tmp = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        tmp.write_all(content).unwrap();
        read_tree_file(tmp.path(), &rf_opts())
    }

    const QUARTETS: &str = "((A:1,B:1):1,(C:1,D:1):1);\n((A:1,C:1):1,(B:1,D:1):1);\n";

    #[test]
    fn test_read_trees_crlf_line_endings() {
        let nexus = b"#NEXUS\r\nBegin trees;\r\n\tTranslate\r\n\t\t1 A,\r\n\t\t2 B,\r\n\t\t3 C,\r\n\t\t4 D\r\n\t\t;\r\ntree STATE_0 = ((1:1,2:1):1,(3:1,4:1):1);\r\ntree STATE_10 = ((1:1,3:1):1,(2:1,4:1):1);\r\nEnd;\r\n";
        let (names, snaps) = read(nexus).unwrap();
        assert_eq!(names, ["STATE_0", "STATE_10"]);
        assert_eq!(snaps.leaf_names, ["A", "B", "C", "D"]);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);

        let newick = b"((A:1,B:1):1,(C:1,D:1):1);\r\n((A:1,C:1):1,(B:1,D:1):1);\r\n";
        let (names, snaps) = read(newick).unwrap();
        assert_eq!(names, ["line1", "line2"]);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    #[test]
    fn test_read_trees_headerless_nexus() {
        let content =
            b"Begin trees;\n\ttree t1 = (A:1,B:1,C:1);\n\ttree t2 = (A:2,B:2,C:2);\nEnd;\n";
        let (names, snaps) = read(content).unwrap();
        assert_eq!(names, ["t1", "t2"]);
        assert_eq!(snaps.len(), 2);
    }

    /// A MrBayes `.t` file: the TRANSLATE block ends on its last entry, so the
    /// tree lines after it are still trees, and they keep their `gen.N` names
    /// although they have no `STATE_`.
    #[test]
    fn test_read_trees_mrbayes() {
        let content = b"#NEXUS\n[ID: 123]\nbegin trees;\n   translate\n       1 A,\n       2 B,\n       3 C,\n       4 D;\n   tree gen.0 = [&U] ((1:1,2:1):1,(3:1,4:1):1);\n   tree gen.100 = [&U] ((1:1,3:1):1,(2:1,4:1):1);\nend;\n";
        let (names, snaps) = read(content).unwrap();
        assert_eq!(names, ["gen.0", "gen.100"]);
        assert_eq!(snaps.leaf_names, ["A", "B", "C", "D"]);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    /// A byte-order mark neither hides a `#NEXUS` header nor a Newick file's
    /// first tree.
    #[test]
    fn test_read_trees_byte_order_mark() {
        let nexus = "\u{feff}#NEXUS\nBegin trees;\n\tTranslate\n\t\t1 A,\n\t\t2 B,\n\t\t3 C,\n\t\t4 D\n\t\t;\ntree STATE_0 = ((1:1,2:1):1,(3:1,4:1):1);\ntree STATE_10 = ((1:1,3:1):1,(2:1,4:1):1);\nEnd;\n";
        let (names, snaps) = read(nexus.as_bytes()).unwrap();
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.leaf_names, ["A", "B", "C", "D"]);

        let newick = format!("\u{feff}{QUARTETS}");
        let (names, snaps) = read(newick.as_bytes()).unwrap();
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    /// A comment line before the first tree once hid that tree: it was glued
    /// to the tree, and the chunk then dropped for not opening with `(`.
    #[test]
    fn test_read_trees_comment_line_before_the_trees() {
        let content = format!("[written by some tool]\n{QUARTETS}");
        let (names, snaps) = read(content.as_bytes()).unwrap();
        assert_eq!(names, ["line2", "line3"]);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    #[test]
    fn test_read_trees_a_single_tree_is_too_few() {
        let result = read(b"((A:1,B:1):1,(C:1,D:1):1);\n");
        assert!(
            matches!(result, Err(ReadTreesError::TooFewTrees { found: 1 })),
            "got {result:?}"
        );
    }

    /// More trees than the read-ahead queue holds still arrive complete and in
    /// order, and give the matrix the same trees give in memory.
    #[test]
    fn test_read_tree_file_streams_many_trees_in_order() {
        let shapes = [
            "((A:1,B:1):1,(C:1,D:1):1,E:1);",
            "((A:1,C:1):1,(B:1,D:1):1,E:1);",
            "((A:1,D:1):1,(B:1,E:1):1,C:1);",
        ];
        let trees: Vec<&str> = (0..3000).map(|k| shapes[k % 3]).collect();
        let (names, snaps) = read_file(trees.join("\n").as_bytes(), ".newick").unwrap();
        assert_eq!(names.len(), trees.len());
        for (k, name) in names.iter().enumerate() {
            assert!(
                name.ends_with(&format!("_line{}", k + 1)),
                "tree {k} is {name}"
            );
        }
        let reference = Snapshots::from_newicks(&trees, false).unwrap();
        assert_eq!(snaps.pairwise_rf(None), reference.pairwise_rf(None));
    }

    /// A read error part-way through fails the whole read rather than
    /// returning the trees before it, on the calling thread and read ahead.
    #[test]
    fn test_read_trees_unreadable_line_fails() {
        let mut content = b"((A:1,B:1):1,(C:1,D:1):1);\n".to_vec();
        content.extend_from_slice(b"((A:1,\xff\xfe:1):1,(C:1,D:1):1);\n");
        content.extend_from_slice(b"((A:1,C:1):1,(B:1,D:1):1);\n");
        for result in [read(&content), read_file(&content, ".newick")] {
            assert!(
                matches!(result, Err(ReadTreesError::Io(_))),
                "got {result:?}"
            );
        }
    }

    #[test]
    fn test_read_trees_bad_tree_fails() {
        let content = b"((A:1,B:1):1,(C:1,D:1):1);\n((A:1,B:1):1,(C:1,E:1):1);\n";
        let result = read(content);
        assert!(
            matches!(result, Err(ReadTreesError::Parse(_))),
            "got {result:?}"
        );
    }

    #[test]
    fn test_read_trees_without_trees_says_which_format() {
        let newick = read(b"[a comment and nothing else]\n");
        assert!(
            matches!(newick, Err(ReadTreesError::NoTrees(TreeFormat::Newick))),
            "got {newick:?}"
        );
        let nexus = read(b"#NEXUS\nBEGIN TAXA;\n\tDIMENSIONS NTAX=2;\nEND;\n");
        assert!(
            matches!(nexus, Err(ReadTreesError::NoTrees(TreeFormat::Nexus))),
            "got {nexus:?}"
        );
    }

    #[test]
    fn test_read_trees_error_messages() {
        let too_few = ReadTreesError::TooFewTrees { found: 1 };
        assert_eq!(
            too_few.to_string(),
            "need at least 2 trees to compute pairwise distances, found 1"
        );
        assert_eq!(
            ReadTreesError::NoTrees(TreeFormat::Nexus).to_string(),
            "no trees found: expected NEXUS `tree <name> = <newick>` lines"
        );
    }

    /// Bindings hand the error across threads and wrap it in their own.
    #[test]
    fn test_read_trees_error_is_send_sync() {
        fn assert_error<E: std::error::Error + Send + Sync + 'static>() {}
        assert_error::<ReadTreesError>();
    }

    // ── read_tree_file ────────────────────────────────────────────────────────

    #[test]
    fn test_read_tree_file_returns_correct_count() {
        let (names, snaps) = read_tree_file(hiv2_path(), &ReadOptions::default()).unwrap();
        assert_eq!(names.len(), 21);
        assert_eq!(snaps.len(), 21);
        assert!(!snaps.leaf_names.is_empty());
        assert!(names.iter().all(|name| name.starts_with("hiv2_STATE_")));
    }

    #[test]
    fn test_read_tree_file_burnin_reduces_count() {
        let opts = ReadOptions {
            burnin_trees: 5,
            ..ReadOptions::default()
        };
        let (names, snaps) = read_tree_file(hiv2_path(), &opts).unwrap();
        assert_eq!(names.len(), 16);
        assert_eq!(snaps.len(), 16);
    }

    #[test]
    fn test_read_tree_file_burnin_leaving_one_tree_is_too_few() {
        let opts = ReadOptions {
            burnin_trees: 20,
            ..ReadOptions::default()
        };
        let result = read_tree_file(hiv2_path(), &opts);
        assert!(
            matches!(result, Err(ReadTreesError::TooFewTrees { found: 1 })),
            "got {result:?}"
        );
    }

    #[test]
    fn test_read_tree_file_nonexistent_is_not_found() {
        let result = read_tree_file("nonexistent.trees", &ReadOptions::default());
        assert!(
            matches!(&result, Err(ReadTreesError::Io(e)) if e.kind() == io::ErrorKind::NotFound),
            "got {result:?}"
        );
    }

    /// The text and the file give the same trees; only the file's stem sets
    /// the names apart.
    #[test]
    fn test_read_tree_file_matches_read_trees() {
        let text = std::fs::read(hiv2_path()).unwrap();
        let (text_names, text_snaps) = read(&text).unwrap();
        let (file_names, file_snaps) = read_tree_file(hiv2_path(), &rf_opts()).unwrap();
        let unprefixed: Vec<_> = file_names
            .iter()
            .map(|name| name.strip_prefix("hiv2_").unwrap())
            .collect();
        assert_eq!(text_names, unprefixed);
        assert_eq!(text_snaps.pairwise_rf(None), file_snaps.pairwise_rf(None));
    }

    #[test]
    fn test_read_tree_file_for_distances_matches_everything() {
        // The CLI reads with `Retain::for_distances`, which skips what a matrix
        // never reads. That must change no distance, and must actually skip it.
        let load = |retain: Retain| {
            let opts = ReadOptions {
                retain,
                ..ReadOptions::default()
            };
            read_tree_file(hiv2_path(), &opts).unwrap().1
        };
        let all = load(Retain::everything());
        let rf = load(Retain::for_distances(false));
        let weighted = load(Retain::for_distances(true));

        assert_eq!(rf.pairwise_rf(None), all.pairwise_rf(None));
        assert_eq!(weighted.pairwise_wrf(None), all.pairwise_wrf(None));
        assert_eq!(weighted.pairwise_kf(None), all.pairwise_kf(None));

        assert_ne!(all.clades.len(), 0);
        assert_eq!((rf.clades.len(), weighted.clades.len()), (0, 0));
        for (r, w) in rf.snapshots.iter().zip(&weighted.snapshots) {
            assert!(r.lengths.is_empty(), "RF path must not store lengths");
            assert_eq!(w.lengths.len(), w.n_splits());
        }
    }

    #[test]
    fn test_read_tree_file_newick_matches_nexus() {
        // hiv2.newick holds the same 21 trees as hiv2.trees, so both files must
        // yield identical RF distances regardless of format.
        let opts = ReadOptions::default();
        let (nexus_names, nexus_snaps) = read_tree_file(hiv2_path(), &opts).unwrap();
        let (newick_names, newick_snaps) = read_tree_file(hiv2_newick_path(), &opts).unwrap();
        assert_eq!(newick_names.len(), nexus_names.len());
        assert_eq!(newick_snaps.len(), nexus_snaps.len());
        assert_eq!(
            newick_snaps.pairwise_rf(None),
            nexus_snaps.pairwise_rf(None)
        );
    }

    // ── load_beast_trees ──────────────────────────────────────────────────────

    /// The deprecated wrapper keeps its contract: errors become an empty
    /// result, and a single tree still loads.
    #[test]
    #[allow(deprecated)]
    fn test_load_beast_trees_keeps_its_contract() {
        let (names, snaps) = load_beast_trees(hiv2_path(), 5, 0, false, false, None);
        assert_eq!((names.len(), snaps.len()), (16, 16));

        let (names, snaps) = load_beast_trees("nonexistent.trees", 0, 0, false, false, None);
        assert!(names.is_empty() && snaps.is_empty());

        let mut tmp = tempfile::Builder::new()
            .suffix(".newick")
            .tempfile()
            .unwrap();
        io::Write::write_all(&mut tmp, b"((A:1,B:1):1,(C:1,D:1):1);\n").unwrap();
        let (names, snaps) = load_beast_trees(tmp.path(), 0, 0, false, false, None);
        assert_eq!((names.len(), snaps.len()), (1, 1));
    }
}
