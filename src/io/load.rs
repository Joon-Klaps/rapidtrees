//! Loading a whole tree file into [`Snapshots`].

use super::reader::{TreeReader, read_ahead};
use crate::snapshot::{Retain, Snapshots};
use std::path::Path;

/// Load and parse all trees from a NEXUS or plain-Newick tree file.
///
/// The format is detected per file: NEXUS trees keep their `tree` header
/// name, plain Newick trees are named `<file stem>_line<n>` after the line
/// they start on. The file is streamed: a chunk of trees is read, parsed and
/// dropped before the next is read, so memory does not grow with the size of
/// the file. Returns `(tree_names, Snapshots)`. On any error, prints to stderr
/// and returns an empty `Snapshots`.
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
    let path = path.as_ref();
    let (translate, reader) =
        match TreeReader::open(path, burnin_trees, burnin_states, use_real_taxa) {
            Ok(opened) => opened,
            Err(e) => {
                eprintln!("Failed to read {path:?}: {e}");
                return (Vec::new(), Snapshots::empty());
            }
        };

    let mut names = Vec::new();
    let mut failed = None;
    // Fused: the builder asks again after a short chunk, and `map_while`
    // would read on past the error that ended it.
    let entries = read_ahead(reader)
        .map_while(|tree| match tree {
            Ok(mut tree) => {
                names.push(std::mem::take(&mut tree.name));
                Some((tree, &translate))
            }
            Err(e) => {
                failed = Some(e);
                None
            }
        })
        .fuse();
    let built = Snapshots::from_newick_iter_opts(entries, rooted, retain);
    match (failed, built) {
        (None, Ok(snaps)) => return (names, snaps),
        (Some(e), _) => eprintln!("Failed to read {path:?}: {e}"),
        (None, Err(e)) => eprintln!("Failed to parse trees in {path:?}: {e}"),
    }
    (Vec::new(), Snapshots::empty())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{hiv2_newick_path, hiv2_path};

    /// Writes `content` to a temporary `*<suffix>` file and loads it for RF,
    /// resolving TRANSLATE.
    fn load_bytes(content: &[u8], suffix: &str) -> (Vec<String>, Snapshots) {
        use std::io::Write;
        let mut tmp = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        tmp.write_all(content).unwrap();
        load_beast_trees(tmp.path(), 0, 0, true, false, Retain::for_distances(false))
    }

    #[test]
    fn test_load_beast_trees_crlf_line_endings() {
        let nexus = b"#NEXUS\r\nBegin trees;\r\n\tTranslate\r\n\t\t1 A,\r\n\t\t2 B,\r\n\t\t3 C,\r\n\t\t4 D\r\n\t\t;\r\ntree STATE_0 = ((1:1,2:1):1,(3:1,4:1):1);\r\ntree STATE_10 = ((1:1,3:1):1,(2:1,4:1):1);\r\nEnd;\r\n";
        let (names, snaps) = load_bytes(nexus, ".trees");
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.leaf_names, ["A", "B", "C", "D"]);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);

        let newick = b"((A:1,B:1):1,(C:1,D:1):1);\r\n((A:1,C:1):1,(B:1,D:1):1);\r\n";
        let (names, snaps) = load_bytes(newick, ".newick");
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    #[test]
    fn test_load_beast_trees_headerless_nexus() {
        let content =
            b"Begin trees;\n\ttree t1 = (A:1,B:1,C:1);\n\ttree t2 = (A:2,B:2,C:2);\nEnd;\n";
        let (names, snaps) = load_bytes(content, ".trees");
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.len(), 2);
    }

    /// A MrBayes `.t` file: the TRANSLATE block ends on its last entry, so the
    /// tree lines after it are still trees.
    #[test]
    fn test_load_beast_trees_mrbayes_translate() {
        let content = b"#NEXUS\n[ID: 123]\nbegin trees;\n   translate\n       1 A,\n       2 B,\n       3 C,\n       4 D;\n   tree gen.0 = [&U] ((1:1,2:1):1,(3:1,4:1):1);\n   tree gen.100 = [&U] ((1:1,3:1):1,(2:1,4:1):1);\nend;\n";
        let (names, snaps) = load_bytes(content, ".t");
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.leaf_names, ["A", "B", "C", "D"]);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    /// A byte-order mark neither hides a `#NEXUS` header nor a Newick file's
    /// first tree.
    #[test]
    fn test_load_beast_trees_byte_order_mark() {
        let nexus = "\u{feff}#NEXUS\nBegin trees;\n\tTranslate\n\t\t1 A,\n\t\t2 B,\n\t\t3 C,\n\t\t4 D\n\t\t;\ntree STATE_0 = ((1:1,2:1):1,(3:1,4:1):1);\ntree STATE_10 = ((1:1,3:1):1,(2:1,4:1):1);\nEnd;\n";
        let (names, snaps) = load_bytes(nexus.as_bytes(), ".trees");
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.leaf_names, ["A", "B", "C", "D"]);

        let newick = "\u{feff}((A:1,B:1):1,(C:1,D:1):1);\n((A:1,C:1):1,(B:1,D:1):1);\n";
        let (names, snaps) = load_bytes(newick.as_bytes(), ".newick");
        assert_eq!(names.len(), 2);
        assert_eq!(snaps.pairwise_rf(None), vec![0, 2, 2, 0]);
    }

    /// More trees than the read-ahead queue holds still arrive complete and in
    /// order, and give the matrix the same trees give in memory.
    #[test]
    fn test_load_beast_trees_streams_many_trees_in_order() {
        let shapes = [
            "((A:1,B:1):1,(C:1,D:1):1,E:1);",
            "((A:1,C:1):1,(B:1,D:1):1,E:1);",
            "((A:1,D:1):1,(B:1,E:1):1,C:1);",
        ];
        let trees: Vec<&str> = (0..3000).map(|k| shapes[k % 3]).collect();
        let (names, snaps) = load_bytes(trees.join("\n").as_bytes(), ".newick");
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

    /// A read error part-way through fails the whole load rather than
    /// returning the trees before it.
    #[test]
    fn test_load_beast_trees_unreadable_line_fails_cleanly() {
        let mut content = b"((A:1,B:1):1,(C:1,D:1):1);\n".to_vec();
        content.extend_from_slice(b"((A:1,\xff\xfe:1):1,(C:1,D:1):1);\n");
        content.extend_from_slice(b"((A:1,C:1):1,(B:1,D:1):1);\n");
        let (names, snaps) = load_bytes(&content, ".newick");
        assert!(names.is_empty());
        assert!(snaps.is_empty());
    }

    #[test]
    fn test_load_beast_trees_bad_tree_fails_cleanly() {
        let content = b"((A:1,B:1):1,(C:1,D:1):1);\n((A:1,B:1):1,(C:1,E:1):1);\n";
        let (names, snaps) = load_bytes(content, ".newick");
        assert!(names.is_empty());
        assert!(snaps.is_empty());
    }

    #[test]
    fn test_load_beast_trees_newick_without_trees_is_empty() {
        let (names, snaps) = load_bytes(b"[a comment and nothing else]\n", ".newick");
        assert!(names.is_empty());
        assert!(snaps.is_empty());
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
            assert_eq!(w.lengths.len(), w.n_splits());
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
