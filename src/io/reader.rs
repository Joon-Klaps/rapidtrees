//! The streaming tree reader: one tree at a time, burn-in applied.

use super::annotations::strip_annotations;
use super::format::TreeFormat;
use super::nexus::{closes_translate, extract_name_state, starts_with_ci, translate_entry};
use std::borrow::Cow;
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::io::{self, BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

/// Buffer between the file and the reader. Every file allocates one, and a
/// larger one reads no faster, not even lines of over a megabyte.
const READ_BUFFER_BYTES: usize = 64 << 10;

/// Bytes of tree text the background reader may hold ahead of the parser:
/// enough to keep it fed, and never fewer than two trees.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
const READ_AHEAD_BYTES: usize = 32 << 20;

/// One tree as read from a tree file.
///
/// It keeps the line it was read in, and its Newick text is a range of that
/// buffer, so reading a tree never copies it.
pub(crate) struct TreeRecord {
    /// Display name, already prefixed with the file stem.
    pub(crate) name: String,
    /// MCMC state it was sampled at; `0` when the format records none.
    pub(crate) state: usize,
    text: String,
    span: Range<usize>,
}

impl TreeRecord {
    /// The tree's Newick string.
    pub(crate) fn newick(&self) -> &str {
        &self.text[self.span.clone()]
    }
}

impl AsRef<str> for TreeRecord {
    fn as_ref(&self) -> &str {
        self.newick()
    }
}

/// Byte offset of `part` in `whole`, which `part` must borrow from.
fn offset_in(whole: &str, part: &str) -> usize {
    part.as_ptr() as usize - whole.as_ptr() as usize
}

/// Reads the trees of a NEXUS or plain-Newick file one at a time.
///
/// Opening reads the header. The first non-empty line decides the format:
/// `#NEXUS`, or a `begin` or `tree` line, is NEXUS, and anything else is
/// Newick. For NEXUS the `TRANSLATE` block comes back complete before the
/// first tree is read. Trees then arrive one by one with burn-in applied, each
/// in the line it was read from, so only the trees in flight are ever in
/// memory, however large the file and however much of it is annotation.
///
/// NEXUS trees are the `name = newick` lines from the first `tree` line to
/// `end;`, named after their header. Newick trees end at `;`, may share a line
/// or wrap across lines, have `[&...]` annotations stripped per line, and are
/// named `<file stem>_line<n>` after the line they start on.
pub(crate) struct TreeReader<R> {
    lines: R,
    /// The line waiting to be processed; empty once the input is exhausted.
    line: String,
    /// Lines read so far, which is the 1-based number of `line`.
    lines_read: usize,
    format: TreeFormat,
    base_name: String,
    /// How the input is named in messages.
    label: String,
    burnin_trees: usize,
    burnin_states: usize,
    /// Trees read so far, burn-in included.
    seen: usize,
    /// Newick only: the text of a tree still waiting for its `;`, and the
    /// line it started on.
    pending: String,
    pending_line: usize,
    /// Trees read but not yet handed out.
    ready: VecDeque<TreeRecord>,
    /// Whether the input held anything besides blank lines.
    has_content: bool,
    done: bool,
}

impl TreeReader<BufReader<fs::File>> {
    /// Open `path` and read its header. Returns the `TRANSLATE` map, empty
    /// unless `use_real_taxa`, and a reader positioned at the first tree.
    pub(crate) fn open(
        path: &Path,
        burnin_trees: usize,
        burnin_states: usize,
        use_real_taxa: bool,
    ) -> io::Result<(HashMap<String, String>, Self)> {
        let base_name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        let lines = BufReader::with_capacity(READ_BUFFER_BYTES, fs::File::open(path)?);
        let label = format!("{path:?}");
        Self::new(
            lines,
            base_name,
            &label,
            None,
            burnin_trees,
            burnin_states,
            use_real_taxa,
        )
    }
}

impl<R: BufRead> TreeReader<R> {
    /// Read the header of `lines`; `format` forces a format instead of
    /// sniffing one.
    fn new(
        lines: R,
        base_name: &str,
        label: &str,
        format: Option<TreeFormat>,
        burnin_trees: usize,
        burnin_states: usize,
        use_real_taxa: bool,
    ) -> io::Result<(HashMap<String, String>, Self)> {
        let mut reader = Self {
            lines,
            line: String::new(),
            lines_read: 0,
            format: TreeFormat::Newick,
            base_name: base_name.to_string(),
            label: label.to_string(),
            burnin_trees,
            burnin_states,
            seen: 0,
            pending: String::new(),
            pending_line: 0,
            ready: VecDeque::new(),
            has_content: false,
            done: false,
        };
        reader.read_line()?;
        // A byte-order mark is not text: left in, it hides `#NEXUS` from the
        // sniff below and a Newick tree's opening `(`.
        if reader.line.starts_with('\u{feff}') {
            reader.line.drain(..'\u{feff}'.len_utf8());
        }
        while !reader.line.is_empty() && reader.line.trim().is_empty() {
            reader.read_line()?;
        }
        let head = reader.line.trim();
        reader.has_content = !head.is_empty();
        let nexus = starts_with_ci(head, "#nexus")
            || starts_with_ci(head, "begin")
            || starts_with_ci(head, "tree ");
        reader.format = format.unwrap_or(if nexus {
            TreeFormat::Nexus
        } else {
            TreeFormat::Newick
        });

        let mut translate = HashMap::new();
        if reader.format == TreeFormat::Nexus {
            reader.read_nexus_header(use_real_taxa, &mut translate)?;
        } else if burnin_states > 0 {
            // Newick trees carry no `STATE_` labels, so a state threshold cannot apply.
            eprintln!(
                "Ignoring burn-in by state for {label}: Newick trees carry no `STATE_` labels"
            );
            reader.burnin_states = 0;
        }
        Ok((translate, reader))
    }

    /// Replace `line` with the next line of input, or empty it at the end.
    fn read_line(&mut self) -> io::Result<()> {
        self.line.clear();
        if self.lines.read_line(&mut self.line)? > 0 {
            self.lines_read += 1;
        }
        Ok(())
    }

    /// Hand `line` over whole, leaving a buffer sized for the next one: the
    /// trees of a run are about the same length, so it rarely has to grow.
    fn take_line(&mut self) -> String {
        let capacity = self.line.len() + self.line.len() / 8;
        std::mem::replace(&mut self.line, String::with_capacity(capacity))
    }

    /// Skip the NEXUS header up to the first `tree` line, reading the
    /// `TRANSLATE` block into `translate` on the way when `use_real_taxa`.
    fn read_nexus_header(
        &mut self,
        use_real_taxa: bool,
        translate: &mut HashMap<String, String>,
    ) -> io::Result<()> {
        while !self.line.is_empty() && !starts_with_ci(self.line.trim(), "tree ") {
            if use_real_taxa && starts_with_ci(self.line.trim(), "translate") {
                self.read_line()?;
                while !self.line.is_empty() {
                    translate.extend(translate_entry(&self.line));
                    if closes_translate(&self.line) {
                        break;
                    }
                    self.read_line()?;
                }
            }
            self.read_line()?;
        }
        Ok(())
    }

    /// Read the current line, queueing any tree it completes, and move on to
    /// the next one.
    fn advance(&mut self) -> io::Result<()> {
        let more = !self.line.is_empty()
            && match self.format {
                TreeFormat::Nexus => self.nexus_line(),
                TreeFormat::Newick => {
                    self.newick_line();
                    true
                }
            };
        if more {
            self.read_line()
        } else {
            self.finish();
            Ok(())
        }
    }

    /// Queue the tree on the current NEXUS line, if it holds one. `false` at
    /// the trees block's `end;`.
    fn nexus_line(&mut self) -> bool {
        let trimmed = self.line.trim();
        if starts_with_ci(trimmed, "end;") {
            return false;
        }
        if let Some((header, body)) = trimmed.split_once(" = ") {
            let (name, state) = extract_name_state(header.trim());
            let body = body.trim();
            let start = offset_in(&self.line, body);
            let span = start..start + body.len();
            let name = format!("{}_{name}", self.base_name);
            let text = self.take_line();
            self.offer(TreeRecord {
                name,
                state,
                text,
                span,
            });
        }
        true
    }

    /// Add the current Newick line to the tree being read, queueing every
    /// tree it completes.
    fn newick_line(&mut self) {
        // A whole tree alone on its line, with nothing stripped from it, is
        // handed over in the line it was read into. Anything else is gathered
        // in `pending` and split at each `;`.
        let alone = {
            let stripped = strip_annotations(self.line.trim());
            let piece = stripped.trim();
            if piece.is_empty() {
                return;
            }
            if self.pending.is_empty() {
                self.pending_line = self.lines_read;
            }
            if self.pending.is_empty()
                && matches!(stripped, Cow::Borrowed(_))
                && piece.find(';') == Some(piece.len() - 1)
            {
                let start = offset_in(&self.line, piece);
                Some(start..start + piece.len())
            } else {
                self.pending.push_str(piece);
                None
            }
        };
        if let Some(span) = alone {
            let text = self.take_line();
            self.offer_newick(text, span);
        }
        while let Some(end) = self.pending.find(';') {
            let rest = self.pending.split_off(end + 1);
            let tree = std::mem::replace(&mut self.pending, rest.trim_start().to_string());
            let span = 0..tree.len();
            self.offer_newick(tree, span);
        }
    }

    /// Queue `text[span]` as a Newick tree read from `pending_line`, if it is
    /// one: a chunk that does not open with `(` is dropped, which keeps a
    /// file of neither format from coming back as one nonsense tree.
    fn offer_newick(&mut self, text: String, span: Range<usize>) {
        if text[span.clone()].starts_with('(') {
            let name = format!("{}_line{}", self.base_name, self.pending_line);
            self.offer(TreeRecord {
                name,
                state: 0,
                text,
                span,
            });
        }
    }

    /// Queue `tree` unless burn-in drops it.
    fn offer(&mut self, tree: TreeRecord) {
        if keep_tree(self.seen, tree.state, self.burnin_trees, self.burnin_states) {
            self.ready.push_back(tree);
        }
        self.seen += 1;
    }

    /// End of input: queue a trailing Newick tree that never met its `;`,
    /// leaving the parser to complain, and warn when there were no trees.
    fn finish(&mut self) {
        self.done = true;
        if self.format == TreeFormat::Newick {
            let tree = std::mem::take(&mut self.pending);
            let span = 0..tree.len();
            self.offer_newick(tree, span);
        }
        if self.seen == 0 && self.has_content {
            let missing = match self.format {
                TreeFormat::Nexus => "no NEXUS `tree ... = ...` lines",
                TreeFormat::Newick => "no `;`-terminated Newick trees",
            };
            eprintln!("No trees found in {}: {missing}", self.label);
        }
    }
}

impl<R: BufRead> Iterator for TreeReader<R> {
    type Item = io::Result<TreeRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.ready.is_empty() && !self.done {
            if let Err(e) = self.advance() {
                self.done = true;
                return Some(Err(e));
            }
        }
        self.ready.pop_front().map(Ok)
    }
}

/// A reader's trees, read on a background thread about [`READ_AHEAD_BYTES`]
/// ahead of the consumer, so that reading the file overlaps parsing it.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
pub(super) fn read_ahead<R>(
    mut reader: TreeReader<R>,
) -> Box<dyn Iterator<Item = io::Result<TreeRecord>>>
where
    R: BufRead + Send + 'static,
{
    let Some(first) = reader.next() else {
        return Box::new(std::iter::empty());
    };
    // The trees of a run are about the same length, so the first sizes the queue.
    let bytes = first.as_ref().map_or(0, |tree| tree.text.capacity());
    let depth = (READ_AHEAD_BYTES / bytes.max(1)).clamp(2, 1024);
    let (tx, rx) = std::sync::mpsc::sync_channel(depth);
    // The thread stops at the end of the file, or as soon as the consumer
    // hangs up and a send fails.
    let mut thread = Some(std::thread::spawn(move || {
        reader.try_for_each(|tree| tx.send(tree))
    }));
    // A panic drops the sender just as the end of the file does, so once the
    // queue is drained the thread is joined, and a panic becomes an error
    // instead of a quietly shortened file.
    let panicked = std::iter::from_fn(move || {
        thread.take()?.join().err()?;
        Some(Err(io::Error::other("the tree reader thread panicked")))
    });
    Box::new(std::iter::once(first).chain(rx).chain(panicked))
}

/// Without the `parallel` feature trees are read inline, and so on wasm32,
/// which has no `std::thread::spawn` even where `parallel` is on.
#[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
pub(super) fn read_ahead<R>(
    reader: TreeReader<R>,
) -> Box<dyn Iterator<Item = io::Result<TreeRecord>>>
where
    R: BufRead + Send + 'static,
{
    Box::new(reader)
}

/// Burn-in predicate shared by both formats: keep tree `idx` (0-based), sampled
/// at MCMC `state`, unless a threshold excludes it.
#[inline]
fn keep_tree(idx: usize, state: usize, burnin_trees: usize, burnin_states: usize) -> bool {
    (burnin_trees == 0 && burnin_states == 0)
        || (burnin_trees > 0 && idx >= burnin_trees)
        || (burnin_states > 0 && state > burnin_states)
}

/// Every kept tree of a file as `(name, newick)`, with its `TRANSLATE` map:
/// the reader, collected. A file that cannot be opened gives nothing.
#[cfg(test)]
fn load_beast_raw<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> (HashMap<String, String>, Vec<(String, String)>) {
    let Ok((translate, reader)) =
        TreeReader::open(path.as_ref(), burnin_trees, burnin_states, use_real_taxa)
    else {
        return (HashMap::new(), Vec::new());
    };
    let trees = reader
        .map_while(Result::ok)
        .map(|tree| (tree.name.clone(), tree.newick().to_string()))
        .collect();
    (translate, trees)
}

/// The trees of a plain Newick text, whatever its first line looks like.
#[cfg(test)]
pub(super) fn collect_newick_trees(content: &str, base_name: &str) -> Vec<TreeRecord> {
    let format = Some(TreeFormat::Newick);
    TreeReader::new(content.as_bytes(), base_name, "input", format, 0, 0, false)
        .map(|(_, reader)| reader.map_while(Result::ok).collect())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{hiv2_newick_path, hiv2_path};

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

    // ── TreeReader ────────────────────────────────────────────────────────────

    /// Read the header and trees of `content` as a file named `f` would be.
    fn reader(content: &str, use_real_taxa: bool) -> (HashMap<String, String>, TreeReader<&[u8]>) {
        TreeReader::new(content.as_bytes(), "f", "input", None, 0, 0, use_real_taxa).unwrap()
    }

    #[test]
    fn test_reader_hands_nexus_trees_over_in_their_line() {
        let (_, mut trees) = reader(
            "#NEXUS\nBegin trees;\ntree STATE_5 = [&R] (A:1,B:1,C:1);\nEnd;\n",
            false,
        );
        let tree = trees.next().unwrap().unwrap();
        assert_eq!(tree.name, "f_STATE_5");
        assert_eq!(tree.state, 5);
        assert_eq!(tree.newick(), "[&R] (A:1,B:1,C:1);");
        // The record keeps the line it was read in: nothing was copied out.
        assert_eq!(tree.text, "tree STATE_5 = [&R] (A:1,B:1,C:1);\n");
        assert!(trees.next().is_none());
    }

    #[test]
    fn test_reader_hands_a_lone_newick_tree_over_in_its_line() {
        let (_, mut trees) = reader("  (A:1,B:1,C:1);  \n", false);
        let tree = trees.next().unwrap().unwrap();
        assert_eq!(tree.name, "f_line1");
        assert_eq!(tree.newick(), "(A:1,B:1,C:1);");
        assert_eq!(tree.text, "  (A:1,B:1,C:1);  \n");
    }

    #[test]
    fn test_reader_translate_block_cut_short_by_end_of_file() {
        let (translate, mut trees) =
            reader("#NEXUS\nBegin trees;\n\tTranslate\n\t\t1 A,\n\t\t2 B", true);
        assert_eq!(translate.len(), 2);
        assert_eq!(translate.get("2").map(String::as_str), Some("B"));
        assert!(trees.next().is_none());
    }

    /// Input that serves `data` and then panics, as a bug in the reader would.
    #[cfg(feature = "parallel")]
    struct PanicsAfter(&'static [u8]);

    #[cfg(feature = "parallel")]
    impl io::Read for PanicsAfter {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            assert!(!self.0.is_empty(), "reader bug");
            let n = buf.len().min(self.0.len());
            buf[..n].copy_from_slice(&self.0[..n]);
            self.0 = &self.0[n..];
            Ok(n)
        }
    }

    /// A panic on the read-ahead thread ends the trees with an error, not with
    /// what looks like the end of the file.
    #[cfg(feature = "parallel")]
    #[test]
    fn test_read_ahead_reports_a_reader_panic() {
        let input = BufReader::new(PanicsAfter(b"(A:1,B:1,C:1);\n(A:2,B:2,C:2);\n"));
        let (_, trees) = TreeReader::new(input, "f", "input", None, 0, 0, false).unwrap();
        let mut trees: Vec<_> = read_ahead(trees).collect();
        assert!(trees.pop().is_some_and(|last| last.is_err()));
        assert!(trees.iter().all(Result::is_ok));
    }

    // ── collect_newick_trees ──────────────────────────────────────────────────

    /// `(line number, newick)` for each tree, which is all these tests assert on.
    fn newick_blocks(content: &str) -> Vec<(String, String)> {
        collect_newick_trees(content, "f")
            .into_iter()
            .map(|tree| (tree.name.clone(), tree.newick().to_string()))
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
}
