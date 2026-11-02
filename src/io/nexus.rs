//! NEXUS primitives: tree headers, tree lines and the `TRANSLATE` block.
//!
//! Readers outside this crate build on `extract_name_state`, `nexus_tree_lines`
//! and `parse_taxon_block`, so treat their signatures as published surface.

use std::collections::HashMap;

/// Split a NEXUS tree header into `(tree_name, state_number)`.
///
/// The name is the first word after `tree`, past the `*` that marks a default
/// tree (`TREE * STATE_0`). The state is the number after `STATE_`, and `0` for
/// a header without one, such as MrBayes's `tree gen.100`, which keeps its name.
///
/// ```
/// # use rapidtrees::io::extract_name_state;
/// let beast = extract_name_state("tree STATE_100 [&lnP=-12.3]");
/// assert_eq!(beast, ("STATE_100".to_string(), 100));
/// let mrbayes = extract_name_state("tree gen.100");
/// assert_eq!(mrbayes, ("gen.100".to_string(), 0));
/// ```
#[must_use]
pub fn extract_name_state(header: &str) -> (String, usize) {
    let (name, state) = name_state(header);
    (name.to_string(), state)
}

/// [`extract_name_state`] with the name borrowed from `header`.
pub(super) fn name_state(header: &str) -> (&str, usize) {
    let name = header
        .split_whitespace()
        .skip(1)
        .find(|&word| word != "*")
        .unwrap_or("");
    // Every byte of a match is ASCII, so the slice after it starts on a char.
    let state: usize = header
        .as_bytes()
        .windows(6)
        .position(|window| window.eq_ignore_ascii_case(b"STATE_"))
        .and_then(|pos| {
            header[pos + 6..]
                .split(|c: char| !c.is_ascii_digit())
                .next()
        })
        .and_then(|digits| digits.parse().ok())
        .unwrap_or(0);
    (name, state)
}

/// Case-insensitive ASCII prefix test that allocates no uppercased copy.
#[inline]
pub(super) fn starts_with_ci(s: &str, prefix: &str) -> bool {
    s.get(..prefix.len())
        .is_some_and(|p| p.eq_ignore_ascii_case(prefix))
}

/// A NEXUS tree line, `tree <name> = <newick>`, as its trimmed header and
/// body, or `None` for any other line. Splits on `" = "`, not `'='`: BEAST
/// headers carry `[&lnP=…]`.
#[inline]
pub(super) fn split_tree_line(line: &str) -> Option<(&str, &str)> {
    let (header, body) = line.split_once(" = ")?;
    Some((header.trim(), body.trim()))
}

/// Every `tree ... = ...` line of a NEXUS trees block, as `(header, body)`.
///
/// For a reader that wants the tree lines themselves, borrowed and unparsed,
/// to pick its own trees from: [`extract_name_state`] reads a header, and
/// [`read_trees`](super::read_trees) is the reader for everything else.
/// Lazy, so [`detect_format`](super::detect_format) can stop at the first hit
/// instead of collecting the whole block just to ask whether one exists.
///
/// ```
/// # use rapidtrees::io::nexus_tree_lines;
/// let nexus = "#NEXUS\nbegin trees;\n  tree STATE_0 = (A,B,C);\nend;\n";
/// let lines: Vec<_> = nexus_tree_lines(nexus).collect();
/// assert_eq!(lines, [("tree STATE_0", "(A,B,C);")]);
/// ```
pub fn nexus_tree_lines(content: &str) -> impl Iterator<Item = (&str, &str)> {
    content
        .lines()
        .map(str::trim)
        .skip_while(|line| !starts_with_ci(line, "tree "))
        .take_while(|line| !starts_with_ci(line, "end;"))
        .filter_map(split_tree_line)
}

/// Read a NEXUS `TRANSLATE` block into a numeric-ID → taxon-name map.
pub fn parse_taxon_block(content: &str) -> HashMap<String, String> {
    let mut translate = HashMap::new();
    let block = content
        .lines()
        .skip_while(|line| !starts_with_ci(line.trim(), "translate"))
        .skip(1);
    for line in block {
        translate.extend(translate_entry(line));
        if closes_translate(line) {
            break;
        }
    }
    translate
}

/// Whether a `TRANSLATE` line ends the block: a line opening with `;`, or the
/// last entry carrying the `;` itself, as MrBayes writes it.
pub(super) fn closes_translate(line: &str) -> bool {
    let line = line.trim();
    line.starts_with(';') || line.ends_with(';')
}

/// One `TRANSLATE` line, `<id> <label>,` or a last `<id> <label>;`, as
/// `(id, label)`, quotes dropped. A line opening with `;` holds none.
pub(super) fn translate_entry(line: &str) -> Option<(String, String)> {
    let line = line.trim();
    if line.starts_with(';') {
        return None;
    }
    let line = line.trim_end_matches(';').trim_end().trim_end_matches(',');
    let mut parts = line.split_whitespace();
    let id = parts.next()?.to_string();
    let label = parts.next()?.trim_matches('\'').to_string();
    Some((id, label))
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// A header without `STATE_` keeps its name: MrBayes and PAUP* write none.
    #[test]
    fn test_extract_name_state_without_state_keeps_the_name() {
        assert_eq!(extract_name_state("tree my_tree"), ("my_tree".into(), 0));
        assert_eq!(extract_name_state("tree gen.100"), ("gen.100".into(), 0));
        assert_eq!(
            extract_name_state("TREE * UNTITLED"),
            ("UNTITLED".into(), 0)
        );
    }

    #[test]
    fn test_extract_name_state_unparsable_state_keeps_the_name() {
        assert_eq!(extract_name_state("tree STATE_x"), ("STATE_x".into(), 0));
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

    /// MrBayes ends the block on the last entry: its `;` is not part of the
    /// label, and the tree lines after it are not entries.
    #[test]
    fn test_parse_taxon_block_semicolon_on_last_entry() {
        let content = "begin trees;\n   translate\n       1 A,\n       2 'B';\n   tree gen.0 = [&U] (1:1,2:1);\nend;\n";
        let map = parse_taxon_block(content);
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("2").map(String::as_str), Some("B"));
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
}
