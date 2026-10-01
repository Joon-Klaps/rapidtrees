//! NEXUS primitives: tree headers, tree lines and the `TRANSLATE` block.
//!
//! `extract_name_state` and `parse_taxon_block` are treetracer-web's API: change
//! their signatures only together with that repo.

use std::collections::HashMap;

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
pub(super) fn starts_with_ci(s: &str, prefix: &str) -> bool {
    s.get(..prefix.len())
        .is_some_and(|p| p.eq_ignore_ascii_case(prefix))
}

/// Every `tree ... = ...` line of a NEXUS trees block, as `(header, body)`.
///
/// Lazy, so [`detect_format`](super::detect_format) can stop at the first hit instead of collecting
/// the whole block just to ask whether one exists.
pub(super) fn nexus_tree_lines(content: &str) -> impl Iterator<Item = (&str, &str)> {
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

/// Read a NEXUS `TRANSLATE` block into a numeric-ID → taxon-name map.
pub fn parse_taxon_block(content: &str) -> HashMap<String, String> {
    content
        .lines()
        .skip_while(|line| !starts_with_ci(line.trim(), "translate"))
        .skip(1)
        .take_while(|line| !line.trim().starts_with(';'))
        .filter_map(translate_entry)
        .collect::<HashMap<_, _>>()
}

/// One `TRANSLATE` line, `<id> <label>,`, as `(id, label)`, quotes dropped.
pub(super) fn translate_entry(line: &str) -> Option<(String, String)> {
    let line = line.trim().trim_end_matches(',');
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
}
