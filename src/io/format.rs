//! Telling NEXUS from plain Newick.

use super::nexus::{nexus_tree_lines, starts_with_ci};

/// Tree-file formats [`load_beast_trees`](super::load_beast_trees) can read.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::reader::collect_newick_trees;
    use crate::io::{hiv2_newick_path, hiv2_path};

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
}
