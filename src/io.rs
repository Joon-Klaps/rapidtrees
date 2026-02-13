use phylotree::tree::Tree;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;

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

pub fn read_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> (HashMap<String, String>, Vec<(String, Tree)>) {
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

    let trees = collect_tree_blocks(&content)
        .into_iter()
        .enumerate()
        //generate tree name & extract state number
        .map(|(idx, tree)| {
            let (name, state) = extract_name_state(tree.header);
            (idx, tree, state, format!("{base_name}_{name}"))
        })
        // Filter out burn-in trees based on count and/or state number if 0 we don't filter
        .filter(|(idx, _tree, state, _name)| {
            (burnin_trees == 0 && burnin_states == 0)
                || (burnin_trees > 0 && *idx >= burnin_trees)
                || (burnin_states > 0 && *state > burnin_states)
        })
        // read in the files
        .filter_map(|(idx, tree, _state, name)| {
            // Strip BEAST annotations from newick string (e.g., [&rate=...])
            // BEAST format: :[&rate=X.XX]length -> :length
            let newick = strip_beast_annotations(&tree.body);
            let mut phylo_tree = match phylotree::tree::Tree::from_newick(&newick) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!(
                        "Failed to parse tree {} at index {}: {}",
                        path.as_ref().display(),
                        idx,
                        e
                    );
                    return None;
                }
            };

            // Rename the leaves with the map
            if use_real_taxa {
                rename_leaf_nodes(&mut phylo_tree, &taxons);
            }

            Some((name, phylo_tree))
        })
        .collect::<Vec<_>>();

    (taxons, trees)
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
