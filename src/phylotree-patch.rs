//! Minimal in-crate replacement for the slice of `phylotree` that `rapidtrees`
//! actually uses, for the `wasm32-unknown-unknown` build.
//!
//! Upstream `phylotree 0.1.3` pulls in `needletail` (C-backed compression),
//! `ptree` (`directories`, no wasm support) and `indicatif` as *unconditional*
//! dependencies, none of which compile for `wasm32-unknown-unknown`. The core
//! only ever needs a Newick parser that yields a tree we can walk, so this
//! module reproduces exactly that surface — nothing more.
//!
//! The public API here intentionally mirrors `phylotree::tree::{Tree, Node,
//! TreeError}` so that `snapshot.rs` and `io.rs` can swap backends with a single
//! `cfg`-gated type alias and keep `Snapshot::from_tree` / `compute_bitsets` /
//! `collect_partitions` byte-for-byte unchanged. A native test cross-checks that
//! this parser produces identical pairwise distance matrices to upstream
//! phylotree (see `snapshot.rs` tests).
//!
//! # Supported Newick
//! Standard Newick: nested `(...)`, optional internal labels, `:branch_length`
//! edges (including scientific notation), single-quoted labels with `''`
//! escaping, `[...]` comments, and a trailing `;`. Underscores are kept literal
//! (not converted to spaces) — the inputs `rapidtrees` handles use numeric IDs
//! or plain alphanumeric taxa, and the cross-check test guards equivalence.

use thiserror::Error;

/// Errors surfaced by the Newick parser and tree accessors.
///
/// Mirrors the subset of `phylotree::tree::TreeError` that `rapidtrees`
/// constructs or matches on (`NodeNotFound` is built in `snapshot.rs`).
#[derive(Debug, Error)]
pub enum TreeError {
    /// A node id was requested that does not exist in the arena.
    #[error("node {0} not found")]
    NodeNotFound(usize),
    /// The Newick string could not be parsed.
    #[error("failed to parse newick: {0}")]
    ParseError(String),
}

/// A single tree node in the arena.
///
/// Field names and types match the `phylotree` node members read by the core:
/// `name`, `children`, and `parent_edge` (the length of the edge above this
/// node, `None` for the root or unspecified edges).
#[derive(Debug, Clone, Default)]
pub struct Node {
    /// Taxon label for leaves; optional internal-node label otherwise.
    pub name: Option<String>,
    /// Arena indices of this node's children (empty for leaves).
    pub children: Vec<usize>,
    /// Length of the edge connecting this node to its parent, if given.
    pub parent_edge: Option<f64>,
}

/// An arena-backed phylogenetic tree parsed from a Newick string.
///
/// Nodes are stored in a flat `Vec`; a node id is its index. The root is always
/// id `0` (the first node created during parsing), but callers should use
/// [`Tree::get_root`] rather than assume that.
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    root: usize,
}

impl Tree {
    /// Parse a Newick string into a [`Tree`].
    ///
    /// BEAST `[&...]` annotations are tolerated as `[...]` comments, but callers
    /// in `rapidtrees` strip them beforehand via `io::strip_beast_annotations`.
    ///
    /// # Errors
    /// Returns [`TreeError::ParseError`] on malformed input (unbalanced
    /// parentheses, missing labels, or invalid branch lengths).
    pub fn from_newick(newick: &str) -> Result<Self, TreeError> {
        let mut parser = Parser::new(newick);
        let root = parser.parse_tree()?;
        Ok(Tree {
            nodes: parser.nodes,
            root,
        })
    }

    /// Borrow the node with id `node_id`.
    ///
    /// # Errors
    /// Returns [`TreeError::NodeNotFound`] if `node_id` is out of range.
    pub fn get(&self, node_id: &usize) -> Result<&Node, TreeError> {
        self.nodes
            .get(*node_id)
            .ok_or(TreeError::NodeNotFound(*node_id))
    }

    /// Mutably borrow the node with id `node_id` (used to apply translate maps).
    ///
    /// # Errors
    /// Returns [`TreeError::NodeNotFound`] if `node_id` is out of range.
    pub fn get_mut(&mut self, node_id: &usize) -> Result<&mut Node, TreeError> {
        self.nodes
            .get_mut(*node_id)
            .ok_or(TreeError::NodeNotFound(*node_id))
    }

    /// Return the id of the root node.
    ///
    /// # Errors
    /// Returns [`TreeError::ParseError`] if the tree is empty.
    pub fn get_root(&self) -> Result<usize, TreeError> {
        if self.nodes.is_empty() {
            return Err(TreeError::ParseError("empty tree".to_string()));
        }
        Ok(self.root)
    }

    /// Return the ids of all leaf nodes (those with no children).
    ///
    /// Order is arena order; callers in `rapidtrees` re-sort by taxon name, so
    /// the ordering here is not significant.
    pub fn get_leaves(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty())
            .map(|(id, _)| id)
            .collect()
    }

    /// Report whether the tree is rooted, matching phylotree's convention:
    /// a tree is rooted when its root has exactly two children.
    ///
    /// # Errors
    /// Returns [`TreeError::ParseError`] if the tree is empty.
    pub fn is_rooted(&self) -> Result<bool, TreeError> {
        let root = self.get_root()?;
        Ok(self.nodes[root].children.len() == 2)
    }
}

/// Recursive-descent Newick parser over a char buffer.
struct Parser {
    chars: Vec<char>,
    pos: usize,
    nodes: Vec<Node>,
}

impl Parser {
    fn new(s: &str) -> Self {
        Parser {
            chars: s.chars().collect(),
            pos: 0,
            nodes: Vec::new(),
        }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn bump(&mut self) -> Option<char> {
        let c = self.peek();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    /// Skip whitespace and `[...]` comment blocks.
    fn skip_trivia(&mut self) {
        loop {
            match self.peek() {
                Some(c) if c.is_whitespace() => {
                    self.bump();
                }
                Some('[') => {
                    while let Some(c) = self.bump() {
                        if c == ']' {
                            break;
                        }
                    }
                }
                _ => break,
            }
        }
    }

    /// Allocate a fresh default node and return its id.
    fn new_node(&mut self) -> usize {
        let id = self.nodes.len();
        self.nodes.push(Node::default());
        id
    }

    /// Parse a complete tree: a subtree, an optional root branch length, `;`.
    fn parse_tree(&mut self) -> Result<usize, TreeError> {
        self.skip_trivia();
        let root = self.parse_subtree()?;
        self.skip_trivia();
        // A root branch length is allowed by the grammar but unused by the core.
        if self.peek() == Some(':') {
            self.bump();
            let _ = self.parse_length()?;
        }
        self.skip_trivia();
        if self.peek() == Some(';') {
            self.bump();
        }
        Ok(root)
    }

    /// Parse a subtree (internal node `(...)` or a leaf) and return its id.
    /// The caller is responsible for the optional `:length` that follows.
    fn parse_subtree(&mut self) -> Result<usize, TreeError> {
        self.skip_trivia();
        let id = self.new_node();

        if self.peek() == Some('(') {
            self.bump(); // consume '('
            loop {
                let child = self.parse_subtree()?;
                self.skip_trivia();
                if self.peek() == Some(':') {
                    self.bump();
                    self.nodes[child].parent_edge = Some(self.parse_length()?);
                }
                self.nodes[id].children.push(child);
                self.skip_trivia();
                match self.peek() {
                    Some(',') => {
                        self.bump();
                    }
                    Some(')') => {
                        self.bump();
                        break;
                    }
                    other => {
                        return Err(TreeError::ParseError(format!(
                            "expected ',' or ')', found {other:?}"
                        )));
                    }
                }
            }
            // Optional internal-node label.
            let label = self.parse_label();
            if !label.is_empty() {
                self.nodes[id].name = Some(label);
            }
        } else {
            let label = self.parse_label();
            if label.is_empty() {
                return Err(TreeError::ParseError(
                    "expected leaf label or '('".to_string(),
                ));
            }
            self.nodes[id].name = Some(label);
        }

        Ok(id)
    }

    /// Parse a node label: either a single-quoted string (with `''` escaping) or
    /// a bare token terminated by Newick punctuation/whitespace.
    fn parse_label(&mut self) -> String {
        self.skip_trivia();
        if self.peek() == Some('\'') {
            self.bump(); // opening quote
            let mut s = String::new();
            while let Some(c) = self.bump() {
                if c == '\'' {
                    if self.peek() == Some('\'') {
                        self.bump();
                        s.push('\'');
                    } else {
                        break;
                    }
                } else {
                    s.push(c);
                }
            }
            s
        } else {
            let mut s = String::new();
            while let Some(c) = self.peek() {
                if matches!(c, '(' | ')' | ',' | ':' | ';' | '[') || c.is_whitespace() {
                    break;
                }
                s.push(c);
                self.bump();
            }
            s
        }
    }

    /// Parse a branch length (decimal, optionally signed, with exponent).
    fn parse_length(&mut self) -> Result<f64, TreeError> {
        self.skip_trivia();
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() || matches!(c, '.' | '-' | '+' | 'e' | 'E') {
                s.push(c);
                self.bump();
            } else {
                break;
            }
        }
        s.parse::<f64>()
            .map_err(|_| TreeError::ParseError(format!("invalid branch length: {s:?}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_simple_tree() {
        let t = Tree::from_newick("((A:1,B:2):3,(C:4,D:5):6);").unwrap();
        assert_eq!(t.get_leaves().len(), 4);
        assert!(t.is_rooted().unwrap());
    }

    #[test]
    fn parses_quoted_labels_and_comments() {
        let t = Tree::from_newick("('A b':1.0,[comment]C:2.0,D:3.0);").unwrap();
        let names: Vec<String> = t
            .get_leaves()
            .iter()
            .filter_map(|id| t.get(id).ok()?.name.clone())
            .collect();
        assert!(names.contains(&"A b".to_string()));
        assert!(names.contains(&"C".to_string()));
        // Trifurcating root → unrooted by phylotree's convention.
        assert!(!t.is_rooted().unwrap());
    }

    #[test]
    fn parses_scientific_branch_lengths() {
        let t = Tree::from_newick("(A:1.5E-3,B:2.0e2);").unwrap();
        let lens: Vec<f64> = t
            .get_leaves()
            .iter()
            .filter_map(|id| t.get(id).ok()?.parent_edge)
            .collect();
        assert!(lens.contains(&0.0015));
        assert!(lens.contains(&200.0));
    }

    #[test]
    fn rejects_unbalanced_parens() {
        assert!(Tree::from_newick("((A:1,B:2);").is_err());
    }

    #[test]
    fn get_mut_allows_renaming_a_leaf() {
        // Mirrors how `io::rename_leaf_nodes` applies a BEAST translate map.
        let mut t = Tree::from_newick("(A:1,B:1);").unwrap();
        let leaf = t.get_leaves()[0];
        t.get_mut(&leaf).unwrap().name = Some("Z".to_string());
        assert_eq!(t.get(&leaf).unwrap().name.as_deref(), Some("Z"));
    }

    // ── Cross-check against upstream phylotree ────────────────────────────────
    //
    // `Snapshot::from_tree` consumes a tree purely as: for every non-root node,
    // the set of leaf names beneath it plus the edge length above it. If this
    // parser and phylotree produce the same such signature for the same Newick,
    // they yield identical `Snapshots` and therefore identical distance
    // matrices. We compare signatures directly (rather than distances) because
    // the `PhyloTree` type alias is fixed per build, so we can't run `from_tree`
    // on both backends in one compilation.

    /// Canonical structural signature: sorted list of
    /// `(sorted leaf names under node, branch length)` for every non-root node.
    /// The macro works on either tree type since both expose the same API.
    macro_rules! signature {
        ($tree:expr) => {{
            use std::collections::{BTreeMap, BTreeSet};
            let tree = &$tree;
            let root = tree.get_root().unwrap();

            // Pre-order traversal collects every node id (parents before children).
            let mut order = vec![root];
            let mut stack = vec![root];
            while let Some(id) = stack.pop() {
                for &c in &tree.get(&id).unwrap().children {
                    order.push(c);
                    stack.push(c);
                }
            }

            // Reverse pre-order = children before parents → bottom-up leaf sets.
            let mut leafset: BTreeMap<usize, BTreeSet<String>> = BTreeMap::new();
            for &id in order.iter().rev() {
                let node = tree.get(&id).unwrap();
                let mut set = BTreeSet::new();
                if node.children.is_empty() {
                    set.insert(node.name.clone().unwrap_or_default());
                } else {
                    for &c in &node.children {
                        for name in &leafset[&c] {
                            set.insert(name.clone());
                        }
                    }
                }
                leafset.insert(id, set);
            }

            let mut sig: Vec<(Vec<String>, String)> = order
                .iter()
                .filter(|&&id| id != root)
                .map(|&id| {
                    let node = tree.get(&id).unwrap();
                    let names: Vec<String> = leafset[&id].iter().cloned().collect();
                    (names, format!("{:.6}", node.parent_edge.unwrap_or(0.0)))
                })
                .collect();
            sig.sort();
            sig
        }};
    }

    fn assert_matches_phylotree(newick: &str) {
        let ours = Tree::from_newick(newick).unwrap();
        let theirs = phylotree::tree::Tree::from_newick(newick).unwrap();
        assert_eq!(
            signature!(ours),
            signature!(theirs),
            "structural signature mismatch for newick: {newick}"
        );
    }

    #[test]
    fn matches_phylotree_on_balanced_rooted_tree() {
        assert_matches_phylotree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
    }

    #[test]
    fn matches_phylotree_on_unrooted_trifurcation() {
        assert_matches_phylotree("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);");
    }

    #[test]
    fn matches_phylotree_on_deep_caterpillar() {
        assert_matches_phylotree(
            "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        );
    }

    #[test]
    fn matches_phylotree_with_missing_branch_lengths() {
        assert_matches_phylotree("((A,B),(C,D));");
    }
}
