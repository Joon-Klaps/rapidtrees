//! Newick text straight to a [`Snapshot`], with no tree in between.
//!
//! phylotree's parser builds a whole `Tree` for every input: a `String` per name
//! and per branch length, and a `Node` with its own child vector per node, which
//! [`Snapshot`] construction then walked once and dropped. On a posterior that
//! round trip was about three quarters of construction time. This reader does
//! the same walk over the text itself. Leaves are numbered as they are read,
//! which is the pre-order that keeps every subtree a contiguous run of
//! `leaf_order`, and a `)` closes a node whose fingerprint is the XOR of its
//! children's. The XOR of every taxon label is known before the first tree, so
//! a split's canonical key is ready the moment its node closes.
//!
//! It accepts what phylotree's parser accepts, and reads it the same way:
//! whitespace outside double quotes is ignored, a double-quoted name keeps its
//! quotes, `[...]` comments (BEAST `[&...]` annotations included) end at the
//! first `]` and are skipped, labels on internal nodes are ignored, and a
//! missing branch length is 0.0.

use super::build::{Part, Snapshot};
use super::fingerprint::Fingerprint;
use rustc_hash::FxHashMap;
use std::borrow::Cow;
use std::collections::HashMap;

/// What every tree of a run shares: how a taxon name becomes a bit and a label.
pub(super) struct RunTables<'a> {
    pub(super) leaf_index: &'a FxHashMap<&'a str, usize>,
    pub(super) labels: &'a [Fingerprint],
    /// The XOR of every label, i.e. the fingerprint of the whole leaf set, so a
    /// split's complement is `fp ^ total`.
    pub(super) total: Fingerprint,
    pub(super) rooted: bool,
    /// Rooted-facts collection distinguishes an omitted edge length from an
    /// explicit zero; ordinary distance paths keep the historical default of
    /// treating an omitted length as zero.
    pub(super) require_explicit_lengths: bool,
}

/// Parse one tree into a [`Snapshot`], checking its leaf set against the run's.
///
/// `index` only feeds the error messages, so they point at the offending tree.
pub(super) fn snapshot(
    newick: &str,
    translate: &HashMap<String, String>,
    index: usize,
    run: &RunTables,
) -> Result<Snapshot, String> {
    let num_leaves = run.labels.len();
    let mismatch = || {
        format!(
            "Tree {index} has a different leaf set than tree 0. All trees must share the same taxa."
        )
    };

    // `seen` is `n` bytes and is what catches a name repeated within one tree:
    // a plain count would let `(A,A,B)` pass against `{A,B,C}`.
    let mut seen = vec![false; num_leaves];
    let shape = Shape {
        rooted: run.rooted,
        num_leaves,
        total: run.total,
        require_explicit_lengths: run.require_explicit_lengths,
    };
    let (leaf_order, parts) = walk(newick, index, shape, |raw| {
        let name = rename(raw, translate).ok_or_else(|| unnamed(index))?;
        let &bit = run.leaf_index.get(name).ok_or_else(mismatch)?;
        if std::mem::replace(&mut seen[bit], true) {
            return Err(format!(
                "Tree {index} repeats the leaf name {name:?}. All leaf names must be unique."
            ));
        }
        Ok((bit as u32, run.labels[bit]))
    })?;

    // Every name was known and none repeated, so a short count is a missing taxon.
    if leaf_order.len() != num_leaves {
        return Err(mismatch());
    }
    Ok(Snapshot {
        parts,
        leaf_order,
        words: num_leaves.div_ceil(64),
    })
}

/// Tree 0's leaf names, renamed through `translate`, in the order they appear.
///
/// Read on their own because they are what the run's tables are built from:
/// no tree, tree 0 included, can be fingerprinted until they exist.
pub(super) fn leaf_names(
    newick: &str,
    translate: &HashMap<String, String>,
) -> Result<Vec<String>, String> {
    let mut names = Vec::new();
    // Rooted and zero-width, so no split is filtered or canonicalised: the
    // parts are thrown away and only the names are wanted.
    let shape = Shape {
        rooted: true,
        num_leaves: 0,
        total: 0,
        require_explicit_lengths: false,
    };
    walk(newick, 0, shape, |raw| {
        names.push(rename(raw, translate).ok_or_else(|| unnamed(0))?.to_owned());
        Ok((0, 0))
    })?;
    Ok(names)
}

/// A leaf's taxon name: its label, or what `translate` maps it to when the tree
/// comes with a translate table. `None` is a label missing from the table, which
/// is reported as an unnamed leaf. `raw` is never empty: a leaf with no label
/// is caught by [`walk`] before it gets here.
fn rename<'a>(raw: &'a str, translate: &'a HashMap<String, String>) -> Option<&'a str> {
    if translate.is_empty() {
        Some(raw)
    } else {
        translate.get(raw).map(String::as_str)
    }
}

fn unnamed(index: usize) -> String {
    format!("Tree {index} has an unnamed leaf. All leaves must be named.")
}

/// How a closed node becomes a [`Part`].
#[derive(Clone, Copy)]
struct Shape {
    rooted: bool,
    num_leaves: usize,
    total: Fingerprint,
    require_explicit_lengths: bool,
}

/// A node whose subtree is complete and whose edge is still being read: a
/// label, comments and a branch length may follow before the next delimiter.
#[derive(Clone, Copy)]
struct Closed {
    fp: Fingerprint,
    first: u32,
    size: u32,
    length: Option<f64>,
}

/// An internal node still reading its children.
struct Open {
    fp: Fingerprint,
    first: u32,
    children: u32,
}

/// Walk one tree, resolving each leaf name through `leaf` to its bit and label,
/// and return `(leaf_order, parts)`.
///
/// Parts arrive in post-order, filtered and canonicalised as
/// [`Snapshot`] documents. Two nodes can only share a canonical key when one
/// sits above the other through unary nodes, or when they are the two children
/// of a bifurcating root, whose splits are each other's complement. The second
/// case is the usual rooted binary tree and is merged directly; the first is rare
/// and falls back to sorting the tree's parts and merging equal keys.
fn walk(
    newick: &str,
    index: usize,
    shape: Shape,
    mut leaf: impl FnMut(&str) -> Result<(u32, Fingerprint), String>,
) -> Result<(Vec<u32>, Vec<Part>), String> {
    let syntax = |what: &str| format!("Failed to parse newick at index {index}: {what}");
    let bytes = newick.as_bytes();

    let mut leaf_order: Vec<u32> = Vec::with_capacity(shape.num_leaves);
    let mut parts: Vec<Part> = Vec::with_capacity(2 * shape.num_leaves);
    let mut open: Vec<Open> = Vec::new();
    let mut closed: Option<Closed> = None;
    // The parts of the root's first two children, and how many children it had.
    let mut root_children: [Option<usize>; 2] = [None, None];
    let mut root_arity = 0;
    let mut unary = false;

    let mut i = 0;
    loop {
        let Some(&byte) = bytes.get(i) else {
            return Err(syntax("the tree does not end in ';'"));
        };
        match byte {
            b if b.is_ascii_whitespace() => i += 1,
            b'(' => {
                if closed.is_some() {
                    return Err(syntax("'(' after a finished subtree"));
                }
                open.push(Open {
                    fp: 0,
                    first: leaf_order.len() as u32,
                    children: 0,
                });
                i += 1;
            }
            b',' | b')' => {
                let at_root = open.len() == 1;
                let Some(parent) = open.last_mut() else {
                    return Err(syntax(&format!("'{}' outside any subtree", byte as char)));
                };
                // Nothing read since the last delimiter is an anonymous leaf.
                let child = closed.take().ok_or_else(|| unnamed(index))?;
                if shape.require_explicit_lengths && child.length.is_none() {
                    return Err(syntax(
                        "a non-root node is missing an explicit branch length",
                    ));
                }
                if shape.require_explicit_lengths
                    && child.length.is_some_and(|length| !length.is_finite())
                {
                    return Err(syntax("a non-root node has a non-finite branch length"));
                }
                parent.fp ^= child.fp;
                parent.children += 1;
                let children = parent.children;
                let part = emit(&mut parts, child, shape);
                if at_root && children <= 2 {
                    root_children[children as usize - 1] = part;
                }

                if byte == b')' {
                    let node = open.pop().expect("a parent was found above");
                    unary |= node.children == 1;
                    if open.is_empty() {
                        root_arity = node.children;
                    }
                    closed = Some(Closed {
                        fp: node.fp,
                        first: node.first,
                        size: leaf_order.len() as u32 - node.first,
                        length: None,
                    });
                }
                i += 1;
            }
            b':' => {
                // A length with nothing before it belongs to an anonymous leaf.
                let node = closed.as_mut().ok_or_else(|| unnamed(index))?;
                // BEAST writes its annotation between the colon and the number,
                // as in `13:[&rate=0.71]16.04`, so comments are skipped here too.
                let start = skip_blank(bytes, i + 1)
                    .ok_or_else(|| syntax("a '[' comment is never closed"))?;
                let (text, next) = token(newick, start);
                // An empty length is a missing one, which phylotree read as 0.0.
                node.length = Some(if text.is_empty() {
                    0.0
                } else {
                    text.parse::<f64>()
                        .map_err(|e| syntax(&format!("invalid branch length {text:?}: {e}")))?
                });
                i = next;
            }
            b'[' => {
                i = skip_blank(bytes, i).ok_or_else(|| syntax("a '[' comment is never closed"))?;
            }
            // A stray ']' is what phylotree tolerated too.
            b']' => i += 1,
            b';' => break,
            _ => {
                let (name, next) = token(newick, i);
                i = next;
                // After a node has closed, a name is an internal label or a
                // support value, which no metric reads.
                if closed.is_none() {
                    let (bit, label) = leaf(&name)?;
                    closed = Some(Closed {
                        fp: label,
                        first: leaf_order.len() as u32,
                        size: 1,
                        length: None,
                    });
                    leaf_order.push(bit);
                }
            }
        }
    }

    if !open.is_empty() {
        return Err(syntax("a '(' is never closed"));
    }
    if closed.is_none() {
        return Err(syntax("the tree is empty"));
    }

    if !shape.rooted {
        if unary {
            parts.sort_unstable_by_key(|p| p.key);
            parts.dedup_by(|later, kept| {
                later.key == kept.key && {
                    kept.length += later.length;
                    true
                }
            });
        } else if root_arity == 2
            && let [Some(a), Some(b)] = root_children
            && parts[a].key == parts[b].key
        {
            parts[a].length += parts[b].length;
            parts.remove(b);
        }
    }
    Ok((leaf_order, parts))
}

/// Append `node`'s part, unless it is the complement of a pendant edge, and
/// return its index.
fn emit(parts: &mut Vec<Part>, node: Closed, shape: Shape) -> Option<usize> {
    // A split of `n - 1` leaves is the other side of a pendant edge, which is
    // kept as the pendant itself. Rooted clades are all kept.
    if !shape.rooted && node.size != 1 && node.size as usize >= shape.num_leaves.saturating_sub(1) {
        return None;
    }
    // A pendant keeps its raw fingerprint: at two taxa `min(h, h ^ total)`
    // would collapse the two pendants onto each other.
    let key = if shape.rooted || node.size == 1 {
        node.fp
    } else {
        node.fp.min(node.fp ^ shape.total)
    };
    parts.push(Part {
        key,
        first: node.first,
        size: node.size,
        length: node.length.unwrap_or(0.0),
    });
    Some(parts.len() - 1)
}

/// Skip whitespace and `[...]` comments from `i`, returning the next byte to
/// read, or `None` for a comment that is never closed.
fn skip_blank(bytes: &[u8], mut i: usize) -> Option<usize> {
    loop {
        match bytes.get(i) {
            Some(b) if b.is_ascii_whitespace() => i += 1,
            Some(b'[') => i += bytes[i..].iter().position(|&b| b == b']')? + 1,
            _ => return Some(i),
        }
    }
}

/// One name or number: everything up to the next delimiter, read as phylotree
/// read it. Unquoted whitespace is dropped, and a double-quoted run is kept
/// whole, quotes included. Returns the text and where the delimiter sits.
fn token(newick: &str, start: usize) -> (Cow<'_, str>, usize) {
    let bytes = newick.as_bytes();
    let mut i = start;
    let mut plain = true;
    while let Some(&b) = bytes.get(i) {
        match b {
            b'(' | b')' | b',' | b':' | b';' | b'[' | b']' => break,
            b'"' => {
                plain = false;
                i += 1;
                i = match bytes[i..].iter().position(|&c| c == b'"') {
                    Some(end) => i + end + 1,
                    None => bytes.len(),
                };
            }
            b if b.is_ascii_whitespace() => {
                plain = false;
                i += 1;
            }
            _ => i += 1,
        }
    }

    let raw = &newick[start..i];
    let trimmed = raw.trim_end_matches(|c: char| c.is_ascii_whitespace());
    if plain
        || !trimmed
            .bytes()
            .any(|b| b == b'"' || b.is_ascii_whitespace())
    {
        return (Cow::Borrowed(trimmed), i);
    }
    // Rare: whitespace inside the token, or a quoted run.
    let mut text = String::with_capacity(raw.len());
    let mut quoted = false;
    for c in raw.chars() {
        if c == '"' {
            quoted = !quoted;
            text.push(c);
        } else if quoted || !c.is_ascii_whitespace() {
            text.push(c);
        }
    }
    (Cow::Owned(text), i)
}
