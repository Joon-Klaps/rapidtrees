//! BEAST `[&...]` annotation stripping.

use std::borrow::Cow;

/// Strip BEAST `[&...]` annotations from a Newick string.
pub fn strip_beast_annotations(newick: &str) -> String {
    strip_annotations(newick).into_owned()
}

/// [`strip_beast_annotations`] without the copy when there is nothing to strip,
/// which is every tree in a plain Newick file.
pub(super) fn strip_annotations(newick: &str) -> Cow<'_, str> {
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
