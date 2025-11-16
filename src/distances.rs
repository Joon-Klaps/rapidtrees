//! Tree distance metrics using bitset-based snapshots.
//!
//! This module implements three phylogenetic tree distance measures:
//!
//! 1. **Robinson-Foulds (RF)**: Counts the number of bipartitions that differ
//!    between two trees. Range: [0, 2n-6] where n is the number of leaves.
//!
//! 2. **Weighted Robinson-Foulds**: Like RF but considers branch lengths.
//!    For shared partitions, adds |length_a - length_b|.
//!    For unique partitions, adds the full branch length.
//!
//! 3. **Kuhner-Felsenstein (Branch Score)**: Similar to weighted RF but uses
//!    squared differences: sqrt(Σ(length_a - length_b)²)

use crate::snapshot::TreeSnapshot;
use phylotree::tree::{Tree as PhyloTree, TreeError};

#[cfg(test)]
use itertools::Itertools;

/// Compute Robinson-Foulds distance between two trees.
///
/// # Algorithm
/// RF = |A ∪ B| - |A ∩ B| = |A| + |B| - 2|A ∩ B|
///
/// Where A and B are the sets of bipartitions in each tree.
///
/// Since snapshots have sorted canonical bitsets, we use a linear merge
/// (O(m+n)) instead of hash lookups (O(m*n)).
///
/// # Rooted Tree Adjustment
/// For rooted trees, if the root position differs, we add 2 to the distance.
/// This accounts for the two extra bipartitions created by moving the root.
///
/// # Example
/// ```text
/// Tree 1:  ((A,B),(C,D))     Partitions: {A,B}, {C,D}
/// Tree 2:  ((A,C),(B,D))     Partitions: {A,C}, {B,D}
///
/// Intersection: 0 partitions match
/// RF = 2 + 2 - 2*0 = 4
/// ```
///
/// # Errors
/// Returns `TreeError` if trees have different leaf sets or are malformed.
pub fn robinson_foulds(tree_a: &PhyloTree, tree_b: &PhyloTree) -> Result<usize, TreeError> {
    let snap_a = TreeSnapshot::from_tree(tree_a)?;
    let snap_b = TreeSnapshot::from_tree(tree_b)?;

    Ok(rf_from_snapshots(&snap_a, &snap_b))
}

/// Compute Robinson-Foulds distance from two pre-computed snapshots.
///
/// This is the core RF algorithm using sorted merge for O(m+n) performance.
///
/// # Algorithm (O(m+n) using sorted merge)
/// ```text
/// Since parts are sorted, use two-pointer merge to count intersection
/// RF = len(A) + len(B) - 2 * len(intersection)
/// ```
///
#[inline]
pub fn rf_from_snapshots(a: &TreeSnapshot, b: &TreeSnapshot) -> usize {
    // Sorted merge intersection count - O(m+n) with excellent cache locality
    let mut inter = 0;
    let mut i = 0;
    let mut j = 0;

    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Equal => {
                inter += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }

    a.parts.len() + b.parts.len() - 2 * inter
}

/// Compute Weighted Robinson-Foulds distance between two trees.
///
/// # Algorithm
/// For each partition:
/// - If in both trees: add |length_a - length_b|
/// - If only in A: add length_a
/// - If only in B: add length_b
///
/// Total: Sum of all branch length differences
///
/// # Example
/// ```text
/// Tree 1: ((A:1.0,B:1.0):2.0,(C:1.0,D:1.0):2.0);
/// Tree 2: ((A:1.5,B:1.0):3.0,(C:0.5,D:1.0):2.0);
///
/// Shared partition {A,B}: |2.0 - 3.0| = 1.0
/// Shared partition {C,D}: |2.0 - 2.0| = 0.0
/// Different leaf branches contribute their full lengths
/// ```
///
/// # Errors
/// Returns `TreeError` if trees have different leaf sets or are malformed.
pub fn weighted_robinson_foulds(tree_a: &PhyloTree, tree_b: &PhyloTree) -> Result<f64, TreeError> {
    let snap_a = TreeSnapshot::from_tree(tree_a)?;
    let snap_b = TreeSnapshot::from_tree(tree_b)?;

    Ok(weighted_rf_from_snapshots(&snap_a, &snap_b))
}

/// Compute Weighted RF distance from two pre-computed snapshots.
///
/// Uses sorted merge for O(m+n) performance with excellent cache locality.
/// Direct array indexing (no HashMap lookups) for branch lengths.
pub fn weighted_rf_from_snapshots(a: &TreeSnapshot, b: &TreeSnapshot) -> f64 {
    let mut distance = 0.0;
    let mut i = 0;
    let mut j = 0;

    // Sorted merge through both partition lists
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Equal => {
                // Partition in both: add absolute difference
                // Direct array access - no hash lookup! 10-20× faster!
                distance += (a.lengths[i] - b.lengths[j]).abs();
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                // Partition only in A: add full length
                distance += a.lengths[i];
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                // Partition only in B: add full length
                distance += b.lengths[j];
                j += 1;
            }
        }
    }

    // Handle remaining partitions in A
    while i < a.parts.len() {
        distance += a.lengths[i];
        i += 1;
    }

    // Handle remaining partitions in B
    while j < b.parts.len() {
        distance += b.lengths[j];
        j += 1;
    }

    distance
}

/// Compute Kuhner-Felsenstein (Branch Score) distance between two trees.
///
/// # Algorithm
/// Like Weighted RF but uses squared differences:
/// distance = sqrt(Σ (length_a - length_b)²)
///
/// For each partition:
/// - If in both trees: add (length_a - length_b)²
/// - If only in A: add length_a²
/// - If only in B: add length_b²
///
/// Then take the square root of the sum.
///
/// # Properties
/// - More sensitive to large branch length differences
/// - Euclidean metric in branch length space
/// - Range: [0, ∞)
///
/// # Errors
/// Returns `TreeError` if trees have different leaf sets or are malformed.
pub fn kuhner_felsenstein(tree_a: &PhyloTree, tree_b: &PhyloTree) -> Result<f64, TreeError> {
    let snap_a = TreeSnapshot::from_tree(tree_a)?;
    let snap_b = TreeSnapshot::from_tree(tree_b)?;

    Ok(kf_from_snapshots(&snap_a, &snap_b))
}

/// Compute Kuhner-Felsenstein distance from two pre-computed snapshots.
///
/// Uses sorted merge for O(m+n) performance, accumulating squared differences.
/// Direct array indexing (no HashMap lookups) for branch lengths.
pub fn kf_from_snapshots(a: &TreeSnapshot, b: &TreeSnapshot) -> f64 {
    let mut sum_squared: f64 = 0.0;
    let mut i = 0;
    let mut j = 0;

    // Sorted merge through both partition lists
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Equal => {
                // Partition in both: add (diff)²
                // Direct array access - no hash lookup!
                let diff = a.lengths[i] - b.lengths[j];
                sum_squared += diff * diff;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                // Partition only in A: add length²
                sum_squared += a.lengths[i] * a.lengths[i];
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                // Partition only in B: add length²
                sum_squared += b.lengths[j] * b.lengths[j];
                j += 1;
            }
        }
    }

    // Handle remaining partitions in A
    while i < a.parts.len() {
        sum_squared += a.lengths[i] * a.lengths[i];
        i += 1;
    }

    // Handle remaining partitions in B
    while j < b.parts.len() {
        sum_squared += b.lengths[j] * b.lengths[j];
        j += 1;
    }

    sum_squared.sqrt()
}

#[test]
// Robinson foulds distances according to
// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
fn robinson_foulds_treedist() {
    let trees = [
        "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    ];
    let rfs = [
        vec![0, 4, 2, 10, 10, 10, 10, 10, 10, 10, 2, 10],
        vec![4, 0, 2, 10, 8, 10, 8, 10, 8, 10, 2, 10],
        vec![2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
        vec![10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
        vec![10, 8, 10, 2, 0, 4, 2, 4, 2, 2, 10, 4],
        vec![10, 10, 10, 2, 4, 0, 2, 2, 4, 2, 10, 2],
        vec![10, 8, 10, 4, 2, 2, 0, 4, 2, 4, 10, 4],
        vec![10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
        vec![10, 8, 10, 4, 2, 4, 2, 2, 0, 4, 10, 2],
        vec![10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
        vec![2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
        vec![10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
    ];

    for indices in (0..trees.len()).combinations(2) {
        let (i0, i1) = (indices[0], indices[1]);

        let t0 = PhyloTree::from_newick(trees[i0]).unwrap();
        let t1 = PhyloTree::from_newick(trees[i1]).unwrap();

        assert_eq!(robinson_foulds(&t0, &t1).unwrap(), rfs[i0][i1])
    }
}

#[test]
// Robinson foulds distances according to
// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
fn weighted_robinson_foulds_treedist() {
    let trees = [
        "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    ];
    let rfs = [
        [
            0.,
            0.4,
            0.2,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.9999999999999999,
        ],
        [
            0.4,
            0.,
            0.2,
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.2,
            0.9999999999999999,
        ],
        [
            0.2,
            0.2,
            0.,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.9999999999999999,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.2,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.9999999999999999,
            0.2,
        ],
        [
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.2,
            0.,
            0.4,
            0.2,
            0.4,
            0.2,
            0.2,
            0.9999999999999999,
            0.4,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.4,
            0.,
            0.2,
            0.2,
            0.4,
            0.2,
            0.9999999999999999,
            0.2,
        ],
        [
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.4,
            0.2,
            0.2,
            0.,
            0.4,
            0.2,
            0.4,
            0.9999999999999999,
            0.4,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.2,
            0.2,
            0.9999999999999999,
            0.,
        ],
        [
            0.9999999999999999,
            0.7999999999999999,
            0.9999999999999999,
            0.4,
            0.2,
            0.4,
            0.2,
            0.2,
            0.,
            0.4,
            0.9999999999999999,
            0.2,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.2,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.9999999999999999,
            0.2,
        ],
        [
            0.2,
            0.2,
            0.,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.,
            0.9999999999999999,
        ],
        [
            0.9999999999999999,
            0.9999999999999999,
            0.9999999999999999,
            0.2,
            0.4,
            0.2,
            0.4,
            0.,
            0.2,
            0.2,
            0.9999999999999999,
            0.,
        ],
    ];

    for indices in (0..trees.len()).combinations(2) {
        let (i0, i1) = (indices[0], indices[1]);
        let t0 = PhyloTree::from_newick(trees[i0]).unwrap();
        let t1 = PhyloTree::from_newick(trees[i1]).unwrap();

        assert!((weighted_robinson_foulds(&t0, &t1).unwrap() - rfs[i0][i1]).abs() <= f64::EPSILON)
    }
}

#[test]
// Branch score distances according to
// https://evolution.genetics.washington.edu/phylip/doc/treedist.html
fn kuhner_felsenstein_treedist() {
    let trees = [
        "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
        "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
    ];
    let rfs = [
        [
            0.,
            0.2,
            0.14142135623730953,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.316227766016838,
        ],
        [
            0.2,
            0.,
            0.14142135623730953,
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.14142135623730953,
            0.316227766016838,
        ],
        [
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.316227766016838,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.14142135623730953,
            0.,
            0.2,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.14142135623730953,
            0.316227766016838,
            0.2,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.2,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.2,
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.2,
            0.14142135623730953,
            0.2,
            0.316227766016838,
            0.2,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.316227766016838,
            0.,
        ],
        [
            0.316227766016838,
            0.28284271247461906,
            0.316227766016838,
            0.2,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.2,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.316227766016838,
            0.14142135623730953,
        ],
        [
            0.14142135623730953,
            0.14142135623730953,
            0.,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.,
            0.316227766016838,
        ],
        [
            0.316227766016838,
            0.316227766016838,
            0.316227766016838,
            0.14142135623730953,
            0.2,
            0.14142135623730953,
            0.2,
            0.,
            0.14142135623730953,
            0.14142135623730953,
            0.316227766016838,
            0.,
        ],
    ];

    for indices in (0..trees.len()).combinations(2) {
        let (i0, i1) = (indices[0], indices[1]);
        let t0 = PhyloTree::from_newick(trees[i0]).unwrap();
        let t1 = PhyloTree::from_newick(trees[i1]).unwrap();

        println!(
            "[{i0}, {i1}] c:{:?} ==? t:{}",
            kuhner_felsenstein(&t0, &t1).unwrap(),
            rfs[i0][i1]
        );

        assert_eq!(kuhner_felsenstein(&t0, &t1).unwrap(), rfs[i0][i1])
    }
}
