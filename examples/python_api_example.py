"""
Example usage of rapidtrees Python API
"""

import rapidtrees

# Example: Compute Robinson-Foulds distances from multiple tree files
def example_rf():
    # List of tree files to compare
    tree_files = [
        "tests/data/hiv1.trees",
        "tests/data/hiv2.trees",
    ]

    # Compute pairwise RF distances
    # - burnin_trees: skip first N trees from each file
    # - burnin_states: skip trees with STATE < N
    # - use_real_taxa: use TRANSLATE block when available
    tree_names, distance_matrix = rapidtrees.pairwise_rf(
        paths=tree_files,
        burnin_trees=10,
        burnin_states=0,
        use_real_taxa=True
    )

    print(f"Computed RF distances for {len(tree_names)} trees")
    print(f"Distance matrix shape: {len(distance_matrix)}x{len(distance_matrix[0])}")

    # Print first few tree names
    print("\nFirst 5 trees:")
    for name in tree_names[:5]:
        print(f"  - {name}")

    # Print a sample of the distance matrix
    print("\nSample distances (first 3x3):")
    for i in range(min(3, len(distance_matrix))):
        row = [distance_matrix[i][j] for j in range(min(3, len(distance_matrix[i])))]
        print(f"  {row}")

    return tree_names, distance_matrix


# Example: Compute Weighted Robinson-Foulds distances
def example_weighted_rf():
    tree_files = ["tests/data/hiv1.trees"]

    tree_names, distance_matrix = rapidtrees.pairwise_weighted_rf(
        paths=tree_files,
        burnin_trees=5,
    )

    print(f"\nWeighted RF distances for {len(tree_names)} trees")
    return tree_names, distance_matrix


# Example: Compute Kuhner-Felsenstein distances
def example_kf():
    tree_files = ["tests/data/hiv1.trees"]

    tree_names, distance_matrix = rapidtrees.pairwise_kf(
        paths=tree_files,
        burnin_trees=5,
    )

    print(f"\nKuhner-Felsenstein distances for {len(tree_names)} trees")
    return tree_names, distance_matrix


if __name__ == "__main__":
    print("=== Robinson-Foulds Distance ===")
    example_rf()

    print("\n=== Weighted Robinson-Foulds Distance ===")
    example_weighted_rf()

    print("\n=== Kuhner-Felsenstein Distance ===")
    example_kf()
