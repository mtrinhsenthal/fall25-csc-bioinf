# Converted for Codon
# Original: biotite.sequence.phylo.nj
# License: 3-Clause BSD License

__name__ = "biotite.sequence.phylo"
__author__ = "Patrick Kunzmann"
__all__ = ["neighbor_joining"]

import numpy as np
from .tree import Tree, TreeNode


MAX_FLOAT: float = np.finfo(np.float64).max


def neighbor_joining(distances: np.ndarray) -> Tree:
    """
    Perform hierarchical clustering using the Neighbor Joining algorithm.

    Unlike UPGMA, this method does not assume a constant evolution rate.
    The resulting tree is unrooted, except that the returned representation
    has a root node with three children.

    Parameters
    ----------
    distances : ndarray, shape (n, n)
        Pairwise distance matrix.

    Returns
    -------
    tree : Tree
        A rooted tree. The `index` attribute in each leaf node corresponds
        to the index of `distances`.

    Raises
    ------
    ValueError
        If the matrix is asymmetric, contains invalid values, or
        has fewer than 4 taxa.
    """

    n = distances.shape[0]
    if distances.shape[1] != n or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if n < 4:
        raise ValueError("At least 4 nodes are required")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    # Initialize leaves
    nodes: list[TreeNode | None] = [TreeNode(index=i) for i in range(n)]
    is_clustered = np.full(n, False, dtype=bool)
    divergence = np.zeros(n, dtype=np.float64)
    corr_distances = np.zeros((n, n), dtype=np.float64)
    distances_v = distances.astype(np.float64, copy=True)
    n_rem_nodes = n

    while True:
        # Calculate divergence for each unclustered node
        for i in range(n):
            if is_clustered[i]:
                continue
            dist_sum = 0.0
            for k in range(n):
                if not is_clustered[k]:
                    dist_sum += distances_v[i, k]
            divergence[i] = dist_sum

        # Compute corrected distance matrix
        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                corr_distances[i, j] = (
                    (n_rem_nodes - 2) * distances_v[i, j]
                    - divergence[i]
                    - divergence[j]
                )

        # Find pair with smallest corrected distance
        dist_min = MAX_FLOAT
        i_min, j_min = -1, -1
        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist = corr_distances[i, j]
                if dist < dist_min:
                    dist_min, i_min, j_min = dist, i, j

        # Exit if all nodes are clustered
        if i_min == -1 or j_min == -1:
            break

        # Compute branch lengths for the new node
        node_dist_i = 0.5 * (
            distances_v[i_min, j_min]
            + (divergence[i_min] - divergence[j_min]) / (n_rem_nodes - 2)
        )
        node_dist_j = 0.5 * (
            distances_v[i_min, j_min]
            + (divergence[j_min] - divergence[i_min]) / (n_rem_nodes - 2)
        )

        if n_rem_nodes > 3:
            # Cluster these two into a new node
            left, right = nodes[i_min], nodes[j_min]
            if left is None or right is None:
                raise RuntimeError("Attempted to cluster missing node")
            nodes[i_min] = TreeNode(
                children=(left, right),
                distances=(node_dist_i, node_dist_j),
            )
            nodes[j_min] = None
            is_clustered[j_min] = True
        else:
            # Final step: combine last three nodes into root
            is_clustered[i_min] = True
            is_clustered[j_min] = True
            k = np.where(~is_clustered)[0][0]
            node_dist_k = 0.5 * (
                distances_v[i_min, k]
                + distances_v[j_min, k]
                - distances_v[i_min, j_min]
            )
            root = TreeNode(
                children=(nodes[i_min], nodes[j_min], nodes[k]),
                distances=(node_dist_i, node_dist_j, node_dist_k),
            )
            return Tree(root)

        # Update distance matrix for new cluster
        for k in range(n):
            if not is_clustered[k] and k != i_min:
                new_dist = 0.5 * (
                    distances_v[i_min, k]
                    + distances_v[j_min, k]
                    - distances_v[i_min, j_min]
                )
                distances_v[i_min, k] = new_dist
                distances_v[k, i_min] = new_dist

        # Update remaining node count
        n_rem_nodes = n - np.count_nonzero(is_clustered)
