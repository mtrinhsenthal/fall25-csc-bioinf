# # Converted for Codon
# # Original source: biotite.sequence.phylo.upgma
# # License: 3-Clause BSD License

# __name__ = "biotite.sequence.phylo"
# __author__ = "Patrick Kunzmann"
# __all__ = ["upgma"]

# import numpy as np
# from .tree import Tree, TreeNode
# from typing import Optional


# MAX_FLOAT: float = np.finfo(np.float64).max


# def upgma(distances: np.ndarray) -> Tree:
#     """
#     Perform hierarchical clustering using the
#     *Unweighted Pair Group Method with Arithmetic Mean* (UPGMA).

#     This algorithm produces leaf nodes with equal distance to the root node.
#     In phylogenetics, this corresponds to a constant evolution rate
#     (molecular clock).

#     Parameters
#     ----------
#     distances : ndarray of shape (n, n)
#         Pairwise distance matrix.

#     Returns
#     -------
#     tree : Tree
#         A rooted binary tree. Each leaf TreeNode's `index` attribute
#         corresponds to a row/column index in `distances`.

#     Raises
#     ------
#     ValueError
#         If the distance matrix is not symmetric,
#         contains NaN or infinity, or has negative values.
#     """

#     n = distances.shape[0]
#     if distances.shape[1] != n or not np.allclose(distances.T, distances):
#         raise ValueError("Distance matrix must be symmetric")
#     if np.isnan(distances).any():
#         raise ValueError("Distance matrix contains NaN values")
#     if (distances >= MAX_FLOAT).any():
#         raise ValueError("Distance matrix contains infinity")
#     if (distances < 0).any():
#         raise ValueError("Distances must be positive")

#     # Initialize leaves
#     # nodes: list[TreeNode | None] = [TreeNode(index=i) for i in range(n)]
#     nodes: list[Optional[TreeNode]] = [TreeNode(index=i) for i in range(n)]
#     is_clustered = np.full(n, False, dtype=bool)
#     cluster_size = np.ones(n, dtype=np.uint32)
#     node_heights = np.zeros(n, dtype=np.float64)
#     distances_v = distances.astype(np.float64, copy=True)

#     while True:
#         # Find minimum distance between unclustered nodes
#         dist_min = MAX_FLOAT
#         i_min, j_min = -1, -1
#         for i in range(n):
#             if is_clustered[i]:
#                 continue
#             for j in range(i):
#                 if is_clustered[j]:
#                     continue
#                 dist = distances_v[i, j]
#                 if dist < dist_min:
#                     dist_min, i_min, j_min = dist, i, j

#         # Exit if all nodes have been clustered
#         if i_min == -1 or j_min == -1:
#             break

#         # Cluster the two closest nodes
#         height = dist_min / 2.0
#         left, right = nodes[i_min], nodes[j_min]
#         if left is None or right is None:
#             raise RuntimeError("Attempted to cluster missing node")

#         nodes[i_min] = TreeNode(
#             children=(left, right),
#             distances=(height - node_heights[i_min], height - node_heights[j_min]),
#         )
#         node_heights[i_min] = height
#         nodes[j_min] = None
#         is_clustered[j_min] = True

#         # Update distance matrix with new cluster means
#         for k in range(n):
#             if not is_clustered[k] and k != i_min:
#                 mean = (
#                     distances_v[i_min, k] * cluster_size[i_min]
#                     + distances_v[j_min, k] * cluster_size[j_min]
#                 ) / (cluster_size[i_min] + cluster_size[j_min])
#                 distances_v[i_min, k] = mean
#                 distan

# upgma.py

from tree import Tree, TreeNode
import numpy as np

MAX_FLOAT: float = np.finfo(np.float64).max

def upgma(distances: np.ndarray) -> Tree:
    """
    Perform hierarchical clustering using UPGMA.
    """
    n: int = distances.shape[0]
    if distances.shape[0] != distances.shape[1] or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    nodes: list[TreeNode] = [TreeNode(index=i) for i in range(n)]
    is_clustered: list[bool] = [False] * n
    cluster_size: list[int] = [1] * n
    node_heights: list[float] = [0.0] * n

    distances_v: np.ndarray = distances.astype(np.float64, copy=True)

    while True:
        dist_min: float = MAX_FLOAT
        i_min: int = -1
        j_min: int = -1

        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist: float = distances_v[i, j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j

        if i_min == -1 or j_min == -1:
            break

        height: float = dist_min / 2
        nodes[i_min] = TreeNode(
            children=(nodes[i_min], nodes[j_min]),
            distance=(height - node_heights[i_min], height - node_heights[j_min])
        )
        node_heights[i_min] = height
        is_clustered[j_min] = True

        for k in range(n):
            if not is_clustered[k] and k != i_min:
                mean: float = (
                    distances_v[i_min, k] * cluster_size[i_min] +
                    distances_v[j_min, k] * cluster_size[j_min]
                ) / (cluster_size[i_min] + cluster_size[j_min])
                distances_v[i_min, k] = mean
                distances_v[k, i_min] = mean

        cluster_size[i_min] += cluster_size[j_min]

    # Find the last non-clustered node (root)
    for i in range(n):
        if not is_clustered[i]:
            return Tree(root=nodes[i])
