# upgma.py
from p_tree import Tree, TreeNode
import numpy as np

MAX_FLOAT = np.finfo(np.float64).max

def upgma(distances: np.ndarray) -> Tree:
    """
    Perform hierarchical clustering using UPGMA (Unweighted Pair Group Method
    with Arithmetic mean).
    """
    n = distances.shape[0]

    # Input checks
    if distances.shape[0] != distances.shape[1] or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    # Keep track of nodes
    # nodes: list[TreeNode | None] = [TreeNode(index=i) for i in range(n)]
    nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n)]
    is_clustered: np.ndarray = np.zeros(n, dtype=bool)
    cluster_size: np.ndarray = np.ones(n, dtype=int)
    node_heights: np.ndarray = np.zeros(n, dtype=np.float64)

    # Working copy of distances
    distances_v: np.ndarray = distances.astype(np.float64, copy=True)

    while True:
        # Find minimum distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist = distances_v[i, j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j

        if i_min == -1 or j_min == -1:
            break  # all nodes clustered

        # Merge the two closest nodes
        height = dist_min / 2
        nodes[i_min] = TreeNode(
            children=(nodes[i_min], nodes[j_min]),
            distances=(height - node_heights[i_min], height - node_heights[j_min])
        )
        node_heights[i_min] = height
        # nodes[j_min] = None
        is_clustered[j_min] = True

        # Update distances
        for k in range(n):
            if not is_clustered[k] and k != i_min:
                mean = (
                    distances_v[i_min, k] * cluster_size[i_min]
                    + distances_v[j_min, k] * cluster_size[j_min]
                ) / (cluster_size[i_min] + cluster_size[j_min])
                distances_v[i_min, k] = mean
                distances_v[k, i_min] = mean

        cluster_size[i_min] += cluster_size[j_min]

    root_idx = -1
    for idx in range(n - 1, -1, -1):
        if not is_clustered[idx]:
            root_idx = idx
            break

    if root_idx == -1:
        raise RuntimeError("UPGMA failed to produce a root node")

    # nodes[root_idx] is guaranteed to be a TreeNode
    return Tree(nodes[root_idx])

    # Return the root node
    # return Tree(nodes[-1])
    # root = nodes[-1]
    # assert root is not None
    # return Tree(root)

    # root_opt = nodes[-1]
    # if root_opt is None:
    #     raise ValueError("UPGMA failed to produce a root node")
    # root: TreeNode = root_opt
    # return Tree(root)

