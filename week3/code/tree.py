from typing import Optional, Tuple, List
import numpy as np

# class TreeError():
#     """An exception that occurs in context of tree topology."""
#     def __init__(self, message="Tree error occurred"):
#         print(message)

@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747)  * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h

class TreeNode:
    _index: int
    _distance: float
    _is_root: bool
    _parent: Optional[TreeNode]
    _children: List[TreeNode]

    def __init__(self, children=None, distances=None, index=None):
        self._is_root = False
        self._distance = 0
        # self._parent = None
        self._parent: Optional[TreeNode] = None
        self._children = []
        # child: TreeNode
        # distance: float
        if index is None:
            # Node is intermediate -> has children
            if children is None or distances is None:
                raise TypeError(
                    "Either reference index (for terminal node) or "
                    "child nodes including the distance "
                    "(for intermediate node) must be set"
                )
            for item in children:
                if not isinstance(item, TreeNode):
                    raise TypeError(
                        f"Expected 'TreeNode', but got '{type(item).__name__}'"
                    )
            for item in distances:
                if not isinstance(item, float) and not isinstance(item, int):
                    raise TypeError(
                        f"Expected 'float' or 'int', "
                        f"but got '{type(item).__name__}'"
                    )
            # if len(children) == 0:
            #     raise TreeError(
            #         "Intermediate nodes must at least contain one child node"
            #     )
            if len(children) != len(distances):
                raise ValueError(
                    "The number of children must equal the number of distances"
                )
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        # raise TreeError(
                        #     "Two child nodes cannot be the same object"
                        # )
                        print("Two child nodes can't be the same object")
            self._index = -1
            # self._children = tuple(children)
            self._children = [i for i in children]
            for child, distance in zip(children, distances):
                child._set_parent(self, distance)
        elif index < 0:
            raise ValueError("Index cannot be negative")
        else:
            # Node is terminal -> has no children
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive"
                )
            self._index = index
            self._children = []

    def _set_parent(self, parent: Optional[TreeNode], distance: float):
        # if self._parent is not None or self._is_root:
        #     raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = distance

    def copy(self):
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            distances = [child.distance for child in self._children]
            children_clones = [child.copy() for child in self._children]
            return TreeNode(children_clones, distances)

    @property
    def index(self):
        return None if self._index == -1 else self._index
    
    @property
    def children(self):
        return self._children
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def distance(self):
        return None if self._parent is None else self._distance


    def is_leaf(self):
        return False if self._index == -1 else True

    def is_root(self):
        return bool(self._is_root)

    def as_root(self):
        # if self._parent is not None:
        #     raise TreeError("Node has parent, cannot be a root node")
        self._is_root = True

    def distance_to(self, node: Optional[TreeNode], topological: bool = False):
        # Sum distances until LCA has been reached
        distance: float = 0.0
        current_node: Optional[TreeNode] = None
        lca: TreeNode = self.lowest_common_ancestor(node)
        # if lca is None:
        #     raise TreeError("The nodes do not have a common ancestor")
        current_node = self
        while current_node is not lca:
            if topological:
                distance += 1
            else:
                distance += current_node._distance
            current_node = current_node._parent
        current_node = node
        while current_node is not lca:
            if topological:
                distance += 1
            else:
                distance += current_node._distance
            current_node = current_node._parent
        return distance
    
    def lowest_common_ancestor(self, node: Optional[TreeNode]):
        lca: Optional[TreeNode] = None
        # Create two paths from the leaves to root
        self_path: List = _create_path_to_root(self)
        other_path: List = _create_path_to_root(node)
        # Reverse Iteration through path (beginning from root)
        # until the paths diverge
        for i in range(-1, -min(len(self_path), len(other_path))-1, -1):
            if self_path[i] is other_path[i]:
                # Same node -> common ancestor
                lca = self_path[i]
            else:
                # Different node -> Not common ancestor
                # -> return last common ancewstor found
                break
        return lca

    def get_leaves(self):
        leaf_list: List = []
        # delegate to 'cdef' method
        # to reduce overhead of recursive function calling
        _get_leaves(self, leaf_list)
        return leaf_list

    def to_newick(self, labels=None, include_distance: bool =True, 
                  round_distance=None):
        if self.is_leaf():
            if labels is not None:
                for label in labels:
                    label = labels[self._index]
                    # Characters that are part of the Newick syntax
                    # are illegal
                    illegal_chars = [",",":",";","(",")"]
                    for char in illegal_chars:
                        if char in label:
                            raise ValueError(
                                f"Label '{label}' contains "
                                f"illegal character '{char}'"
                            )
            else:
                label = str(self._index)
            if include_distance:
                if round_distance is None:
                    return f"{label}:{self._distance}"
                else:
                    return f"{label}:{self._distance:.{round_distance}f}"
            else:
                return f"{label}"
        else:
            # Build string in a recursive way
            child_strings = [child.to_newick(
                labels, include_distance, round_distance
            ) for child in self._children]
            if include_distance:
                if round_distance is None:
                    return f"({','.join(child_strings)}):{self._distance}"
                else:
                    return (
                        f"({','.join(child_strings)}):"
                        f"{self._distance:.{round_distance}f}"
                    )
            else:
                return f"({','.join(child_strings)})"

    @staticmethod
    def from_newick(newick: str, labels: List =None):
        subnewick_start_i: int = -1
        subnewick_stop_i: int  = -1
        level: int = 0
        comma_pos: List
        children: List
        distances: List
        pos: int
        next_pos: int
        
        # Ignore any whitespace
        newick = "".join(newick.split())

        # Find brackets belonging to sub-newick
        # e.g. (A:0.1,B:0.2):0.5
        #      ^           ^
        for i in range(len(newick)):
            char = newick[i]
            if char == "(":
                subnewick_start_i = i
                break
            if char == ")":
                raise InvalidFileError("Bracket closed before it was opened")
        for i in reversed(range(len(newick))):
            char = newick[i]
            if char == ")":
                subnewick_stop_i = i+1
                break
            if char == "(":
                raise InvalidFileError("Bracket was opened but not closed")
        
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            # No brackets -> no sub-newwick -> Leaf node
            label_and_distance = newick
            try:
                label, distance = label_and_distance.split(":")
                distance = float(distance)
            except ValueError:
                # No colon -> No distance is provided
                distance = 0
                label = label_and_distance
            index = int(label) if labels is None else labels.index(label)
            return TreeNode(index=index), distance
        
        else:
            # Intermediate node
            if subnewick_stop_i == len(newick):
                # Node with neither distance nor label
                label = None
                distance = 0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                try:
                    label, distance = label_and_distance.split(":")
                    distance = float(distance)
                except ValueError:
                    # No colon -> No distance is provided
                    distance = 0
                    label = label_and_distance
                # Label of intermediate nodes is discarded 
                distance = float(distance)
            
            subnewick = newick[subnewick_start_i+1 : subnewick_stop_i-1]
            if len(subnewick) == 0:
                raise InvalidFileError(
                    "Intermediate node must at least have one child"
                )
            # Parse childs
            # Split subnewick at ',' if ',' is at current level
            # (not in a subsubnewick)
            comma_pos = []
            for i, char in enumerate(subnewick):
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == ",":
                    if level == 0:
                        comma_pos.append(i)
                if level < 0:
                    raise InvalidFileError(
                        "Bracket closed before it was opened"
                    )
        
            children = []
            distances = []
            # Recursive tree construction
            for i, pos in enumerate(comma_pos):
                if i == 0:
                    # (A,B),(C,D),(E,F)
                    # -----
                    child, dist = TreeNode.from_newick(
                        subnewick[:pos], labels=labels
                    )
                else:
                    # (A,B),(C,D),(E,F)
                    #       -----
                    prev_pos = comma_pos[i-1]
                    child, dist = TreeNode.from_newick(
                        subnewick[prev_pos+1 : pos], labels=labels
                    )
                children.append(child)
                distances.append(dist)
            # Node after last comma
            # (A,B),(C,D),(E,F)
            #             -----
            if len(comma_pos) != 0:
                child, dist = TreeNode.from_newick(
                    subnewick[comma_pos[-1]+1:], labels=labels
                )
            else:
                # Single child node:
                child, dist = TreeNode.from_newick(
                    subnewick, labels=labels
                )
            children.append(child)
            distances.append(dist)
            return TreeNode(children, distances), distance

    def __str__(self):
        return self.to_newick()
    
    def __eq__(self, item):
        if not isinstance(item, TreeNode):
            return False
        node: TreeNode = item
        if self._distance != node._distance:
            return False
        if self._index !=-1:
            if self._index != node._index:
                return False
        else:
            if set(self._children) != set(node._children):
                return False
        return True

    def __hash__(self):
        # Order of children is not important -> set
        children_set = set(self._children) \
                       if self._children is not None else None
        return hash((self._index, children_set, self._distance))


class Tree():
    
    def __init__(self, root: TreeNode):
        root.as_root()
        self._root = root
        
        leaves_unsorted: List = self._root.get_leaves()
        leaf_count: int = len(leaves_unsorted)
        indices: np.ndarray = np.array(
            [leaf.index for leaf in leaves_unsorted]
        )
        self._leaves = [None] * leaf_count
        for i in range(len(indices)):
            index = indices[i]
            if index >= leaf_count or index < 0:
                # raise TreeError("The tree's indices are out of range")
                print("The tree's indices are out of range")
            self._leaves[index] = leaves_unsorted[i]
    
    def __copy_create__(self):
        return Tree(self._root.copy())
    
    @property
    def root(self):
        return self._root
    
    @property
    def leaves(self):
        # return copy.copy(self._leaves)
        return self._leaves[:]  # shallow copy

    def get_distance(self, index1, index2, topological: bool=False):
        return self._leaves[index1].distance_to(
            self._leaves[index2], topological
        )
    
    def to_newick(self, labels=None, include_distance: bool=True, 
                  round_distance=None):
        return self._root.to_newick(
            labels, include_distance, round_distance
        ) + ";"
    
    @staticmethod
    def from_newick(newick: str, labels: List=None):
        newick = newick.strip()
        if len(newick) == 0:
            raise InvalidFileError("Newick string is empty")
        # Remove terminal colon as required by 'TreeNode.from_newick()'
        if newick[-1] == ";":
            newick = newick[:-1]
        root, distance = TreeNode.from_newick(newick, labels)
        return Tree(root)

    def __str__(self):
        return self.to_newick()
    
    def __len__(self):
        return len(self._leaves)
    
    def __eq__(self, item):
        if not isinstance(item, Tree):
            return False
        return self._root == item._root
    
    def __hash__(self):
        return hash(self._root)


# cdef functions
def _get_leaves(node: TreeNode, leaf_list: List[TreeNode]):
    child: TreeNode
    if node._index == -1:
        # Intermediate node -> Recursive calls
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        # Node itself is leaf node -> add node -> terminate
        leaf_list.append(node)

def _get_leaf_count(node: TreeNode) -> int:
    child: TreeNode
    count: int = 0
    if node._index == -1:
        # Intermediate node -> Recursive calls
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        # Leaf node -> return count of itself = 1
        return 1


def _create_path_to_root(node: Optional[TreeNode]) -> List:
    """
    Create a list of nodes representing the path from this node to the
    specified node
    """
    path: List = []
    current_node: Optional[TreeNode] = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path

def as_binary(tree_or_node):
    """
    as_binary(tree_or_node)

    Convert a tree into a binary tree.

    In general a :class:`TreeNode` can have more or less than two
    children.
    However guide trees usually expect each intermediate node to have
    exactly two child nodes.
    This function creates a binary :class:`Tree` (or :class:`TreeNode`)
    for the given :class:`Tree` (or :class:`TreeNode`):
    Intermediate nodes that have only a single child are deleted and its
    parent node is directly connected to its child node.
    Intermediate nodes that have more than two childs are divided into
    multiple nodes (distances are preserved).
    
    Parameters
    ----------
    tree_or_node : Tree or TreeNode
        The tree or node to be converted into a binary tree or node.
    
    Returns
    -------
    binary_tree_or_node : Tree or TreeNode
        The converted tree or node.
    """
    if isinstance(tree_or_node, Tree):
        node, _ = _as_binary(tree_or_node.root)
        return Tree(node)
    elif isinstance(tree_or_node, TreeNode):
        node, _ = _as_binary(tree_or_node)
        return _as_binary(node)
    else:
        raise TypeError(
            f"Expected 'Tree' or 'TreeNode', not {type(tree_or_node).__name__}"
        )

def _as_binary(node: TreeNode):
    child: TreeNode
    current_div_node: TreeNode
    children: tuple
    rem_children: list
    distances: list
    distance: float

    children = node.children
    if children is None:
        # Leaf node
        return TreeNode(index=node.index), node.distance
    elif len(children) == 1:
        # Intermediate node with one child
        # -> Omit node and directly connect its child to its parent
        # The distances are added
        #
        #      |--            |--   
        #      |              |   
        # --|--|--   ->   ----|--  
        #      |              |   
        #      |--            |-- 
        #
        child, distance = _as_binary(node.children[0])
        if node.is_root():
            # Child is new root -> No distance to parent
            return child, None
        else:
            return child, node.distance + distance
    elif len(children) > 2:
        # Intermediate node with more than two childs
        # -> Create a new node having two childs:
        #    - One of the childs of the original node
        #    - The original node with one child less (distance = 0)
        # Repeat until all children are put into binary nodes
        #
        #   |--          |--
        #   |          --|  |--
        # --|--   ->     |--|
        #   |               |--
        #   |--
        #
        # The remaining children
        rem_children, distances = [list(tup) for tup in zip(
            *[_as_binary(child) for child in children]
        )]
        current_div_node = None
        while len(rem_children) > 0:
            if current_div_node is None:
                # The bottom-most node is created
                #-> Gets two of the remaining childs
                current_div_node = TreeNode(
                    rem_children[:2],
                    distances[:2]
                )
                # Pop the two utilized remaining childs from the list
                rem_children.pop(0)
                rem_children.pop(0)
                distances.pop(0)
                distances.pop(0)
            else:
                # A node is created that gets one remaining child
                # and the intermediate node from the last step
                current_div_node = TreeNode(
                    (current_div_node, rem_children[0]),
                    (0, distances[0]) 
                )
                # Pop the utilized remaining child from the list
                rem_children.pop(0)
                distances.pop(0)
        return current_div_node, node.distance
    else:
        # Intermediate node with exactly two childs
        # -> Keep node unchanged
        binary_children, distances = [list(tup) for tup in zip(
            *[_as_binary(child) for child in children]
        )]
        return TreeNode(binary_children, distances), node.distance










# class TreeNode:
#     def __init__(
#         self,
#         children: List[TreeNode],
#         distances: Optional[Tuple] = None,
#         index: Optional[int] = None
#     ):
#         self._is_root: bool = False
#         self._distance: float = 0.0
#         self._parent: Optional[TreeNode] = None

#         if index is None:
#             # Intermediate node
#             if children is None or distances is None:
#                 raise TypeError(
#                     "Either reference index (for terminal node) or "
#                     "child nodes including the distance (for intermediate node) must be set"
#                 )
#             for item in children:
#                 if not isinstance(item, TreeNode):
#                     raise TypeError(f"Expected 'TreeNode', but got '{type(item).__name__}'")
#             for item in distances:
#                 if not isinstance(item, (float, int)):
#                     raise TypeError(f"Expected 'float' or 'int', but got '{type(item).__name__}'")
#             # if len(children) == 0:
#             #     raise TreeError("Intermediate nodes must at least contain one child node")
#             if len(children) != len(distances):
#                 raise ValueError("The number of children must equal the number of distances")
#             # for i in range(len(children)):
#             #     for j in range(len(children)):
#             #         if i != j and children[i] is children[j]:
#             #             raise TreeError("Two child nodes cannot be the same object")
#             self._index = -1
#             self._children = tuple(children)
#             for child, distance in zip(children, distances):
#                 child._set_parent(self, distance)
#         else:
#             # Terminal node
#             if children is not None or distances is not None:
#                 raise TypeError("Reference index and child nodes are mutually exclusive")
#             if index < 0:
#                 raise ValueError("Index cannot be negative")
#             self._index = index
#             self._children = None

#     def _set_parent(self, parent: "TreeNode", distance: float):
#         self._parent = parent
#         self._distance = distance

#     @property
#     def index(self) -> Optional[int]:
#         return None if self._index == -1 else self._index

#     @property
#     def children(self) -> Optional[Tuple]:
#         return self._children

#     @property
#     def parent(self) -> Optional[TreeNode]:
#         return self._parent

#     @property
#     def distance(self) -> Optional[float]:
#         return None if self._parent is None else self._distance

#     def __str__(self) -> str:
#         return self.to_newick()

#     def __eq__(self, other) -> bool:
#         if not isinstance(other, TreeNode):
#             return False
#         node = other
#         if self._distance != node._distance:
#             return False
#         if self._index != -1:
#             return self._index == node._index
#         else:
#             return set(self._children) == set(node._children)

#     def __hash__(self) -> int:
#         children_set = set(self._children) if self._children else None
#         return hash((self._index, children_set, self._distance))


# class Tree:
#     def __init__(self, root: TreeNode):
#         root._is_root = True
#         self._root: TreeNode = root

#         leaves_unsorted: List[TreeNode] = self._root.get_leaves()
#         leaf_count = len(leaves_unsorted)
#         indices: np.ndarray = np.array([leaf.index for leaf in leaves_unsorted])

#         self._leaves: List[Optional[TreeNode]] = [None] * leaf_count
#         for i in range(len(indices)):
#             index = indices[i]
#             # if index >= leaf_count or index < 0:
#             #     raise TreeError("The tree's indices are out of range")
#             self._leaves[index] = leaves_unsorted[i]

#     def get_distance(self, index1: int, index2: int, topological: bool = False) -> float:
#         return self._leaves[index1].distance_to(self._leaves[index2], topological)

#     def __str__(self) -> str:
#         return self.to_newick()

#     def __len__(self) -> int:
#         return len(self._leaves)

#     def __eq__(self, other) -> bool:
#         if not isinstance(other, Tree):
#             return False
#         return self._root == other._root

#     def __hash__(self) -> int:
#         return hash(self._root)
