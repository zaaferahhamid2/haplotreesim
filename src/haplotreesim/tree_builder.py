"""
Tree Structure Generation for HaploTreeSim.

This module builds clone tree topologies.
"""

import numpy as np
from typing import List, Dict, Optional


class TreeBuilder:
    """
    Generates clone tree structures.
    
    For Week 6, we implement simple tree topologies.
    Future versions will support more complex structures.
    """
    
    def __init__(self, rng: np.random.Generator, num_clones: int):
        """
        Initialize the tree builder.
        
        Args:
            rng: NumPy random number generator
            num_clones: Total number of clones (including root)
        """
        self.rng = rng
        self.num_clones = num_clones
    
    def build_tree(self, tree_type: str = "star") -> Dict[int, Optional[int]]:
        """
        Build a tree structure.
        
        Args:
            tree_type: Type of tree topology:
                - "star": All clones are children of root (Week 6 default)
                - "linear": Linear chain of clones
                - "balanced": Balanced binary tree
                - "random": Random tree structure
        
        Returns:
            Dictionary mapping clone_id -> parent_id (root has parent=None)
        """
        if tree_type == "star":
            return self._build_star_tree()
        elif tree_type == "linear":
            return self._build_linear_tree()
        elif tree_type == "balanced":
            return self._build_balanced_tree()
        elif tree_type == "random":
            return self._build_random_tree()
        else:
            raise ValueError(f"Unknown tree type: {tree_type}")
    
    def _build_star_tree(self) -> Dict[int, Optional[int]]:
        """
        Build a star tree: all non-root clones are children of root.
        
        Returns:
            Parent mapping
        """
        tree = {0: None}  # Root has no parent
        
        for clone_id in range(1, self.num_clones):
            tree[clone_id] = 0  # All children of root
        
        return tree
    
    def _build_linear_tree(self) -> Dict[int, Optional[int]]:
        """
        Build a linear tree: 0 -> 1 -> 2 -> 3 -> ...
        
        Returns:
            Parent mapping
        """
        tree = {0: None}  # Root has no parent
        
        for clone_id in range(1, self.num_clones):
            tree[clone_id] = clone_id - 1  # Parent is previous clone
        
        return tree
    
    def _build_balanced_tree(self) -> Dict[int, Optional[int]]:
        """
        Build a balanced binary tree.
        
        Returns:
            Parent mapping
        """
        tree = {0: None}  # Root has no parent
        
        for clone_id in range(1, self.num_clones):
            parent_id = (clone_id - 1) // 2  # Binary tree parent formula
            tree[clone_id] = parent_id
        
        return tree
    
    def _build_random_tree(self) -> Dict[int, Optional[int]]:
        """
        Build a random tree by randomly selecting parents.
        
        Ensures tree is valid (no cycles, single root).
        
        Returns:
            Parent mapping
        """
        tree = {0: None}  # Root has no parent
        
        for clone_id in range(1, self.num_clones):
            # Choose parent from existing clones
            possible_parents = list(range(clone_id))
            parent_id = self.rng.choice(possible_parents)
            tree[clone_id] = parent_id
        
        return tree
    
    def get_children(self, tree: Dict[int, Optional[int]]) -> Dict[int, List[int]]:
        """
        Convert parent mapping to children mapping.
        
        Args:
            tree: Parent mapping (clone_id -> parent_id)
        
        Returns:
            Children mapping (parent_id -> [child_ids])
        """
        children = {i: [] for i in range(self.num_clones)}
        
        for clone_id, parent_id in tree.items():
            if parent_id is not None:
                children[parent_id].append(clone_id)
        
        return children
    
    def get_tree_depth(self, tree: Dict[int, Optional[int]]) -> int:
        """
        Calculate the maximum depth of the tree.
        
        Args:
            tree: Parent mapping
        
        Returns:
            Maximum depth (root is depth 0)
        """
        depths = {0: 0}  # Root has depth 0
        
        # Calculate depth for each clone
        for clone_id in range(1, self.num_clones):
            parent_id = tree[clone_id]
            depths[clone_id] = depths[parent_id] + 1
        
        return max(depths.values())
