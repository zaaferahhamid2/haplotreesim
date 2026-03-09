"""
Beta-Splitting Tree Generator for HaploTreeSim.

Implements the Beta-splitting model from Section 3.2 of the updated paper:
- Binary tree with K leaves (extant clones)
- Clone proportions via Beta distribution
- Normalized branch lengths
- Events tied to branch length
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


@dataclass
class TreeNode:
    """
    Node in the Beta-splitting tree.
    
    Attributes:
        node_id: Unique identifier
        parent_id: Parent node ID (-1 for root)
        interval: Tuple (a, b) defining mass interval [a, b)
        perc: Percentage/mass of this node
        edge_length: Normalized branch length (τ_v in paper)
        depth: Depth in tree (root = 1)
        is_leaf: Whether this is a leaf node
    """
    node_id: int
    parent_id: int
    interval: Tuple[float, float]
    perc: float
    edge_length: float
    depth: int
    is_leaf: bool


class BetaSplittingTreeBuilder:
    """
    Generates clone trees using Beta-splitting model.
    
    Implements Section 3.2: Binary tree construction via iterative
    Beta-distributed splits of leaf nodes.
    """
    
    def __init__(
        self,
        rng: np.random.Generator,
        num_clones: int,
        alpha_tree: float = 0.5,
        beta_tree: float = 0.3
    ):
        """
        Initialize the Beta-splitting tree builder.
        
        Args:
            rng: NumPy random number generator
            num_clones: Number of leaves (K in paper)
            alpha_tree: Beta distribution α parameter (α+1 used in Beta)
            beta_tree: Beta distribution β parameter (β+1 used in Beta)
        """
        self.rng = rng
        self.num_clones = num_clones
        self.alpha_tree = alpha_tree
        self.beta_tree = beta_tree
    
    def build_tree(self) -> Tuple[List[TreeNode], Dict[int, Optional[int]]]:
        """
        Build a Beta-splitting tree with K leaves.
        
        Returns:
            Tuple of:
                - List of all tree nodes (internal + leaves)
                - Parent mapping (node_id -> parent_id)
        """
        # Total nodes in binary tree: 2K - 1
        total_nodes = 2 * self.num_clones - 1
        
        # Sample branch lengths (one per edge)
        # Equation (12): τ̃_v ~ Exponential(1)
        raw_branch_lengths = self.rng.exponential(1.0, total_nodes)
        
        # Equation (13): Normalize to sum to 1
        branch_lengths = raw_branch_lengths / raw_branch_lengths.sum()
        
        # Sample Beta splits for K-1 internal nodes
        # Equation (7): ξ_j ~ Beta(α+1, β+1)
        beta_splits = self.rng.beta(
            self.alpha_tree + 1,
            self.beta_tree + 1,
            self.num_clones - 1
        )
        
        # Sample uniform values to choose which leaf to split
        uniform_splits = self.rng.uniform(0.0, 1.0, self.num_clones - 1)
        
        # Initialize tree with root
        nodes = []
        
        # Root node (node 0)
        root = TreeNode(
            node_id=0,
            parent_id=-1,
            interval=(0.0, 1.0),
            perc=1.0,
            edge_length=branch_lengths[0],
            depth=1,
            is_leaf=False
        )
        nodes.append(root)
        
        # Create initial two children of root
        left_child = TreeNode(
            node_id=1,
            parent_id=0,
            interval=(0.0, beta_splits[0]),
            perc=beta_splits[0],
            edge_length=branch_lengths[1],
            depth=2,
            is_leaf=True
        )
        
        right_child = TreeNode(
            node_id=2,
            parent_id=0,
            interval=(beta_splits[0], 1.0),
            perc=1.0 - beta_splits[0],
            edge_length=branch_lengths[2],
            depth=2,
            is_leaf=True
        )
        
        nodes.extend([left_child, right_child])
        
        # Iteratively split leaves until we have K leaves
        num_leaves = 2
        node_counter = 2
        split_idx = 1
        
        while num_leaves < self.num_clones:
            # Find which leaf to split based on uniform draw
            u = uniform_splits[split_idx]
            
            # Find leaf whose interval contains u
            leaf_to_split = None
            leaf_idx = None
            
            for i, node in enumerate(nodes):
                if node.is_leaf:
                    a, b = node.interval
                    if a <= u < b:
                        leaf_to_split = node
                        leaf_idx = i
                        break
            
            if leaf_to_split is None:
                # Shouldn't happen, but fallback to first leaf
                for i, node in enumerate(nodes):
                    if node.is_leaf:
                        leaf_to_split = node
                        leaf_idx = i
                        break
            
            # Convert leaf to internal node
            nodes[leaf_idx].is_leaf = False
            
            # Calculate split point within parent's interval
            # Equation (7): m_ℓ = a_ℓ + ξ_j(b_ℓ - a_ℓ)
            a, b = leaf_to_split.interval
            split_point = a + beta_splits[split_idx] * (b - a)
            
            # Create two new leaf children
            # Equations (8-9)
            node_counter += 1
            left_new = TreeNode(
                node_id=node_counter,
                parent_id=leaf_to_split.node_id,
                interval=(a, split_point),
                perc=leaf_to_split.perc * beta_splits[split_idx],
                edge_length=branch_lengths[node_counter],
                depth=leaf_to_split.depth + 1,
                is_leaf=True
            )
            
            node_counter += 1
            right_new = TreeNode(
                node_id=node_counter,
                parent_id=leaf_to_split.node_id,
                interval=(split_point, b),
                perc=leaf_to_split.perc * (1.0 - beta_splits[split_idx]),
                edge_length=branch_lengths[node_counter],
                depth=leaf_to_split.depth + 1,
                is_leaf=True
            )
            
            nodes.extend([left_new, right_new])
            
            num_leaves += 1
            split_idx += 1
        
        # Build parent mapping
        parent_map = {}
        for node in nodes:
            parent_map[node.node_id] = node.parent_id if node.parent_id != -1 else None
        
        return nodes, parent_map
    
    def get_leaf_nodes(self, nodes: List[TreeNode]) -> List[TreeNode]:
        """Get only the leaf nodes (extant clones)."""
        return [n for n in nodes if n.is_leaf]
    
    def get_clone_proportions(self, nodes: List[TreeNode]) -> np.ndarray:
        """
        Get clone mixing proportions π from leaf nodes.
        
        Equation (11): π_k = perc(k) for leaves k ∈ {1,...,K}
        
        Returns:
            Array of length K with clone proportions (sum to 1)
        """
        leaves = self.get_leaf_nodes(nodes)
        # Sort by node_id to ensure consistent ordering
        leaves_sorted = sorted(leaves, key=lambda n: n.node_id)
        proportions = np.array([leaf.perc for leaf in leaves_sorted])
        
        # Verify they sum to 1 (within numerical precision)
        assert abs(proportions.sum() - 1.0) < 1e-10, \
            f"Clone proportions must sum to 1, got {proportions.sum()}"
        
        return proportions
