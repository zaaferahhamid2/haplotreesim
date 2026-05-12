"""
Tree Metrics (Week 13)
"""

import numpy as np
from typing import List, Dict, Tuple, Set, Optional


class TreeNode:
    def __init__(self, node_id: int, parent_id: Optional[int] = None):
        self.node_id = node_id
        self.parent_id = parent_id
        self.children = []
    
    def add_child(self, child):
        self.children.append(child)


def build_tree_from_edges(edges: List[Tuple[int, int]]) -> Dict[int, TreeNode]:
    nodes = {}
    all_nodes = set()
    for parent, child in edges:
        all_nodes.add(parent)
        all_nodes.add(child)
    
    for node_id in all_nodes:
        nodes[node_id] = TreeNode(node_id)
    
    for parent_id, child_id in edges:
        nodes[child_id].parent_id = parent_id
        nodes[parent_id].add_child(nodes[child_id])
    
    return nodes


def get_leaf_descendants(node: TreeNode, nodes: Dict[int, TreeNode]) -> Set[int]:
    if len(node.children) == 0:
        return {node.node_id}
    descendants = set()
    for child in node.children:
        descendants.update(get_leaf_descendants(child, nodes))
    return descendants


def get_bipartitions(nodes: Dict[int, TreeNode]) -> Set[frozenset]:
    root = None
    for node in nodes.values():
        if node.parent_id is None:
            root = node
            break
    
    if root is None:
        return set()
    
    leaves = {nid for nid, node in nodes.items() if len(node.children) == 0}
    bipartitions = set()
    
    for node in nodes.values():
        if len(node.children) == 0 or node.parent_id is None:
            continue
        descendants = get_leaf_descendants(node, nodes)
        if 0 < len(descendants) < len(leaves):
            bipartitions.add(frozenset(descendants))
    
    return bipartitions


def compute_robinson_foulds_distance(
    true_edges: List[Tuple[int, int]],
    pred_edges: List[Tuple[int, int]],
    normalize: bool = True
) -> float:
    true_nodes = build_tree_from_edges(true_edges)
    pred_nodes = build_tree_from_edges(pred_edges)
    
    true_bip = get_bipartitions(true_nodes)
    pred_bip = get_bipartitions(pred_nodes)
    
    rf_distance = len(true_bip.symmetric_difference(pred_bip))
    
    if normalize:
        true_leaves = {nid for nid, n in true_nodes.items() if len(n.children) == 0}
        n_leaves = len(true_leaves)
        if n_leaves <= 3:
            return 0.0
        max_rf = 2 * (n_leaves - 3)
        return rf_distance / max_rf
    
    return float(rf_distance)


def get_lca(node1_id: int, node2_id: int, nodes: Dict[int, TreeNode]) -> int:
    ancestors1 = set()
    current = node1_id
    while current is not None:
        ancestors1.add(current)
        current = nodes[current].parent_id
    
    current = node2_id
    while current is not None:
        if current in ancestors1:
            return current
        current = nodes[current].parent_id
    return None


def get_node_depth(node_id: int, nodes: Dict[int, TreeNode]) -> int:
    depth = 0
    current = node_id
    while nodes[current].parent_id is not None:
        depth += 1
        current = nodes[current].parent_id
    return depth


def compute_ancestor_descendant_accuracy(
    true_edges: List[Tuple[int, int]],
    pred_edges: List[Tuple[int, int]],
    leaf_mapping: Optional[Dict[int, int]] = None
) -> Dict[str, float]:
    true_nodes = build_tree_from_edges(true_edges)
    pred_nodes = build_tree_from_edges(pred_edges)
    
    true_leaves = [nid for nid, n in true_nodes.items() if len(n.children) == 0]
    pred_leaves = [nid for nid, n in pred_nodes.items() if len(n.children) == 0]
    
    if leaf_mapping:
        pred_leaves = [leaf_mapping.get(p, p) for p in pred_leaves]
    
    common = sorted(set(true_leaves) & set(pred_leaves))
    
    if len(common) < 2:
        return {'accuracy': 1.0, 'correct': 0, 'total': 0}
    
    correct, total = 0, 0
    
    for i in range(len(common)):
        for j in range(i + 1, len(common)):
            leaf1, leaf2 = common[i], common[j]
            
            true_lca = get_lca(leaf1, leaf2, true_nodes)
            pred_lca = get_lca(leaf1, leaf2, pred_nodes)
            
            true_rel = (get_node_depth(leaf1, true_nodes) - get_node_depth(true_lca, true_nodes),
                       get_node_depth(leaf2, true_nodes) - get_node_depth(true_lca, true_nodes))
            pred_rel = (get_node_depth(leaf1, pred_nodes) - get_node_depth(pred_lca, pred_nodes),
                       get_node_depth(leaf2, pred_nodes) - get_node_depth(pred_lca, pred_nodes))
            
            if true_rel == pred_rel:
                correct += 1
            total += 1
    
    return {'accuracy': float(correct / total if total > 0 else 1.0), 
            'correct': int(correct), 'total': int(total)}


def compute_event_matching(
    true_events: List[Dict],
    pred_events: List[Dict],
    tolerance: int = 1
) -> Dict[str, float]:
    if not true_events and not pred_events:
        return {'precision': 1.0, 'recall': 1.0, 'f1': 1.0, 'tp': 0, 'fp': 0, 'fn': 0}
    
    matched_true, matched_pred = set(), set()
    
    for i, pred in enumerate(pred_events):
        for j, true in enumerate(true_events):
            if j in matched_true:
                continue
            
            overlaps = (pred['start'] - tolerance <= true['end'] and 
                       true['start'] - tolerance <= pred['end'])
            if not overlaps:
                continue
            
            pred_hap = pred.get('haplotype', 'A')
            true_hap = true.get('haplotype', 'A')
            hap_match = pred_hap == true_hap or (pred_hap in ['A','B'] and true_hap in ['A','B'])
            
            if hap_match:
                matched_true.add(j)
                matched_pred.add(i)
                break
    
    tp = len(matched_pred)
    fp = len(pred_events) - tp
    fn = len(true_events) - len(matched_true)
    
    prec = tp / (tp + fp) if tp + fp > 0 else 0.0
    rec = tp / (tp + fn) if tp + fn > 0 else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0.0
    
    return {'precision': float(prec), 'recall': float(rec), 'f1': float(f1),
            'tp': int(tp), 'fp': int(fp), 'fn': int(fn)}


def compute_all_tree_metrics(
    true_edges: List[Tuple[int, int]],
    pred_edges: List[Tuple[int, int]],
    true_events: Optional[List[Dict]] = None,
    pred_events: Optional[List[Dict]] = None,
    leaf_mapping: Optional[Dict[int, int]] = None
) -> Dict[str, float]:
    results = {}
    results['rf_distance'] = compute_robinson_foulds_distance(true_edges, pred_edges)
    
    ad = compute_ancestor_descendant_accuracy(true_edges, pred_edges, leaf_mapping)
    results['ancestor_descendant_accuracy'] = ad['accuracy']
    
    if true_events and pred_events:
        em = compute_event_matching(true_events, pred_events)
        results.update({f'event_{k}': v for k, v in em.items() if k in ['precision','recall','f1']})
    
    return results
