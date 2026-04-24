"""
HaploTreeSim simulator core.

This module implements the main simulator class that generates synthetic
scDNA-seq data with haplotype-specific CNAs and clone tree ground truth.
"""

import numpy as np
import json
import pickle
from typing import List, Tuple, Dict, Optional
from .data_models import (
    Bin, Segment, HaplotypeBlock, CNAEvent, Clone, Cell,
    SimulationConfig, Haplotype
)


class HaploTreeSimulator:
    """
    Main simulator class for HaploTreeSim.
    
    This class orchestrates the simulation process:
    1. Initialize genome (bins, segments, haplotype blocks)
    2. Generate clone tree and haplotype-specific CN profiles
    3. Sample cells from clones
    4. Generate observations (read counts and allele counts)
    
    Attributes:
        config: Simulation configuration parameters
        rng: NumPy random number generator
        bins: List of genomic bins
        segments: List of segments
        haplotype_blocks: List of haplotype blocks
        clones: List of clones (tree nodes)
        cells: List of simulated cells
    """
    
    def __init__(self, config: SimulationConfig):
        """
        Initialize the simulator.
        
        Args:
            config: Simulation configuration
        """
        self.config = config
        self.rng = np.random.default_rng(config.random_seed)
        
        # Data structures (to be populated)
        self.bins: List[Bin] = []
        self.segments: List[Segment] = []
        
        # WGD tracking
        self.wgd_occurred: bool = False
        self.wgd_node: Optional[int] = None  # Which node has WGD
        self.haplotype_blocks: List[HaplotypeBlock] = []
        self.clones: List[Clone] = []
        self.cells: List[Cell] = []
        
        # Mappings
        self.bin_to_segment: Dict[int, int] = {}  # bin index -> segment index
        self.segment_to_block: Dict[int, int] = {}  # segment index -> block index
    

    def export_tree_structure(self, filepath: str):
        """
        Export tree structure for reuse across simulations.
        
        Saves:
        - Node IDs, parent relationships, branch lengths
        - Clone proportions
        - Tree topology
        
        Args:
            filepath: Path to save tree structure (.pkl or .json)
        """
        tree_data = {
            'num_clones': self.config.num_clones,
            'alpha_tree': self.config.alpha_tree,
            'beta_tree': self.config.beta_tree,
            'nodes': [],
            'clone_proportions': self._clone_proportions.tolist() if hasattr(self, '_clone_proportions') else None,
            'leaf_clone_indices': self._leaf_clone_indices if hasattr(self, '_leaf_clone_indices') else None
        }
        
        # Export tree nodes
        for node in self._tree_nodes:
            tree_data['nodes'].append({
                'node_id': int(node.node_id),
                'parent_id': int(node.parent_id) if node.parent_id is not None else None,
                'edge_length': float(node.edge_length),
                'is_leaf': bool(node.is_leaf)
            })
        
        if filepath.endswith('.json'):
            with open(filepath, 'w') as f:
                json.dump(tree_data, f, indent=2)
        else:
            with open(filepath, 'wb') as f:
                pickle.dump(tree_data, f)
        
        print(f"✓ Tree structure exported to: {filepath}")
    
    def import_tree_structure(self, filepath: str):
        """
        Import previously saved tree structure.
        
        Args:
            filepath: Path to tree structure file (.pkl or .json)
        """
        if filepath.endswith('.json'):
            with open(filepath, 'r') as f:
                tree_data = json.load(f)
        else:
            with open(filepath, 'rb') as f:
                tree_data = pickle.load(f)
        
        # Recreate tree nodes
        from .beta_tree_builder import TreeNode
        
        self._tree_nodes = []
        for node_data in tree_data['nodes']:
            node = TreeNode(
                node_id=node_data['node_id'],
                parent_id=node_data['parent_id'],
                edge_length=node_data['edge_length'],
                is_leaf=node_data['is_leaf']
            )
            self._tree_nodes.append(node)
        
        self._clone_proportions = np.array(tree_data['clone_proportions'])
        self._leaf_clone_indices = tree_data['leaf_clone_indices']
        
        print(f"✓ Tree structure imported from: {filepath}")
        print(f"  Nodes: {len(self._tree_nodes)}, Leaves: {len([n for n in self._tree_nodes if n.is_leaf])}")


    def _initialize_clones_from_imported_tree(self):
        """
        Initialize clones from imported tree structure.
        Called after import_tree_structure() to continue simulation.
        """
        from .event_generator import EventGenerator
        from .event_applier import EventApplier
        
        # Initialize event infrastructure
        event_generator = EventGenerator(
            rng=self.rng,
            num_bins=self.config.num_bins,
            bin_length=self.config.bin_width,
            chromosome=self.config.chromosome,
            lambda_events=self.config.lambda_events,
            lambda_amplitude=self.config.lambda_amplitude,
            prob_wgd=self.config.prob_wgd,
            gain_prob=self.config.gain_prob,
            prob_focal=self.config.prob_focal,
            prob_arm_given_broad=self.config.prob_arm_given_broad,
            focal_length_min=self.config.focal_length_min,
            focal_length_max_fraction=self.config.focal_length_max_fraction,
            prob_haplotype_A=self.config.prob_haplotype_A,
        )
        
        event_applier = EventApplier(max_copy_number=self.config.max_copy_number)
        
        # Sample WGD placement
        wgd_node_id = self._sample_wgd_placement()
        if wgd_node_id is not None:
            self.wgd_occurred = True
            self.wgd_node = wgd_node_id
            print(f"  WGD will be placed on edge to node {wgd_node_id}")
        
        # Create root clone
        cn_A_init, cn_B_init = self.config.get_root_cn()
        root_cn_A = np.full(self.config.num_bins, cn_A_init, dtype=int)
        root_cn_B = np.full(self.config.num_bins, cn_B_init, dtype=int)
        
        root_clone = Clone(
            index=0,
            parent_index=None,
            cn_profile_A=root_cn_A,
            cn_profile_B=root_cn_B,
            events=[],
            is_root=True
        )
        
        node_id_to_clone_idx = {}
        self.clones = [root_clone]
        node_id_to_clone_idx[0] = 0
        
        # Process nodes in order
        nodes_by_id = {node.node_id: node for node in self._tree_nodes}
        
        for node in sorted(self._tree_nodes, key=lambda n: n.node_id):
            if node.node_id == 0:
                continue
            
            parent_node = nodes_by_id[node.parent_id]
            parent_clone = self.clones[node_id_to_clone_idx[parent_node.node_id]]
            
            # Generate events
            num_events = self.rng.poisson(self.config.lambda_events * node.edge_length)
            
            events = []
            for _ in range(num_events):
                event = event_generator._generate_single_event()
                events.append(event)
            
            # Check for WGD
            if self.wgd_occurred and node.node_id == self.wgd_node:
                from .data_models import CNAEvent, Haplotype
                wgd = CNAEvent(
                    start_bin=0,
                    end_bin=self.config.num_bins - 1,
                    haplotype=Haplotype.WGD,
                    amplitude=0
                )
                events = [wgd] + events
            
            # Apply events
            cn_A, cn_B = event_applier.apply_events(parent_clone, events)
            
            # Create clone
            child_clone = Clone(
                index=len(self.clones),
                parent_index=node_id_to_clone_idx[parent_node.node_id],
                cn_profile_A=cn_A,
                cn_profile_B=cn_B,
                events=events,
                is_root=False
            )
            
            self.clones.append(child_clone)
            node_id_to_clone_idx[node.node_id] = child_clone.index
        
        # Compute ploidy
        for clone in self.clones:
            clone.ploidy = float(np.mean(clone.total_cn()))
        
        # Detect segments
        self._detect_segments_from_clones()
        
        print(f"  Created {len(self.clones)} clones from imported tree")
        print(f"  Total CNA events: {sum(len(c.events) for c in self.clones)}")

    def run(self) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Run the full simulation pipeline.
        
        Returns:
            Tuple of:
                - read_counts: Array of shape (N, B) with read counts
                - allele_counts: Tuple of (alternate, reference, total) arrays,
                                each of shape (N, S)
        """
        print("Initializing genome...")
        self._initialize_genome()
        
        print(f"Generating clone tree with {self.config.num_clones} clones...")
        self._generate_clone_tree()

        # Compute ploidy for each clone (Section 2: auxiliary truth)
        for clone in self.clones:
            clone.ploidy = float(np.mean(clone.total_cn()))
        
        print("Detecting segments from CNA breakpoints...")
        self._detect_segments_from_clones()
        
        print(f"Sampling {self.config.num_cells} cells...")
        self._sample_cells()
        
        print("Generating read-depth observations...")
        read_counts = self._generate_read_counts()
        
        print("Generating allelic observations...")
        allele_counts = self._generate_allele_counts()
        
        print("Simulation complete!")
        return read_counts, allele_counts
    
    def _initialize_genome(self):
        """
        Initialize genome representation (bins, segments, haplotype blocks).
        
        For Week 5 deliverable: creates uniform bins with a single segment
        covering the entire genome.
        """
        # Create bins
        self.bins = []
        for i in range(self.config.num_bins):
            start = i * self.config.bin_length
            end = start + self.config.bin_length
            bin_obj = Bin(
                index=i,
                chromosome="chr1",  # Simplified: single chromosome
                start=start,
                end=end,
                length=self.config.bin_length
            )
            self.bins.append(bin_obj)
        
        # Create a single segment covering all bins (for now)
        # In future weeks, segments will be defined by CNA breakpoints
        segment = Segment(
            index=0,
            bin_indices=set(range(self.config.num_bins)),
            start_bin=0,
            end_bin=self.config.num_bins - 1,
            length=self.config.num_bins * self.config.bin_length,
            haplotype_block=0
        )
        self.segments = [segment]
        
        # Create a single haplotype block
        hap_block = HaplotypeBlock(
            index=0,
            segment_indices=[0],
            orientation=1,
            alternate_haplotype=Haplotype.A
        )
        self.haplotype_blocks = [hap_block]
        
        # Build mappings
        for bin_idx in range(self.config.num_bins):
            self.bin_to_segment[bin_idx] = 0
        self.segment_to_block[0] = 0
        
        print(f"  Created {len(self.bins)} bins, {len(self.segments)} segments, "
              f"{len(self.haplotype_blocks)} haplotype blocks")
    
    def _generate_clone_tree(self):
        """
        Generate clone tree with haplotype-specific copy number profiles.
        
        Updated: Now uses Beta-splitting tree model (Section 3.2).
        """
        from .event_generator import EventGenerator
        from .event_applier import EventApplier
        from .beta_tree_builder import BetaSplittingTreeBuilder
        
        print(f"  Building Beta-splitting tree (α={self.config.alpha_tree}, β={self.config.beta_tree})...")
        
        # Build Beta-splitting tree
        tree_builder = BetaSplittingTreeBuilder(
            rng=self.rng,
            num_clones=self.config.num_clones,
            alpha_tree=self.config.alpha_tree,
            beta_tree=self.config.beta_tree
        )
        
        tree_nodes, parent_map = tree_builder.build_tree()
        clone_proportions = tree_builder.get_clone_proportions(tree_nodes)
        
        # Store for later use in cell sampling
        self._clone_proportions = clone_proportions
        self._tree_nodes = tree_nodes
        
        # Get root copy number
        cn_A_init, cn_B_init = self.config.get_root_cn()
        
        # Create root clone
        root_cn_A = np.full(self.config.num_bins, cn_A_init, dtype=int)
        root_cn_B = np.full(self.config.num_bins, cn_B_init, dtype=int)
        
        root_clone = Clone(
            index=0,
            parent_index=None,
            cn_profile_A=root_cn_A,
            cn_profile_B=root_cn_B,
            events=[],
            is_root=True
        )
        
        # Map tree nodes to clones (only leaves become clones)
        node_id_to_clone_idx = {}
        self.clones = [root_clone]
        node_id_to_clone_idx[0] = 0
        
        # Initialize event generator and applier
        event_generator = EventGenerator(
            rng=self.rng,
            num_bins=self.config.num_bins,
            bin_length=self.config.bin_width,
            chromosome=self.config.chromosome,
            lambda_events=self.config.lambda_events,
            lambda_amplitude=self.config.lambda_amplitude,
            prob_wgd=self.config.prob_wgd,
            gain_prob=self.config.gain_prob,
            prob_focal=self.config.prob_focal,
            prob_arm_given_broad=self.config.prob_arm_given_broad,
            focal_length_min=self.config.focal_length_min,
            focal_length_max_fraction=self.config.focal_length_max_fraction,
            prob_haplotype_A=self.config.prob_haplotype_A,
        )
        
        event_applier = EventApplier(max_copy_number=self.config.max_copy_number)
        
        # Process tree nodes in breadth-first order
        # Build clones for all nodes (not just leaves)
        nodes_by_id = {node.node_id: node for node in tree_nodes}
        
        # Sample WGD placement BEFORE processing nodes (Equations 18-19)
        wgd_node_id = self._sample_wgd_placement()
        if wgd_node_id is not None:
            self.wgd_occurred = True
            self.wgd_node = wgd_node_id  # Store node_id
            print(f"  WGD will be placed on edge to node {wgd_node_id}")
        
        # Process nodes in order of node_id (ensures parents before children)
        for node in sorted(tree_nodes, key=lambda n: n.node_id):
            if node.node_id == 0:
                continue  # Root already created                print(f"    node.node_id = {node.node_id}")
            
            parent_node = nodes_by_id[node.parent_id]
            parent_clone = self.clones[node_id_to_clone_idx[parent_node.node_id]]
            
            # Generate events for this edge
            # Equation (21): M_v ~ Poisson(λ_E * τ_v) 
            # Events tied to branch length (no +1 offset)
            num_events = self.rng.poisson(self.config.lambda_events * node.edge_length)
            
            # Generate focal/arm/chr events
            events = []
            for _ in range(num_events):
                event = event_generator._generate_single_event()
                events.append(event)
            
            # Check if THIS node is the WGD node (Equation 19)
            # If so, prepend WGD event (applied BEFORE other events)
            if self.wgd_occurred and node.node_id == self.wgd_node:
                from .data_models import CNAEvent, Haplotype
                wgd = CNAEvent(
                    start_bin=0,
                    end_bin=self.config.num_bins - 1,
                    haplotype=Haplotype.WGD,
                    amplitude=0
                )
                events = [wgd] + events  # WGD applied FIRST
            
            # Apply events to get child CN profile
            cn_A, cn_B = event_applier.apply_events(parent_clone, events)
            
            # Create clone for this node
            child_clone = Clone(
                index=len(self.clones),
                parent_index=node_id_to_clone_idx[parent_node.node_id],
                cn_profile_A=cn_A,
                cn_profile_B=cn_B,
                events=events,
                is_root=False
            )
            
            self.clones.append(child_clone)
            node_id_to_clone_idx[node.node_id] = child_clone.index
        
        # Map leaf nodes to final clone indices for sampling
        leaf_nodes = [n for n in tree_nodes if n.is_leaf]
        self._leaf_clone_indices = [node_id_to_clone_idx[n.node_id] for n in sorted(leaf_nodes, key=lambda n: n.node_id)]
        
        total_events = sum(len(clone.events) for clone in self.clones)
        num_leaves = len(leaf_nodes)
        
        print(f"  Created {len(self.clones)} clones ({num_leaves} leaves)")
        print(f"  Total CNA events: {total_events}")
        print(f"  Clone proportions: {clone_proportions}")


    def _sample_wgd_placement(self):
        """
        Sample WGD placement according to Equations 18-19.
        
        Returns:
            Node ID where WGD occurs, or None if no WGD
        """
        # Step 1: Check if WGD occurs (probability p_WGD)
        if self.rng.random() >= self.config.prob_wgd:
            return None
        
        # Step 2: Find early edges (depth <= d_WGD)
        early_edges = []
        for node in self._tree_nodes:
            if node.node_id == 0:  # Skip root
                continue
            
            # Compute depth (number of edges from root)
            depth = self._get_node_depth(node)
            
            if depth <= self.config.d_wgd:
                early_edges.append(node)
        
        if not early_edges:
            return None
        
        # Step 3: Sample edge proportional to branch length (Equation 18)
        branch_lengths = np.array([node.edge_length for node in early_edges])
        probs = branch_lengths / branch_lengths.sum()
        
        selected_idx = self.rng.choice(len(early_edges), p=probs)
        wgd_node = early_edges[selected_idx]
        
        return wgd_node.node_id
    
    def _get_node_depth(self, node):
        """
        Compute depth of a node (number of edges from root).
        """
        depth = 0
        current = node
        nodes_by_id = {n.node_id: n for n in self._tree_nodes}
        
        # Walk up to root (parent_id = -1 or None means root)
        while current.parent_id not in (-1, None):
            depth += 1
            current = nodes_by_id[current.parent_id]
        
        return depth

    def _sample_cells(self):
        """
        Sample cells from clones with optional normal/doublet contamination.
        
        For Week 5: uniform sampling from clones, no contamination.
        """
        # Use clone frequencies from Beta-splitting tree
        # Only sample from leaf clones (extant clones)
        if hasattr(self, '_clone_proportions'):
            clone_frequencies = self._clone_proportions
        else:
            # Fallback to uniform (shouldn't happen)
            clone_frequencies = np.ones(self.config.num_clones) / self.config.num_clones
        
        # Sample library sizes and allelic coverage factors
        library_sizes = self._sample_library_sizes(self.config.num_cells)
        allelic_coverages = self._sample_allelic_coverages(self.config.num_cells)
        
        # Assign cells to leaf clones only
        if hasattr(self, '_leaf_clone_indices'):
            # Sample from leaf clones using Beta-splitting proportions
            leaf_indices = self.rng.choice(
                len(self._leaf_clone_indices),
                size=self.config.num_cells,
                p=clone_frequencies
            )
            # Map to actual clone indices
            clone_assignments = np.array([self._leaf_clone_indices[i] for i in leaf_indices])
        else:
            # Fallback
            clone_assignments = self.rng.choice(
                self.config.num_clones,
                size=self.config.num_cells,
                p=clone_frequencies
            )
        
        # Create cell objects with doublet sampling
        self.cells = []
        num_doublets = 0
        
        for n in range(self.config.num_cells):
            # Check if this cell is a doublet
            is_doublet = self.rng.random() < self.config.prob_doublet
            
            if is_doublet:
                # Sample two clones independently
                k1 = self.rng.choice(len(self._leaf_clone_indices), p=self._clone_proportions)
                k2 = self.rng.choice(len(self._leaf_clone_indices), p=self._clone_proportions)
                
                # Map to actual clone indices
                clone1 = self._leaf_clone_indices[k1]
                clone2 = self._leaf_clone_indices[k2]
                
                doublet_pair = (clone1, clone2)
                num_doublets += 1
                
                # For clone assignment, use first clone
                clone_assign = clone1
            else:
                doublet_pair = None
                clone_assign = clone_assignments[n]
            
            cell = Cell(
                index=n,
                clone_assignment=clone_assign,
                library_size=library_sizes[n],
                allelic_coverage=allelic_coverages[n],
                is_normal=False,
                is_doublet=is_doublet,
                doublet_clones=doublet_pair
            )
            self.cells.append(cell)
        
        if num_doublets > 0:
            print(f"  Sampled {num_doublets} doublets ({100*num_doublets/self.config.num_cells:.1f}%)")
        
        print(f"  Sampled {len(self.cells)} cells")
        print(f"  Clone distribution: {np.bincount(clone_assignments)}")
    
    def _sample_library_sizes(self, num_cells: int) -> np.ndarray:
        """
        Sample library size factors (α_n) for cells.
        
        Uses log-normal distribution to ensure positive values.
        """
        mean = self.config.mean_library_size
        cv = self.config.library_size_cv
        
        # Log-normal parameters
        sigma = np.sqrt(np.log(1 + cv**2))
        mu = np.log(mean) - 0.5 * sigma**2
        
        return self.rng.lognormal(mu, sigma, size=num_cells)
    
    def _sample_allelic_coverages(self, num_cells: int) -> np.ndarray:
        """
        Sample allelic coverage factors (β_n) for cells.
        
        Uses log-normal distribution to ensure positive values.
        """
        mean = self.config.mean_allelic_coverage
        cv = self.config.allelic_coverage_cv
        
        # Log-normal parameters
        sigma = np.sqrt(np.log(1 + cv**2))
        mu = np.log(mean) - 0.5 * sigma**2
        
        return self.rng.lognormal(mu, sigma, size=num_cells)
    
    
    def _sample_bin_biases(self) -> np.ndarray:
        """
        Sample bin bias factors (κ_b) for GC/mappability effects.
        
        Uses log-normal distribution centered at 1.0.
        
        Returns:
            Array of length B with bias factors (mean = 1.0)
        """
        if not self.config.use_gc_bias:
            # No bias
            return np.ones(self.config.num_bins)
        
        # Sample from log-normal and normalize to have mean 1
        biases = self.rng.lognormal(0, 0.3, size=self.config.num_bins)
        biases = biases / biases.mean()  # Normalize so mean = 1
        
        return biases

    def _generate_read_counts(self) -> np.ndarray:
        """
        Generate read-depth counts (x_{n,b}) via negative binomial model.
        
        Returns:
            Array of shape (N, B) with read counts
        """
        N = self.config.num_cells
        B = self.config.num_bins
        read_counts = np.zeros((N, B), dtype=int)
        
        for n, cell in enumerate(self.cells):
            clone = self.clones[cell.clone_assignment]
            
            for b in range(B):
                # Total copy number at this bin
                tcn = clone.cn_profile_A[b] + clone.cn_profile_B[b]
                
                # Bin bias (κ_b)
                if not hasattr(self, 'bin_biases'):
                    self.bin_biases = self._sample_bin_biases()
                kappa_b = self.bin_biases[b]
                
                # Mean read count: μ_{n,b} = α_n * κ_b * (TCN / 2)
                mu = cell.library_size * kappa_b * (tcn / 2.0)
                
                # Negative binomial: variance = μ + μ²/θ
                # We use the (n, p) parameterization where:
                #   n = θ, p = θ/(θ + μ)
                if mu > 0:
                    n_param = self.config.theta_x
                    p_param = n_param / (n_param + mu)
                    read_counts[n, b] = self.rng.negative_binomial(n_param, p_param)
                else:
                    read_counts[n, b] = 0
            
            # Store in cell object
            cell.read_counts = read_counts[n, :]
        
        return read_counts
    
    def _generate_allele_counts(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate allele counts (a_{n,s}, r_{n,s}) via Beta-Binomial model.
        
        Returns:
            Tuple of (alternate, reference, total) arrays, each of shape (N, S)
        """
        N = self.config.num_cells
        S = len(self.segments)
        
        alternate_counts = np.zeros((N, S), dtype=int)
        reference_counts = np.zeros((N, S), dtype=int)
        total_counts = np.zeros((N, S), dtype=int)
        
        for n, cell in enumerate(self.cells):
            for s, segment in enumerate(self.segments):
                # Sample number of heterozygous SNPs in this segment
                M_s = self.rng.poisson(self.config.snp_density * segment.length)
                
                if M_s == 0:
                    continue
                
                # Handle doublets differently (Equation 15)
                if cell.is_doublet:
                    k1, k2 = cell.doublet_clones
                    clone1 = self.clones[k1]
                    clone2 = self.clones[k2]
                    
                    # Get CN for both clones
                    cn_A1, cn_B1, tcn1 = clone1.get_segment_cn(segment)
                    cn_A2, cn_B2, tcn2 = clone2.get_segment_cn(segment)
                    
                    hap_block = self.haplotype_blocks[segment.haplotype_block]
                    
                    # Get alternate CN for each clone
                    if hap_block.alternate_haplotype == Haplotype.A:
                        cn_alt1 = cn_A1
                        cn_alt2 = cn_A2
                    else:
                        cn_alt1 = cn_B1
                        cn_alt2 = cn_B2
                    
                    # CN-weighted mixture (Equation 15)
                    epsilon = 1e-10
                    p_doublet = (cn_alt1 + cn_alt2) / (tcn1 + tcn2 + epsilon)
                    
                    # Total depth for doublet (sum of both clones)
                    tcn_doublet = tcn1 + tcn2
                    if tcn_doublet > 0:
                        t_ns = self.rng.poisson(cell.allelic_coverage * M_s * (tcn_doublet / 2.0))
                    else:
                        t_ns = 0
                    
                    if t_ns == 0:
                        continue
                    
                    # Apply phase orientation
                    if hap_block.orientation == -1:
                        p_doublet = 1.0 - p_doublet
                    
                    p_alt = p_doublet
                    
                else:
                    # Singlet cell
                    clone = self.clones[cell.clone_assignment]
                    
                    # Get average copy numbers over segment
                    cn_A, cn_B, tcn = clone.get_segment_cn(segment)
                    
                    # Total allelic depth: Poisson(β_n * M_s * TCN / 2)
                    if tcn > 0:
                        t_ns = self.rng.poisson(cell.allelic_coverage * M_s * (tcn / 2.0))
                    else:
                        t_ns = 0
                    
                    if t_ns == 0:
                        continue
                    
                    # Expected alternate fraction
                    hap_block = self.haplotype_blocks[segment.haplotype_block]
                    if hap_block.alternate_haplotype == Haplotype.A:
                        p_alt = cn_A / (tcn + 1e-10)
                    else:
                        p_alt = cn_B / (tcn + 1e-10)
                    
                    # Apply phase orientation
                    if hap_block.orientation == -1:
                        p_alt = 1.0 - p_alt
                
                # Apply LOH floor to prevent degenerate Beta (Equation 50)
                p_alt = np.clip(p_alt, self.config.epsilon_floor, 1.0 - self.config.epsilon_floor)
                
                # Add sequencing error contamination (Equation 52)
                p_alt = (1.0 - self.config.epsilon_seq) * p_alt + self.config.epsilon_seq * 0.5
                
                # Beta-Binomial: sample q ~ Beta(p*ν, (1-p)*ν), then a ~ Binomial(t, q)
                alpha_beta = p_alt * self.config.nu_a
                beta_beta = (1.0 - p_alt) * self.config.nu_a
                
                q = self.rng.beta(alpha_beta, beta_beta)
                a_ns = self.rng.binomial(t_ns, q)
                r_ns = t_ns - a_ns
                
                alternate_counts[n, s] = a_ns
                reference_counts[n, s] = r_ns
                total_counts[n, s] = t_ns
            
            # Store in cell object
            cell.allele_counts = (
                alternate_counts[n, :],
                reference_counts[n, :],
                total_counts[n, :]
            )
        
        return alternate_counts, reference_counts, total_counts
    

    def _detect_segments_from_clones(self):
        """
        Detect segment boundaries from clone CN profiles.
        
        This is called AFTER clone generation, so segments reflect
        the actual CNA breakpoints.
        """
        from .segment_detector import SegmentDetector
        
        detector = SegmentDetector(num_bins=self.config.num_bins)
        self.segments = detector.detect_segments(self.clones)
        
        # Update bin_to_segment mapping
        self.bin_to_segment = {}
        for segment in self.segments:
            for bin_idx in segment.bin_indices:
                self.bin_to_segment[bin_idx] = segment.index
        
        # Create haplotype blocks with phase-switch process
        self.haplotype_blocks = []
        self._create_haplotype_blocks_with_phase_switches()


    def _create_haplotype_blocks_with_phase_switches(self):
        """
        Create haplotype blocks at chromosome-arm level.
        
        Section 3.1: "We define haplotype blocks at the chromosome-arm level:
        each chromosome arm constitutes one block"
        
        For single chromosome: H = 2 (p-arm and q-arm)
        Phase orientation sampled once per arm with probability p_switch.
        """
        from .chromosome_data import get_arm_boundaries
        
        if len(self.segments) == 0:
            return
        
        # Get chromosome arm boundaries
        p_end, q_start = get_arm_boundaries(self.config.chromosome)
        p_arm_end_bin = p_end // self.config.bin_width
        q_arm_start_bin = q_start // self.config.bin_width
        
        # Create two blocks: one for p-arm, one for q-arm
        self.haplotype_blocks = []
        
        for arm_idx in range(2):  # 0=p-arm, 1=q-arm
            # Sample orientation (Equation 51)
            # η_h = -1 with probability p_switch, else +1
            if self.rng.random() < self.config.prob_phase_switch:
                orientation = -1
            else:
                orientation = +1
            
            # Default: alternate haplotype is A
            alternate_hap = Haplotype.A
            
            hap_block = HaplotypeBlock(
                index=arm_idx,
                segment_indices=[],
                orientation=orientation,
                alternate_haplotype=alternate_hap
            )
            self.haplotype_blocks.append(hap_block)
        
        # Assign each segment to its arm block
        for segment in self.segments:
            # Check if segment is in p-arm or q-arm based on its position
            if segment.end_bin < p_arm_end_bin:
                block_idx = 0  # p-arm
            else:
                block_idx = 1  # q-arm
            
            segment.haplotype_block = block_idx
            self.haplotype_blocks[block_idx].segment_indices.append(segment.index)
            self.segment_to_block[segment.index] = block_idx

    def get_ground_truth(self) -> Dict:
        """
        Extract ground truth data for evaluation.
        
        Returns:
            Dictionary containing:
                - clone_tree: Adjacency list representation of tree
                - clone_assignments: Array of shape (N,) with cell->clone mapping
                - cn_profiles_A: Array of shape (K, B) with haplotype A CNs
                - cn_profiles_B: Array of shape (K, B) with haplotype B CNs
                - segments: List of segment objects
                - events: Dict mapping clone_idx -> list of events
        """
        # Clone tree as adjacency list
        clone_tree = {}
        for clone in self.clones:
            if not clone.is_root:
                parent = clone.parent_index
                if parent not in clone_tree:
                    clone_tree[parent] = []
                clone_tree[parent].append(clone.index)
        
        # Clone assignments
        clone_assignments = np.array([cell.clone_assignment for cell in self.cells])
        
        # CN profiles
        cn_profiles_A = np.vstack([clone.cn_profile_A for clone in self.clones])
        cn_profiles_B = np.vstack([clone.cn_profile_B for clone in self.clones])
        
        # Events
        events = {clone.index: clone.events for clone in self.clones}
        
        return {
            "clone_tree": clone_tree,
            "clone_assignments": clone_assignments,
            "cn_profiles_A": cn_profiles_A,
            "cn_profiles_B": cn_profiles_B,
            "segments": self.segments,
            "events": events,
        }
