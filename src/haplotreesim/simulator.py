"""
HaploTreeSim simulator core.

This module implements the main simulator class that generates synthetic
scDNA-seq data with haplotype-specific CNAs and clone tree ground truth.
"""

import numpy as np
import json
import pickle
from typing import List, Tuple, Dict, Optional
from scipy import sparse
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
        self.chromosome_bin_offsets: Dict[str, Tuple[int, int]] = {}
        self.chromosome_boundary_bins: List[int] = []
        self.chromosome_intervals: List[Tuple[str, int, int]] = []
        self.chromosome_arm_intervals: List[Tuple[str, str, int, int]] = []
        
        # Mappings
        self.bin_to_segment: Dict[int, int] = {}  # bin index -> segment index
        self.bin_to_block: Dict[int, int] = {}  # bin index -> haplotype block index
        self.segment_to_block: Dict[int, int] = {}  # segment index -> block index
        self.bin_snp_counts: Optional[np.ndarray] = None  # M_b for allelic simulation
        self.snp_bin_offsets: Optional[np.ndarray] = None
        self.snp_bins: Optional[np.ndarray] = None
        self.snp_positions: Optional[np.ndarray] = None
        self.snp_chromosomes: Optional[np.ndarray] = None
        self.snp_bin_local_indices: Optional[np.ndarray] = None
        self.snp_alt_counts = sparse.csr_matrix((0, config.num_cells), dtype=np.int32)
        self.snp_ref_counts = sparse.csr_matrix((0, config.num_cells), dtype=np.int32)

    def _create_event_generator(self):
        """Create an event generator using the simulator's current genome layout."""
        from .event_generator import EventGenerator

        return EventGenerator(
            rng=self.rng,
            num_bins=self.config.num_bins,
            bin_length=self.config.bin_width,
            chromosome=self.config.chromosome,
            chromosomes=self.config.chromosomes,
            chromosome_intervals=self.chromosome_intervals,
            arm_intervals=self.chromosome_arm_intervals,
            lambda_events=self.config.lambda_events,
            lambda_amplitude=self.config.lambda_amplitude,
            prob_wgd=self.config.prob_wgd,
            gain_prob=self.config.gain_prob,
            prob_focal=self.config.prob_focal,
            prob_arm_given_broad=self.config.prob_arm_given_broad,
            focal_length_min=self.config.focal_length_min,
            focal_length_max_fraction=self.config.focal_length_max_fraction,
            prob_haplotype_A=self.config.prob_haplotype_A,
            max_copy_number=self.config.max_copy_number,
            max_event_attempts=self.config.max_event_attempts,
        )
    

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
                'interval': [float(node.interval[0]), float(node.interval[1])],
                'perc': float(node.perc),
                'edge_length': float(node.edge_length),
                'depth': int(node.depth),
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

        if tree_data.get('num_clones') is not None:
            # A fixed tree defines its extant leaf count.
            self.config.num_clones = int(tree_data['num_clones'])
        
        # Simply store the tree data - we'll recreate nodes later
        self._imported_tree_data = tree_data
        self._clone_proportions = np.array(tree_data['clone_proportions'])
        self._leaf_clone_indices = tree_data['leaf_clone_indices']
        
        # Create minimal tree node objects for compatibility
        class SimpleTreeNode:
            def __init__(self, node_id, parent_id, edge_length, is_leaf):
                self.node_id = node_id
                self.parent_id = parent_id
                self.edge_length = edge_length
                self.is_leaf = is_leaf
        
        self._tree_nodes = []
        for node_data in tree_data['nodes']:
            node = SimpleTreeNode(
                node_id=node_data['node_id'],
                parent_id=node_data['parent_id'],
                interval=tuple(node_data.get('interval', (0.0, 0.0))),
                perc=node_data.get('perc', 0.0),
                edge_length=node_data['edge_length'],
                depth=node_data.get('depth', 0),
                is_leaf=node_data['is_leaf']
            )
            self._tree_nodes.append(node)
        
        print(f"✓ Tree structure imported from: {filepath}")
        print(f"  Nodes: {len(self._tree_nodes)}, Leaves: {len([n for n in self._tree_nodes if n.is_leaf])}")

    def _generate_ordered_edge_events(
        self,
        event_generator,
        parent_clone: Clone,
        node,
        include_wgd: bool
    ) -> List[CNAEvent]:
        """
        Generate an ordered, state-dependent event list for one tree edge.
        """
        events: List[CNAEvent] = []
        state_A = parent_clone.cn_profile_A.copy()
        state_B = parent_clone.cn_profile_B.copy()

        if include_wgd:
            wgd = CNAEvent(
                start_bin=0,
                end_bin=self.config.num_bins - 1,
                haplotype=Haplotype.WGD,
                amplitude=0,
                event_time=0.0,
                scale_class="wgd"
            )
            events.append(wgd)
            state_A = np.clip(state_A * 2, 0, self.config.max_copy_number)
            state_B = np.clip(state_B * 2, 0, self.config.max_copy_number)

        num_proposals = self.rng.poisson(self.config.lambda_events * node.edge_length)
        events.extend(
            event_generator.generate_valid_events(
                num_proposals=num_proposals,
                edge_length=node.edge_length,
                cn_A=state_A,
                cn_B=state_B
            )
        )

        return events


    def _initialize_clones_from_imported_tree(self):
        """
        Initialize clones from imported tree structure.
        Called after import_tree_structure() to continue simulation.
        """
        from .event_applier import EventApplier
        
        # Initialize event infrastructure
        event_generator = self._create_event_generator()
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
            
            events = self._generate_ordered_edge_events(
                event_generator=event_generator,
                parent_clone=parent_clone,
                node=node,
                include_wgd=self.wgd_occurred and node.node_id == self.wgd_node
            )
            
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

    def run(
        self,
        tree_structure_path: Optional[str] = None
    ) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Run the full simulation pipeline.

        Args:
            tree_structure_path: Optional exported tree structure to reuse.
        
        Returns:
            Tuple of:
                - read_counts: Array of shape (N, B) with read counts
                - allele_counts: Tuple of (alternate, reference, total) arrays,
                                each of shape (N, B)
        """
        print("Initializing genome...")
        self._initialize_genome()

        self.wgd_occurred = False
        self.wgd_node = None

        if tree_structure_path is None:
            print(f"Generating clone tree with {self.config.num_clones} clones...")
            self._generate_clone_tree()

            # Compute ploidy for each clone (Section 2: auxiliary truth)
            for clone in self.clones:
                clone.ploidy = float(np.mean(clone.total_cn()))

            print("Detecting segments from CNA breakpoints...")
            self._detect_segments_from_clones()
        else:
            print(f"Importing clone tree from {tree_structure_path}...")
            self.import_tree_structure(tree_structure_path)
            self._initialize_clones_from_imported_tree()
        
        print(f"Sampling {self.config.num_cells} cells...")
        self._sample_cells()
        
        print("Generating read-depth observations...")
        read_counts = self._generate_read_counts()
        
        print("Generating allelic observations...")
        allele_counts = self._generate_allele_counts()
        
        print("Simulation complete!")
        return read_counts, allele_counts

    def run_with_tree_structure(
        self,
        filepath: str
    ) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Run the simulator using a previously exported fixed tree."""
        return self.run(tree_structure_path=filepath)
    
    def _initialize_genome(self):
        """
        Initialize genome representation (bins, segments, haplotype blocks).
        """
        from .chromosome_data import create_bins_for_chromosome, get_chromosome_length

        self.bins = []
        self.segments = []
        self.chromosome_bin_offsets = {}
        self.chromosome_boundary_bins = []
        self.chromosome_intervals = []

        computed_counts = {
            chromosome: create_bins_for_chromosome(chromosome, self.config.bin_width)
            for chromosome in self.config.chromosomes
        }
        computed_total = sum(computed_counts.values())
        if len(self.config.chromosomes) == 1 and self.config.num_bins != computed_total:
            chromosome_counts = {self.config.chromosomes[0]: self.config.num_bins}
            use_reference_lengths = False
        elif self.config.num_bins == computed_total:
            chromosome_counts = computed_counts
            use_reference_lengths = True
        else:
            raise ValueError(
                "Explicit num_bins with multiple chromosomes is ambiguous; "
                "omit num_bins so it can be computed from chromosomes and bin_width."
            )

        global_bin = 0
        for chromosome in self.config.chromosomes:
            chrom_first_bin = global_bin
            chrom_length = get_chromosome_length(chromosome)
            bin_count = chromosome_counts[chromosome]

            if chrom_first_bin > 0:
                self.chromosome_boundary_bins.append(chrom_first_bin)

            for local_bin in range(bin_count):
                start = local_bin * self.config.bin_length
                if use_reference_lengths:
                    end = min(start + self.config.bin_length, chrom_length)
                else:
                    end = start + self.config.bin_length

                if end <= start:
                    continue

                bin_obj = Bin(
                    index=global_bin,
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    length=end - start
                )
                self.bins.append(bin_obj)
                global_bin += 1

            chrom_last_bin = global_bin - 1
            self.chromosome_bin_offsets[chromosome] = (chrom_first_bin, chrom_last_bin)
            self.chromosome_intervals.append((chromosome, chrom_first_bin, chrom_last_bin))

            segment = Segment(
                index=len(self.segments),
                bin_indices=set(range(chrom_first_bin, chrom_last_bin + 1)),
                start_bin=chrom_first_bin,
                end_bin=chrom_last_bin,
                length=sum(bin_obj.length for bin_obj in self.bins[chrom_first_bin:chrom_last_bin + 1]),
                haplotype_block=0
            )
            self.segments.append(segment)

        self.chromosome_arm_intervals = self._build_chromosome_arm_intervals()
        
        # Create a single haplotype block
        hap_block = HaplotypeBlock(
            index=0,
            segment_indices=[segment.index for segment in self.segments],
            orientation=1,
            alternate_haplotype=Haplotype.A
        )
        self.haplotype_blocks = [hap_block]
        
        # Build mappings
        self.bin_to_segment = {}
        self.bin_to_block = {}
        self.segment_to_block = {}
        for segment in self.segments:
            for bin_idx in segment.bin_indices:
                self.bin_to_segment[bin_idx] = segment.index
                self.bin_to_block[bin_idx] = 0
            self.segment_to_block[segment.index] = 0
        
        print(f"  Created {len(self.bins)} bins, {len(self.segments)} segments, "
              f"{len(self.haplotype_blocks)} haplotype blocks")

    def _build_chromosome_arm_intervals(self) -> List[Tuple[str, str, int, int]]:
        """Build chromosome-arm intervals in global bin coordinates."""
        from .chromosome_data import get_arm_boundaries

        arm_intervals: List[Tuple[str, str, int, int]] = []

        for chromosome, (chrom_start, chrom_end) in self.chromosome_bin_offsets.items():
            p_end, _ = get_arm_boundaries(chromosome)
            p_bins = [
                bin_obj.index
                for bin_obj in self.bins[chrom_start:chrom_end + 1]
                if (bin_obj.start + bin_obj.end) / 2 <= p_end
            ]
            q_bins = [
                bin_obj.index
                for bin_obj in self.bins[chrom_start:chrom_end + 1]
                if (bin_obj.start + bin_obj.end) / 2 > p_end
            ]

            if p_bins:
                arm_intervals.append((chromosome, "p", min(p_bins), max(p_bins)))
            if q_bins:
                arm_intervals.append((chromosome, "q", min(q_bins), max(q_bins)))

        return arm_intervals
    
    def _generate_clone_tree(self):
        """
        Generate clone tree with haplotype-specific copy number profiles.
        
        Updated: Now uses Beta-splitting tree model (Section 3.2).
        """
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
        event_generator = self._create_event_generator()
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
            
            events = self._generate_ordered_edge_events(
                event_generator=event_generator,
                parent_clone=parent_clone,
                node=node,
                include_wgd=self.wgd_occurred and node.node_id == self.wgd_node
            )
            
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
        
        # Create cell objects with contamination sampling
        self.cells = []
        num_normals = 0
        num_doublets = 0
        
        for n in range(self.config.num_cells):
            is_normal = self.rng.random() < self.config.prob_normal
            is_doublet = False

            if is_normal:
                # Normal cells retain the root diploid genotype.
                doublet_pair = None
                clone_assign = 0
                num_normals += 1
            elif self.rng.random() < self.config.prob_doublet:
                # Sample two clones independently
                is_doublet = True
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
                is_normal=is_normal,
                is_doublet=is_doublet,
                doublet_clones=doublet_pair
            )
            self.cells.append(cell)
        
        if num_normals > 0:
            print(f"  Sampled {num_normals} normal cells ({100*num_normals/self.config.num_cells:.1f}%)")
        if num_doublets > 0:
            print(f"  Sampled {num_doublets} doublets ({100*num_doublets/self.config.num_cells:.1f}%)")
        
        print(f"  Sampled {len(self.cells)} cells")
        sampled_assignments = np.array([cell.clone_assignment for cell in self.cells])
        print(f"  Clone distribution: {np.bincount(sampled_assignments)}")
    
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
            for b in range(B):
                # Bin bias (κ_b)
                if not hasattr(self, 'bin_biases'):
                    self.bin_biases = self._sample_bin_biases()
                kappa_b = self.bin_biases[b]

                if cell.is_doublet:
                    k1, k2 = cell.doublet_clones
                    read_counts[n, b] = (
                        self._sample_clone_read_count(cell, self.clones[k1], b, kappa_b)
                        + self._sample_clone_read_count(cell, self.clones[k2], b, kappa_b)
                    )
                else:
                    clone = self.clones[cell.clone_assignment]
                    read_counts[n, b] = self._sample_clone_read_count(cell, clone, b, kappa_b)
            
            # Store in cell object
            cell.read_counts = read_counts[n, :]
        
        return read_counts

    def _sample_clone_read_count(self, cell: Cell, clone: Clone, bin_index: int, bin_bias: float) -> int:
        """Sample one clone contribution to a cell-bin read-depth count."""
        tcn = clone.cn_profile_A[bin_index] + clone.cn_profile_B[bin_index]
        mu = cell.library_size * bin_bias * (tcn / 2.0)

        if mu <= 0:
            return 0

        # NumPy uses the NB (n, p) parameterization for mean mu and size theta.
        n_param = self.config.theta_x
        p_param = n_param / (n_param + mu)
        return int(self.rng.negative_binomial(n_param, p_param))
    
    def _generate_segment_allele_counts_legacy(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Deprecated segment-level allele-count generator retained for reference.
        
        Returns:
            Tuple of deprecated segment-level arrays.
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
                
                # Deprecated segment-level emission path.
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
    

    def _sample_snp_alternate_counts(self, total_depths: np.ndarray, p_alt: float) -> np.ndarray:
        """
        Sample alternate reads for one cell-bin's SNP-level observations.

        The default Beta-Binomial concentration is per-SNP when
        allelic_concentration_scale is "snps"; summing many SNPs then naturally
        increases evidence at the bin level. A sensitivity mode scales by each
        SNP's observed allelic depth, and a debug mode removes the Beta layer.
        """
        p_alt = np.clip(p_alt, self.config.epsilon_floor, 1.0 - self.config.epsilon_floor)
        p_alt = (1.0 - self.config.epsilon_seq) * p_alt + self.config.epsilon_seq * 0.5
        total_depths = np.asarray(total_depths, dtype=int)
        alternate = np.zeros(total_depths.shape, dtype=np.int32)
        positive = total_depths > 0
        if not np.any(positive):
            return alternate

        if self.config.disable_allelic_beta:
            q_alt = np.full(int(np.sum(positive)), p_alt, dtype=float)
        else:
            if self.config.allelic_concentration_scale == "depth":
                evidence = total_depths[positive].astype(float)
            else:
                evidence = np.ones(int(np.sum(positive)), dtype=float)
            nu_eff = self.config.nu_a * evidence
            alpha_beta = p_alt * nu_eff
            beta_beta = (1.0 - p_alt) * nu_eff
            q_alt = self.rng.beta(alpha_beta, beta_beta)

        alternate[positive] = self.rng.binomial(total_depths[positive], q_alt)
        return alternate

    def _generate_snp_metadata(self) -> None:
        """Generate usable heterozygous SNP positions and bin mappings."""
        B = self.config.num_bins
        self.bin_snp_counts = self.rng.poisson(
            [self.config.snp_density * bin_obj.length for bin_obj in self.bins]
        ).astype(np.int32)
        self.snp_bin_offsets = np.empty(B + 1, dtype=np.int64)
        self.snp_bin_offsets[0] = 0
        np.cumsum(self.bin_snp_counts, out=self.snp_bin_offsets[1:])
        total_snps = int(self.snp_bin_offsets[-1])

        self.snp_bins = np.repeat(np.arange(B, dtype=np.int32), self.bin_snp_counts)
        self.snp_positions = np.empty(total_snps, dtype=np.int64)
        self.snp_bin_local_indices = np.empty(total_snps, dtype=np.int32)
        self.snp_chromosomes = np.empty(total_snps, dtype=object)

        for b, bin_obj in enumerate(self.bins):
            start = int(self.snp_bin_offsets[b])
            end = int(self.snp_bin_offsets[b + 1])
            snp_count = end - start
            if snp_count == 0:
                continue

            replace = snp_count > bin_obj.length
            offsets = self.rng.choice(bin_obj.length, size=snp_count, replace=replace)
            self.snp_positions[start:end] = bin_obj.start + 1 + np.sort(offsets)
            self.snp_bin_local_indices[start:end] = np.arange(snp_count, dtype=np.int32)
            self.snp_chromosomes[start:end] = bin_obj.chromosome

    def _cell_bin_allelic_parameters(self, cell: Cell, bin_index: int) -> tuple[float, float]:
        """Return expected alternate fraction and total copy number for one cell-bin."""
        hap_block = self.haplotype_blocks[self.bin_to_block.get(bin_index, 0)]

        if cell.is_doublet:
            k1, k2 = cell.doublet_clones
            clone1 = self.clones[k1]
            clone2 = self.clones[k2]
            cn_A1 = clone1.cn_profile_A[bin_index]
            cn_B1 = clone1.cn_profile_B[bin_index]
            cn_A2 = clone2.cn_profile_A[bin_index]
            cn_B2 = clone2.cn_profile_B[bin_index]
            tcn = cn_A1 + cn_B1 + cn_A2 + cn_B2
            if tcn == 0:
                return 0.5, 0.0
            cn_alt = cn_A1 + cn_A2 if hap_block.alternate_haplotype == Haplotype.A else cn_B1 + cn_B2
        else:
            clone = self.clones[cell.clone_assignment]
            cn_A = clone.cn_profile_A[bin_index]
            cn_B = clone.cn_profile_B[bin_index]
            tcn = cn_A + cn_B
            if tcn == 0:
                return 0.5, 0.0
            cn_alt = cn_A if hap_block.alternate_haplotype == Haplotype.A else cn_B

        p_alt = cn_alt / (tcn + 1e-10)
        if hap_block.orientation == -1:
            p_alt = 1.0 - p_alt
        return float(p_alt), float(tcn)

    def _generate_allele_counts(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate bin-level allele counts (a_{n,b}, r_{n,b}).

        Returns:
            Tuple of (alternate_counts, reference_counts, total_counts), each
            with shape (N, B), where B is the number of bins.
        """
        N = self.config.num_cells
        B = self.config.num_bins

        alternate_counts = np.zeros((N, B), dtype=np.int32)
        reference_counts = np.zeros((N, B), dtype=np.int32)
        total_counts = np.zeros((N, B), dtype=np.int32)
        self._generate_snp_metadata()
        total_snps = int(self.snp_bin_offsets[-1])
        alt_rows = []
        alt_cols = []
        alt_data = []
        ref_rows = []
        ref_cols = []
        ref_data = []

        for n, cell in enumerate(self.cells):
            for b, snp_count in enumerate(self.bin_snp_counts):
                if snp_count == 0:
                    continue

                p_alt, tcn = self._cell_bin_allelic_parameters(cell, b)
                if tcn == 0:
                    continue

                depth_mean = cell.allelic_coverage * snp_count * (tcn / 2.0)
                total_depth = int(self.rng.poisson(depth_mean))
                if total_depth == 0:
                    continue

                snp_depths = self.rng.multinomial(total_depth, np.full(int(snp_count), 1.0 / snp_count))
                alt_depths = self._sample_snp_alternate_counts(snp_depths, p_alt)
                ref_depths = snp_depths - alt_depths

                snp_start = int(self.snp_bin_offsets[b])
                alt_nonzero = alt_depths > 0
                if np.any(alt_nonzero):
                    snp_rows = snp_start + np.where(alt_nonzero)[0]
                    alt_rows.append(snp_rows.astype(np.int64, copy=False))
                    alt_cols.append(np.full(len(snp_rows), n, dtype=np.int32))
                    alt_data.append(alt_depths[alt_nonzero].astype(np.int32, copy=False))

                ref_nonzero = ref_depths > 0
                if np.any(ref_nonzero):
                    snp_rows = snp_start + np.where(ref_nonzero)[0]
                    ref_rows.append(snp_rows.astype(np.int64, copy=False))
                    ref_cols.append(np.full(len(snp_rows), n, dtype=np.int32))
                    ref_data.append(ref_depths[ref_nonzero].astype(np.int32, copy=False))

                alt_depth = int(alt_depths.sum())
                ref_depth = int(ref_depths.sum())
                alternate_counts[n, b] = alt_depth
                reference_counts[n, b] = ref_depth
                total_counts[n, b] = total_depth

            cell.allele_counts = (
                alternate_counts[n, :],
                reference_counts[n, :],
                total_counts[n, :]
            )

        if alt_data:
            self.snp_alt_counts = sparse.coo_matrix(
                (
                    np.concatenate(alt_data),
                    (np.concatenate(alt_rows), np.concatenate(alt_cols)),
                ),
                shape=(total_snps, N),
                dtype=np.int32,
            ).tocsr()
        else:
            self.snp_alt_counts = sparse.csr_matrix((total_snps, N), dtype=np.int32)

        if ref_data:
            self.snp_ref_counts = sparse.coo_matrix(
                (
                    np.concatenate(ref_data),
                    (np.concatenate(ref_rows), np.concatenate(ref_cols)),
                ),
                shape=(total_snps, N),
                dtype=np.int32,
            ).tocsr()
        else:
            self.snp_ref_counts = sparse.csr_matrix((total_snps, N), dtype=np.int32)
        return alternate_counts, reference_counts, total_counts


    def _detect_segments_from_clones(self):
        """
        Detect segment boundaries from clone CN profiles.
        
        This is called AFTER clone generation, so segments reflect
        the actual CNA breakpoints.
        """
        from .segment_detector import SegmentDetector
        
        detector = SegmentDetector(
            num_bins=self.config.num_bins,
            bin_width=self.config.bin_width,
            mandatory_breakpoints=self.chromosome_boundary_bins,
            bin_lengths=[bin_obj.length for bin_obj in self.bins],
        )
        leaf_clones = [self.clones[i] for i in self._leaf_clone_indices]
        self.segments = detector.detect_segments(leaf_clones)
        
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
        if len(self.segments) == 0:
            return

        self.haplotype_blocks = []
        self.bin_to_block = {}
        self.segment_to_block = {}

        for chromosome, arm, start_bin, end_bin in self.chromosome_arm_intervals:
            if self.rng.random() < self.config.prob_phase_switch:
                orientation = -1
            else:
                orientation = +1

            hap_block = HaplotypeBlock(
                index=len(self.haplotype_blocks),
                segment_indices=[],
                orientation=orientation,
                alternate_haplotype=Haplotype.A
            )
            self.haplotype_blocks.append(hap_block)

            for bin_idx in range(start_bin, end_bin + 1):
                self.bin_to_block[bin_idx] = hap_block.index
        
        for segment in self.segments:
            block_idx = self.bin_to_block.get(segment.start_bin, 0)
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
            "bin_snp_counts": getattr(self, "bin_snp_counts", None),
            "snp_bins": getattr(self, "snp_bins", None),
            "snp_positions": getattr(self, "snp_positions", None),
            "snp_chromosomes": getattr(self, "snp_chromosomes", None),
            "snp_bin_local_indices": getattr(self, "snp_bin_local_indices", None),
            "events": events,
        }
