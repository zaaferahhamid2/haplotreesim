"""
Example configuration for HaploTreeSim.

This is a minimal configuration that produces diploid output (Week 5 deliverable).
"""

# Minimal diploid configuration - no CNAs yet
minimal_diploid = {
    # Genome
    "num_bins": 1000,
    "bin_length": 50000,  # 50kb bins
    
    # Clones (all diploid for now)
    "num_clones": 3,
    "max_copy_number": 8,
    "root_type": "diploid",
    
    # Events (disabled for Week 5)
    "lambda_events": 0.0,  # No events yet
    "lambda_amplitude": 0.0,
    "prob_wgd": 0.0,
    "prob_mirror": 0.0,
    
    # Cells
    "num_cells": 100,
    "prob_normal": 0.0,
    "prob_doublet": 0.0,
    
    # Read-depth model
    "mean_library_size": 50.0,
    "library_size_cv": 0.3,
    "theta_x": 10.0,
    
    # Allelic model
    "snp_density": 0.001,  # 1 SNP per kb
    "mean_allelic_coverage": 30.0,
    "allelic_coverage_cv": 0.3,
    "nu_a": 20.0,
    "prob_phase_switch": 0.0,
    
    # Random seed
    "random_seed": 42,
}


# Low-pass shallow coverage configuration (0.01x)
low_pass_config = {
    # Genome
    "num_bins": 5000,
    "bin_length": 50000,
    
    # Clones
    "num_clones": 5,
    "max_copy_number": 8,
    "root_type": "diploid",
    
    # Events
    "lambda_events": 0.0,  # Will be enabled in future weeks
    "lambda_amplitude": 0.0,
    "prob_wgd": 0.0,
    "prob_mirror": 0.0,
    
    # Cells
    "num_cells": 200,
    "prob_normal": 0.0,
    "prob_doublet": 0.0,
    
    # Read-depth model (adjusted for 0.01x coverage)
    "mean_library_size": 10.0,  # Low coverage
    "library_size_cv": 0.4,
    "theta_x": 5.0,
    
    # Allelic model (sparse due to low coverage)
    "snp_density": 0.001,
    "mean_allelic_coverage": 5.0,  # Very sparse
    "allelic_coverage_cv": 0.4,
    "nu_a": 15.0,
    "prob_phase_switch": 0.01,
    
    # Random seed
    "random_seed": 123,
}


# WGD root configuration
wgd_root_config = {
    # Genome
    "num_bins": 1000,
    "bin_length": 50000,
    
    # Clones (WGD root)
    "num_clones": 3,
    "max_copy_number": 8,
    "root_type": "wgd",  # All clones start at (2,2)
    
    # Events
    "lambda_events": 0.0,
    "lambda_amplitude": 0.0,
    "prob_wgd": 0.0,
    "prob_mirror": 0.0,
    
    # Cells
    "num_cells": 100,
    "prob_normal": 0.0,
    "prob_doublet": 0.0,
    
    # Read-depth model
    "mean_library_size": 50.0,
    "library_size_cv": 0.3,
    "theta_x": 10.0,
    
    # Allelic model
    "snp_density": 0.001,
    "mean_allelic_coverage": 30.0,
    "allelic_coverage_cv": 0.3,
    "nu_a": 20.0,
    "prob_phase_switch": 0.0,
    
    # Random seed
    "random_seed": 456,
}
