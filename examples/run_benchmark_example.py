"""
Example: Running the benchmark pipeline

Demonstrates Week 14 deliverable
"""

from haplotreesim import SimulationConfig
from haplotreesim.benchmark import BenchmarkConfig, run_benchmark

print("="*70)
print("BENCHMARK PIPELINE EXAMPLE")
print("="*70)

# Configure simulation
sim_config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100,
    lambda_events=2.0,
    prob_wgd=0.3,
    random_seed=42
)

# Configure benchmark
benchmark_config = BenchmarkConfig(
    sim_config=sim_config,
    output_dir="benchmark_results",
    methods=['dummy'],  # Start with oracle baseline
    breakpoint_tolerance=2,
    event_tolerance=1
)

# Run benchmark
results = run_benchmark(benchmark_config)

print("\n" + "="*70)
print("BENCHMARK COMPLETE!")
print("="*70)
print(f"\nResults saved to: benchmark_results/")
print(f"  - benchmark_results.json (full results)")
print(f"  - ground_truth.json (ground truth data)")
