"""
Benchmark Pipeline (Week 14)
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
import time

from .simulator import HaploTreeSimulator
from .data_models import SimulationConfig
from .metrics import (
    compute_all_hscn_metrics,
    compute_breakpoint_metrics,
    compute_clone_assignment_metrics,
    compute_all_tree_metrics
)


def _json_serialize(obj):
    """Recursively convert numpy/custom types to JSON-serializable Python types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, dict):
        return {str(k): _json_serialize(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_json_serialize(item) for item in obj]
    elif hasattr(obj, '__dataclass_fields__'):
        return _json_serialize(asdict(obj))
    elif hasattr(obj, 'start_bin'):
        return {
            'start_bin': int(obj.start_bin),
            'end_bin': int(obj.end_bin),
            'haplotype': str(obj.haplotype),
            'amplitude': int(obj.amplitude)
        }
    else:
        return obj


@dataclass
class BenchmarkConfig:
    """Configuration for benchmark run"""
    sim_config: SimulationConfig
    output_dir: str = "benchmark_results"
    methods: List[str] = None
    breakpoint_tolerance: int = 2
    event_tolerance: int = 1

    def __post_init__(self):
        if self.methods is None:
            self.methods = ['dummy']


class BaseMethod:
    """Base class for inference methods"""

    def __init__(self, name: str):
        self.name = name

    def run(self, read_counts, allele_counts, ground_truth):
        raise NotImplementedError


class DummyMethod(BaseMethod):
    """Dummy baseline that returns ground truth (oracle)"""

    def __init__(self):
        super().__init__("dummy_oracle")

    def run(self, read_counts, allele_counts, ground_truth):
        return {
            'cn_A': ground_truth['cn_A'],
            'cn_B': ground_truth['cn_B'],
            'clone_labels': ground_truth['clone_labels'],
            'breakpoints': ground_truth['breakpoints'],
            'tree_edges': [],
            'events': []
        }


class BenchmarkRunner:
    """Main benchmark orchestrator"""

    def __init__(self, config: BenchmarkConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.methods = self._initialize_methods()

    def _initialize_methods(self):
        methods = {}
        for method_name in self.config.methods:
            if method_name == 'dummy':
                methods['dummy'] = DummyMethod()
            else:
                raise NotImplementedError(f"Method {method_name} not yet implemented")
        return methods

    def run(self):
        results = {
            'config': asdict(self.config.sim_config),
            'methods': {},
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }

        print("="*70)
        print("HAPLOTREESIM BENCHMARK PIPELINE")
        print("="*70)

        print("\n1. Generating simulation...")
        sim_data, ground_truth = self._generate_simulation()
        print(f"   Generated {sim_data['read_counts'].shape[0]} cells, "
              f"{len(ground_truth['segments'])} segments")

        self._save_ground_truth(ground_truth)

        for method_name, method in self.methods.items():
            print(f"\n2. Running method: {method_name}")

            method_start = time.time()
            predictions = method.run(
                sim_data['read_counts'],
                sim_data['allele_counts'],
                ground_truth
            )
            method_time = time.time() - method_start

            print(f"   Completed in {method_time:.2f}s")

            print(f"3. Computing metrics for {method_name}...")
            metrics = self._compute_metrics(predictions, ground_truth)

            results['methods'][method_name] = {
                'predictions': self._serialize_predictions(predictions),
                'metrics': _json_serialize(metrics),
                'runtime_seconds': method_time
            }

            self._print_method_summary(method_name, metrics)

        results_file = self.output_dir / "benchmark_results.json"
        with open(results_file, 'w') as f:
            json.dump(_json_serialize(results), f, indent=2)

        print(f"\n✓ Results saved to: {results_file}")
        print("="*70)

        return results

    def _generate_simulation(self):
        sim = HaploTreeSimulator(self.config.sim_config)
        read_counts, allele_counts = sim.run()

        ground_truth = sim.get_ground_truth()

        n_cells = len(ground_truth['clone_assignments'])
        n_segments = len(ground_truth['segments'])

        cn_A = np.zeros((n_cells, n_segments), dtype=int)
        cn_B = np.zeros((n_cells, n_segments), dtype=int)

        for cell_idx, clone_idx in enumerate(ground_truth['clone_assignments']):
            clone = sim.clones[clone_idx]
            for seg_idx, segment in enumerate(ground_truth['segments']):
                start = segment.start_bin
                end = segment.end_bin
                cn_A[cell_idx, seg_idx] = int(np.mean(clone.cn_profile_A[start:end+1]))
                cn_B[cell_idx, seg_idx] = int(np.mean(clone.cn_profile_B[start:end+1]))

        ground_truth['cn_A'] = cn_A
        ground_truth['cn_B'] = cn_B
        ground_truth['clone_labels'] = np.array(ground_truth['clone_assignments'])

        breakpoints = [seg.start_bin for seg in ground_truth['segments'][1:]]
        ground_truth['breakpoints'] = np.array(breakpoints)

        return {
            'read_counts': read_counts,
            'allele_counts': allele_counts
        }, ground_truth

    def _compute_metrics(self, predictions, ground_truth):
        metrics = {}

        hscn = compute_all_hscn_metrics(
            ground_truth['cn_A'],
            ground_truth['cn_B'],
            predictions['cn_A'],
            predictions['cn_B'],
            ground_truth['clone_labels']
        )
        metrics['hscn'] = hscn

        if 'breakpoints' in predictions and len(predictions['breakpoints']) > 0:
            bp = compute_breakpoint_metrics(
                ground_truth['breakpoints'],
                predictions['breakpoints'],
                tolerance=self.config.breakpoint_tolerance
            )
            metrics['breakpoints'] = bp

        if 'clone_labels' in predictions:
            clone = compute_clone_assignment_metrics(
                ground_truth['clone_labels'],
                predictions['clone_labels']
            )
            metrics['clone_assignment'] = clone

        if 'tree_edges' in predictions and len(predictions['tree_edges']) > 0:
            tree = compute_all_tree_metrics(
                ground_truth.get('tree_edges', []),
                predictions['tree_edges'],
                ground_truth.get('events', []),
                predictions.get('events', [])
            )
            metrics['tree'] = tree

        return metrics

    def _serialize_predictions(self, predictions):
        return _json_serialize(predictions)

    def _save_ground_truth(self, ground_truth):
        gt_file = self.output_dir / "ground_truth.json"

        gt_serialized = {}
        for key, value in ground_truth.items():
            if isinstance(value, np.ndarray):
                gt_serialized[key] = value.tolist()
            elif key == 'segments':
                gt_serialized[key] = [
                    {
                        'start_bin': seg.start_bin,
                        'end_bin': seg.end_bin,
                        'index': seg.index
                    }
                    for seg in value
                ]
            elif key == 'events' and isinstance(value, dict):
                gt_serialized[key] = {
                    str(k): [
                        {
                            'start_bin': int(e.start_bin),
                            'end_bin': int(e.end_bin),
                            'haplotype': str(e.haplotype),
                            'amplitude': int(e.amplitude)
                        }
                        for e in events
                    ]
                    for k, events in value.items()
                }
            elif key == 'tree':
                continue
            else:
                gt_serialized[key] = value

        with open(gt_file, 'w') as f:
            json.dump(_json_serialize(gt_serialized), f, indent=2)

    def _print_method_summary(self, method_name, metrics):
        print(f"\n   Results for {method_name}:")

        if 'hscn' in metrics:
            hscn = metrics['hscn']
            print(f"     HSCN Error: {hscn['hscn_error']:.4f}")
            print(f"     LOH F1: {hscn['loh_f1']:.4f}")
            print(f"     MSR: {hscn['msr']:.4f}")

        if 'breakpoints' in metrics:
            bp = metrics['breakpoints']
            print(f"     Breakpoint F1: {bp['f1']:.4f}")

        if 'clone_assignment' in metrics:
            clone = metrics['clone_assignment']
            print(f"     ARI: {clone['ari']:.4f}")
            print(f"     NMI: {clone['nmi']:.4f}")


def run_benchmark(config: BenchmarkConfig):
    runner = BenchmarkRunner(config)
    return runner.run()


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description="HaploTreeSim Benchmark Pipeline")
    parser.add_argument("--config", type=str, required=True, help="Path to JSON config file")
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        cfg = json.load(f)

    sim_cfg = SimulationConfig(**cfg.get("sim_config", {}))
    bench_cfg = BenchmarkConfig(
        sim_config=sim_cfg,
        output_dir=cfg.get("output_dir", "benchmark_results"),
        methods=cfg.get("methods", ["dummy"]),
        breakpoint_tolerance=cfg.get("breakpoint_tolerance", 2),
        event_tolerance=cfg.get("event_tolerance", 1),
    )

    run_benchmark(bench_cfg)
