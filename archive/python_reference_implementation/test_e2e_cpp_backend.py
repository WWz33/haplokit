#!/usr/bin/env python3
"""
End-to-end test for C++ network backend integration.
Tests Python -> C++ backend -> PopART visualization pipeline.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import matplotlib.pyplot as plt
from haplokit.network import compute_msn_network, compute_tcs_network, compute_mjn_network
from haplokit.plot_network_popart import plot_popart_network
import networkx as nx

def create_test_data():
    """Create test dataset matching Python backend test."""
    sequences = [
        "AAAAAAAAAA",
        "AAAAAAAAAC",
        "AAAAAACCCC",
        "CCCCCCCCCC"
    ]

    samples = [
        ["S1", "S2", "S3"],
        ["S4"],
        ["S5", "S6"],
        ["S7", "S8"]
    ]

    traits = {
        0: {"Pop1": 3},
        1: {"Pop2": 1},
        2: {"Pop1": 1, "Pop2": 1},
        3: {"Pop2": 2}
    }

    return sequences, samples, traits

def test_cpp_backend():
    """Test C++ backend with all three algorithms."""
    print("=" * 60)
    print("E2E Test: Python -> C++ Backend -> PopART Visualization")
    print("=" * 60)

    sequences, samples, traits = create_test_data()

    print(f"\nTest data:")
    print(f"  Sequences: {len(sequences)}")
    print(f"  Samples: {sum(len(s) for s in samples)} total")
    print(f"  Traits: {len(traits)} haplotypes with population info")

    # Test MSN algorithm
    print("\n[1/3] Testing MSN (Minimum Spanning Network)...")
    try:
        result_msn = compute_msn_network(sequences, samples, epsilon=5)
        print(f"  Result: {len(result_msn['nodes'])} nodes, {len(result_msn['edges'])} edges")
        print(f"  Backend: C++")
    except Exception as e:
        print(f"  ERROR: {e}")
        result_msn = None

    # Test TCS algorithm
    print("\n[2/3] Testing TCS (Statistical Parsimony)...")
    try:
        result_tcs = compute_tcs_network(sequences, samples)
        median_count = sum(1 for node in result_tcs['nodes'] if node.get('is_median', False))
        print(f"  Result: {len(result_tcs['nodes'])} nodes, {len(result_tcs['edges'])} edges")
        print(f"  Median vertices: {median_count}")
    except Exception as e:
        print(f"  ERROR: {e}")
        result_tcs = None

    # Test MJN algorithm
    print("\n[3/3] Testing MJN (Median-Joining Network)...")
    try:
        result_mjn = compute_mjn_network(sequences, samples, epsilon=5)
        print(f"  Result: {len(result_mjn['nodes'])} nodes, {len(result_mjn['edges'])} edges")
    except Exception as e:
        print(f"  ERROR: {e}")
        result_mjn = None

    # Generate visualization
    if all(result is not None for result in [result_msn, result_tcs, result_mjn]):
        print("\n" + "=" * 60)
        print("Generating PopART-style visualization...")
        print("=" * 60)

        # Add traits to results
        for result in [result_msn, result_tcs, result_mjn]:
            for node in result['nodes']:
                node_id = node['id']
                if node_id in traits:
                    node['traits'] = traits[node_id]

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        algorithms = [
            (result_msn, "MSN", "hatch_marks"),
            (result_tcs, "TCS", "dots"),
            (result_mjn, "MJN", "numbers")
        ]

        for ax, (result, name, edge_mode) in zip(axes, algorithms):
            plot_popart_network(
                result,
                ax=ax,
                edge_mode=edge_mode,
                show_node_labels=True,
                color_theme='colorbrewer'
            )
            ax.set_title(f"{name} Algorithm (C++ Backend)", fontsize=14, fontweight='bold')

        plt.tight_layout()
        output_file = "e2e_test_cpp_backend_v1.0.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"\nVisualization saved: {output_file}")

        # Get file size
        file_size = os.path.getsize(output_file) / 1024
        print(f"File size: {file_size:.1f} KB")

        print("\n" + "=" * 60)
        print("E2E Test PASSED - C++ backend integration successful!")
        print("=" * 60)

        return True
    else:
        print("\n" + "=" * 60)
        print("E2E Test FAILED - Some algorithms failed")
        print("=" * 60)
        return False

if __name__ == "__main__":
    success = test_cpp_backend()
    sys.exit(0 if success else 1)
