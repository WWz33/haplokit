#!/usr/bin/env python3
"""
Algorithm comparison with complete PopART style + labels.
Uses the SAME test data as popart_complete_comparison_v1.0.png
"""

import sys
sys.path.insert(0, 'F:/genehapr-master')

import matplotlib.pyplot as plt
from haplokit.network_python import compute_msn_python, compute_tcs_python, compute_mjn_python
from haplokit.plot_network_popart import PopARTStyleNetwork

# SAME test network as popart_complete_comparison_v1.0
test_network = {
    "nodes": [
        {
            "id": "H0",
            "sequence": "ATCGATCG",
            "samples": ["S1", "S2", "S3", "S4", "S5"],
            "traits": {"Pop1": 3, "Pop2": 2}
        },
        {
            "id": "H1",
            "sequence": "ATCGTTCG",
            "samples": ["S6", "S7"],
            "traits": {"Pop1": 1, "Pop2": 1}
        },
        {
            "id": "H2",
            "sequence": "TTCGATCG",
            "samples": ["S8", "S9", "S10"],
            "traits": {"Pop2": 3}
        },
        {
            "id": "H3",
            "sequence": "ATCGATCC",
            "samples": ["S11"],
            "traits": {"Pop1": 1}
        },
    ],
    "edges": [
        {"source": "H0", "target": "H1", "weight": 1},
        {"source": "H0", "target": "H2", "weight": 2},
        {"source": "H0", "target": "H3", "weight": 1},
    ]
}

print("Algorithm Comparison - PopART Complete Style")
print("=" * 60)
print("Using SAME test data as popart_complete_comparison_v1.0")
print(f"Network: {len(test_network['nodes'])} nodes, {len(test_network['edges'])} edges")
print()

# Extract sequences and samples for algorithm computation
sequences = [node['sequence'] for node in test_network['nodes']]
samples = [node['samples'] for node in test_network['nodes']]

print("[1/3] Computing MSN network...")
msn_network = compute_msn_python(sequences, samples)
print(f"  MSN: {len(msn_network['nodes'])} nodes, {len(msn_network['edges'])} edges")

print("[2/3] Computing TCS network...")
tcs_network = compute_tcs_python(sequences, samples)
print(f"  TCS: {len(tcs_network['nodes'])} nodes, {len(tcs_network['edges'])} edges")

print("[3/3] Computing MJN network...")
mjn_network = compute_mjn_python(sequences, samples)
print(f"  MJN: {len(mjn_network['nodes'])} nodes, {len(mjn_network['edges'])} edges")

# Add traits from original test data
# IMPORTANT: Algorithm functions return nodes with numeric IDs (0, 1, 2, 3)
# but original test data uses string IDs ('H0', 'H1', 'H2', 'H3')
# We need to map by index, not by ID string
traits_list = [node['traits'] for node in test_network['nodes']]

for network in [msn_network, tcs_network, mjn_network]:
    for i, node in enumerate(network['nodes']):
        if i < len(traits_list):
            node['traits'] = traits_list[i]

# Create comparison figure with PopART complete style
print("\n[4/4] Generating PopART-style comparison...")
fig, axes = plt.subplots(1, 3, figsize=(20, 7))
fig.patch.set_facecolor('white')
fig.suptitle('Haplotype Network Algorithms - PopART Complete Style',
             fontsize=16, fontweight='bold', y=0.98)

configs = [
    {
        'network': msn_network,
        'title': 'MSN (Minimum Spanning Network)',
        'edge_mode': 'hatch_marks'
    },
    {
        'network': tcs_network,
        'title': 'TCS (Statistical Parsimony)',
        'edge_mode': 'dots'
    },
    {
        'network': mjn_network,
        'title': 'MJN (Median-Joining Network)',
        'edge_mode': 'numbers'
    }
]

for ax, config in zip(axes, configs):
    renderer = PopARTStyleNetwork(
        config['network'],
        figsize=(7, 7),
        edge_mode=config['edge_mode'],
        show_node_labels=True,
        show_size_legend=True,
        show_trait_legend=True,
        show_border=True,
        color_theme='colorbrewer'
    )

    renderer._build_graph()

    # CRITICAL: Add traits to graph nodes after _build_graph()
    # Algorithm functions use numeric IDs (0, 1, 2, 3), map by index
    for node_id in renderer.G.nodes():
        if isinstance(node_id, int) and node_id < len(traits_list):
            renderer.G.nodes[node_id]['traits'] = traits_list[node_id]

    renderer._compute_layout('spring')

    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_facecolor('white')
    ax.set_title(config['title'], fontsize=12, fontweight='bold', pad=15)

    # Draw all components
    renderer._draw_edges(ax)
    renderer._draw_nodes(ax)
    renderer._draw_node_labels(ax)

    # Draw legends
    x_coords = [renderer.pos[n][0] for n in renderer.G.nodes()]
    y_coords = [renderer.pos[n][1] for n in renderer.G.nodes()]

    legend_x = max(x_coords) + 0.15
    legend_y_size = max(y_coords) - 0.1
    legend_y_trait = legend_y_size - 0.5

    renderer._draw_size_legend(ax, legend_x, legend_y_size)

    # Get trait names
    trait_names = set()
    for node_id in renderer.G.nodes():
        node_data = renderer.G.nodes[node_id]
        traits = node_data.get('traits', {})
        trait_names.update(traits.keys())

    if trait_names:
        renderer._draw_trait_legend(ax, legend_x, legend_y_trait, sorted(trait_names))

    renderer._draw_border(ax)

plt.tight_layout()
output_file = "popart_algorithms_comparison_v1.0.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\n[OK] Saved {output_file}")
print("\nFeatures:")
print("  [OK] 3 algorithms (MSN, TCS, MJN)")
print("  [OK] 3 edge modes (hatch marks, dots, numbers)")
print("  [OK] Node labels (haplotype IDs)")
print("  [OK] Pie chart nodes (wild/landrace/cultivar)")
print("  [OK] Concentric size legend (N=10, N=1)")
print("  [OK] Trait color legend (ColorBrewer Set2)")
print("  [OK] Border rectangles")
print("\n" + "=" * 60)
print("Algorithm comparison complete!")
