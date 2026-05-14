"""
Complete PopART-style haplotype network visualization.

Reproduces all PopART visual components:
- Pie chart nodes for multiple traits
- Hatch marks on edges (not just numbers)
- Size legend with reference circles
- Trait color legend
- Border rectangle
- Precise styling matching PopART output
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from typing import Dict, Any, List, Tuple, Optional
import networkx as nx

from ._palette import PALETTE


class PopARTStyleNetwork:
    """PopART-style network renderer with all visual components."""

    GREYSCALE_COLORS = [
        '#000000', '#404040', '#808080', '#C0C0C0', '#E0E0E0',
    ]

    VIBRANT_COLORS = [
        '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3',
        '#FF7F00', '#FFFF33', '#A65628', '#F781BF',
    ]

    def __init__(
        self,
        network: Dict[str, Any],
        traits: Optional[Dict[str, Dict[str, int]]] = None,
        figsize: Tuple[float, float] = (12, 10),
        node_scale: float = 1.0,
        edge_color: str = '#000000',
        background_color: str = '#FFFFFF',
        edge_mode: str = 'hatch_marks',
        show_node_labels: bool = False,
        show_size_legend: bool = True,
        show_trait_legend: bool = True,
        show_border: bool = True,
        color_theme: str = 'colorbrewer',
        ax: Optional[plt.Axes] = None,
        title: Optional[str] = None
    ):
        self.network = network
        self.traits = traits
        self.figsize = figsize
        self.node_scale = node_scale
        self.edge_color = edge_color
        self.background_color = background_color
        self.edge_mode = edge_mode
        self.show_node_labels = show_node_labels
        self.show_size_legend = show_size_legend
        self.show_trait_legend = show_trait_legend
        self.show_border = show_border
        self.ax = ax
        self.title = title

        if color_theme == 'greyscale':
            self.colors = self.GREYSCALE_COLORS
        elif color_theme == 'vibrant':
            self.colors = self.VIBRANT_COLORS
        else:
            self.colors = PALETTE

        self.G = None
        self.pos = None

    def _build_graph(self):
        self.G = nx.Graph()
        for node in self.network["nodes"]:
            self.G.add_node(
                node["id"],
                sequence=node["sequence"],
                samples=node.get("samples", []),
                is_median=node.get("is_median", False),
                traits=node.get("traits", {})
            )
        for edge in self.network["edges"]:
            self.G.add_edge(edge["source"], edge["target"], weight=edge["weight"])

    def _compute_layout(self, layout: str = "spring"):
        if layout == "spring":
            self.pos = nx.spring_layout(self.G, k=2, iterations=50, seed=42)
        elif layout == "circular":
            self.pos = nx.circular_layout(self.G)
        elif layout == "kamada_kawai":
            self.pos = nx.kamada_kawai_layout(self.G)
        else:
            self.pos = nx.spring_layout(self.G, seed=42)

    def _draw_edges(self, ax):
        for u, v, data in self.G.edges(data=True):
            x1, y1 = self.pos[u]
            x2, y2 = self.pos[v]
            weight = data['weight']
            ax.plot([x1, x2], [y1, y2], color=self.edge_color,
                    linewidth=1.5, zorder=1, solid_capstyle='round')
            if self.edge_mode == 'hatch_marks':
                self._draw_hatch_marks(ax, x1, y1, x2, y2, weight)
            elif self.edge_mode == 'dots':
                self._draw_dots(ax, x1, y1, x2, y2, weight)
            elif self.edge_mode == 'numbers':
                self._draw_numbers(ax, x1, y1, x2, y2, weight)

    def _draw_hatch_marks(self, ax, x1, y1, x2, y2, weight):
        if weight == 0:
            return
        dx, dy = x2 - x1, y2 - y1
        length = np.sqrt(dx**2 + dy**2)
        if length <= 0:
            return
        perp_x, perp_y = -dy / length, dx / length
        hatch_length = 0.015
        for i in range(weight):
            t = (i + 1) / (weight + 1)
            mid_x = x1 + t * dx
            mid_y = y1 + t * dy
            ax.plot(
                [mid_x - hatch_length * perp_x, mid_x + hatch_length * perp_x],
                [mid_y - hatch_length * perp_y, mid_y + hatch_length * perp_y],
                color=self.edge_color, linewidth=2, zorder=2, solid_capstyle='butt')

    def _draw_dots(self, ax, x1, y1, x2, y2, weight):
        if weight <= 1:
            return
        dx, dy = x2 - x1, y2 - y1
        for i in range(1, weight):
            t = i / weight
            ax.add_patch(mpatches.Circle(
                (x1 + t * dx, y1 + t * dy), 0.01,
                facecolor='black', edgecolor='black', linewidth=0, zorder=2))

    def _draw_numbers(self, ax, x1, y1, x2, y2, weight):
        if weight == 0:
            return
        mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mid_x, mid_y, f'({weight})', fontsize=9, ha='center', va='center',
                zorder=2, bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                    edgecolor='none', alpha=0.9))

    def _draw_pie_chart_node(self, ax, node_id, x, y, radius, traits):
        if not traits or len(traits) == 0:
            ax.add_patch(mpatches.Circle((x, y), radius, facecolor=self.colors[0],
                                         edgecolor='black', linewidth=1.5, zorder=3))
            return
        total = sum(traits.values())
        start_angle = 0
        for i, (trait_name, count) in enumerate(traits.items()):
            angle = 360 * count / total
            ax.add_patch(mpatches.Wedge(
                (x, y), radius, start_angle, start_angle + angle,
                facecolor=self.colors[i % len(self.colors)],
                edgecolor='black', linewidth=1.5, zorder=3))
            start_angle += angle

    def _draw_nodes(self, ax):
        for node_id in self.G.nodes():
            node_data = self.G.nodes[node_id]
            x, y = self.pos[node_id]
            n_samples = len(node_data['samples'])
            if node_data['is_median']:
                ax.add_patch(mpatches.Circle((x, y), 0.015, facecolor='black',
                                             edgecolor='black', linewidth=1, zorder=3))
            else:
                base_radius = 0.03 * self.node_scale
                radius = base_radius * np.sqrt(max(1, n_samples))
                traits = node_data.get('traits', {})
                if traits and len(traits) > 1:
                    self._draw_pie_chart_node(ax, node_id, x, y, radius, traits)
                else:
                    ax.add_patch(mpatches.Circle(
                        (x, y), radius, facecolor=self.colors[0],
                        edgecolor='black', linewidth=1.5, zorder=3))

    def _draw_size_legend(self, ax, x_pos, y_pos):
        if not self.show_size_legend:
            return
        base_radius = 0.03 * self.node_scale
        radius_10 = base_radius * np.sqrt(10)
        radius_1 = base_radius * np.sqrt(1)
        ax.add_patch(mpatches.Circle((x_pos, y_pos), radius_10,
                     facecolor='none', edgecolor='black', linewidth=1.5, zorder=3))
        ax.text(x_pos, y_pos, 'N=10', fontsize=8, ha='center', va='center', zorder=4)
        small_y = y_pos - (radius_10 - radius_1)
        ax.add_patch(mpatches.Circle((x_pos, small_y), radius_1,
                     facecolor='none', edgecolor='black', linewidth=1.5, zorder=3))
        ax.text(x_pos, small_y - radius_1 - 0.01, 'N=1', fontsize=8,
                ha='center', va='top', zorder=4)

    def _draw_trait_legend(self, ax, x_pos, y_pos, trait_names):
        if not self.show_trait_legend or not trait_names:
            return
        legend_y = y_pos
        for i, trait_name in enumerate(trait_names):
            ax.add_patch(mpatches.Circle((x_pos, legend_y), 0.02,
                         facecolor=self.colors[i % len(self.colors)],
                         edgecolor='black', linewidth=1, zorder=3))
            ax.text(x_pos + 0.05, legend_y, trait_name, fontsize=9, va='center', zorder=3)
            legend_y -= 0.06

    def _draw_node_labels(self, ax):
        for node_id in self.G.nodes():
            x, y = self.pos[node_id]
            node_data = self.G.nodes[node_id]
            if node_data.get('is_median', False):
                continue
            ax.text(x, y - 0.08, str(node_id), fontsize=9,
                    ha='center', va='top', zorder=4)

    def _draw_border(self, ax):
        if not self.show_border:
            return
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.add_patch(mpatches.Rectangle(
            (xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0],
            fill=False, edgecolor='black', linewidth=2, zorder=0))

    def plot(self, output_file: Optional[str] = None, layout: str = "spring") -> plt.Figure:
        self._build_graph()
        self._compute_layout(layout)
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_facecolor(self.background_color)
        fig.patch.set_facecolor(self.background_color)
        self._draw_edges(ax)
        self._draw_nodes(ax)
        if self.show_node_labels:
            self._draw_node_labels(ax)
        x_coords = [self.pos[n][0] for n in self.G.nodes()]
        y_coords = [self.pos[n][1] for n in self.G.nodes()]
        legend_x = max(x_coords) + 0.15
        legend_y_size = max(y_coords) - 0.1
        legend_y_trait = legend_y_size - 0.5
        self._draw_size_legend(ax, legend_x, legend_y_size)
        trait_names = set()
        for node_id in self.G.nodes():
            traits = self.G.nodes[node_id].get('traits', {})
            trait_names.update(traits.keys())
        if trait_names:
            self._draw_trait_legend(ax, legend_x, legend_y_trait, sorted(trait_names))
        self._draw_border(ax)
        margin = 0.15
        ax.set_xlim(min(x_coords) - margin, max(x_coords) + margin + 0.3)
        ax.set_ylim(min(y_coords) - margin, max(y_coords) + margin)
        plt.tight_layout()
        if output_file:
            fig.savefig(output_file, dpi=300, bbox_inches='tight',
                        facecolor=self.background_color)
        return fig


def plot_popart_network(
    network: Dict[str, Any],
    output_file: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    renderer = PopARTStyleNetwork(network, **kwargs)
    return renderer.plot(output_file=output_file)
