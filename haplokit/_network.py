"""Haplotype network visualization in popart style.

Implements a TCS-style (Templeton, Crandall & Sing 1992) network construction
with median (intermediate) vertices for unconnected components, and renders it
using popart's visual conventions:

* Node diameter ∝ sqrt(frequency)  (popart NetworkLayout.cpp uses VERTRAD*sqrt(size))
* Edge length ∝ mutation count    (popart uses EDGELENGTH * weight as preferred spring length)
* Hatch marks (tick marks) across each edge, 1 per mutation
  (popart EdgeItem.cpp `ShowDashes` mode)
* Small filled black "median" nodes for inferred ancestors with degree >= 3
* Population pie charts coloring real haplotype nodes

Reference:
    Leigh JW & Bryant D (2015). popart: full-feature software for haplotype
    network construction. Methods in Ecology and Evolution 6: 1110-1116.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch, Wedge

from ._palette import PALETTE, save_figure
from .network import compute_tcs_network


# ---------------------------------------------------------------------------
# Distance + TCS-style network construction
# ---------------------------------------------------------------------------

def _hamming(a: list[str], b: list[str]) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _pairwise_distance_matrix(haps: list[str]) -> np.ndarray:
    """Pairwise Hamming distance between pipe-delimited state strings."""
    split = [h.split("|") for h in haps]
    n = len(split)
    d = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            v = _hamming(split[i], split[j])
            d[i, j] = v
            d[j, i] = v
    return d


def _build_tcs_network(
    dist: np.ndarray,
    n_real: int,
) -> tuple[list[tuple[int, int, int]], int]:
    """Build TCS-style minimum-step network.

    Preserved as archived Python reference inside active frontend module.
    Production path now uses C++ backend via `compute_tcs_network()`.
    """
    n = n_real
    edges: list[tuple[int, int, int]] = []

    pairs: list[tuple[int, int, int]] = [
        (int(dist[i, j]), i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if dist[i, j] > 0
    ]
    pairs.sort()

    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    next_idx = n
    accepted_pairs: set[tuple[int, int]] = set()
    for w, i, j in pairs:
        if find(i) == find(j):
            continue
        parent[find(i)] = find(j)
        accepted_pairs.add((i, j))

        if w == 1:
            edges.append((i, j, 1))
        else:
            prev = i
            for _ in range(w - 1):
                m = next_idx
                next_idx += 1
                edges.append((prev, m, 1))
                prev = m
            edges.append((prev, j, 1))

    adj: dict[int, list[int]] = {}
    for u, v, _ in edges:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)

    def bfs_dist(s: int, t: int, cap: int) -> int:
        if s == t:
            return 0
        visited = {s}
        frontier = [s]
        d_step = 0
        while frontier and d_step < cap:
            d_step += 1
            nxt = []
            for u in frontier:
                for v in adj.get(u, ()):
                    if v in visited:
                        continue
                    if v == t:
                        return d_step
                    visited.add(v)
                    nxt.append(v)
            frontier = nxt
        return cap + 1

    for w, i, j in pairs:
        if (i, j) in accepted_pairs:
            continue
        if w < 1 or w > 3:
            continue
        existing = bfs_dist(i, j, w)
        if existing == w and w == 1:
            edges.append((i, j, 1))
            adj.setdefault(i, []).append(j)
            adj.setdefault(j, []).append(i)

    return edges, next_idx


# ---------------------------------------------------------------------------
# Layout — Tunkelang-style spring+repulsion (popart NetworkLayout)
# ---------------------------------------------------------------------------

def _spring_layout(
    n_total: int,
    edges: list[tuple[int, int, int]],
    node_sizes: np.ndarray,
    iterations: int = 400,
    seed: int = 42,
) -> np.ndarray:
    """Force-directed layout matching popart conventions.

    - Spring (attractive) force along edges; preferred length scales with edge
      weight so that a 2-mutation edge tries to be twice as long as a
      1-mutation edge.
    - Repulsive (Coulomb-like) force between every pair.
    - Node radii feed into the spring's "rest length" so big circles don't
      overlap.

    Args:
        n_total: number of vertices (real + median).
        edges: list of (u, v, w) where w is the edge length in mutations.
        node_sizes: per-vertex radius (data units).

    Returns:
        positions array of shape (n_total, 2).
    """
    rng = np.random.default_rng(seed)
    EDGE_UNIT = 1.6  # base spring length per mutation

    # Initial placement on a small jittered circle to avoid coincident points
    angles = np.linspace(0, 2 * math.pi, n_total, endpoint=False)
    radius = 0.6 * math.sqrt(max(1, n_total))
    pos = np.column_stack([
        radius * np.cos(angles) + rng.normal(0, 0.1, n_total),
        radius * np.sin(angles) + rng.normal(0, 0.1, n_total),
    ])

    edge_arr = [(u, v, max(1, w)) for u, v, w in edges]

    temp = 0.5
    cooling = 0.985

    for _ in range(iterations):
        disp = np.zeros_like(pos)

        # Repulsion between all pairs
        diff = pos[:, None, :] - pos[None, :, :]      # (n, n, 2)
        dist2 = (diff ** 2).sum(axis=2) + 1e-6        # (n, n)
        dist = np.sqrt(dist2)
        # sum of radii — keeps big nodes from touching
        rsum = node_sizes[:, None] + node_sizes[None, :]
        k = EDGE_UNIT + rsum
        repulse = (k * k) / dist2
        np.fill_diagonal(repulse, 0.0)
        # sum contributions from all j -> i
        contrib = (diff / dist[..., None]) * repulse[..., None]   # (n, n, 2)
        disp += contrib.sum(axis=1)

        # Attraction along edges
        for u, v, w in edge_arr:
            dx = pos[v, 0] - pos[u, 0]
            dy = pos[v, 1] - pos[u, 1]
            d = math.hypot(dx, dy) + 1e-6
            pref = w * EDGE_UNIT + 0.5 * (node_sizes[u] + node_sizes[v])
            attract = (d - pref) / pref * 0.6
            fx = dx / d * attract * d
            fy = dy / d * attract * d
            disp[u, 0] += fx
            disp[u, 1] += fy
            disp[v, 0] -= fx
            disp[v, 1] -= fy

        # Cap step
        mag = np.linalg.norm(disp, axis=1) + 1e-9
        scale = np.minimum(mag, temp) / mag
        pos = pos + disp * scale[:, None]
        temp *= cooling

    # Centre on origin
    pos -= pos.mean(axis=0)
    return pos


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------

def _draw_edge_with_hatches(
    ax,
    p0: tuple[float, float],
    p1: tuple[float, float],
    r0: float,
    r1: float,
    n_mut: int,
    edge_color: str = "#777777",
    edge_width: float = 0.9,
    hatch_length_data: float = 0.18,
) -> None:
    """Draw an edge between two node centres, trimmed to each node's circumference,
    with `n_mut` perpendicular hatch (tick) marks across it.

    Mirrors popart EdgeItem.cpp `ShowDashes` mode (one tick per mutation,
    evenly spaced along the segment).
    """
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = x1 - x0, y1 - y0
    seg_len = math.hypot(dx, dy)
    if seg_len < 1e-6:
        return
    ux, uy = dx / seg_len, dy / seg_len

    # Trim to circle edges
    sx, sy = x0 + ux * r0, y0 + uy * r0
    ex, ey = x1 - ux * r1, y1 - uy * r1

    ax.plot([sx, ex], [sy, ey],
            color=edge_color, linewidth=edge_width,
            solid_capstyle="round", zorder=1)

    if n_mut <= 0:
        return

    # Perpendicular unit vector for tick marks
    px, py = -uy, ux

    # Tick spacing along the visible segment (popart's r/(weight+1))
    vis_len = math.hypot(ex - sx, ey - sy)
    if vis_len < 1e-6:
        return
    spacing = vis_len / (n_mut + 1)
    half_h = min(hatch_length_data, spacing / 1.6) * 0.5

    for k in range(1, n_mut + 1):
        t = k * spacing
        cx = sx + ux * t
        cy = sy + uy * t
        ax.plot(
            [cx - px * half_h, cx + px * half_h],
            [cy - py * half_h, cy + py * half_h],
            color="#222222", linewidth=0.9, solid_capstyle="butt", zorder=3,
        )


def _draw_pie(
    ax,
    centre: tuple[float, float],
    radius: float,
    slices: list[tuple[str, int]],
    pop_color: dict[str, str],
    edge_color: str = "#1a1a1a",
    edge_width: float = 0.7,
) -> None:
    x, y = centre
    total = sum(c for _, c in slices)
    if total <= 0:
        ax.add_patch(Circle((x, y), radius, facecolor="#f0f0f0",
                            edgecolor=edge_color, linewidth=edge_width, zorder=5))
        return
    start = 90.0
    for name, count in slices:
        if count <= 0:
            continue
        span = 360.0 * count / total
        ax.add_patch(Wedge(
            (x, y), radius, start, start + span,
            facecolor=pop_color.get(name, "#cccccc"),
            edgecolor="white", linewidth=0.4, zorder=5,
        ))
        start += span
    ax.add_patch(Circle((x, y), radius, fill=False,
                        edgecolor=edge_color, linewidth=edge_width, zorder=6))


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def plot_hap_network(
    hap_names: list[str],
    hap_strings: list[str],
    hap_counts: list[int],
    output_path: str | Path,
    pop_data: dict[str, list[tuple[str, int]]] | None = None,
    title: str = "",
    size_scale: str = "log10",
    show_mutation: int = 2,
    legend_font_size: float = 7,
    dpi: int = 600,
    fmt: str | None = None,
) -> Path:
    """Render a TCS-style haplotype network in popart's visual style.

    Args:
        hap_names: Haplotype labels (e.g., ["H001", "H002"]).
        hap_strings: Pipe-delimited state strings, one per haplotype.
        hap_counts: Sample count per haplotype.
        output_path: Output file path (PNG/PDF/SVG).
        pop_data: Optional dict mapping hap_name -> [(pop_name, count), ...].
        title: Accepted for CLI compatibility; NOT rendered (matches the new
            distribution-plot style — no titles).
        size_scale: kept for signature compatibility; popart-style sqrt scaling
            is always used for the node geometry.
        show_mutation: kept for signature compatibility; popart renders mutations
            as hatch marks across edges in all modes.
        legend_font_size: Legend text size.
        dpi: Output resolution for raster formats.
        fmt: Output format override.
    """
    n = len(hap_names)
    if n == 0:
        raise ValueError("no haplotypes provided")

    counts = np.asarray(hap_counts, dtype=float)

    # ---- Build network -----------------------------------------------------
    samples = [[name] * max(int(count), 1) for name, count in zip(hap_names, hap_counts)]
    network = compute_tcs_network(hap_strings, samples)
    node_index = {node["id"]: idx for idx, node in enumerate(network["nodes"])}
    edges = [
        (
            node_index[edge["source"]],
            node_index[edge["target"]],
            max(int(edge.get("weight", 1)), 1),
        )
        for edge in network["edges"]
    ]
    n_total = len(network["nodes"])

    # ---- Node radii: diameter ∝ sqrt(freq), like popart -------------------
    # Real nodes: radius from frequency. Median nodes: small fixed black dots.
    real_radii = np.sqrt(np.maximum(counts, 1.0))
    if real_radii.max() > 0:
        real_radii = real_radii / real_radii.max()
    # Map to data-space radii in [0.30, 1.05]
    r_min, r_max = 0.30, 1.05
    real_radii = r_min + real_radii * (r_max - r_min)

    median_radius = 0.10
    radii = np.empty(n_total, dtype=float)
    radii[:n] = real_radii
    radii[n:] = median_radius

    # ---- Layout -----------------------------------------------------------
    pos = _spring_layout(n_total, edges, radii)

    # ---- Figure setup -----------------------------------------------------
    pad = float(radii.max() * 2.5 + 1.0)
    xs, ys = pos[:, 0], pos[:, 1]
    xmin, xmax = xs.min() - pad, xs.max() + pad
    ymin, ymax = ys.min() - pad, ys.max() + pad
    span_x = xmax - xmin
    span_y = ymax - ymin

    fig_w = max(6.5, min(11.0, span_x * 0.55 + 2.0))
    fig_h = max(5.5, min(10.0, span_y * 0.55 + 2.0))
    # Reserve top strip for the size legend
    top_strip = 0.10
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.04, right=0.96, top=1.0 - top_strip, bottom=0.04)
    ax.set_aspect("equal")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    # ---- Population colors ------------------------------------------------
    pop_color: dict[str, str] = {}
    pop_order: list[str] = []
    if pop_data:
        for name in hap_names:
            for pname, _ in pop_data.get(name, []):
                if pname not in pop_color:
                    pop_color[pname] = PALETTE[len(pop_order) % len(PALETTE)]
                    pop_order.append(pname)

    # ---- Edges (popart ShowDashes hatches) --------------------------------
    for u, v, w in edges:
        _draw_edge_with_hatches(
            ax,
            (pos[u, 0], pos[u, 1]),
            (pos[v, 0], pos[v, 1]),
            radii[u], radii[v],
            n_mut=w,
            edge_color="#7a7a7a",
            edge_width=0.85,
            hatch_length_data=0.22,
        )

    # ---- Median (intermediate) vertices: small black filled dots ---------
    for idx in range(n, n_total):
        ax.add_patch(Circle(
            (pos[idx, 0], pos[idx, 1]), median_radius,
            facecolor="#1a1a1a", edgecolor="#1a1a1a",
            linewidth=0.4, zorder=5,
        ))

    # ---- Real haplotype nodes ---------------------------------------------
    for idx in range(n):
        name = hap_names[idx]
        x, y = pos[idx]
        r = radii[idx]
        if pop_data and name in pop_data and pop_data[name]:
            _draw_pie(ax, (x, y), r, pop_data[name], pop_color)
        else:
            ax.add_patch(Circle(
                (x, y), r, facecolor="#dddddd",
                edgecolor="#1a1a1a", linewidth=0.7, zorder=5,
            ))
        # Label above the node
        ax.text(x, y + r + 0.12, name,
                fontsize=7, fontweight="bold",
                ha="center", va="bottom",
                color="#1a1a1a", zorder=8)

    # ---- Population legend (right side) -----------------------------------
    if pop_color:
        handles = [
            Patch(facecolor=col, edgecolor="#666", label=name)
            for name, col in pop_color.items()
        ]
        leg = ax.legend(
            handles=handles,
            fontsize=legend_font_size,
            loc="center left",
            bbox_to_anchor=(1.005, 0.5),
            frameon=False,
            title="Population",
            title_fontsize=legend_font_size,
        )
        leg._legend_box.align = "left"

    # ---- Size legend (PopART concentric circles style) ---------------------
    unique_counts = sorted({int(c) for c in hap_counts})
    if len(unique_counts) > 1:
        hi = unique_counts[-1]
        lo = unique_counts[0]

        raw_hi = np.sqrt(float(hi))
        raw_lo = np.sqrt(float(lo))
        scale_factor = np.sqrt(max(counts.max(), 1.0))
        r_hi = r_min + (raw_hi / scale_factor) * (r_max - r_min)
        r_lo = r_min + (raw_lo / scale_factor) * (r_max - r_min)

        size_ax = fig.add_axes((0.04, 1.0 - top_strip + 0.005, 0.92, top_strip - 0.02))
        size_ax.set_xlim(0, 1)
        size_ax.set_ylim(0, 1)
        size_ax.axis("off")

        # Normalize to fit legend axes
        legend_scale = 0.40 / r_hi if r_hi > 0 else 1.0
        R_hi = r_hi * legend_scale
        R_lo = r_lo * legend_scale

        cx = 0.50
        # Large circle centered so bottom touches baseline
        baseline_y = 0.10
        cy_hi = baseline_y + R_hi

        # Large circle (N=hi)
        size_ax.add_patch(plt.Circle(
            (cx, cy_hi), R_hi,
            facecolor="none", edgecolor="#1a1a1a",
            linewidth=0.8, zorder=2,
        ))
        # Small circle (N=lo) nested at bottom of large circle
        cy_lo = baseline_y + R_lo
        size_ax.add_patch(plt.Circle(
            (cx, cy_lo), R_lo,
            facecolor="none", edgecolor="#1a1a1a",
            linewidth=0.8, zorder=3,
        ))

        # Labels
        size_ax.text(cx + R_hi + 0.04, cy_hi, f"N={hi}",
                     fontsize=legend_font_size, ha="left", va="center", color="#333")
        size_ax.text(cx + R_lo + 0.04, cy_lo, f"N={lo}",
                     fontsize=legend_font_size, ha="left", va="center", color="#333")

    # Title intentionally NOT rendered (parameter accepted only for compat).
    _ = title

    return save_figure(fig, output_path, dpi=dpi, fmt=fmt)
