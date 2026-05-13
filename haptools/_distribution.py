from __future__ import annotations

import json
import math
from pathlib import Path
from urllib.request import urlopen

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Wedge, Patch
from matplotlib.collections import PatchCollection

from ._palette import PALETTE, save_figure

# GeoJSON sources keyed by database name
_GEO_SOURCES: dict[str, str] = {
    "china": "https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json",
}


def _load_geojson(source: str) -> list[dict]:
    """Load GeoJSON features from URL or local file."""
    if source.startswith(("http://", "https://", "ftp://")):
        data = json.loads(urlopen(source, timeout=20).read())
    else:
        data = json.loads(Path(source).read_text(encoding="utf-8"))
    return data.get("features", [])


def _geom_to_polygons(geom: dict) -> list[list[list[float]]]:
    """Extract coordinate rings from a GeoJSON geometry as flat polygon paths."""
    rings: list[list[list[float]]] = []
    gtype = geom["type"]

    if gtype == "Polygon":
        rings.extend(geom["coordinates"])
    elif gtype == "MultiPolygon":
        for poly in geom["coordinates"]:
            rings.extend(poly)
    elif gtype in ("LineString", "MultiLineString"):
        lines = [geom["coordinates"]] if gtype == "LineString" else geom["coordinates"]
        for line in lines:
            rings.append(line)
    return rings


def plot_hap_distribution(
    samples: list[dict],
    hap_names: list[str],
    output_path: str | Path,
    title: str = "",
    hap_colors: list[str] | None = None,
    database: str = "china",
    geo_source: str | None = None,
    symbol_lim: tuple[float, float] = (0.3, 1.2),
    show_labels: bool = True,
    label_font_size: float = 5.5,
    legend_font_size: float = 7,
    dpi: int = 600,
    fmt: str | None = None,
) -> Path:
    """Haplotype geographic distribution with map + pie charts.

    Reads a GeoJSON boundary file for the base map layer, then overlays
    pie charts at each sample location showing haplotype proportions.
    Symbol size scaled by sqrt(total count), matching R's hapDistribution.

    Args:
        samples: List of dicts with keys 'lon', 'lat', 'hap'.
        hap_names: Haplotype names to display.
        output_path: Output file path (suffix sets format).
        title: Figure title.
        hap_colors: Colors per haplotype (auto from PALETTE if None).
        database: Map source key ('china') or 'none' for no map.
        geo_source: Direct URL/path to GeoJSON (overrides database).
        symbol_lim: (min, max) radius scaling range.
        show_labels: Show total-count labels at each location.
        label_font_size: Font size for location labels.
        legend_font_size: Font size for legend.
        dpi: Output resolution.
    """
    if not samples:
        raise ValueError("no samples provided")

    n_haps = len(hap_names)
    if hap_colors is None:
        hap_colors = [PALETTE[i % len(PALETTE)] for i in range(n_haps)]
    color_map = dict(zip(hap_names, hap_colors))

    # Aggregate samples by (lon, lat)
    loc_data: dict[tuple[float, float], dict[str, int]] = {}
    for s in samples:
        key = (round(s["lon"], 2), round(s["lat"], 2))
        if key not in loc_data:
            loc_data[key] = {h: 0 for h in hap_names}
        hap = s.get("hap", "")
        if hap in loc_data[key]:
            loc_data[key][hap] += 1

    lons = [k[0] for k in loc_data]
    lats = [k[1] for k in loc_data]
    lon_min, lon_max = min(lons), max(lons)
    lat_min, lat_max = min(lats), max(lats)
    lon_span = lon_max - lon_min if len(set(lons)) > 1 else 10
    lat_span = lat_max - lat_min if len(set(lats)) > 1 else 10

    pad_lon = max(lon_span * 0.12, 2)
    pad_lat = max(lat_span * 0.12, 2)

    fig_w = max(6, lon_span * 0.08 + 2)
    fig_h = max(5, lat_span * 0.08 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_aspect("equal")

    # ── Draw map from GeoJSON ──
    src = geo_source or _GEO_SOURCES.get(database)
    if src:
        features = _load_geojson(src)
        map_patches: list[Polygon] = []
        for feat in features:
            geom = feat.get("geometry")
            if not geom:
                continue
            for ring in _geom_to_polygons(geom):
                if len(ring) < 3:
                    continue
                map_patches.append(Polygon(ring, closed=True))
        if map_patches:
            pc = PatchCollection(
                map_patches,
                facecolor="#f5f5f0", edgecolor="#9e9e9e",
                linewidth=0.3, zorder=1,
            )
            ax.add_collection(pc)

    # Set extent after map is drawn
    ax.set_xlim(lon_min - pad_lon, lon_max + pad_lon)
    ax.set_ylim(lat_min - pad_lat, lat_max + pad_lat)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["left"].set_color("#b0b0b0")
    ax.spines["bottom"].set_color("#b0b0b0")
    ax.tick_params(labelsize=6, colors="#767676", length=3, width=0.5)
    ax.set_xlabel("Longitude", fontsize=7, color="#4d4d4d")
    ax.set_ylabel("Latitude", fontsize=7, color="#4d4d4d")

    # ── Sqrt frequency scaling (R symbol.lim logic) ──
    totals = [sum(d.values()) for d in loc_data.values()]
    sq = [math.sqrt(t) for t in totals]
    max_sq = max(sq) if sq else 1
    min_sq = min(v for v in sq if v > 0) if any(v > 0 for v in sq) else 0
    dif = max_sq - min_sq

    base = min(lon_span, lat_span) * 0.04

    for (lon, lat), counts in loc_data.items():
        total = sum(counts.values())
        if total == 0:
            continue

        s = math.sqrt(total)
        if dif > 0:
            norm = (s - min_sq) / dif
            radius = norm * (symbol_lim[1] - symbol_lim[0]) + symbol_lim[0]
        else:
            radius = (symbol_lim[0] + symbol_lim[1]) / 2
        scaled_r = radius * base

        # Pie slices
        start = 90
        for hap in hap_names:
            cnt = counts.get(hap, 0)
            if cnt == 0:
                continue
            span = 360 * cnt / total
            ax.add_patch(Wedge(
                (lon, lat), scaled_r,
                start, start + span,
                facecolor=color_map[hap],
                edgecolor="white", linewidth=0.3, zorder=5,
            ))
            start += span

        if show_labels:
            ax.text(
                lon, lat, str(total),
                ha="center", va="center",
                fontsize=label_font_size, fontweight="bold",
                color="#272727", zorder=6,
            )

    # ── Legend ──
    handles = [
        Patch(facecolor=color_map[h], edgecolor="#999", label=h)
        for h in hap_names
    ]
    ax.legend(handles=handles, fontsize=legend_font_size, loc="upper left", framealpha=0.9)

    if title:
        fig.suptitle(title, fontsize=9, fontweight="bold", y=0.96, color="#272727")

    fig.tight_layout(pad=1.5)
    return save_figure(fig, output_path, dpi=dpi, fmt=fmt)
