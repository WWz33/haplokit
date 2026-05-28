from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Sequence

from ._palette import PALETTE, _ensure_mpl, save_figure
from ._phenotype import (
    PhenotypeError,
    PhenotypeRecord,
    group_values,
    load_phenotype_dataset,
    pairwise_statistics,
    significance_label,
    sort_haplotype_labels,
)


def plot_hap_phenotype_box(
    haplotypes,
    phenotypes: str | Path | None = None,
    trait: str | None = None,
    output_path: str | Path | None = None,
    *,
    min_hap_size: int = 5,
    method: str = "welch",
    comparisons: Sequence[tuple[str, str]] | None = None,
    hap_delimiter: str = "auto",
    phenotype_delimiter: str = "auto",
    title: str | None = None,
    fmt: str | None = None,
) -> Path:
    """Plot phenotype distributions by haplotype as a boxplot with jittered samples."""
    if output_path is None:
        raise PhenotypeError("output_path is required")
    if trait is None:
        raise PhenotypeError("trait is required")

    if phenotypes is None:
        records = tuple(haplotypes)
        if not records or not isinstance(records[0], PhenotypeRecord):
            raise PhenotypeError("phenotypes path is required unless haplotypes is a PhenotypeRecord sequence")
        hap_order = sort_haplotype_labels({record.haplotype for record in records})
    else:
        dataset = load_phenotype_dataset(
            haplotypes,
            phenotypes,
            traits=[trait],
            hap_delimiter=hap_delimiter,
            phenotype_delimiter=phenotype_delimiter,
        )
        records = dataset.records
        hap_order = dataset.haplotypes

    grouped = group_values(records, trait, min_hap_size=min_hap_size, haplotypes=hap_order)
    if not grouped:
        raise PhenotypeError(f"no haplotypes retained for trait {trait!r} with min_hap_size={min_hap_size}")

    _ensure_mpl()
    import matplotlib.pyplot as plt

    labels = list(grouped)
    data = [grouped[label] for label in labels]
    fig_width = max(5.0, 1.15 * len(labels) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 4.8))

    boxplot_kwargs = {
        "patch_artist": True,
        "widths": 0.58,
        "showfliers": False,
        "medianprops": {"color": "#222222", "linewidth": 1.2},
        "whiskerprops": {"color": "#555555", "linewidth": 1.0},
        "capprops": {"color": "#555555", "linewidth": 1.0},
        "boxprops": {"edgecolor": "#444444", "linewidth": 1.0},
    }
    try:
        box = ax.boxplot(data, tick_labels=labels, **boxplot_kwargs)
    except TypeError:
        box = ax.boxplot(data, labels=labels, **boxplot_kwargs)
    for index, patch in enumerate(box["boxes"]):
        patch.set_facecolor(PALETTE[index % len(PALETTE)])
        patch.set_alpha(0.55)

    rng = random.Random(1729)
    for x_pos, values in enumerate(data, start=1):
        jitter = [x_pos + rng.uniform(-0.12, 0.12) for _ in values]
        ax.scatter(jitter, values, s=18, color="#222222", alpha=0.72, linewidths=0, zorder=3)

    ax.set_xlabel("Haplotype")
    ax.set_ylabel(trait)
    ax.set_title(title or f"{trait} by haplotype")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, alpha=0.8)
    ax.tick_params(axis="x", rotation=30 if max(len(label) for label in labels) > 7 else 0)

    if comparisons:
        _annotate_comparisons(ax, records, trait, labels, comparisons, min_hap_size, method)

    return save_figure(fig, output_path, fmt=fmt)


def _annotate_comparisons(
    ax,
    records: Sequence[PhenotypeRecord],
    trait: str,
    labels: Sequence[str],
    comparisons: Sequence[tuple[str, str]],
    min_hap_size: int,
    method: str,
) -> None:
    rows = pairwise_statistics(records, traits=[trait], min_hap_size=min_hap_size, method=method)
    row_map = {frozenset((str(row["group1"]), str(row["group2"]))): row for row in rows}
    label_to_x = {label: idx for idx, label in enumerate(labels, start=1)}

    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1.0)
    base = y_max + span * 0.06
    step = span * 0.1
    drawn = 0

    for group1, group2 in comparisons:
        if group1 not in label_to_x or group2 not in label_to_x:
            continue
        row = row_map.get(frozenset((group1, group2)))
        if row is None:
            continue
        p_adjusted = row.get("p_adjusted", math.nan)
        label = significance_label(float(p_adjusted)) if isinstance(p_adjusted, (float, int)) else "NA"
        x1 = label_to_x[group1]
        x2 = label_to_x[group2]
        y = base + drawn * step
        _draw_bracket(ax, x1, x2, y, step * 0.28, label)
        drawn += 1

    if drawn:
        ax.set_ylim(top=base + drawn * step)


def _draw_bracket(ax, x1: float, x2: float, y: float, height: float, label: str) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color="#333333", linewidth=0.9)
    ax.text((x1 + x2) / 2, y + height, label, ha="center", va="bottom", fontsize=9)
