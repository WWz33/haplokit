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
    population_pairwise_statistics,
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
    population_file: str | Path | None = None,
    population_delimiter: str = "auto",
    populations: Sequence[str] | None = None,
    title: str | None = None,
    fmt: str | None = None,
    figsize: tuple[float, float] | None = None,
) -> Path:
    """Plot phenotype distributions by haplotype as boxplots with jittered samples."""
    if output_path is None:
        raise PhenotypeError("output_path is required")
    if trait is None:
        raise PhenotypeError("trait is required")

    if phenotypes is None:
        records = tuple(haplotypes)
        if not records or not isinstance(records[0], PhenotypeRecord):
            raise PhenotypeError("phenotypes path is required unless haplotypes is a PhenotypeRecord sequence")
        hap_order = sort_haplotype_labels({record.haplotype for record in records})
        population_order = _resolve_population_strata(records, populations)
    else:
        dataset = load_phenotype_dataset(
            haplotypes,
            phenotypes,
            traits=[trait],
            hap_delimiter=hap_delimiter,
            phenotype_delimiter=phenotype_delimiter,
            population_file=population_file,
            population_delimiter=population_delimiter,
        )
        records = dataset.records
        hap_order = dataset.haplotypes
        population_order = list(populations) if populations is not None else list(dataset.populations or [None])

    grouped_panels = [
        (population, group_values(records, trait, min_hap_size=min_hap_size, haplotypes=hap_order, population=population))
        for population in population_order
    ]
    grouped_panels = [(population, grouped) for population, grouped in grouped_panels if grouped]
    if not grouped_panels:
        raise PhenotypeError(f"no haplotypes retained for trait {trait!r} with min_hap_size={min_hap_size}")

    _ensure_mpl()
    if any(population for population, _ in grouped_panels):
        fig = _plot_population_grouped_box(
            tuple(records),
            trait,
            grouped_panels,
            hap_order,
            min_hap_size,
            method,
            comparisons,
            title,
            figsize,
        )
    else:
        fig = _plot_haplotype_box(
            tuple(records),
            trait,
            grouped_panels[0][1],
            min_hap_size,
            method,
            comparisons,
            title,
            figsize,
        )

    fig.tight_layout()
    return save_figure(fig, output_path, fmt=fmt)


def _plot_haplotype_box(
    records: Sequence[PhenotypeRecord],
    trait: str,
    grouped: dict[str, list[float]],
    min_hap_size: int,
    method: str,
    comparisons: Sequence[tuple[str, str]] | None,
    title: str | None,
    figsize: tuple[float, float] | None,
):
    import matplotlib.pyplot as plt

    labels = list(grouped)
    data = [grouped[label] for label in labels]
    fig_width = max(5.0, 1.15 * len(labels) + 2.0)
    fig, ax = plt.subplots(figsize=figsize or (fig_width, 4.8))
    boxplot_kwargs = _boxplot_kwargs()
    try:
        box = ax.boxplot(data, tick_labels=labels, **boxplot_kwargs)
    except TypeError:
        box = ax.boxplot(data, labels=labels, **boxplot_kwargs)
    _color_boxes(box, [PALETTE[index % len(PALETTE)] for index in range(len(labels))])
    _scatter_points(ax, data, list(range(1, len(labels) + 1)), width=0.58)

    ax.set_xlabel("Haplotype")
    ax.set_ylabel(trait)
    ax.set_title(title or f"{trait} by haplotype")
    _style_box_axes(ax)
    ax.tick_params(axis="x", rotation=30 if max(len(label) for label in labels) > 7 else 0)

    if comparisons:
        _annotate_haplotype_comparisons(ax, records, trait, labels, comparisons, min_hap_size, method)
    return fig


def _plot_population_grouped_box(
    records: Sequence[PhenotypeRecord],
    trait: str,
    grouped_panels: Sequence[tuple[str | None, dict[str, list[float]]]],
    hap_order: Sequence[str],
    min_hap_size: int,
    method: str,
    comparisons: Sequence[tuple[str, str]] | None,
    title: str | None,
    figsize: tuple[float, float] | None,
):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    population_panels = [(str(population), grouped) for population, grouped in grouped_panels if population]
    retained_haps = [hap for hap in hap_order if any(hap in grouped for _, grouped in population_panels)]
    if not retained_haps:
        raise PhenotypeError(f"no haplotypes retained for trait {trait!r} with min_hap_size={min_hap_size}")

    if len(retained_haps) == 1:
        return _plot_single_haplotype_by_population(
            records,
            trait,
            population_panels,
            retained_haps[0],
            min_hap_size,
            method,
            title,
            figsize,
        )

    population_labels = [population for population, _ in population_panels]
    slot_width = 0.72 / max(len(retained_haps), 1)
    box_width = min(0.22, slot_width * 0.78)
    offsets = [(index - (len(retained_haps) - 1) / 2) * slot_width for index in range(len(retained_haps))]
    position_map: dict[tuple[str, str], float] = {}
    data: list[list[float]] = []
    positions: list[float] = []
    colors: list[str] = []

    for population_index, (population, grouped) in enumerate(population_panels, start=1):
        for hap_index, haplotype in enumerate(retained_haps):
            values = grouped.get(haplotype)
            if not values:
                continue
            position = population_index + offsets[hap_index]
            position_map[(population, haplotype)] = position
            data.append(values)
            positions.append(position)
            colors.append(PALETTE[hap_index % len(PALETTE)])

    fig_width = max(6.0, len(population_labels) * (1.15 + 0.2 * len(retained_haps)) + 2.3)
    fig, ax = plt.subplots(figsize=figsize or (fig_width, 4.9))
    box = ax.boxplot(data, positions=positions, widths=box_width, **_boxplot_kwargs())
    _color_boxes(box, colors)
    _scatter_points(ax, data, positions, width=box_width)

    ax.set_xticks(range(1, len(population_labels) + 1))
    ax.set_xticklabels(population_labels)
    ax.set_xlim(0.45, len(population_labels) + 0.55)
    ax.set_xlabel("Population")
    ax.set_ylabel(trait)
    ax.set_title(title or f"{trait} by population and haplotype")
    _style_box_axes(ax)
    ax.legend(
        handles=[
            Patch(facecolor=PALETTE[index % len(PALETTE)], edgecolor="#444444", label=haplotype, alpha=0.55)
            for index, haplotype in enumerate(retained_haps)
        ],
        fontsize=5,
        loc="upper left",
        frameon=False,
    )

    annotations = _within_population_annotations(
        records,
        trait,
        population_labels,
        min_hap_size,
        method,
        position_map,
        comparisons,
    )
    annotations.extend(
        _between_population_annotations(
            records,
            trait,
            retained_haps,
            population_labels,
            min_hap_size,
            method,
            position_map,
        )
    )
    _draw_stat_annotations_inside(ax, annotations)
    return fig


def _plot_single_haplotype_by_population(
    records: Sequence[PhenotypeRecord],
    trait: str,
    population_panels: Sequence[tuple[str, dict[str, list[float]]]],
    haplotype: str,
    min_hap_size: int,
    method: str,
    title: str | None,
    figsize: tuple[float, float] | None,
):
    import matplotlib.pyplot as plt

    population_labels = [population for population, grouped in population_panels if haplotype in grouped]
    data = [grouped[haplotype] for _, grouped in population_panels if haplotype in grouped]
    positions = list(range(1, len(population_labels) + 1))
    position_map = {(population, haplotype): position for population, position in zip(population_labels, positions)}

    fig_width = max(5.0, 1.05 * len(population_labels) + 2.0)
    fig, ax = plt.subplots(figsize=figsize or (fig_width, 4.8))
    box = ax.boxplot(data, positions=positions, widths=0.58, **_boxplot_kwargs())
    _color_boxes(box, [PALETTE[0] for _ in data])
    _scatter_points(ax, data, positions, width=0.58)

    ax.set_xticks(positions)
    ax.set_xticklabels(population_labels)
    ax.set_xlabel("Population")
    ax.set_ylabel(trait)
    ax.set_title(title or f"{trait} by population ({haplotype})")
    _style_box_axes(ax)

    annotations = _between_population_annotations(
        records,
        trait,
        [haplotype],
        population_labels,
        min_hap_size,
        method,
        position_map,
    )
    _draw_stat_annotations_inside(ax, annotations)
    return fig


def _boxplot_kwargs() -> dict[str, object]:
    return {
        "patch_artist": True,
        "showfliers": False,
        "medianprops": {"color": "#222222", "linewidth": 1.2},
        "whiskerprops": {"color": "#555555", "linewidth": 1.0},
        "capprops": {"color": "#555555", "linewidth": 1.0},
        "boxprops": {"edgecolor": "#444444", "linewidth": 1.0},
    }


def _style_box_axes(ax) -> None:
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#333333")
        spine.set_linewidth(1.0)
    ax.tick_params(colors="#222222", width=0.9, length=4)


def _color_boxes(box, colors: Sequence[str]) -> None:
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)


def _scatter_points(ax, data: Sequence[Sequence[float]], positions: Sequence[float], width: float) -> None:
    rng = random.Random(1729)
    jitter_width = max(width * 0.22, 0.035)
    for x_pos, values in zip(positions, data):
        jitter = [x_pos + rng.uniform(-jitter_width, jitter_width) for _ in values]
        ax.scatter(jitter, values, s=18, color="#222222", alpha=0.72, linewidths=0, zorder=3)


def _resolve_population_strata(
    records: Sequence[PhenotypeRecord],
    populations: Sequence[str] | None,
) -> list[str | None]:
    if populations is not None:
        return list(populations)
    discovered = sort_haplotype_labels({record.population for record in records if record.population})
    return discovered if discovered else [None]


def _annotate_haplotype_comparisons(
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
    annotations: list[tuple[float, float, str]] = []

    for group1, group2 in comparisons:
        if group1 not in label_to_x or group2 not in label_to_x:
            continue
        row = row_map.get(frozenset((group1, group2)))
        if row is None:
            continue
        annotations.append((label_to_x[group1], label_to_x[group2], _row_significance(row)))
    _draw_stat_annotations_inside(ax, annotations)


def _within_population_annotations(
    records: Sequence[PhenotypeRecord],
    trait: str,
    populations: Sequence[str],
    min_hap_size: int,
    method: str,
    position_map: dict[tuple[str, str], float],
    comparisons: Sequence[tuple[str, str]] | None,
) -> list[tuple[float, float, str]]:
    allowed = {frozenset((group1, group2)) for group1, group2 in comparisons} if comparisons else None
    annotations: list[tuple[float, float, str]] = []
    for population in populations:
        rows = pairwise_statistics(
            records,
            traits=[trait],
            min_hap_size=min_hap_size,
            method=method,
            populations=[population],
        )
        for row in rows:
            group1 = str(row["group1"])
            group2 = str(row["group2"])
            if allowed is not None and frozenset((group1, group2)) not in allowed:
                continue
            x1 = position_map.get((population, group1))
            x2 = position_map.get((population, group2))
            if x1 is not None and x2 is not None:
                annotations.append((x1, x2, _row_significance(row)))
    return annotations


def _between_population_annotations(
    records: Sequence[PhenotypeRecord],
    trait: str,
    haplotypes: Sequence[str],
    populations: Sequence[str],
    min_hap_size: int,
    method: str,
    position_map: dict[tuple[str, str], float],
) -> list[tuple[float, float, str]]:
    rows = population_pairwise_statistics(
        records,
        traits=[trait],
        min_hap_size=min_hap_size,
        method=method,
        haplotypes=haplotypes,
        populations=populations,
    )
    annotations: list[tuple[float, float, str]] = []
    for row in rows:
        haplotype = str(row["haplotype"])
        population1 = str(row["group1"])
        population2 = str(row["group2"])
        x1 = position_map.get((population1, haplotype))
        x2 = position_map.get((population2, haplotype))
        if x1 is not None and x2 is not None:
            annotations.append((x1, x2, _row_significance(row)))
    return annotations


def _row_significance(row: dict[str, object]) -> str:
    p_adjusted = row.get("p_adjusted", math.nan)
    return significance_label(float(p_adjusted)) if isinstance(p_adjusted, (float, int)) else "NA"


def _draw_stat_annotations_inside(ax, annotations: Sequence[tuple[float, float, str]]) -> None:
    if not annotations:
        return
    sorted_annotations = sorted(annotations, key=lambda item: abs(item[1] - item[0]), reverse=True)
    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1.0)
    extra = min(1.0, max(0.18, 0.045 * len(sorted_annotations)))
    ax.set_ylim(y_min, y_max + span * extra)
    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1.0)
    start = y_max - span * 0.045
    step = span * 0.038
    height = span * 0.012

    for index, (x1, x2, label) in enumerate(sorted_annotations):
        y = start - index * step
        if y <= y_min + span * 0.04:
            break
        _draw_bracket(ax, x1, x2, y, height, label)


def _draw_bracket(ax, x1: float, x2: float, y: float, height: float, label: str) -> None:
    if x1 > x2:
        x1, x2 = x2, x1
    ax.plot([x1, x1, x2, x2], [y - height, y, y, y - height], color="#333333", linewidth=0.9)
    ax.text((x1 + x2) / 2, y - height * 0.25, label, ha="center", va="top", fontsize=9)
