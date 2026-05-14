from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

from ._palette import (
    GOLDEN, SIDEBAR_BG, POP_BG, FUNC_COLORS, PALETTE,
    is_dark, allele_palette, make_legend_handles, save_figure,
)
from ._transform import (
    _unique_alleles, _indel_footnotes, transform_for_display,
    read_hap_summary_tsv,
)
from ._gff import parse_gff_features, classify_positions


def plot_hap_table(
    rows: list[list[str]],
    output_path: str | Path,
    pop_data: dict[str, str] | None = None,
    gff_path: str | Path | None = None,
    title: str = "",
    font_size: float = 7,
    title_font_size: float = 9,
    dpi: int = 600,
    fmt: str | None = None,
) -> Path:
    """Publication-quality haplotype summary table.

    Nature-style: Arial, editable text, luminance-aware label colors,
    600 DPI, SVG+PDF+TIFF export.
    CHR → title; header = POS + ALLELE (golden ratio, white bg).
    When gff_path: thin colored strip above table for functional category,
    SnpEff-style legend only (no allele legend).
    """
    if not rows:
        raise ValueError("input table is empty")

    rows, region_title = transform_for_display(rows, pop_data)

    n_rows = len(rows)
    n_cols = max(len(r) for r in rows)
    padded = [r + [""] * (n_cols - len(r)) for r in rows]

    # ── Header vs data rows ──
    header_set = {"POS", "ALLELE"}
    is_header = [r[0] in header_set for r in padded]

    # ── Find variant/split boundary from ALLELE row ──
    allele_row = next((r for r in padded if r[0] == "ALLELE"), None)
    var_end = n_cols
    pop_names_detected: list[str] = []
    if allele_row:
        for c in range(1, len(allele_row)):
            v = allele_row[c].strip()
            if v in {"n/N", "Accession"} or (v and v not in {"", "NA", "."}
                                               and "/" not in v
                                               and not any(ch.isdigit() for ch in v)
                                               and v not in PALETTE):
                var_end = c
                break
        for c in range(var_end, len(allele_row)):
            v = allele_row[c].strip()
            if v and v != "n/N" and v != "Accession":
                pop_names_detected.append(v)

    n_pop = len(pop_names_detected)

    alleles = _unique_alleles(padded, var_end)
    colors = allele_palette(alleles)

    padded, footnote = _indel_footnotes(padded, var_end)
    n_cols = max(len(r) for r in padded)
    padded = [r + [""] * (n_cols - len(r)) for r in padded]

    # ── GFF functional classification ──
    pos_row = next((r for r in padded if r[0] == "POS"), None)
    func_map: dict[int, str] = {}  # col index → category
    has_gff = False
    if gff_path and pos_row:
        chrom = region_title.split(":")[0] if ":" in region_title else ""
        positions: list[int] = []
        for c in range(1, var_end):
            v = pos_row[c].strip()
            if v.isdigit():
                positions.append(int(v))
        if positions and chrom:
            features = parse_gff_features(gff_path)
            func_cats = classify_positions(positions, chrom, features)
            for i, cat in enumerate(func_cats):
                func_map[i + 1] = cat
            has_gff = bool(func_map)

    # ── Row heights: header uses golden ratio ──
    ch_data = 0.42
    ch_header = ch_data * GOLDEN
    strip_h = 0.12  # thin functional category strip
    strip_total = strip_h if has_gff else 0.0

    row_heights = [ch_header if is_header[r] else ch_data for r in range(n_rows)]
    row_y = [strip_total] * (n_rows + 1)
    for r in range(n_rows):
        row_y[r + 1] = row_y[r] + row_heights[r]
    total_h = row_y[n_rows]

    cw = 0.72
    total_w = n_cols * cw
    fig_w = max(5, total_w + 1.0)
    fig_h = max(3, total_h + 1.6)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, total_h)
    ax.invert_yaxis()
    ax.axis("off")

    # ── Draw functional category strip above POS row ──
    if has_gff:
        for c, cat in func_map.items():
            x0 = c * cw
            fc = FUNC_COLORS.get(cat, "#cccccc")
            ax.add_patch(Rectangle(
                (x0, 0), cw, strip_h,
                facecolor=fc, edgecolor="white", linewidth=0.4,
            ))

    # ── Draw table cells ──
    for r, row in enumerate(padded):
        y0 = row_y[r]
        rh = row_heights[r]
        for c, val in enumerate(row):
            x0 = c * cw

            # ── background ──
            is_pop_col = n_pop > 0 and c >= var_end and c < var_end + n_pop
            is_acc_col = n_pop > 0 and c == var_end + n_pop
            is_stat_col = n_pop == 0 and c >= var_end

            if is_header[r]:
                bg = "#ffffff"
            elif is_pop_col:
                bg = POP_BG
            elif is_acc_col or is_stat_col:
                bg = SIDEBAR_BG
            elif val in colors:
                bg = colors[val]
            elif val in {"", "NA", "."}:
                bg = "#ffffff"
            else:
                bg = "#ffffff"

            ax.add_patch(Rectangle(
                (x0, y0), cw, rh,
                facecolor=bg, edgecolor="white", linewidth=0.6,
            ))

            # ── text style ──
            weight = "normal"
            size = font_size
            text_color = "#272727"

            if is_header[r]:
                weight = "bold"
            elif c == 0:
                # Haplotype name: bold, 1.5x larger
                weight = "bold"
                size = font_size * 1.5
                text_color = "#272727"
            elif val in colors:
                text_color = "white" if is_dark(bg) else "#272727"
            elif val in {"", "NA", "."}:
                text_color = "#b0b0b0"

            ax.text(
                x0 + cw / 2, y0 + rh / 2, val,
                ha="center", va="center",
                fontsize=size, fontweight=weight, color=text_color,
            )

    # ── separator between header and data ──
    allele_row_idx = next((i for i, r in enumerate(padded) if r[0] == "ALLELE"), None)
    if allele_row_idx is not None:
        sep_y = row_y[allele_row_idx + 1]
        ax.plot([0, total_w], [sep_y, sep_y], color="#767676", linewidth=1.0, zorder=5)

    # ── title ──
    display_title = title or region_title
    if title and region_title:
        display_title = f"{region_title} — {title}"
    if display_title:
        fig.suptitle(
            display_title, fontsize=title_font_size,
            fontweight="bold", y=0.96, ha="center", color="#272727",
        )

    # ── SnpEff-style legend (only when --gff) ──
    if has_gff and func_map:
        seen_cats: dict[str, str] = {}
        for cat in func_map.values():
            if cat not in seen_cats:
                seen_cats[cat] = FUNC_COLORS.get(cat, "#cccccc")
        legend_handles = make_legend_handles(seen_cats)
        fig.legend(
            handles=legend_handles, fontsize=font_size - 0.5,
            loc="lower center", ncol=min(len(legend_handles), 6),
            framealpha=0.9, handlelength=1.2, handleheight=0.8,
            columnspacing=1.0,
        )

    if footnote:
        fig.text(
            0.03, 0.01, footnote,
            fontsize=font_size - 2, color="#767676", style="italic",
        )

    fig.subplots_adjust(top=0.90, bottom=0.12, left=0.03, right=0.97)
    return save_figure(fig, output_path, dpi=dpi, fmt=fmt)
