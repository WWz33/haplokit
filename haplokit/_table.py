from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

from ._palette import (
    GOLDEN, SIDEBAR_BG, POP_BG, PALETTE,
    is_dark, allele_palette, save_figure,
)
from ._transform import (
    _unique_alleles, _indel_footnotes, transform_for_display,
    read_hap_summary_tsv,
)
from ._genemodel import (
    CDS_COLOR, UTR5_COLOR, UTR3_COLOR, BACKBONE,
    load_genes, select_view, draw_gene_track,
)

# guide line / SNP tick colors
_GUIDE = "#9aa0a6"
_TICK_GENIC = "#c0392b"
_TICK_INTERGENIC = "#9aa0a6"

# Per-theme table geometry. The gene-model track, header/Hap styling, separators
# and title are shared; themes differ only in cell proportions and row heights.
#   detailed — square-ish cells (0.72 x 0.42 ≈ 1.71:1), header shrunk by GOLDEN.
#   compact  — flat wide cells (254:63 ≈ 4.03:1) flush with no white edge,
#              header ≈0.85x data height, tighter spacing, narrower label column.
_THEMES = {
    "detailed": {"ch_data": 0.42, "cw": 0.72, "hdr_ratio": GOLDEN,
                 "edge": "white", "lw": 0.6, "hspace": 0.06, "legend_gap": 0.012, "model_reserve": 1.3,
                 "label_cw_ratio": 1.0, "pos_row_shrink": 1.0},
    "compact":  {"ch_data": 0.32, "cw": 0.32 * (254 / 63), "hdr_ratio": 0.85,
                 "edge": "none", "lw": 0.0, "hspace": 0.02, "legend_gap": 0.001, "model_reserve": 1.0,
                 "label_cw_ratio": 0.7, "pos_row_shrink": 0.6},
}


def plot_hap_table(
    rows: list[list[str]],
    output_path: str | Path,
    pop_data: dict[str, str] | None = None,
    gff_path: str | Path | None = None,
    title: str = "",
    font_size: float = 9,
    title_font_size: float = 9,
    dpi: int = 600,
    fmt: str | None = None,
    figsize: tuple[float, float] | None = None,
    theme: str = "detailed",
) -> Path:
    """Publication-quality haplotype summary table.

    Nature-style: Arial, editable text, luminance-aware label colors,
    600 DPI, SVG+PDF+TIFF export. CHR → title; header = POS + ALLELE.

    When ``gff_path`` is given and a gene can be anchored on the variant
    positions, a pyGenomeTracks-style gene model is drawn on top of the table
    (backbone/introns, CDS, UTR, strand arrow) with SNP ticks and guide lines
    to the matching allele columns, plus a compact CDS/UTR/intron legend. With
    no GFF (or no gene found) only the table is drawn.

    ``theme`` selects the table look: ``detailed`` (default) or ``compact``.
    """
    if not rows:
        raise ValueError("input table is empty")
    if theme not in _THEMES:
        raise ValueError(f"unknown theme {theme!r}; expected one of {sorted(_THEMES)}")
    geom = _THEMES[theme]

    rows, region_title = transform_for_display(rows, pop_data)

    n_rows = len(rows)
    n_cols = max(len(r) for r in rows)
    padded = [r + [""] * (n_cols - len(r)) for r in rows]

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

    # ── Variant positions (allele data columns 1..var_end-1) ──
    pos_row = next((r for r in padded if r[0] == "POS"), None)
    positions: list[int] = []
    if pos_row:
        for c in range(1, var_end):
            v = pos_row[c].strip()
            positions.append(int(v) if v.isdigit() else None)

    # ── Gene model: load genes + select a view only when a GFF is supplied ──
    genes: list = []
    region_lo = region_hi = None
    chrom = region_title.split(":")[0] if ":" in region_title else ""
    if gff_path and chrom:
        snp_positions = [p for p in positions if p is not None]
        if snp_positions:
            genes = load_genes(gff_path, chrom)
            region_lo, region_hi, genes = select_view(genes, snp_positions)
    draw_model = bool(genes) and region_lo is not None and region_hi is not None

    # ── Row heights ──
    ch_data = geom["ch_data"]
    ch_header = ch_data * geom["hdr_ratio"]
    row_heights = [ch_header if is_header[r] else ch_data for r in range(n_rows)]

    # Shrink POS row height (compact theme)
    pos_row_idx = next((i for i, r in enumerate(padded) if r[0] == "POS"), None)
    if pos_row_idx is not None:
        row_heights[pos_row_idx] = ch_header * geom["pos_row_shrink"]

    row_y = [0.0] * (n_rows + 1)
    for r in range(n_rows):
        row_y[r + 1] = row_y[r] + row_heights[r]
    total_h = row_y[n_rows]

    # ── Column widths: narrower first column (compact theme) ──
    cw = geom["cw"]
    label_cw = cw * geom["label_cw_ratio"]
    total_w = label_cw + (n_cols - 1) * cw
    fig_w = max(5, total_w + 1.0)
    fig_h = max(3, total_h + 1.6) + (geom["model_reserve"] if draw_model else 0.0)

    # ── Figure: gene-model band on top of the table when drawn ──
    if draw_model:
        fig = plt.figure(figsize=figsize or (fig_w, fig_h))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.0, total_h], hspace=geom["hspace"])
        ax_model = fig.add_subplot(gs[0])
        ax = fig.add_subplot(gs[1])
    else:
        fig, ax = plt.subplots(figsize=figsize or (fig_w, fig_h))
        ax_model = None
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, total_h)
    ax.invert_yaxis()
    ax.axis("off")

    # ── Draw table cells ──
    for r, row in enumerate(padded):
        y0 = row_y[r]
        rh = row_heights[r]
        x0 = 0.0
        for c, val in enumerate(row):
            # First column uses narrower width (compact theme)
            col_w = label_cw if c == 0 else cw

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
            else:
                bg = "#ffffff"

            ax.add_patch(Rectangle(
                (x0, y0), col_w, rh,
                facecolor=bg, edgecolor=geom["edge"], linewidth=geom["lw"],
            ))

            # text style — bold everywhere except Hap labels (c==0, non-header)
            if c == 0 and not is_header[r]:
                weight = "normal"  # Hap01, Hap02... not bold
            else:
                weight = "bold"    # POS, ALLELE, all data cells, n/N col

            text_color = "#000000"  # pure black for table text
            if not is_header[r] and c != 0:
                if val in colors:
                    text_color = "white" if is_dark(bg) else "#000000"
                elif val in {"", "NA", "."}:
                    text_color = "#b0b0b0"

            ax.text(
                x0 + col_w / 2, y0 + rh / 2, val,
                ha="center", va="center",
                fontsize=font_size, fontweight=weight, color=text_color,
            )
            x0 += col_w

    # ── Gene model track + SNP ticks + guide lines ──
    if draw_model:
        yc, h_cds = draw_gene_track(ax_model, genes, region_lo, region_hi)
        # align the model to span ONLY the allele data columns (1..var_end-1)
        pos_t = ax.get_position()
        pos_m = ax_model.get_position()
        data_cols_start_frac = label_cw / total_w
        data_cols_width_frac = ((var_end - 1) * cw) / total_w
        new_x = pos_t.x0 + data_cols_start_frac * pos_t.width
        new_w = data_cols_width_frac * pos_t.width
        ax_model.set_position([new_x, pos_m.y0, new_w, pos_m.height])

        overshoot = 0.08
        tick_bottom = yc - h_cds / 2 - overshoot
        tick_top = yc + h_cds / 2 + overshoot

        def _genic(p: int) -> bool:
            return any(g.lo <= p <= g.hi for g in genes)

        for i, c in enumerate(range(1, var_end)):
            if i >= len(positions) or positions[i] is None:
                continue
            gpos = positions[i]
            # Adjust column center for variable-width first column
            col_center = label_cw + i * cw + cw / 2
            inter = not _genic(gpos)
            tick_col = _TICK_INTERGENIC if inter else _TICK_GENIC
            con = ConnectionPatch(
                xyA=(col_center, 0), coordsA=ax.transData,
                xyB=(gpos, tick_bottom), coordsB=ax_model.transData,
                color=_GUIDE, lw=0.7, alpha=0.7, zorder=0,
            )
            fig.add_artist(con)
            ax_model.plot([gpos, gpos], [tick_bottom, tick_top],
                          color=tick_col, lw=1.2, alpha=0.9, zorder=6,
                          linestyle="--" if inter else "-")
            ax_model.plot([gpos], [tick_top], marker="v", markersize=3.5,
                          color=tick_col, zorder=7)

        # CDS / 5' UTR / 3' UTR legend, centered just above the model
        pos_m = ax_model.get_position()
        ly = min(pos_m.y1 + geom["legend_gap"], 0.985)
        sw, gap = 0.018, 0.006
        items = [(CDS_COLOR, "CDS"), (UTR5_COLOR, "5' UTR"), (UTR3_COLOR, "3' UTR")]
        widths = [sw + gap + 0.014 * len(lab) + 0.03 for _, lab in items]
        lx = pos_m.x0 + pos_m.width / 2 - sum(widths) / 2
        for (col, lab), w in zip(items, widths):
            fig.patches.append(Rectangle(
                (lx, ly), sw, 0.016, facecolor=col, edgecolor="white",
                lw=0.4, transform=fig.transFigure, figure=fig, zorder=10))
            fig.text(lx + sw + gap, ly + 0.008, lab, ha="left", va="center",
                     fontsize=font_size - 0.5, color="#000000")
            lx += w

    # ── title ──
    display_title = title or region_title
    if title and region_title:
        display_title = f"{region_title} — {title}"
    if display_title:
        fig.suptitle(
            display_title, fontsize=title_font_size,
            fontweight="bold", y=0.98, ha="center", color="#000000",
        )

    if footnote:
        fig.text(
            0.03, 0.01, footnote,
            fontsize=9, color="#767676", style="italic",
        )

    if not draw_model:
        fig.subplots_adjust(top=0.90, bottom=0.12, left=0.03, right=0.97)
    return save_figure(fig, output_path, dpi=dpi, fmt=fmt)
