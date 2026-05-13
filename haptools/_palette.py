from __future__ import annotations

from matplotlib.patches import Patch

# ── Mandatory Nature-style rcParams ─────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",       # text as <text> nodes
    "pdf.fonttype": 42,           # editable TrueType in PDF
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

# ColorBrewer Set2 qualitative — print-friendly, colorblind-safe
PALETTE = [
    "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
    "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#666666",
]

# Genomics scaffolding greys
META_BG = "#e8e8e8"
ALLELE_BG = "#d9d9d9"
SIDEBAR_BG = "#f2f2f2"
POP_BG = "#eef0f4"

# Functional category colors (SnpEff-inspired, genomics palette)
FUNC_COLORS = {
    "CDS": "#3775BA",
    "5'UTR": "#7BAFD4",
    "3'UTR": "#7BAFD4",
    "UTR": "#7BAFD4",
    "exon": "#B4C0E4",
    "intron": "#D8D8D8",
    "upstream": "#A6D854",
    "downstream": "#A6D854",
    "intergenic": "#F2F2F2",
}

GOLDEN = 0.618


def is_dark(hex_color: str, threshold: float = 0.5) -> bool:
    """Return True if hex color is dark (use white text)."""
    c = hex_color.lstrip("#")
    r, g, b = int(c[0:2], 16) / 255, int(c[2:4], 16) / 255, int(c[4:6], 16) / 255
    return (0.299 * r + 0.587 * g + 0.114 * b) < threshold


def allele_palette(alleles: list[str]) -> dict[str, str]:
    """Map alleles to ColorBrewer Set2 colors."""
    return {a: PALETTE[i % len(PALETTE)] for i, a in enumerate(alleles)}


def make_legend_handles(mapping: dict[str, str]) -> list[Patch]:
    """Create legend Patch handles from a name→color mapping."""
    return [Patch(facecolor=c, edgecolor="#999", label=n) for n, c in mapping.items()]


_FMT_MAP = {".svg": "svg", ".pdf": "pdf", ".png": "png", ".tiff": "tiff"}


def save_figure(fig, output_path: str | Path, dpi: int = 600, fmt: str | None = None) -> Path:
    """Save figure in requested format (default png from path suffix or png).

    Args:
        fig: matplotlib Figure.
        output_path: Output file path; suffix determines format if fmt is None.
        dpi: Raster resolution (png/tiff).
        fmt: Override format: 'svg', 'pdf', 'png', 'tiff'.
    """
    from pathlib import Path
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if fmt is None:
        fmt = _FMT_MAP.get(out.suffix.lower(), "png")
    target = out.with_suffix(f".{fmt}")

    save_kw: dict = {"bbox_inches": "tight"}
    if fmt in ("png", "tiff"):
        save_kw["dpi"] = dpi
    fig.savefig(target, **save_kw)
    plt.close(fig)
    return target.resolve()
