"""Gene-model track for the haplotype table figure.

Parses a GFF3/GTF annotation into per-gene loci (each carrying its longest
transcript's CDS/exon/UTR features), selects a sensible genomic view around the
variant positions, and draws a pyGenomeTracks-style model on a matplotlib axis:
a continuous backbone (introns / intergenic context), thick CDS boxes, thin
UTR/exon boxes, and a terminal CDS arrow encoding strand (flybase style).

This is the structural complement to ``_gff.classify_positions`` (which only
labels individual SNP positions by overlapping feature). The two are kept
separate because the model needs transcript-grouped geometry and strand, while
classification needs only a flat per-position lookup.
"""
from __future__ import annotations

import csv
from pathlib import Path

from matplotlib.patches import Rectangle, Polygon

# pyGenomeTracks-ish palette (structural, not the allele palette)
# ColorBrewer Paired: deep blue CDS, light blue 5' UTR, light green 3' UTR
CDS_COLOR = "#1F78B4"
UTR5_COLOR = "#A6CEE3"
UTR3_COLOR = "#B2DF8A"
BACKBONE = "#1a1a1a"

FLANK = 2000  # bp of context drawn on each side of a single hit gene


def _attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


class Gene:
    """One gene locus and the features of its longest (max-CDS) transcript."""

    def __init__(self, gid: str, chrom: str, lo: int, hi: int, strand: str):
        self.id = gid
        self.chrom = chrom
        self.lo = lo
        self.hi = hi
        self.strand = strand
        self.feats: list[tuple[str, int, int]] = []  # (kind, start, end)


def load_genes(gff_path: str | Path, chrom: str) -> list[Gene]:
    """Parse all genes on `chrom`; attach each one's longest-transcript features."""
    genes: dict[str, Gene] = {}
    tx_parent: dict[str, str] = {}  # transcript ID → parent gene ID
    tx_feats: dict[str, list[tuple[str, int, int]]] = {}  # transcript ID → features

    with Path(gff_path).open(encoding="utf-8", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 9 or row[0].startswith("#") or row[0] != chrom:
                continue
            ftype = row[2].lower()
            try:
                start, end = int(row[3]), int(row[4])
            except ValueError:
                continue
            strand = row[6]
            attrs = _attrs(row[8])
            if ftype == "gene":
                gid = attrs.get("ID", "")
                genes[gid] = Gene(gid, chrom, start, end, strand)
            elif ftype in {"mrna", "transcript"}:
                tx_parent[attrs.get("ID", "")] = attrs.get("Parent", "")
            else:
                # Distinguish 5' UTR, 3' UTR, CDS, exon
                if "utr" in ftype:
                    if "five" in ftype or "5" in ftype:
                        kind = "5UTR"
                    elif "three" in ftype or "3" in ftype:
                        kind = "3UTR"
                    else:
                        kind = "UTR"  # generic fallback
                elif ftype == "cds":
                    kind = "CDS"
                elif ftype == "exon":
                    kind = "exon"
                else:
                    kind = None
                if kind is None:
                    continue
                tx_feats.setdefault(attrs.get("Parent", ""), []).append((kind, start, end))

    # pick the longest (max total CDS) transcript per gene
    by_gene_tx: dict[str, dict[str, list]] = {}
    for tx, parent in tx_parent.items():
        by_gene_tx.setdefault(parent, {})[tx] = tx_feats.get(tx, [])

    for gid, gene in genes.items():
        txs = by_gene_tx.get(gid, {})
        if not txs:
            continue
        best = max(txs, key=lambda t: sum(e - s + 1 for k, s, e in txs[t] if k == "CDS"))
        gene.feats = sorted(txs[best], key=lambda f: f[1])

    return sorted(genes.values(), key=lambda g: g.lo)


def select_view(genes: list[Gene], snps: list[int]):
    """Apply the view-selection rule. Returns (view_lo, view_hi, hit_genes).

    * SNPs inside ONE gene → that gene + 2 kb flank, flank truncated at any
      neighbouring gene so we never bleed into another locus.
    * SNPs spanning 2+ genes → true coords, no flank, union of hit genes + SNPs.

    Returns (None, None, []) when there are no SNPs or no genes to anchor on.
    """
    if not snps or not genes:
        return None, None, []
    smin, smax = min(snps), max(snps)
    # genes that actually contain at least one SNP
    hit = [g for g in genes if any(g.lo <= p <= g.hi for p in snps)]

    if len(hit) <= 1:
        g = hit[0] if hit else min(genes, key=lambda x: abs(x.lo - smin))
        view_lo, view_hi = g.lo - FLANK, g.hi + FLANK
        for other in genes:
            if other.id == g.id:
                continue
            if other.hi < g.lo:  # neighbour upstream
                view_lo = max(view_lo, other.hi + 1)
            if other.lo > g.hi:  # neighbour downstream
                view_hi = min(view_hi, other.lo - 1)
        view_lo = min(view_lo, smin)  # never crop a SNP out of frame
        view_hi = max(view_hi, smax)
        return view_lo, view_hi, hit

    view_lo = min(min(g.lo for g in hit), smin)
    view_hi = max(max(g.hi for g in hit), smax)
    return view_lo, view_hi, hit


def draw_one_gene(ax, gene: Gene, region_lo: int, region_hi: int, label: bool = True):
    """Draw one gene model (flybase style) in genomic x-coords on `ax`."""
    yc = 0.5
    h_cds = 0.16
    h_utr = h_cds / 2

    feats = [(k, max(s, region_lo), min(e, region_hi))
             for k, s, e in gene.feats if min(e, region_hi) >= max(s, region_lo)]
    exons = [f for f in feats if f[0] == "exon"]
    cds = sorted([f for f in feats if f[0] == "CDS"], key=lambda f: f[1])
    utrs_5 = [f for f in feats if f[0] == "5UTR"]
    utrs_3 = [f for f in feats if f[0] == "3UTR"]
    utrs_generic = [f for f in feats if f[0] == "UTR"]

    g_lo = max(gene.lo, region_lo)
    g_hi = min(gene.hi, region_hi)

    # thin boxes: prefer explicit exons as base layer
    for _, s, e in (exons if exons else (utrs_5 + utrs_3 + utrs_generic)):
        ax.add_patch(Rectangle((s, yc - h_utr / 2), e - s, h_utr,
                               facecolor="#D3D3D3", edgecolor="white", lw=0.6, zorder=2))

    # 5' UTR in light blue
    for _, s, e in utrs_5:
        ax.add_patch(Rectangle((s, yc - h_utr / 2), e - s, h_utr,
                               facecolor=UTR5_COLOR, edgecolor="white", lw=0.6, zorder=3))

    # 3' UTR in light green
    for _, s, e in utrs_3:
        ax.add_patch(Rectangle((s, yc - h_utr / 2), e - s, h_utr,
                               facecolor=UTR3_COLOR, edgecolor="white", lw=0.6, zorder=3))

    # Generic UTR (fallback if no 5'/3' distinction)
    for _, s, e in utrs_generic:
        ax.add_patch(Rectangle((s, yc - h_utr / 2), e - s, h_utr,
                               facecolor=UTR5_COLOR, edgecolor="white", lw=0.6, zorder=3))

    arrow_w = (region_hi - region_lo) * 0.025
    for i, (_, s, e) in enumerate(cds):
        is_terminal = (gene.strand == "+" and i == len(cds) - 1) or \
                      (gene.strand == "-" and i == 0)
        clipped = (e >= region_hi) if gene.strand == "+" else (s <= region_lo)
        # terminal CDS becomes an arrow only at the gene's TRUE end inside the view
        true_end = (gene.hi <= region_hi) if gene.strand == "+" else (gene.lo >= region_lo)
        if is_terminal and true_end and not clipped and (e - s) > arrow_w:
            if gene.strand == "+":
                verts = [(s, yc - h_cds / 2), (s, yc + h_cds / 2),
                         (e - arrow_w, yc + h_cds / 2), (e, yc),
                         (e - arrow_w, yc - h_cds / 2)]
            else:
                verts = [(e, yc - h_cds / 2), (e, yc + h_cds / 2),
                         (s + arrow_w, yc + h_cds / 2), (s, yc),
                         (s + arrow_w, yc - h_cds / 2)]
            ax.add_patch(Polygon(verts, closed=True, facecolor=CDS_COLOR,
                                 edgecolor="white", lw=0.6, zorder=4))
        else:
            ax.add_patch(Rectangle((s, yc - h_cds / 2), e - s, h_cds,
                                   facecolor=CDS_COLOR, edgecolor="white", lw=0.6, zorder=4))

    if label:
        ax.text((g_lo + g_hi) / 2, yc + h_cds / 2 + 0.14, gene.id,
                ha="center", va="bottom", fontsize=7, color="#000000")


def draw_gene_track(ax, genes: list[Gene], region_lo: int, region_hi: int):
    """Draw all hit genes on one shared track. Returns (yc, h_cds)."""
    yc, h_cds = 0.5, 0.16
    xpad = (region_hi - region_lo) * 0.04
    ax.set_xlim(region_lo - xpad, region_hi + xpad)
    ax.set_ylim(0, 1)
    ax.axis("off")
    # one continuous backbone across the FULL view, so flanks (single-gene) and
    # intergenic gaps (multi-gene) still read as on-scale genome
    ax.plot([region_lo, region_hi], [yc, yc], color=BACKBONE, lw=1.6, zorder=1)
    multi = len(genes) > 1
    for g in genes:
        draw_one_gene(ax, g, region_lo, region_hi, label=multi)
    return yc, h_cds
