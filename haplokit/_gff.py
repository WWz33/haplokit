from __future__ import annotations

import csv
from pathlib import Path


def parse_gff_features(path: str | Path) -> list[dict]:
    """Parse GFF3 file, return list of feature dicts with type, chrom, start, end."""
    features: list[dict] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 9 or row[0].startswith("#"):
                continue
            try:
                start = int(row[3])
                end = int(row[4])
            except ValueError:
                continue
            features.append({
                "seqid": row[0],
                "type": row[2],
                "start": start,
                "end": end,
                "strand": row[6],
            })
    return features


# Priority-ordered feature types for classification
_CLASSIFY_ORDER = [
    "CDS",
    "five_prime_UTR", "3primeUTR", "three_prime_UTR",
    "exon",
    "intron",
    "mRNA", "gene",
]


def classify_positions(
    positions: list[int],
    chrom: str,
    features: list[dict],
) -> list[str]:
    """Classify each variant position by its most specific overlapping GFF feature.

    Returns a list of functional category strings, one per position.
    Categories: CDS, 5'UTR, 3'UTR, exon, intron, intergenic.
    """
    # Filter features on same chromosome
    chrom_features = [f for f in features if f["seqid"] == chrom]

    # Build gene boundaries for intron/intergenic classification
    gene_regions: list[tuple[int, int]] = []
    for f in chrom_features:
        if f["type"] == "gene":
            gene_regions.append((f["start"], f["end"]))

    # Build exon intervals for intron detection
    exon_intervals: list[tuple[int, int]] = []
    for f in chrom_features:
        if f["type"] == "exon":
            exon_intervals.append((f["start"], f["end"]))

    results: list[str] = []
    for pos in positions:
        category = "intergenic"

        # Check if within any gene
        in_gene = any(gstart <= pos <= gend for gstart, gend in gene_regions)
        if in_gene:
            category = "intron"  # default within gene but not in exon

        # Find most specific overlapping feature (priority order)
        for ftype in _CLASSIFY_ORDER:
            for f in chrom_features:
                if f["type"] == ftype and f["start"] <= pos <= f["end"]:
                    # Map to canonical category names
                    if ftype == "CDS":
                        category = "CDS"
                    elif ftype in ("five_prime_UTR",):
                        category = "5'UTR"
                    elif ftype in ("three_prime_UTR", "3primeUTR"):
                        category = "3'UTR"
                    elif ftype == "exon" and category == "intron":
                        category = "exon"
                    # gene/mRNA don't override more specific types
                    break
            if category in ("CDS", "5'UTR", "3'UTR", "exon"):
                break

        results.append(category)

    return results
