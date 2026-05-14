"""Publication-quality haplotype visualization (Nature-style).

Public API re-exported from submodules.
"""

from ._transform import read_hap_summary_tsv, read_popgroup
from ._table import plot_hap_table
from ._distribution import plot_hap_distribution
from ._network import plot_hap_network
from ._gff import parse_gff_features, classify_positions
from ._palette import PALETTE, FUNC_COLORS

__all__ = [
    "read_hap_summary_tsv",
    "read_popgroup",
    "plot_hap_table",
    "plot_hap_distribution",
    "plot_hap_network",
    "parse_gff_features",
    "classify_positions",
    "PALETTE",
    "FUNC_COLORS",
]
