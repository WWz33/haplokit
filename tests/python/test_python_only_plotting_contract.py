from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit.plot import plot_hap_distribution, plot_hap_table, read_hap_summary_tsv
from haplokit._transform import transform_for_display


def test_haplokit_python_surface_does_not_reference_r_runtime() -> None:
    forbidden_tokens = ("Rscript", "rpy2", "plotHapTable(")
    for file in (ROOT / "haplokit").glob("*.py"):
        content = file.read_text(encoding="utf-8")
        for token in forbidden_tokens:
            assert token not in content, f"{file} unexpectedly references {token}"


def test_python_plotter_renders_pdf_and_svg_from_hap_summary_tsv(tmp_path: Path) -> None:
    tsv = tmp_path / "hap_summary.tsv"
    tsv.write_text(
        "\n".join(
            [
                "CHR\tscaffold_1\tscaffold_1\tHaplotypes: \t2",
                "POS\t4300\t4345\tIndividuals: \t2",
                "INFO\t.\t.\tVariants: \t2",
                "ALLELE\tA/T\tG/C\tAccession\tfreq",
                "H001\tA\tC\tS1\t1",
                "H002\tT\tG\tS2\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    rows = read_hap_summary_tsv(tsv)
    pdf_out = tmp_path / "plot.pdf"
    svg_out = tmp_path / "plot.svg"
    pdf_rendered = plot_hap_table(rows, pdf_out, figsize=(6, 4))
    svg_rendered = plot_hap_table(rows, svg_out, figsize=(6, 4))
    assert pdf_rendered.exists()
    assert pdf_rendered.suffix == ".pdf"
    assert svg_rendered.exists()
    assert svg_rendered.suffix == ".svg"


def test_transform_ignores_persisted_population_columns_for_display() -> None:
    rows = [
        ["CHR", "scaffold_1", "scaffold_1", "Haplotypes: ", "2"],
        ["POS", "4300", "4345", "Individuals: ", "3"],
        ["INFO", ".", ".", "Variants: ", "2"],
        ["ALLELE", "A/T", "G/C", "PopA_n", "PopA_Accession", "PopB_n", "PopB_Accession", "Accession", "freq"],
        ["Hap01", "A", "C", "999/999", "ignored", "999/999", "ignored", "S1;S3", "2"],
        ["Hap02", "T", "G", "999/999", "ignored", "999/999", "ignored", "S2", "1"],
    ]

    transformed, title = transform_for_display(rows, {"S1": "PopA", "S2": "PopA", "S3": "PopB"})

    assert title == "scaffold_1:4300-4345"
    assert transformed == [
        ["POS", "4300", "4345"],
        ["ALLELE", "A/T", "G/C", "PopA", "PopB", "n/N"],
        ["Hap01", "A", "C", "1/2", "1/1", "2/3"],
        ["Hap02", "T", "G", "1/2", "0/1", "1/3"],
    ]


def test_distribution_plot_accepts_map_facecolor(tmp_path: Path) -> None:
    samples = [
        {"lon": -99.13, "lat": 19.43, "hap": "Hap01"},
        {"lon": 116.40, "lat": 39.90, "hap": "Hap02"},
        {"lon": 2.35, "lat": 48.86, "hap": "Hap01"},
    ]
    out = tmp_path / "global_hap_map.pdf"

    rendered = plot_hap_distribution(
        samples,
        ["Hap01", "Hap02"],
        out,
        database="none",
        map_facecolor="#e6f2ff",
        fmt="pdf",
        figsize=(6, 4),
        show_labels=True,
    )

    assert rendered.exists()
    assert rendered.suffix == ".pdf"
