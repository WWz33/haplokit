from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit._phenotype import (
    PhenotypeRecord,
    load_phenotype_dataset,
    pairwise_statistics,
    population_pairwise_statistics,
    read_population_groups,
    read_sample_haplotypes,
    summarize_groups,
)
import haplokit._phenotype_plot as phenotype_plot
from haplokit.cli import build_parser, main
from haplokit.plot import plot_hap_phenotype_box


def _write_hapresult(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "CHR\tchr1\tHaplotypes: \t3",
                "POS\t10\tIndividuals: \t7",
                "INFO\t.\tVariants: \t1",
                "ALLELE\tA/G\tAccession",
                "Hap01\tA\tS1",
                "Hap01\tA\tS2",
                "Hap01\tA\tS3",
                "Hap02\tG\tS4",
                "Hap02\tG\tS5",
                "Hap02\tG\tS6",
                "Hap03\tG\tS7",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_phenotype(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "sample,yield,height",
                "S1,1,10",
                "S2,2,11",
                "S3,3,12",
                "S4,7,20",
                "S5,8,21",
                "S6,9,22",
                "S7,100,30",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_population(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "sample\tpopulation",
                "S1\tPopA",
                "S2\tPopA",
                "S3\tPopB",
                "S4\tPopA",
                "S5\tPopA",
                "S6\tPopB",
                "S7\tPopB",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_balanced_population_box_inputs(haplotypes: Path, phenotype: Path, population: Path) -> None:
    haplotypes.write_text(
        "\n".join(
            [
                "sample\thaplotype",
                "A1\tHap01",
                "A2\tHap01",
                "A3\tHap02",
                "A4\tHap02",
                "B1\tHap01",
                "B2\tHap01",
                "B3\tHap02",
                "B4\tHap02",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    phenotype.write_text(
        "\n".join(
            [
                "sample,yield",
                "A1,1.0",
                "A2,1.2",
                "A3,2.0",
                "A4,2.2",
                "B1,3.0",
                "B2,3.2",
                "B3,4.0",
                "B4,4.2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    population.write_text(
        "\n".join(
            [
                "sample\tpopulation",
                "A1\tPopA",
                "A2\tPopA",
                "A3\tPopA",
                "A4\tPopA",
                "B1\tPopB",
                "B2\tPopB",
                "B3\tPopB",
                "B4\tPopB",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_read_sample_haplotypes_accepts_hapresult_and_simple_table(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    _write_hapresult(hapresult)
    records = read_sample_haplotypes(hapresult)
    assert records[0].sample == "S1"
    assert records[0].haplotype == "Hap01"
    assert len(records) == 7

    simple = tmp_path / "sample_hap.csv"
    simple.write_text("samples,haplotypes\nA,H1\nB,H2\n", encoding="utf-8")
    simple_records = read_sample_haplotypes(simple)
    assert [(item.sample, item.haplotype) for item in simple_records] == [("A", "H1"), ("B", "H2")]

    summary = tmp_path / "hap_summary.tsv"
    summary.write_text(
        "\n".join(
            [
                "CHR\tchr1\tHaplotypes: \t2",
                "POS\t10\tIndividuals: \t3",
                "INFO\t.\tVariants: \t1",
                "ALLELE\tA/G\tAccession\tfreq",
                "Hap01\tA\tS1;S2\t2",
                "Hap02\tG\tS3\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    summary_records = read_sample_haplotypes(summary)
    assert [(item.sample, item.haplotype) for item in summary_records] == [
        ("S1", "Hap01"),
        ("S2", "Hap01"),
        ("S3", "Hap02"),
    ]


def test_phenotype_statistics_filters_small_haplotypes_and_writes_pairs(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    _write_hapresult(hapresult)
    _write_phenotype(phenotype)

    dataset = load_phenotype_dataset(hapresult, phenotype, traits=["yield"])
    rows = pairwise_statistics(dataset.records, traits=dataset.traits, min_hap_size=3, method="welch")
    summary = summarize_groups(dataset.records, traits=dataset.traits, min_hap_size=3)

    assert dataset.matched_sample_count == 7
    assert len(rows) == 1
    assert rows[0]["group1"] == "Hap01"
    assert rows[0]["group2"] == "Hap02"
    assert rows[0]["count1"] == 3
    assert rows[0]["count2"] == 3
    assert rows[0]["mean1"] == 2.0
    assert rows[0]["mean2"] == 8.0
    assert {item["haplotype"] for item in summary} == {"Hap01", "Hap02"}


def test_phenotype_statistics_ignores_missing_values_and_reports_effective_n(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    _write_hapresult(hapresult)
    phenotype.write_text(
        "\n".join(
            [
                "sample,yield",
                "S1,1",
                "S2,NA",
                "S3,3",
                "S4,7",
                "S5,8",
                "S6,9",
                "S7,.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    dataset = load_phenotype_dataset(hapresult, phenotype, traits=["yield"])
    rows = pairwise_statistics(dataset.records, traits=dataset.traits, min_hap_size=2, method="welch")
    summary = summarize_groups(dataset.records, traits=dataset.traits, min_hap_size=2)

    assert dataset.matched_sample_count == 7
    assert len(dataset.records) == 5
    assert rows[0]["count1"] == 2
    assert rows[0]["count2"] == 3
    assert rows[0]["effective_n"] == 5
    assert {item["effective_n"] for item in summary} == {5}


def test_phenotype_cli_reports_effective_n_for_missing_values(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    stats_out = tmp_path / "stats.tsv"
    summary_out = tmp_path / "summary.tsv"
    _write_hapresult(hapresult)
    phenotype.write_text(
        "\n".join(
            [
                "sample,yield",
                "S1,1",
                "S2,NA",
                "S3,3",
                "S4,7",
                "S5,8",
                "S6,9",
                "S7,.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "phenotype",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--min-hap-size",
            "2",
            "--output",
            str(stats_out),
            "--summary-output",
            str(summary_out),
        ]
    )

    assert exit_code == 0
    stats_lines = stats_out.read_text(encoding="utf-8").splitlines()
    stat_header = stats_lines[0].split("\t")
    stat_values = stats_lines[1].split("\t")
    summary_lines = summary_out.read_text(encoding="utf-8").splitlines()
    summary_header = summary_lines[0].split("\t")
    summary_values = summary_lines[1].split("\t")

    assert stat_header[-1] == "effective_n"
    assert stat_values[-1] == "5"
    assert summary_header[-1] == "effective_n"
    assert summary_values[-1] == "5"


def test_phenotype_statistics_stratify_by_population_group(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    population = tmp_path / "popgroup.tsv"
    _write_hapresult(hapresult)
    _write_phenotype(phenotype)
    _write_population(population)

    pop_map = read_population_groups(population)
    dataset = load_phenotype_dataset(hapresult, phenotype, traits=["yield"], population_file=population)
    rows = pairwise_statistics(
        dataset.records,
        traits=dataset.traits,
        min_hap_size=2,
        method="welch",
        populations=dataset.populations,
    )
    summary = summarize_groups(dataset.records, traits=dataset.traits, min_hap_size=2, populations=dataset.populations)

    assert pop_map["S1"] == "PopA"
    assert dataset.populations == ("PopA", "PopB")
    assert len(rows) == 1
    assert rows[0]["population"] == "PopA"
    assert rows[0]["group1"] == "Hap01"
    assert rows[0]["group2"] == "Hap02"
    assert {item["population"] for item in summary} == {"PopA"}


def test_phenotype_population_groups_include_unknown_for_missing_samples(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    population = tmp_path / "partial_popgroup.tsv"
    _write_hapresult(hapresult)
    _write_phenotype(phenotype)
    population.write_text(
        "\n".join(
            [
                "sample\tpopulation",
                "S1\tPopA",
                "S2\tPopA",
                "S4\tPopA",
                "S5\tPopA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    dataset = load_phenotype_dataset(hapresult, phenotype, traits=["yield"], population_file=population)
    summary = summarize_groups(dataset.records, traits=dataset.traits, min_hap_size=1, populations=dataset.populations)

    assert dataset.populations == ("PopA", "Unknown")
    assert {item["population"] for item in summary} == {"PopA", "Unknown"}


def test_population_pairwise_statistics_compares_populations_within_haplotype() -> None:
    records = (
        PhenotypeRecord("A1", "Hap01", "yield", 1.0, "PopA"),
        PhenotypeRecord("A2", "Hap01", "yield", 1.2, "PopA"),
        PhenotypeRecord("B1", "Hap01", "yield", 3.0, "PopB"),
        PhenotypeRecord("B2", "Hap01", "yield", 3.2, "PopB"),
        PhenotypeRecord("A3", "Hap02", "yield", 5.0, "PopA"),
    )

    rows = population_pairwise_statistics(
        records,
        traits=["yield"],
        haplotypes=["Hap01"],
        populations=["PopA", "PopB"],
        min_hap_size=2,
    )

    assert len(rows) == 1
    assert rows[0]["haplotype"] == "Hap01"
    assert rows[0]["group1"] == "PopA"
    assert rows[0]["group2"] == "PopB"
    assert rows[0]["count1"] == 2
    assert rows[0]["count2"] == 2
    assert rows[0]["effective_n"] == 4


def test_phenotype_cli_stats_and_plot_box_write_outputs(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    stats_out = tmp_path / "stats.tsv"
    summary_out = tmp_path / "summary.tsv"
    box_out = tmp_path / "yield_box.pdf"
    _write_hapresult(hapresult)
    _write_phenotype(phenotype)

    parser = build_parser()
    args = parser.parse_args(
        [
            "phenotype",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--min-hap-size",
            "3",
            "--output",
            str(stats_out),
        ]
    )
    assert args.command == "phenotype"
    assert args.plot_box is False
    assert args.phenotypes == str(phenotype)

    alias_args = parser.parse_args(
        [
            "phenotype",
            "--hapresult",
            str(hapresult),
            "--phenotype",
            str(phenotype),
        ]
    )
    assert alias_args.phenotypes == str(phenotype)

    exit_code = main(
        [
            "pheno",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--min-hap-size",
            "3",
            "--output",
            str(stats_out),
            "--summary-output",
            str(summary_out),
        ]
    )
    assert exit_code == 0
    assert stats_out.exists()
    assert "Hap01\tHap02" in stats_out.read_text(encoding="utf-8")
    assert summary_out.exists()

    box_exit = main(
        [
            "phenotype",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--min-hap-size",
            "3",
            "--output",
            str(stats_out),
            "--plot-box",
            "--box-output",
            str(box_out),
            "--plot-format",
            "pdf",
            "--comparison",
            "Hap01,Hap02",
        ]
    )
    assert box_exit == 0
    assert box_out.exists()
    assert box_out.suffix == ".pdf"


def test_phenotype_cli_rejects_invalid_comparison() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "phenotype",
                "--hapresult",
                "hapresult.tsv",
                "--phenotypes",
                "phenotype.csv",
                "--trait",
                "yield",
                "--comparison",
                "Hap01",
            ]
        )


def test_phenotype_cli_accepts_short_options(capsys: pytest.CaptureFixture[str]) -> None:
    parser = build_parser()
    stat_args = parser.parse_args(
        [
            "pheno",
            "-H",
            "hapresult.tsv",
            "-P",
            "phenotype.csv",
            "-t",
            "yield",
            "-o",
            "stats.tsv",
            "-s",
            "summary.tsv",
            "-m",
            "3",
            "-M",
            "student",
            "-a",
            "none",
            "-p",
            "popgroup.tsv",
            "-d",
            "tab",
            "-D",
            "comma",
            "-G",
            "tab",
        ]
    )
    assert stat_args.command == "pheno"
    assert stat_args.plot_box is False
    assert stat_args.haplotypes == "hapresult.tsv"
    assert stat_args.phenotypes == "phenotype.csv"
    assert stat_args.summary_output == "summary.tsv"
    assert stat_args.min_hap_size == 3
    assert stat_args.method == "student"
    assert stat_args.adjust == "none"
    assert stat_args.population_file == "popgroup.tsv"
    assert stat_args.hap_delimiter == "tab"
    assert stat_args.phenotype_delimiter == "comma"
    assert stat_args.population_delimiter == "tab"
    assert stat_args.comparison is None

    box_args = parser.parse_args(
        [
            "phenotype",
            "-H",
            "hapresult.tsv",
            "-P",
            "phenotype.csv",
            "-t",
            "yield",
            "-p",
            "popgroup.tsv",
            "-B",
            "-o",
            "stats.tsv",
            "-b",
            "box.svg",
            "-F",
            "svg",
            "-z",
            "6x4",
            "-m",
            "3",
            "-M",
            "mannwhitney",
            "-c",
            "Hap01,Hap02",
            "-d",
            "tab",
            "-D",
            "comma",
            "-G",
            "tab",
            "-T",
            "Yield by haplotype",
        ]
    )
    assert box_args.plot_box is True
    assert box_args.box_output == "box.svg"
    assert box_args.population_file == "popgroup.tsv"
    assert box_args.population_delimiter == "tab"
    assert box_args.figsize == (6.0, 4.0)
    assert box_args.comparison == [("Hap01", "Hap02")]
    assert box_args.title == "Yield by haplotype"

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "phenotype",
                "-H",
                "hapresult.tsv",
                "-P",
                "phenotype.csv",
                "-t",
                "yield",
                "--figsize",
                "6,0",
            ]
        )

    with pytest.raises(SystemExit):
        parser.parse_args(["phenotype", "--help"])
    out = capsys.readouterr().out
    assert "Haplotype/phenotype input options" in out
    assert "Phenotype test options" in out
    assert "phenotype table; first column is sample ID" in out
    assert "-m MIN_HAP_SIZE" in out
    assert "--min-hap-size MIN_HAP_SIZE" in out


def test_plot_hap_phenotype_box_is_exported(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    _write_hapresult(hapresult)
    _write_phenotype(phenotype)

    rendered = plot_hap_phenotype_box(
        hapresult,
        phenotype,
        "yield",
        tmp_path / "yield.svg",
        min_hap_size=3,
        fmt="svg",
        figsize=(6, 4),
    )

    assert rendered.exists()
    assert rendered.suffix == ".svg"


def test_plot_hap_phenotype_box_groups_population_boxes(tmp_path: Path) -> None:
    hapresult = tmp_path / "haplotypes.tsv"
    phenotype = tmp_path / "phenotype.csv"
    population = tmp_path / "popgroup.tsv"
    _write_balanced_population_box_inputs(hapresult, phenotype, population)

    rendered = plot_hap_phenotype_box(
        hapresult,
        phenotype,
        "yield",
        tmp_path / "yield_by_population.svg",
        min_hap_size=2,
        population_file=population,
        fmt="svg",
    )

    assert rendered.exists()
    assert rendered.suffix == ".svg"
    rendered_text = rendered.read_text(encoding="utf-8")
    assert "PopA" in rendered_text
    assert "PopB" in rendered_text
    assert "Hap01" in rendered_text
    assert "Hap02" in rendered_text
    assert "*" in rendered_text


def test_plot_hap_phenotype_box_filters_records_per_population_panel(tmp_path: Path) -> None:
    records = (
        PhenotypeRecord("A1", "Hap01", "yield", 1.0, "PopA"),
        PhenotypeRecord("A2", "Hap01", "yield", 1.2, "PopA"),
        PhenotypeRecord("A3", "Hap02", "yield", 2.0, "PopA"),
        PhenotypeRecord("A4", "Hap02", "yield", 2.2, "PopA"),
        PhenotypeRecord("B1", "Hap03", "yield", 3.0, "PopB"),
        PhenotypeRecord("B2", "Hap03", "yield", 3.2, "PopB"),
        PhenotypeRecord("B3", "Hap04", "yield", 4.0, "PopB"),
        PhenotypeRecord("B4", "Hap04", "yield", 4.2, "PopB"),
    )

    rendered = plot_hap_phenotype_box(
        records,
        trait="yield",
        output_path=tmp_path / "population_filtered.svg",
        min_hap_size=2,
        fmt="svg",
    )

    rendered_text = rendered.read_text(encoding="utf-8")
    assert "PopA" in rendered_text
    assert "PopB" in rendered_text
    assert rendered_text.count("Hap01") == 1
    assert rendered_text.count("Hap02") == 1
    assert rendered_text.count("Hap03") == 1
    assert rendered_text.count("Hap04") == 1


def test_haplotype_box_emits_default_pairwise_annotations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (
        PhenotypeRecord("S1", "Hap01", "yield", 1.0),
        PhenotypeRecord("S2", "Hap01", "yield", 1.2),
        PhenotypeRecord("S3", "Hap02", "yield", 2.0),
        PhenotypeRecord("S4", "Hap02", "yield", 2.2),
        PhenotypeRecord("S5", "Hap03", "yield", 3.0),
        PhenotypeRecord("S6", "Hap03", "yield", 3.2),
    )
    captured: list[tuple[float, float, str]] = []

    def capture_annotations(_ax, annotations):
        captured.extend(annotations)

    monkeypatch.setattr(phenotype_plot, "_draw_stat_annotations_inside", capture_annotations)

    plot_hap_phenotype_box(
        records,
        trait="yield",
        output_path=tmp_path / "hap_pairwise.svg",
        min_hap_size=2,
        fmt="svg",
    )

    distances = sorted(round(abs(x2 - x1), 2) for x1, x2, _ in captured)
    assert len(captured) == 3
    assert distances == [1.0, 1.0, 2.0]


def test_haplotype_box_explicit_comparisons_filter_annotations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (
        PhenotypeRecord("S1", "Hap01", "yield", 1.0),
        PhenotypeRecord("S2", "Hap01", "yield", 1.2),
        PhenotypeRecord("S3", "Hap02", "yield", 2.0),
        PhenotypeRecord("S4", "Hap02", "yield", 2.2),
        PhenotypeRecord("S5", "Hap03", "yield", 3.0),
        PhenotypeRecord("S6", "Hap03", "yield", 3.2),
    )
    captured: list[tuple[float, float, str]] = []

    def capture_annotations(_ax, annotations):
        captured.extend(annotations)

    monkeypatch.setattr(phenotype_plot, "_draw_stat_annotations_inside", capture_annotations)

    plot_hap_phenotype_box(
        records,
        trait="yield",
        output_path=tmp_path / "hap_filtered.svg",
        min_hap_size=2,
        comparisons=[("Hap01", "Hap03")],
        fmt="svg",
    )

    assert len(captured) == 1
    assert round(abs(captured[0][1] - captured[0][0]), 2) == 2.0


def test_population_grouped_box_emits_within_and_between_annotations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (
        PhenotypeRecord("A1", "Hap01", "yield", 1.0, "PopA"),
        PhenotypeRecord("A2", "Hap01", "yield", 1.2, "PopA"),
        PhenotypeRecord("A3", "Hap02", "yield", 2.0, "PopA"),
        PhenotypeRecord("A4", "Hap02", "yield", 2.2, "PopA"),
        PhenotypeRecord("B1", "Hap01", "yield", 3.0, "PopB"),
        PhenotypeRecord("B2", "Hap01", "yield", 3.2, "PopB"),
        PhenotypeRecord("B3", "Hap02", "yield", 4.0, "PopB"),
        PhenotypeRecord("B4", "Hap02", "yield", 4.2, "PopB"),
    )
    captured: list[tuple[float, float, str]] = []

    def capture_annotations(_ax, annotations):
        captured.extend(annotations)

    monkeypatch.setattr(phenotype_plot, "_draw_stat_annotations_inside", capture_annotations)

    plot_hap_phenotype_box(
        records,
        trait="yield",
        output_path=tmp_path / "grouped.svg",
        min_hap_size=2,
        fmt="svg",
    )

    distances = sorted(round(abs(x2 - x1), 2) for x1, x2, _ in captured)
    assert len(captured) == 4
    assert distances[:2] == [0.36, 0.36]
    assert distances[2:] == [1.0, 1.0]


def test_single_haplotype_population_box_emits_only_between_annotations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (
        PhenotypeRecord("A1", "Hap01", "yield", 1.0, "PopA"),
        PhenotypeRecord("A2", "Hap01", "yield", 1.2, "PopA"),
        PhenotypeRecord("B1", "Hap01", "yield", 3.0, "PopB"),
        PhenotypeRecord("B2", "Hap01", "yield", 3.2, "PopB"),
    )
    captured: list[tuple[float, float, str]] = []

    def capture_annotations(_ax, annotations):
        captured.extend(annotations)

    monkeypatch.setattr(phenotype_plot, "_draw_stat_annotations_inside", capture_annotations)

    plot_hap_phenotype_box(
        records,
        trait="yield",
        output_path=tmp_path / "single_hap.svg",
        min_hap_size=2,
        fmt="svg",
    )

    assert len(captured) == 1
    assert round(abs(captured[0][1] - captured[0][0]), 2) == 1.0


def test_phenotype_cli_plot_box_accepts_population_group(tmp_path: Path) -> None:
    hapresult = tmp_path / "haplotypes.tsv"
    phenotype = tmp_path / "phenotype.csv"
    population = tmp_path / "popgroup.tsv"
    stats_out = tmp_path / "box_stats.tsv"
    box_out = tmp_path / "yield_by_population.svg"
    _write_balanced_population_box_inputs(hapresult, phenotype, population)

    exit_code = main(
        [
            "phenotype",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--population",
            str(population),
            "--min-hap-size",
            "2",
            "--output",
            str(stats_out),
            "--plot-box",
            "--box-output",
            str(box_out),
            "--plot-format",
            "svg",
        ]
    )

    assert exit_code == 0
    assert stats_out.exists()
    assert box_out.exists()
    rendered_text = box_out.read_text(encoding="utf-8")
    assert "PopA" in rendered_text
    assert "PopB" in rendered_text
    assert "*" in rendered_text


def test_phenotype_cli_stats_accept_population_group(tmp_path: Path) -> None:
    hapresult = tmp_path / "hapresult.tsv"
    phenotype = tmp_path / "phenotype.csv"
    population = tmp_path / "popgroup.tsv"
    stats_out = tmp_path / "stratified_stats.tsv"
    _write_hapresult(hapresult)
    _write_phenotype(phenotype)
    _write_population(population)

    exit_code = main(
        [
            "phenotype",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--population",
            str(population),
            "--min-hap-size",
            "2",
            "--output",
            str(stats_out),
        ]
    )

    assert exit_code == 0
    lines = stats_out.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[:4] == ["trait", "population", "group1", "group2"]
    assert lines[1].split("\t")[:4] == ["yield", "PopA", "Hap01", "Hap02"]
