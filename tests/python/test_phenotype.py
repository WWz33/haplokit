from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit._phenotype import (
    load_phenotype_dataset,
    pairwise_statistics,
    read_population_groups,
    read_sample_haplotypes,
    summarize_groups,
)
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


def test_phenotype_cli_stat_and_box_write_outputs(tmp_path: Path) -> None:
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
            "stat",
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
    assert args.phenotype_command == "stat"
    assert args.phenotypes == str(phenotype)

    alias_args = parser.parse_args(
        [
            "phenotype",
            "stat",
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
            "stat",
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
            "box",
            "--hapresult",
            str(hapresult),
            "--phenotypes",
            str(phenotype),
            "--trait",
            "yield",
            "--min-hap-size",
            "3",
            "--output",
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
                "box",
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
            "stat",
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
    assert stat_args.phenotype_command == "stat"
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

    box_args = parser.parse_args(
        [
            "phenotype",
            "box",
            "-H",
            "hapresult.tsv",
            "-P",
            "phenotype.csv",
            "-t",
            "yield",
            "-o",
            "box.svg",
            "-F",
            "svg",
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
            "-T",
            "Yield by haplotype",
        ]
    )
    assert box_args.phenotype_command == "box"
    assert box_args.comparison == [("Hap01", "Hap02")]
    assert box_args.title == "Yield by haplotype"

    with pytest.raises(SystemExit):
        parser.parse_args(["phenotype", "stat", "--help"])
    out = capsys.readouterr().out
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
    )

    assert rendered.exists()
    assert rendered.suffix == ".svg"


def test_phenotype_cli_stat_accepts_population_group(tmp_path: Path) -> None:
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
            "stat",
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
