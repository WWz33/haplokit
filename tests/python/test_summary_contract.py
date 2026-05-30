from __future__ import annotations

import sys
from types import SimpleNamespace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit.cli import Selector, _tsv_paths_for_selector, _write_selector_result_txt, _write_selector_summary_txt
from haplokit.summary_contract import (
    build_hap_label_map,
    hap_samples,
    hap_states,
    hap_summary_states,
    with_population_breakdown,
)


def test_build_hap_label_map_prefers_contract_id() -> None:
    summary_row = {
        "haplotypes": [
            {"hap": "0/0|1/1", "id": "Hap01"},
            {"hap": "1/1|0/0"},
        ]
    }

    assert build_hap_label_map(summary_row) == {
        "0/0|1/1": "Hap01",
        "1/1|0/0": "Hap02",
    }


def test_hap_states_converts_exact_genotype_states_to_display_alleles() -> None:
    sites = [
        {"allele": "A/T"},
        {"allele": "G/C,A"},
    ]

    assert hap_states("0/0|2/2", sites, "exact") == ["A", "A"]
    assert hap_states("A001", sites, "max-diff") == ["A001"]


def test_hap_summary_states_and_samples_fall_back_to_detail_rows() -> None:
    hap = {"hap": "0/0|1/1"}
    detail_row = {
        "accessions": [
            {"sample": "S1", "hap": "0/0|1/1"},
            {"sample": "S2", "hap": "1/1|0/0"},
        ]
    }
    sites = [{"allele": "A/T"}, {"allele": "G/C"}]

    assert hap_summary_states(hap, sites, "exact") == ["A", "C"]
    assert hap_samples(hap, detail_row) == ["S1"]


def test_summary_contract_uses_backend_states_and_samples_when_present() -> None:
    hap = {
        "hap": "0/0|1/1",
        "states": ["A", "C"],
        "samples": ["S1", "S2"],
    }

    assert hap_summary_states(hap, [], "max-diff") == ["A", "C"]
    assert hap_samples(hap, {"accessions": []}) == ["S1", "S2"]


def test_with_population_breakdown_adds_frequency_labels(tmp_path) -> None:
    pop_file = tmp_path / "pop.tsv"
    pop_file.write_text("S1\tPopA\nS2\tPopA\nS3\tPopB\n", encoding="utf-8")
    haplotypes = [
        {"id": "Hap01", "samples": ["S1", "S3"]},
        {"id": "Hap02", "samples": ["S2"]},
    ]

    enriched = with_population_breakdown(haplotypes, str(pop_file))

    assert enriched[0]["populations"] == [
        {
            "population": "PopA",
            "count": 1,
            "total": 2,
            "frequency": 0.5,
            "frequency_label": "1/2",
        },
        {
            "population": "PopB",
            "count": 1,
            "total": 1,
            "frequency": 1.0,
            "frequency_label": "1/1",
        },
    ]
    assert enriched[1]["populations"][0]["frequency_label"] == "1/2"


def test_single_selector_default_tsv_paths_use_region_slug(tmp_path) -> None:
    selector = Selector(
        payload={"type": "region", "chrom": "scaffold_1", "start": 4300, "end": 5000},
        region="scaffold_1:4300-5000",
    )
    args = SimpleNamespace(output_file=str(tmp_path / "out"))

    summary_path, result_path = _tsv_paths_for_selector(args, selector, 0, 1)

    assert summary_path.name == "hap_summary_scaffold_1_4300_5000.tsv"
    assert result_path.name == "hapresult_scaffold_1_4300_5000.tsv"


def test_single_selector_explicit_prefix_keeps_legacy_tsv_names(tmp_path) -> None:
    selector = Selector(
        payload={"type": "region", "chrom": "scaffold_1", "start": 4300, "end": 5000},
        region="scaffold_1:4300-5000",
    )
    args = SimpleNamespace(output_file=str(tmp_path / "custom.tsv"))

    summary_path, result_path = _tsv_paths_for_selector(args, selector, 0, 1)

    assert summary_path.name == "custom.hap_summary.tsv"
    assert result_path.name == "custom.hapresult.tsv"


def test_population_tsv_writers_keep_accession_tails(tmp_path) -> None:
    selector = Selector(
        payload={"type": "region", "chrom": "chr1", "start": 10, "end": 10},
        region="chr1:10-10",
    )
    summary_row = {
        "annotation": {"mode": "none"},
        "sites": [{"chrom": "chr1", "pos": 10, "allele": "A/T"}],
        "haplotype_count": 1,
        "sample_count": 3,
        "variant_count": 1,
        "grouping_method": "exact",
        "haplotypes": [
            {
                "id": "Hap01",
                "hap": "0/0",
                "states": ["A"],
                "samples": ["S1", "S2"],
                "count": 2,
                "populations": [
                    {
                        "population": "PopA",
                        "frequency_label": "1/2",
                        "samples": ["S1"],
                        "total": 2,
                    },
                    {
                        "population": "PopB",
                        "frequency_label": "1/1",
                        "samples": ["S2"],
                        "total": 1,
                    },
                ],
            }
        ],
    }
    detail_row = {
        "accessions": [
            {"hap": "0/0", "sample": "S1", "population": "PopA"},
            {"hap": "0/0", "sample": "S2", "population": "PopB"},
        ]
    }

    summary_path = tmp_path / "hap_summary.tsv"
    result_path = tmp_path / "hapresult.tsv"
    _write_selector_summary_txt(selector, summary_row, detail_row, summary_path)
    _write_selector_result_txt(selector, summary_row, detail_row, result_path)

    summary_rows = [line.split("\t") for line in summary_path.read_text(encoding="utf-8").splitlines()]
    assert summary_rows[3] == ["ALLELE", "A/T", "PopA_n", "PopA_Accession", "PopB_n", "PopB_Accession", "Accession", "freq"]
    assert summary_rows[4][-2:] == ["S1;S2", "2"]
    assert summary_rows[4][2:6] == ["1/2", "S1", "1/1", "S2"]

    result_rows = [line.split("\t") for line in result_path.read_text(encoding="utf-8").splitlines()]
    assert result_rows[3][-2:] == ["Population", "Accession"]
    assert result_rows[4][-2:] == ["PopA", "S1"]
