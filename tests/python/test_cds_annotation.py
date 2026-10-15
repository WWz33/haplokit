from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit.cli import main


DATA_DIR = ROOT / "data"
GFF = DATA_DIR / "annotation.gff"
VCF = DATA_DIR / "var.sorted.vcf.gz"


def _read_tsv(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _cds_spans() -> dict[str, list[tuple[int, int]]]:
    spans: dict[str, list[tuple[int, int]]] = {}
    for line in GFF.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        transcript = re.search(r"Parent=([^;]+)", fields[8]).group(1)
        spans.setdefault(transcript, []).append((int(fields[3]), int(fields[4])))
    for transcript in spans:
        spans[transcript].sort()
    return spans


def _variant_refs() -> dict[int, str]:
    import gzip

    refs: dict[int, str] = {}
    with gzip.open(VCF, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            refs[int(fields[1])] = fields[3]
    return refs


def _write_cds_fixture(path: Path) -> Path:
    """Synthesize a CDS FASTA consistent with the test GFF and VCF.

    The genome backbone is all 'A'; single-base REF alleles from the VCF are
    written back so the annotation's REF cross-check passes. Two codons are
    forced: ATG at the CDS start (both transcripts start at 4332) and TAC at
    5211-5213 so that the g.5213 C>G variant produces a stop codon.
    """
    genome = ["A"] * 8000
    for pos, ref in _variant_refs().items():
        if len(ref) == 1:
            genome[pos - 1] = ref
    genome[4332 - 1 : 4334] = ["A", "T", "G"]
    genome[5211 - 1 : 5213] = ["T", "A", "C"]
    backbone = "".join(genome)

    with path.open("w", encoding="utf-8") as handle:
        for transcript, spans in _cds_spans().items():
            sequence = "".join(backbone[start - 1 : end] for start, end in spans)
            handle.write(f">{transcript}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")
    return path


@pytest.fixture()
def cds_fasta(tmp_path: Path) -> Path:
    return _write_cds_fixture(tmp_path / "cds.fa")


def test_cds_annotation_writes_positions_and_amino_acid_changes(tmp_path: Path, cds_fasta: Path) -> None:
    out_dir = tmp_path / "out"
    exit_code = main(
        [
            "view",
            str(VCF),
            "-r",
            "scaffold_1:4300-7910",
            "--gff",
            str(GFF),
            "--cds",
            str(cds_fasta),
            "--output-file",
            str(out_dir),
        ]
    )
    assert exit_code == 0

    table = _read_tsv(out_dir / "cds_anno.tsv")
    assert table[0] == ["CHR", "scaffold_1", "Variants: ", "9"]
    assert table[1][:3] == ["POS", "4300", "4345"]
    assert table[2][0] == "CDS_POS"
    assert table[3][0] == "CDS_VAR"

    sites = {row[0]: row for row in (table[1], table[2], table[3])}
    positions = dict(zip(sites["POS"][1:], range(len(sites["POS"]) - 1)))

    def cell(row: str, pos: str) -> str:
        return sites[row][1 + positions[pos]]

    # Outside the CDS and multi-base REF/ALT are left as '-'.
    assert cell("CDS_POS", "4300") == "-"
    assert cell("CDS_POS", "4910") == "-"
    assert cell("CDS_POS", "4950") == "-"
    assert cell("CDS_POS", "5037") == "-"
    assert cell("CDS_POS", "5209") == "-"

    # Both transcripts cover g.4345; one entry per transcript, in ID order.
    assert cell("CDS_POS", "4345") == "test1G0387.1:c.14T>A;test1G0387.2:c.14T>A"
    assert (
        cell("CDS_VAR", "4345")
        == "test1G0387.1:p.Ile5Lys(missense_variant);test1G0387.2:p.Ile5Lys(missense_variant)"
    )

    # Two ALT alleles on one site become two entries separated by ';'.
    assert cell("CDS_POS", "4879") == "test1G0387.1:c.548T>A;test1G0387.1:c.548T>G"
    assert (
        cell("CDS_VAR", "4879")
        == "test1G0387.1:p.Ile183Lys(missense_variant);test1G0387.1:p.Ile183Arg(missense_variant)"
    )

    # Forced TAC codon makes the C>G change a stop_gained.
    assert cell("CDS_POS", "5213") == "test1G0387.1:c.882C>G"
    assert cell("CDS_VAR", "5213") == "test1G0387.1:p.Tyr294Ter(stop_gained)"

    # g.6856 sits in transcript .1 only.
    assert cell("CDS_POS", "6856") == "test1G0387.1:c.2525A>G"
    assert cell("CDS_VAR", "6856") == "test1G0387.1:p.Lys842Arg(missense_variant)"


def test_cds_without_gff_is_rejected(tmp_path: Path, cds_fasta: Path) -> None:
    with pytest.raises(RuntimeError, match="--cds requires --gff3"):
        main(
            [
                "view",
                str(VCF),
                "-r",
                "scaffold_1:4300-7910",
                "--cds",
                str(cds_fasta),
                "--output-file",
                str(tmp_path / "out"),
            ]
        )


def test_cds_annotation_absent_without_flag(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    exit_code = main(
        [
            "view",
            str(VCF),
            "-r",
            "scaffold_1:4300-7910",
            "--gff",
            str(GFF),
            "--output-file",
            str(out_dir),
        ]
    )
    assert exit_code == 0
    assert not (out_dir / "cds_anno.tsv").exists()