from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit._backend import find_haplokit_cpp
from haplokit._transform import transform_for_display
from haplokit.cli import main


DATA_DIR = ROOT / "data"


def _cpp_bin() -> Path:
    return find_haplokit_cpp(repo_root=ROOT, package_dir=ROOT / "haplokit", auto_build=False)


def _network_bin() -> Path:
    return ROOT / "build-wsl" / "haplokit_network_backend"


def _run_view(extra: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(_cpp_bin()), "view-json", str(DATA_DIR / "var.sorted.vcf.gz"), *extra],
        capture_output=True,
        text=True,
        check=False,
    )


def _run_network(algorithm: str, payload: dict) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(_network_bin()), algorithm],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.fixture()
def network_available() -> Path:
    backend = _network_bin()
    if not backend.exists():
        pytest.skip(f"network backend not built at {backend}")
    return backend


# ---------------------------------------------------------------------------
# Region parsing: inverted and out-of-range coordinates are errors, not empty
# results.
# ---------------------------------------------------------------------------

def test_inverted_region_is_rejected() -> None:
    completed = _run_view(["scaffold_1:5209-4300", "--by", "region", "--output", "summary"])
    assert completed.returncode != 0
    assert "start must be <= end" in completed.stderr
    assert completed.stdout.strip() == ""


def test_non_numeric_region_coordinate_is_rejected() -> None:
    completed = _run_view(["scaffold_1:abc", "--by", "region", "--output", "summary"])
    assert completed.returncode != 0
    assert "not a number" in completed.stderr


def test_zero_based_region_start_is_rejected() -> None:
    completed = _run_view(["scaffold_1:0-100", "--by", "region", "--output", "summary"])
    assert completed.returncode != 0
    assert "start must be >= 1" in completed.stderr


def test_valid_region_still_succeeds() -> None:
    completed = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary"])
    assert completed.returncode == 0
    payload = json.loads(completed.stdout)
    assert payload["variant_count"] == 7
    assert payload["haplotype_count"] == 8


# ---------------------------------------------------------------------------
# Missing index: a clear error instead of SIGSEGV (rc 139 / signal 11).
# ---------------------------------------------------------------------------

def test_missing_vcf_index_reports_error_not_crash(tmp_path: Path) -> None:
    unindexed = tmp_path / "unindexed.vcf.gz"
    unindexed.write_bytes((DATA_DIR / "var.sorted.vcf.gz").read_bytes())

    completed = subprocess.run(
        [str(_cpp_bin()), "view-json", str(unindexed), "scaffold_1:4300",
         "--by", "site", "--output", "summary"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode >= 0, f"process died from signal {-completed.returncode}"
    assert completed.returncode != 0
    assert "failed to add indexed VCF/BCF reader" in completed.stderr


# ---------------------------------------------------------------------------
# Frequency denominators count only samples that received a haplotype.
# ---------------------------------------------------------------------------

def test_frequency_denominator_matches_analyzed_samples() -> None:
    completed = _run_view(["scaffold_1:4300-5000", "--by", "region", "--output", "summary"])
    assert completed.returncode == 0
    payload = json.loads(completed.stdout)

    analyzed = payload["analyzed_sample_count"]
    assert 0 < analyzed <= payload["sample_count"]
    assert sum(hap["count"] for hap in payload["haplotypes"]) == analyzed
    for hap in payload["haplotypes"]:
        assert hap["total"] == analyzed
        assert hap["frequency_label"] == f"{hap['count']}/{analyzed}"


def test_population_totals_exclude_samples_outside_the_selection(tmp_path: Path) -> None:
    samples = tmp_path / "samples.list"
    samples.write_text("C1\nC2\nC3\n", encoding="utf-8")

    completed = _run_view([
        "scaffold_1:4300-5209", "--by", "region", "--output", "summary",
        "--samples-file", str(samples),
        "--population-file", str(DATA_DIR / "popgroup.txt"),
    ])
    assert completed.returncode == 0
    payload = json.loads(completed.stdout)

    assert payload["sample_count"] == 3
    first = payload["haplotypes"][0]
    assert first["populations"], "expected population breakdown rows"
    for population in first["populations"]:
        count, total = population["frequency_label"].split("/")
        assert int(total) <= 3, (
            f"population {population['population']} total {total} exceeds the 3 selected samples"
        )
        assert int(count) <= int(total)


# ---------------------------------------------------------------------------
# --max-diff: per-site states, a real pattern, and 0 means exact grouping.
# ---------------------------------------------------------------------------

def test_max_diff_keeps_per_site_states_and_real_pattern() -> None:
    completed = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary",
                           "--max-diff", "0.5"])
    assert completed.returncode == 0
    payload = json.loads(completed.stdout)

    assert payload["grouping_mode"] == "approx-region"
    site_count = len(payload["sites"])
    assert site_count == 7
    for hap in payload["haplotypes"]:
        assert len(hap["states"]) == site_count, (
            f"{hap['id']} has {len(hap['states'])} state cells for {site_count} sites"
        )
        assert hap["hap"].startswith("A"), "cluster key should stay the A### label"
        assert hap["pattern"].count("|") == site_count - 1, (
            f"{hap['id']} pattern {hap['pattern']!r} is not a per-site pattern"
        )


def test_max_diff_zero_equals_exact_grouping() -> None:
    exact = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary"])
    zero = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary",
                      "--max-diff", "0"])
    assert exact.returncode == 0 and zero.returncode == 0

    exact_payload = json.loads(exact.stdout)
    zero_payload = json.loads(zero.stdout)

    assert zero_payload["grouping_mode"] == "strict-region"
    assert zero_payload["grouping_method"] == "exact"
    assert [h["hap"] for h in zero_payload["haplotypes"]] == \
           [h["hap"] for h in exact_payload["haplotypes"]]
    assert [h["count"] for h in zero_payload["haplotypes"]] == \
           [h["count"] for h in exact_payload["haplotypes"]]


def test_max_diff_summary_tsv_rows_match_header_width(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    exit_code = main([
        "view", str(DATA_DIR / "var.sorted.vcf.gz"),
        "-r", "scaffold_1:4300-5209", "-m", "0.5",
        "--output-file", str(out_dir),
    ])
    assert exit_code == 0

    summary = out_dir / "hap_summary_scaffold_1_4300_5209.tsv"
    rows = [line.split("\t") for line in summary.read_text(encoding="utf-8").splitlines()]
    header = rows[3]
    for row in rows[4:]:
        assert len(row) == len(header), (
            f"row {row[0]} has {len(row)} fields against a {len(header)}-field header"
        )


# ---------------------------------------------------------------------------
# Table rendering must not depend on the haplotype prefix.
# ---------------------------------------------------------------------------

def test_table_transform_handles_custom_hap_prefix() -> None:
    rows = [
        ["CHR", "chr1", "Haplotypes: ", "1"],
        ["POS", "100", "Individuals: ", "2"],
        ["INFO", ".", "Variants: ", "1"],
        ["ALLELE", "A/T", "Accession", "freq"],
        ["A01", "A", "S1;S2", "2"],
    ]

    transformed, _title = transform_for_display(rows, None)

    assert transformed[1] == ["ALLELE", "A/T", "n/N"]
    assert transformed[2] == ["A01", "A", "2/2"], "sample list leaked into the count column"


def test_table_transform_abbreviates_indels_with_custom_prefix() -> None:
    from haplokit._transform import _indel_footnotes

    rows = [
        ["CHR", "chr1", "Haplotypes: ", "1"],
        ["POS", "100", "Individuals: ", "2"],
        ["INFO", ".", "Variants: ", "1"],
        ["ALLELE", "GCCTA/A", "Accession", "freq"],
        ["A01", "GCCTA", "S1;S2", "2"],
    ]

    # var_end marks the first non-allele column (Accession at index 2).
    abbreviated, footnote = _indel_footnotes(rows, 2)

    assert abbreviated[4][1] == "i1", "long indel was not abbreviated for a custom prefix"
    assert abbreviated[3][1] == "i1/A", "ALLELE row must be abbreviated as well"
    assert "i1" in footnote


# ---------------------------------------------------------------------------
# Network: per-site distance, MJN distinct from MSN, and clean failure on
# malformed input.
# ---------------------------------------------------------------------------

def test_network_distance_counts_sites_not_characters(network_available: Path) -> None:
    completed = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary"])
    payload = json.loads(completed.stdout)
    sequences = [hap["pattern"] for hap in payload["haplotypes"]]
    samples = [hap["samples"] for hap in payload["haplotypes"]]

    net = _run_network("tcs", {"sequences": sequences, "samples": samples, "epsilon": 0})
    assert net.returncode == 0
    network = json.loads(net.stdout)

    weights = {(edge["source"], edge["target"]): edge["weight"] for edge in network["edges"]}
    # Hap01 differs from Hap02 at one site and from Hap03 at two.
    assert weights[(0, 1)] == 1, "one differing site should weigh 1, not 2"
    assert weights[(0, 2)] == 2


def test_tcs_leaves_no_isolated_median_vertices(network_available: Path) -> None:
    completed = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary"])
    payload = json.loads(completed.stdout)
    sequences = [hap["pattern"] for hap in payload["haplotypes"]]
    samples = [hap["samples"] for hap in payload["haplotypes"]]

    net = _run_network("tcs", {"sequences": sequences, "samples": samples, "epsilon": 0})
    network = json.loads(net.stdout)

    degree = {node["id"]: 0 for node in network["nodes"]}
    for edge in network["edges"]:
        degree[edge["source"]] += 1
        degree[edge["target"]] += 1

    isolated_medians = [
        node["id"] for node in network["nodes"]
        if node.get("is_median") and degree[node["id"]] == 0
    ]
    assert isolated_medians == [], f"collapsed medians left isolated: {isolated_medians}"


def test_tcs_intermediates_stay_on_site_boundaries(network_available: Path) -> None:
    sequences = ["0/0|0/0|0/0|0/0|0/0", "1/1|1/1|1/1|1/1|1/1"]
    net = _run_network("tcs", {"sequences": sequences, "samples": [["s1"], ["s2"]], "epsilon": 0})
    assert net.returncode == 0
    network = json.loads(net.stdout)

    for node in network["nodes"]:
        if not node.get("is_median"):
            continue
        for token in node["sequence"].split("|"):
            assert token in {"0/0", "1/1"}, (
                f"intermediate token {token!r} is not a valid genotype"
            )


def test_mjn_with_epsilon_adds_medians_unlike_msn(network_available: Path) -> None:
    completed = _run_view(["scaffold_1:4300-5209", "--by", "region", "--output", "summary"])
    payload = json.loads(completed.stdout)
    sequences = [hap["pattern"] for hap in payload["haplotypes"]]
    samples = [hap["samples"] for hap in payload["haplotypes"]]

    msn = json.loads(_run_network("msn", {"sequences": sequences, "samples": samples, "epsilon": 2}).stdout)
    mjn = json.loads(_run_network("mjn", {"sequences": sequences, "samples": samples, "epsilon": 2}).stdout)

    msn_medians = sum(1 for node in msn["nodes"] if node.get("is_median"))
    mjn_medians = sum(1 for node in mjn["nodes"] if node.get("is_median"))
    assert msn_medians == 0
    assert mjn_medians > 0, "MJN with epsilon 2 produced no median vertices"

    degree = {node["id"]: 0 for node in mjn["nodes"]}
    for edge in mjn["edges"]:
        degree[edge["source"]] += 1
        degree[edge["target"]] += 1
    isolated = [n["id"] for n in mjn["nodes"] if n.get("is_median") and degree[n["id"]] == 0]
    assert isolated == [], f"MJN left isolated median vertices: {isolated}"


def test_network_backend_rejects_malformed_json(network_available: Path) -> None:
    completed = subprocess.run(
        [str(_network_bin()), "msn"],
        input="not json",
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode != 0
    assert "JSON parse error" in completed.stderr
    assert completed.stdout.strip() == ""


# ---------------------------------------------------------------------------
# BED coordinates are 0-based half-open.
# ---------------------------------------------------------------------------

def test_bed_left_edge_is_half_open(tmp_path: Path) -> None:
    vcf = tmp_path / "boundary.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t99\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\n"
        "chr1\t101\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\n",
        encoding="utf-8",
    )

    import pysam
    compressed = tmp_path / "boundary.vcf.gz"
    pysam.tabix_compress(str(vcf), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)

    bed = tmp_path / "windows.bed"
    bed.write_text("chr1\t99\t100\twin\n", encoding="utf-8")

    completed = subprocess.run(
        [str(_cpp_bin()), "view-bed-jsonl", str(compressed), str(bed), "--output", "summary"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    payload = json.loads(completed.stdout.splitlines()[0])
    # BED chr1 99 100 covers exactly 1-based base 100.
    assert [site["pos"] for site in payload["sites"]] == [100]


def test_bed_inverted_window_is_rejected(tmp_path: Path) -> None:
    bed = tmp_path / "bad.bed"
    bed.write_text("scaffold_1\t5000\t4300\tbad\n", encoding="utf-8")

    completed = subprocess.run(
        [str(_cpp_bin()), "view-bed-jsonl", str(DATA_DIR / "var.sorted.vcf.gz"), str(bed),
         "--output", "summary"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode != 0
    assert "start must be <= end" in completed.stderr
