from __future__ import annotations

import csv
import math
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence


class PhenotypeError(ValueError):
    """Raised when phenotype input cannot be parsed or analyzed."""


@dataclass(frozen=True)
class SampleHaplotype:
    sample: str
    haplotype: str


@dataclass(frozen=True)
class PhenotypeRecord:
    sample: str
    haplotype: str
    trait: str
    value: float
    population: str | None = None


@dataclass(frozen=True)
class PhenotypeDataset:
    records: tuple[PhenotypeRecord, ...]
    traits: tuple[str, ...]
    haplotypes: tuple[str, ...]
    populations: tuple[str, ...]
    sample_count: int
    matched_sample_count: int


STAT_COLUMNS = [
    "trait",
    "population",
    "group1",
    "group2",
    "count1",
    "count2",
    "mean1",
    "mean2",
    "std1",
    "std2",
    "anova_f",
    "anova_p",
    "method",
    "pairwise_stat",
    "p_value",
    "p_adjusted",
    "significance",
    "reject",
    "effective_n",
]

SUMMARY_COLUMNS = [
    "trait",
    "population",
    "haplotype",
    "count",
    "mean",
    "std",
    "min",
    "max",
    "effective_n",
]

_SAMPLE_ALIASES = {"sample", "samples", "accession", "accessions", "id", "ids", "individual", "individuals"}
_HAP_ALIASES = {"hap", "haps", "haplotype", "haplotypes", "group", "groups", "hap_group", "hap_groups"}
_NA_VALUES = {"", "na", "nan", "null", "none", "."}


def read_sample_haplotypes(path: str | Path, delimiter: str = "auto") -> list[SampleHaplotype]:
    """Read sample-to-haplotype assignments from haplokit or two-column tables."""
    rows = _read_csv_rows(path, delimiter)
    if not rows:
        raise PhenotypeError(f"{path} is empty")

    allele_index = next((idx for idx, row in enumerate(rows) if row and row[0] == "ALLELE"), None)
    if allele_index is not None:
        return _read_hapresult_rows(rows, path)

    header = [cell.strip() for cell in rows[0]]
    lowered = [_normalize_header(cell) for cell in header]
    sample_idx = _find_alias_index(lowered, _SAMPLE_ALIASES)
    hap_idx = _find_alias_index(lowered, _HAP_ALIASES)
    data_rows = rows[1:]

    if sample_idx is None or hap_idx is None:
        if len(rows[0]) < 2:
            raise PhenotypeError("sample-haplotype input must contain at least two columns")
        sample_idx = 0
        hap_idx = 1
        data_rows = rows

    records = [
        SampleHaplotype(sample=row[sample_idx].strip(), haplotype=row[hap_idx].strip())
        for row in data_rows
        if len(row) > max(sample_idx, hap_idx) and row[sample_idx].strip() and row[hap_idx].strip()
    ]
    _reject_duplicate_samples(records, "sample-haplotype")
    return records


def read_phenotype_table(
    path: str | Path,
    traits: Sequence[str] | None = None,
    delimiter: str = "auto",
) -> tuple[dict[str, dict[str, str]], tuple[str, ...]]:
    """Read raw phenotype values keyed by sample ID."""
    rows = _read_csv_rows(path, delimiter)
    if not rows:
        raise PhenotypeError(f"{path} is empty")
    header = [cell.strip() for cell in rows[0]]
    if len(header) < 2:
        raise PhenotypeError("phenotype table must contain a sample column and at least one trait")

    sample_idx = 0
    trait_columns = [name for idx, name in enumerate(header) if idx != sample_idx and name]
    selected_traits = _resolve_traits(trait_columns, traits)

    phenotypes: dict[str, dict[str, str]] = {}
    for line_no, row in enumerate(rows[1:], start=2):
        if not row or len(row) <= sample_idx:
            continue
        sample = row[sample_idx].strip()
        if not sample:
            continue
        if sample in phenotypes:
            raise PhenotypeError(f"duplicate sample ID in phenotype table at line {line_no}: {sample}")
        values: dict[str, str] = {}
        for trait in selected_traits:
            column_idx = header.index(trait)
            values[trait] = row[column_idx].strip() if column_idx < len(row) else ""
        phenotypes[sample] = values

    if not phenotypes:
        raise PhenotypeError("phenotype table did not contain any sample rows")
    return phenotypes, tuple(selected_traits)


def read_population_groups(path: str | Path, delimiter: str = "auto") -> dict[str, str]:
    """Read sample-to-population assignments from a two-column CSV/TSV file."""
    rows = _read_csv_rows(path, delimiter)
    if not rows:
        raise PhenotypeError(f"{path} is empty")

    lowered = [_normalize_header(cell) for cell in rows[0]]
    sample_idx = _find_alias_index(lowered, _SAMPLE_ALIASES)
    group_idx = _find_alias_index(lowered, {"population", "pop", "group", "pop_group", "population_group"})
    data_rows = rows[1:]
    if sample_idx is None or group_idx is None:
        if len(rows[0]) < 2:
            raise PhenotypeError("population group input must contain at least two columns")
        sample_idx = 0
        group_idx = 1
        data_rows = rows

    groups: dict[str, str] = {}
    for line_no, row in enumerate(data_rows, start=2):
        if len(row) <= max(sample_idx, group_idx):
            continue
        sample = row[sample_idx].strip()
        population = row[group_idx].strip()
        if not sample or not population:
            continue
        if sample in groups:
            raise PhenotypeError(f"duplicate sample ID in population group table at line {line_no}: {sample}")
        groups[sample] = population
    if not groups:
        raise PhenotypeError("population group table did not contain any sample rows")
    return groups


def load_phenotype_dataset(
    haplotypes: str | Path | Sequence[SampleHaplotype],
    phenotypes: str | Path,
    *,
    traits: Sequence[str] | None = None,
    hap_delimiter: str = "auto",
    phenotype_delimiter: str = "auto",
    population_file: str | Path | None = None,
    population_delimiter: str = "auto",
) -> PhenotypeDataset:
    """Join sample haplotypes with phenotype values in long-table form."""
    hap_records = (
        list(haplotypes)
        if not isinstance(haplotypes, (str, Path))
        else read_sample_haplotypes(haplotypes, delimiter=hap_delimiter)
    )
    if not hap_records:
        raise PhenotypeError("no sample-haplotype assignments were found")

    phenotype_rows, selected_traits = read_phenotype_table(phenotypes, traits=traits, delimiter=phenotype_delimiter)
    population_map = read_population_groups(population_file, delimiter=population_delimiter) if population_file else {}
    hap_by_sample = {record.sample: record.haplotype for record in hap_records}
    matched_samples = sorted(set(hap_by_sample) & set(phenotype_rows))
    if not matched_samples:
        raise PhenotypeError("no overlapping sample IDs between haplotype and phenotype inputs")

    records: list[PhenotypeRecord] = []
    for sample in matched_samples:
        haplotype = hap_by_sample[sample]
        population = population_map.get(sample, "Unknown") if population_map else None
        values = phenotype_rows[sample]
        for trait in selected_traits:
            parsed = _parse_numeric(values.get(trait, ""))
            if parsed is None:
                continue
            records.append(
                PhenotypeRecord(
                    sample=sample,
                    haplotype=haplotype,
                    trait=trait,
                    value=parsed,
                    population=population,
                )
            )

    if not records:
        raise PhenotypeError("no numeric phenotype values were found for overlapping samples")

    return PhenotypeDataset(
        records=tuple(records),
        traits=tuple(selected_traits),
        haplotypes=tuple(sort_haplotype_labels({record.haplotype for record in hap_records})),
        populations=tuple(
            dict.fromkeys([*population_map.values(), *(record.population for record in records if record.population)])
        )
        if population_map
        else (),
        sample_count=len(hap_records),
        matched_sample_count=len(matched_samples),
    )


def group_values(
    records: Iterable[PhenotypeRecord],
    trait: str,
    *,
    min_hap_size: int = 1,
    haplotypes: Sequence[str] | None = None,
    population: str | None = None,
) -> dict[str, list[float]]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for record in records:
        if record.trait == trait and (population is None or record.population == population):
            grouped[record.haplotype].append(record.value)

    labels = list(haplotypes) if haplotypes is not None else sort_haplotype_labels(grouped)
    return {label: grouped[label] for label in labels if len(grouped.get(label, [])) >= min_hap_size}


def summarize_groups(
    records: Iterable[PhenotypeRecord],
    *,
    traits: Sequence[str] | None = None,
    min_hap_size: int = 1,
    populations: Sequence[str] | None = None,
) -> list[dict[str, object]]:
    record_list = list(records)
    selected_traits = list(traits) if traits is not None else _traits_from_records(record_list)
    selected_populations = _population_strata(record_list, populations)
    rows: list[dict[str, object]] = []
    for trait in selected_traits:
        for population in selected_populations:
            grouped = group_values(
                record_list,
                trait,
                min_hap_size=min_hap_size,
                population=population,
            )
            effective_n = sum(len(values) for values in grouped.values())
            for haplotype, values in grouped.items():
                rows.append(
                    {
                        "trait": trait,
                        "population": population or "ALL",
                        "haplotype": haplotype,
                        "count": len(values),
                        "mean": _mean(values),
                        "std": _std(values),
                        "min": min(values),
                        "max": max(values),
                        "effective_n": effective_n,
                    }
                )
    return rows


def pairwise_statistics(
    records: Iterable[PhenotypeRecord],
    *,
    traits: Sequence[str] | None = None,
    min_hap_size: int = 5,
    method: str = "welch",
    adjust: str = "bonferroni",
    alpha: float = 0.05,
    populations: Sequence[str] | None = None,
) -> list[dict[str, object]]:
    """Compute per-trait ANOVA plus pairwise haplotype comparisons."""
    stats = _require_scipy()
    record_list = list(records)
    selected_traits = list(traits) if traits is not None else _traits_from_records(record_list)
    selected_populations = _population_strata(record_list, populations)
    rows: list[dict[str, object]] = []

    for trait in selected_traits:
        for population in selected_populations:
            grouped = group_values(record_list, trait, min_hap_size=min_hap_size, population=population)
            if len(grouped) < 2:
                continue

            labels = list(grouped)
            values = [grouped[label] for label in labels]
            effective_n = sum(len(group) for group in values)
            anova_f, anova_p = _anova(stats, values)
            pair_rows = _pairwise_rows(stats, labels, values, method)
            pair_count = len(pair_rows)

            effective_adjust = "none" if method == "tukey" else adjust
            for group1, group2, pair_stat, p_value in pair_rows:
                vals1 = grouped[group1]
                vals2 = grouped[group2]
                p_adjusted = _adjust_p_value(p_value, pair_count, effective_adjust)
                rows.append(
                    {
                        "trait": trait,
                        "population": population or "ALL",
                        "group1": group1,
                        "group2": group2,
                        "count1": len(vals1),
                        "count2": len(vals2),
                        "mean1": _mean(vals1),
                        "mean2": _mean(vals2),
                        "std1": _std(vals1),
                        "std2": _std(vals2),
                        "anova_f": anova_f,
                        "anova_p": anova_p,
                        "method": method,
                        "pairwise_stat": pair_stat,
                        "p_value": p_value,
                        "p_adjusted": p_adjusted,
                        "significance": significance_label(p_adjusted),
                        "reject": _finite(p_adjusted) and p_adjusted < alpha,
                        "effective_n": effective_n,
                    }
                )
    return rows


def population_pairwise_statistics(
    records: Iterable[PhenotypeRecord],
    *,
    traits: Sequence[str] | None = None,
    min_hap_size: int = 5,
    method: str = "welch",
    adjust: str = "bonferroni",
    alpha: float = 0.05,
    haplotypes: Sequence[str] | None = None,
    populations: Sequence[str] | None = None,
) -> list[dict[str, object]]:
    """Compute per-haplotype phenotype comparisons across populations."""
    stats = _require_scipy()
    record_list = list(records)
    selected_traits = list(traits) if traits is not None else _traits_from_records(record_list)
    selected_haplotypes = list(haplotypes) if haplotypes is not None else sort_haplotype_labels(
        {record.haplotype for record in record_list}
    )
    selected_populations = _population_strata(record_list, populations)
    rows: list[dict[str, object]] = []

    for trait in selected_traits:
        for haplotype in selected_haplotypes:
            grouped: dict[str, list[float]] = defaultdict(list)
            for record in record_list:
                if record.trait == trait and record.haplotype == haplotype and record.population:
                    grouped[record.population].append(record.value)
            grouped = {
                population: grouped[population]
                for population in selected_populations
                if population and len(grouped.get(population, [])) >= min_hap_size
            }
            if len(grouped) < 2:
                continue

            labels = list(grouped)
            values = [grouped[label] for label in labels]
            effective_n = sum(len(group) for group in values)
            anova_f, anova_p = _anova(stats, values)
            pair_rows = _pairwise_rows(stats, labels, values, method)
            pair_count = len(pair_rows)
            effective_adjust = "none" if method == "tukey" else adjust
            for group1, group2, pair_stat, p_value in pair_rows:
                vals1 = grouped[group1]
                vals2 = grouped[group2]
                p_adjusted = _adjust_p_value(p_value, pair_count, effective_adjust)
                rows.append(
                    {
                        "trait": trait,
                        "haplotype": haplotype,
                        "group1": group1,
                        "group2": group2,
                        "count1": len(vals1),
                        "count2": len(vals2),
                        "mean1": _mean(vals1),
                        "mean2": _mean(vals2),
                        "std1": _std(vals1),
                        "std2": _std(vals2),
                        "anova_f": anova_f,
                        "anova_p": anova_p,
                        "method": method,
                        "pairwise_stat": pair_stat,
                        "p_value": p_value,
                        "p_adjusted": p_adjusted,
                        "significance": significance_label(p_adjusted),
                        "reject": _finite(p_adjusted) and p_adjusted < alpha,
                        "effective_n": effective_n,
                    }
                )
    return rows


def write_stat_tsv(rows: Sequence[dict[str, object]], output_path: str | Path) -> Path:
    return _write_tsv(STAT_COLUMNS, rows, output_path)


def write_summary_tsv(rows: Sequence[dict[str, object]], output_path: str | Path) -> Path:
    return _write_tsv(SUMMARY_COLUMNS, rows, output_path)


def significance_label(p_value: float) -> str:
    if not _finite(p_value):
        return "NA"
    if p_value < 0.0001:
        return "****"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def sort_haplotype_labels(labels: Iterable[str]) -> list[str]:
    def key(label: str) -> tuple[object, ...]:
        parts = re.split(r"(\d+)", str(label))
        return tuple(int(part) if part.isdecimal() else part.lower() for part in parts)

    return sorted((str(label) for label in labels), key=key)


def _read_hapresult_rows(rows: list[list[str]], path: str | Path) -> list[SampleHaplotype]:
    header_idx = next(idx for idx, row in enumerate(rows) if row and row[0] == "ALLELE")
    header = rows[header_idx]
    records: list[SampleHaplotype] = []

    if header[-1] == "Accession":
        sample_idx = len(header) - 1
        for row in rows[header_idx + 1 :]:
            if len(row) <= sample_idx:
                continue
            haplotype = row[0].strip()
            for sample in _split_samples(row[sample_idx]):
                records.append(SampleHaplotype(sample=sample, haplotype=haplotype))
    elif len(header) >= 2 and header[-2] == "Accession":
        sample_idx = len(header) - 2
        for row in rows[header_idx + 1 :]:
            if len(row) <= sample_idx:
                continue
            haplotype = row[0].strip()
            for sample in _split_samples(row[sample_idx]):
                records.append(SampleHaplotype(sample=sample, haplotype=haplotype))
    else:
        raise PhenotypeError(f"{path} does not look like a haplokit hapresult.tsv file")

    if not records:
        raise PhenotypeError(f"{path} did not contain sample-haplotype rows")
    _reject_duplicate_samples(records, "haplokit hapresult")
    return records


def _read_csv_rows(path: str | Path, delimiter: str) -> list[list[str]]:
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    selected = _infer_delimiter(lines, delimiter)
    return [[cell.strip() for cell in row] for row in csv.reader(lines, delimiter=selected)]


def _infer_delimiter(lines: Sequence[str], delimiter: str) -> str:
    if delimiter == "tab":
        return "\t"
    if delimiter == "comma":
        return ","
    if delimiter != "auto":
        raise PhenotypeError("delimiter must be one of: auto, tab, comma")
    first = lines[0]
    return "\t" if first.count("\t") >= first.count(",") else ","


def _normalize_header(value: str) -> str:
    return value.strip().lower().replace(" ", "_").replace("-", "_")


def _find_alias_index(headers: Sequence[str], aliases: set[str]) -> int | None:
    return next((idx for idx, header in enumerate(headers) if header in aliases), None)


def _resolve_traits(columns: Sequence[str], requested: Sequence[str] | None) -> list[str]:
    if requested is None:
        return list(columns)

    selected: list[str] = []
    for trait in requested:
        value = str(trait)
        if value in columns:
            selected.append(value)
            continue
        if value.isdecimal():
            index = int(value) - 1
            if 0 <= index < len(columns):
                selected.append(columns[index])
                continue
        raise PhenotypeError(f"phenotype trait not found: {trait}")
    return list(dict.fromkeys(selected))


def _parse_numeric(value: str | None) -> float | None:
    if value is None or value.strip().lower() in _NA_VALUES:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _split_samples(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def _reject_duplicate_samples(records: Sequence[SampleHaplotype], source: str) -> None:
    counts = Counter(record.sample for record in records)
    duplicates = [sample for sample, count in counts.items() if count > 1]
    if duplicates:
        preview = ", ".join(duplicates[:5])
        raise PhenotypeError(f"duplicate sample ID in {source}: {preview}")


def _traits_from_records(records: Sequence[PhenotypeRecord]) -> list[str]:
    return list(dict.fromkeys(record.trait for record in records))


def _population_strata(records: Sequence[PhenotypeRecord], populations: Sequence[str] | None) -> list[str | None]:
    if populations is not None:
        return list(populations)
    discovered = sort_haplotype_labels({record.population for record in records if record.population})
    return discovered if discovered else [None]


def _require_scipy():
    try:
        from scipy import stats
    except ImportError as exc:
        raise PhenotypeError("phenotype statistics require scipy; install haplokit with scipy available") from exc
    return stats


def _anova(stats, values: Sequence[Sequence[float]]) -> tuple[float, float]:
    try:
        result = stats.f_oneway(*values)
    except Exception:
        return math.nan, math.nan
    return float(result.statistic), float(result.pvalue)


def _pairwise_rows(stats, labels: Sequence[str], values: Sequence[Sequence[float]], method: str):
    if method == "tukey":
        if not hasattr(stats, "tukey_hsd"):
            raise PhenotypeError("method 'tukey' requires scipy.stats.tukey_hsd")
        result = stats.tukey_hsd(*values)
        rows = []
        for i, j in combinations(range(len(labels)), 2):
            rows.append((labels[i], labels[j], _array_float(result.statistic, i, j), _array_float(result.pvalue, i, j)))
        return rows

    rows = []
    for i, j in combinations(range(len(labels)), 2):
        vals1 = values[i]
        vals2 = values[j]
        if method == "welch":
            result = stats.ttest_ind(vals1, vals2, equal_var=False, nan_policy="omit")
        elif method == "student":
            result = stats.ttest_ind(vals1, vals2, equal_var=True, nan_policy="omit")
        elif method == "mannwhitney":
            result = stats.mannwhitneyu(vals1, vals2, alternative="two-sided")
        else:
            raise PhenotypeError(f"unsupported phenotype test method: {method}")
        rows.append((labels[i], labels[j], float(result.statistic), float(result.pvalue)))
    return rows


def _array_float(values, i: int, j: int) -> float:
    try:
        return float(values[i, j])
    except TypeError:
        return float(values[i][j])


def _adjust_p_value(p_value: float, pair_count: int, adjust: str) -> float:
    if not _finite(p_value):
        return math.nan
    if adjust == "none":
        return p_value
    if adjust == "bonferroni":
        return min(1.0, p_value * pair_count)
    raise PhenotypeError(f"unsupported p-value adjustment: {adjust}")


def _mean(values: Sequence[float]) -> float:
    return float(statistics.fmean(values))


def _std(values: Sequence[float]) -> float:
    return float(statistics.stdev(values)) if len(values) > 1 else 0.0


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def _write_tsv(columns: Sequence[str], rows: Sequence[dict[str, object]], output_path: str | Path) -> Path:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_format_cell(row.get(column, "")) for column in columns])
    return out.resolve()


def _format_cell(value: object) -> str:
    if value is True:
        return "TRUE"
    if value is False:
        return "FALSE"
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else "NA"
    return str(value)
