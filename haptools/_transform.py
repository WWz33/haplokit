from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


def read_hap_summary_tsv(path: str | Path) -> list[list[str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t")]


def read_popgroup(path: str | Path) -> dict[str, str]:
    """Read tab-separated ID POP file. Returns {sample_id: pop_name}."""
    pop: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and row[0].strip():
                pop[row[0].strip()] = row[1].strip()
    return pop


def _unique_alleles(rows: list[list[str]], var_end: int) -> list[str]:
    meta = {"POS", "ALLELE"}
    alleles: set[str] = set()
    for row in rows:
        if not row or row[0] in meta:
            continue
        for cell in row[1:var_end]:
            v = cell.strip()
            if v and v not in {"", "NA", ".", "/"}:
                for token in v.replace("/", ",").split(","):
                    t = token.strip()
                    if t and t not in {"", "NA", "."}:
                        alleles.add(t)
    return sorted(alleles)


def _indel_footnotes(rows: list[list[str]], var_end: int) -> tuple[list[list[str]], str]:
    """Replace long indel alleles with i1/i2/... and return footnote string."""
    meta = {"POS", "ALLELE"}
    allele_row = next((r for r in rows if r and r[0] == "ALLELE"), None)
    if allele_row is None:
        return rows, ""

    threshold = 2
    indel_alleles: set[str] = set()
    for cell in allele_row[1:var_end]:
        for token in cell.replace("/", ",").split(","):
            t = token.strip()
            if len(t) > threshold:
                indel_alleles.add(t)

    if not indel_alleles:
        return rows, ""

    sorted_indels = sorted(indel_alleles)
    notes = {a: f"i{i+1}" for i, a in enumerate(sorted_indels)}
    footnote = "; ".join(f"{v}:{k}" for k, v in notes.items())

    new_rows: list[list[str]] = []
    for row in rows:
        new_row = list(row)
        if row[0] in meta or row[0].startswith("H"):
            for c in range(1, min(var_end, len(new_row))):
                cell = new_row[c]
                for k, v in notes.items():
                    cell = cell.replace(k, v)
                new_row[c] = cell
        new_rows.append(new_row)

    return new_rows, footnote


def transform_for_display(
    rows: list[list[str]],
    pop_data: dict[str, str] | None = None,
) -> tuple[list[list[str]], str]:
    """Transform raw hap_summary rows for publication display.

    Returns (transformed_rows, region_title).

    - Remove CHR row (→ title), remove INFO row
    - Remove freq column
    - Change Accession from sample IDs to count/total format
    - If pop_data provided, add population breakdown columns before n/N
    """
    if not rows:
        return rows, ""

    # ── Extract region info from CHR/POS rows before removal ──
    chrom = ""
    positions: list[str] = []
    total_ind = ""
    _acc_col_early: int | None = None
    for row in rows:
        if row and row[0] == "ALLELE":
            for i, cell in enumerate(row):
                if cell.strip() == "Accession":
                    _acc_col_early = i
                    break
            break
    for row in rows:
        if not row:
            continue
        if row[0] == "CHR":
            chrom = row[1].strip() if len(row) > 1 else ""
        elif row[0] == "POS":
            for i, cell in enumerate(row):
                cs = cell.strip()
                if cs.startswith("Individuals"):
                    total_ind = row[i + 1].strip() if i + 1 < len(row) else ""
                elif cs.isdigit() and i < (_acc_col_early or len(row)):
                    positions.append(cs)

    region_title = ""
    if chrom and positions:
        region_title = f"{chrom}:{min(positions)}-{max(positions)}"

    # ── Locate Accession column ──
    acc_col: int | None = None
    for row in rows:
        if row and row[0] == "ALLELE":
            for i, cell in enumerate(row):
                if cell.strip() == "Accession":
                    acc_col = i
                    break
            break

    # ── Population info ──
    pop_names: list[str] = []
    pop_totals: dict[str, int] = {}
    if pop_data:
        pop_counts = Counter(pop_data.values())
        pop_names = sorted(pop_counts.keys())
        pop_totals = dict(pop_counts)

    # ── Find freq column ──
    freq_col: int | None = None
    for row in rows:
        if row and row[0] == "ALLELE":
            if len(row) > 1 and row[-1].strip() == "freq":
                freq_col = len(row) - 1
            break

    # ── Transform rows ──
    new_rows: list[list[str]] = []
    for row in rows:
        if not row:
            continue
        if row[0] in {"CHR", "INFO"}:
            continue

        new_row = list(row)

        if freq_col is not None and len(new_row) > freq_col:
            new_row = new_row[:freq_col]

        if new_row[0] == "POS" and acc_col is not None:
            new_row = new_row[:acc_col]

        if new_row[0] == "ALLELE" and acc_col is not None:
            pop_headers = pop_names + ["n/N"]
            new_row = new_row[:acc_col] + pop_headers

        elif new_row[0].startswith("H") and acc_col is not None:
            target_col = acc_col
            if target_col < len(new_row):
                raw = new_row[target_col].strip()
                sample_ids = [s.strip() for s in raw.split(";") if s.strip()]
                count = len(sample_ids)

                pop_cols: list[str] = []
                if pop_data and pop_names:
                    for pname in pop_names:
                        pc = sum(1 for s in sample_ids if pop_data.get(s) == pname)
                        pop_cols.append(f"{pc}/{pop_totals.get(pname, 0)}")

                acc_text = f"{count}/{total_ind}" if total_ind else str(count)
                pop_cols.append(acc_text)
                new_row = new_row[:target_col] + pop_cols

        new_rows.append(new_row)

    return new_rows, region_title
