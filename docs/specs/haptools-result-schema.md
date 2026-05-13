# haptools Result Schema

## Scope

This schema defines the default v1 `haptools view` output contract in `tsv` mode.

## Default Output Files

- `hapresult*.tsv`
- `hap_summary*.tsv`

For a single selector, default names are:

- `hapresult.tsv`
- `hap_summary.tsv`

For BED batch selectors, file names are selector-scoped with region slug suffixes.

## TSV Row Contract

Both files follow the same first-four-row structure:

1. `CHR`
2. `POS`
3. `INFO`
4. `ALLELE`

Rules:

- First column is the `Hap` identity axis (`CHR/POS/INFO/ALLELE/H###`).
- `hapresult` last column is `Accession`.
- `hap_summary` last two columns are `Accession` and `freq`.
- Data rows after `ALLELE` are hap rows (`H001`, `H002`, ...).

## Selector Semantics

- `-r chr:pos` is one-site mode.
- `-r chr:start-end` is region mode.
- `-R bed` is region mode per BED row.

## Optional Artifacts

- `--plot` writes selector-scoped PDF files.
- `--gff/--gff3` writes `gff_ann_summary.tsv` and places annotation text in `INFO` row cells.

## Compatibility Mode

`--output-format jsonl` remains available for machine-readable metadata compatibility, but is not the default output contract.
