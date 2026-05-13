# haptools view CLI

## Purpose

`haptools view` is the primary v1 command surface for strict hap grouping over indexed VCF/BCF inputs.

The command is designed to preserve a `bcftools`-like selector vocabulary while exposing hap-specific grouping and output controls.

## Shared Selectors

- `-r, --region chr:start-end`
- `-r, --region chr:pos`
- `-R, --regions-file regions.bed`
- `-S, --samples-file sample.list`

Rules:

- Exactly one of `-r` or `-R` is required.
- `-r chr:pos` is always interpreted as one-site selection.
- `-r chr:start-end` is always interpreted as region selection.
- BED rows remain independent units throughout the command.

## Grouping Controls

- `--by auto|region|site`
- `--by region`
- `--by site`

Rules:

- Default is `--by auto`.
- `--by auto` infers mode from selector shape: `chr:pos -> site`, `chr:start-end -> region`.
- `--by region` and `--by site` are compatibility overrides and must not conflict with `-r` selector semantics.
- v1 does not expose an approximate grouping mode in the CLI contract.

## Preprocessing Controls

- `--impute`

Rules:

- `--impute` is a boolean switch.
- In v1 it means: replace missing values with the reference state before grouping.
- This preprocessing choice must be recorded in output metadata.

## Annotation Controls

- `--gff path/to/file.gff3` (alias: `--gff3`)

Rules:

- Annotation is interval-level only in v1.
- The payload may include overlap status, nearest feature, distance, strand, and feature type.
- This is not a full variant-effect predictor.

## Output Controls

- `--output summary|detail`
- `--output-format tsv|jsonl`
- `--output-file path`
- `--plot`

Rules:

- Default output format is `tsv`.
- In `tsv` mode, output is hapresult/hap_summary style and writes both files.
- `jsonl` remains available as a compatibility mode.
- `--plot` is a `view` flag, not a separate subcommand.

## Examples

```bash
haptools view in.vcf.gz -r chr1:1000-2000
haptools view in.vcf.gz -R genes.bed -S sample.list
haptools view in.vcf.gz -r chr1:1450
haptools view in.vcf.gz -r chr1:1000-2000 --impute
haptools view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3
haptools view in.vcf.gz -r chr1:1000-2000 --plot
```
