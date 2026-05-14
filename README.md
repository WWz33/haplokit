# haplokit

CLI haplotype viewer with bcftools-like selectors, C++ backend, and Python plotting.

<!-- README-I18N:START -->

**English** | [汉语](./README.zh-CN.md)

<!-- README-I18N:END -->

## Installation

```bash
pip install haplokit
```

> Source build requires Linux/WSL, Python 3.10+, C++17 toolchain, CMake 3.22+ — see [Contributing](#contributing).

## Quick Start

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

Output:

- `out/hapresult.tsv` — per-sample haplotype detail
- `out/hap_summary.tsv` — haplotype count summary

## Usage Scenarios

### 1. Region query — strict haplotype grouping

Identify all distinct haplotypes in a genomic region.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
```

Produces `hapresult.tsv` + `hap_summary.tsv` in `out/`. Each haplotype row shows the exact allele pattern; samples with any heterozygous or missing call are excluded.

### 2. Single-site query

Analyze haplotype at one variant position.

```bash
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

`--by` auto-resolves to `site` for `chr:pos` selectors.

### 3. Gene annotation + figure

Overlay gene structure on the haplotype table.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
```

`genes.gff3` format (standard GFF3):

```text
chr1	.	gene	1000	3000	.	+	.	ID=gene1;Name=GeneA
chr1	.	CDS	1200	1500	.	+	0	ID=cds1;Parent=gene1
```

Adds SnpEff-style functional category strip (CDS, UTR, exon, intron, intergenic) above variant positions. Writes figure (`out/*.png`) + `gff_ann_summary.tsv`.

<img src="plottable.png" alt="Haplotype summary table" width="800">

Figure components:

- **Title**: region + overlapping gene name (when `--gff` provided)
- **Function strip** (`--gff` only): colored bar classifying each variant by functional category
- **POS / ALLELE rows**: variant positions and alternate alleles
- **Haplotype rows** (H001, H002, ...): allele per position; empty = reference
- **Population columns** (`--population`): sample counts per haplotype per group
- **n/N**: haplotype frequency
- **Legend** (`--gff` only): functional category colors
- **Indel footnotes**: multi-allele indels annotated with superscript markers

### 4. Population grouping

Compare haplotype distributions across populations.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

`popgroup.txt` (tab-separated: `sample<TAB>population`):

```text
C1	wild
C2	wild
C13	landrace
```

Adds population columns to the table and figure.

### 5. Geographic distribution map

Map haplotype composition at sampling locations.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo sample_geo.txt --plot --output-file out
```

`sample_geo.txt` (tab-separated: `ID<TAB>longitude<TAB>latitude`):

```text
ID	longitude	latitude
C1	116.40	39.90
C4	121.47	31.23
```

<img src="plotmap.png" alt="Haplotype geographic distribution" width="600">

Figure components:

- **Pie charts**: haplotype composition per location; size ∝ √(sample count)
- **Count labels**: total samples at center of each pie
- **Legend**: haplotype color key
- **Base map**: GeoJSON province boundaries (China)

### 6. BED batch processing

Process multiple regions in one run.

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

`regions.bed` (≥3 tab-separated columns):

```text
chr1	1000	2000
chr2	5000	6000
```

Each BED row is processed independently. Output files are suffixed by region slug (`_chr1_1000_2000`).

### 7. Approximate grouping

Cluster similar haplotypes within a tolerance.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff` (0–1): haplotypes differing at ≤ 20% of positions merge into one group. Grouping mode changes from `strict-region` to `approx-region`.

### 8. Sample subset + imputation

Restrict analysis to specific samples; fill missing calls as reference.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list` (one sample ID per line):

```text
C1
C5
C16
```

`--impute` treats missing GT as `0/0`, increasing sample retention.

## Full Parameters

```
haplokit view <input_vcf> (-r <region> | -R <regions.bed>) [options]
```

`<input_vcf>` must be an indexed VCF/BCF (`.vcf.gz` + `.tbi`, or BCF index).

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-r, --region` | string | — | `chr:start-end` or `chr:pos` |
| `-R, --regions-file` | path | — | BED file (≥3 tab-separated columns) |
| `-S, --samples-file` | path | — | One sample ID per line |
| `--by` | `auto\|region\|site` | `auto` | Grouping mode; auto infers from selector shape |
| `--impute` | flag | off | Impute missing GT as reference |
| `-g, --gff` | path | — | GFF3/GTF for gene annotation |
| `-p, --population` | path | — | Tab-separated sample → population map |
| `--output` | `summary\|detail` | `summary` | JSONL mode only; TSV always writes both |
| `--output-format` | `tsv\|jsonl` | `tsv` | Output format |
| `--output-file` | path | — | Output directory, prefix, or JSONL file |
| `--plot` | flag | off | Generate haplotype table figure |
| `--plot-format` | `png\|pdf\|svg\|tiff` | `png` | Figure format |
| `--max-diff` | float [0,1] | — | Approximate grouping threshold |
| `--geo` | path | — | Sample geographic coordinates for map |

Selector rules: `-r` and `-R` are mutually exclusive and one is required. `--by site` only valid with `-r chr:pos`.

## Backend

C++ backend (`haplokit_cpp`) handles VCF reading and haplotype grouping. Discovery order:

1. `HAPLOKIT_CPP_BIN` env var
2. Packaged binary: `haplokit/_bin/haplokit_cpp`
3. Repo build: `build-wsl/haplokit_cpp` → `build/haplokit_cpp`
4. Fallback: auto-run `cmake` build

Vendored libraries:

- **[htslib](https://github.com/samtools/htslib)** — VCF/BCF reading with indexed random access
- **[gffsub](https://github.com/WWz33/gffsub)** — GFF3/GTF parsing with overlap/nearest-gene queries

## Contributing

```bash
cmake -S . -B build-wsl && cmake --build build-wsl -j12
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python
ctest --test-dir build-wsl --output-on-failure
```

## Acknowledgements

Inspired by geneHapR:

> Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9

## License

GPL-3.0-or-later
