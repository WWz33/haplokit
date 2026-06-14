# haplokit

Command-line haplotype analysis for indexed VCF/BCF data, with a C++ data plane and Python statistics/plotting layer.

<!-- README-I18N:START -->

**English** | [Chinese](./README.zh-CN.md)

<!-- README-I18N:END -->

## Table of Contents

- [Capabilities](#capabilities)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Haplotype Workflows](#haplotype-workflows)
  - [Interval or Single-Site Haplotype Calling](#interval-or-single-site-haplotype-calling)
  - [Gene Annotation and Haplotype Table Figure](#gene-annotation-and-haplotype-table-figure)
  - [Population Grouping](#population-grouping)
  - [Geographic Distribution](#geographic-distribution)
  - [Haplotype Network](#haplotype-network)
- [Phenotype Statistics](#phenotype-statistics)
- [Other Workflows](#other-workflows)
  - [BED Batch Processing](#bed-batch-processing)
  - [Approximate Grouping](#approximate-grouping)
  - [Sample Subset and Missing-Call Imputation](#sample-subset-and-missing-call-imputation)
- [Command Reference](#command-reference)
- [Input/Output Formats](#inputoutput-formats)
- [Citation](#citation)
- [Support](#support)
- [License](#license)

## Capabilities

| Module | Purpose | Typical output |
| --- | --- | --- |
| `view` | Extract haplotypes from a genomic interval, single site, gene ID, gene list, or BED file | `hapresult.tsv`, `hap_summary.tsv` |
| Annotation | Resolve gene selectors and annotate variant positions with gene model context | `gff_ann_summary.tsv`, annotated haplotype table figure |
| Population summary | Count haplotypes by population group | population columns in tables and figures |
| Geographic map | Draw haplotype composition at sampling locations | map figure with pie charts and count scale |
| Network | Build PopART-style haplotype networks with MSN, TCS, or MJN | network figure with population pies and mutation ticks |
| `phenotype` | Join haplotypes with numeric traits, run pairwise tests, and draw boxplots | `phenotype_stats.tsv`, phenotype summary TSV, boxplot |

## Installation

```bash
pip install haplokit
```

<details>
<summary><b>From Source (Advanced)</b></summary>

#### Requirements

- Linux/WSL
- Python 3.10+
- C++17 compiler
- CMake 3.22+
- `make`
- Native libraries for htslib

#### Conda/mamba

```bash
mamba install -c conda-forge compilers make cmake libcurl zlib bzip2 xz
python -m pip install --no-cache-dir haplokit
```

#### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install -y build-essential make cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
python -m pip install --no-cache-dir haplokit
```

#### From source checkout

```bash
pip install .
```

#### Development install

```bash
pip install -e .
```

#### Custom backend path

If the C++ backend is built outside the package, point the CLI to it:

```bash
export HAPLOKIT_CPP_BIN=/path/to/haplokit_cpp
```

#### Troubleshooting linker errors

| Error | Install |
| --- | --- |
| `cannot find -lcurl` | `libcurl` / `libcurl4-openssl-dev` |
| `cannot find -lbz2` | `bzip2` / `libbz2-dev` |
| `cannot find -llzma` | `xz` / `liblzma-dev` |
| `cannot find -lz` | `zlib` / `zlib1g-dev` |

</details>

## Quick Start

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

This extracts haplotypes from the region `scaffold_1:4300-5000` and writes results to the `out/` directory:

| File | Meaning |
| --- | --- |
| `out/hapresult.tsv` | Haplotype allele pattern and sample accessions |
| `out/hap_summary.tsv` | Haplotype counts and frequencies |

## Haplotype Workflows

### Interval or Single-Site Haplotype Calling

```bash
# Region mode: group by full allele pattern across the interval
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out

# Site mode: group by alleles at a single position
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

<details>
<summary><b>Details</b></summary>

Interval selectors group by the full allele pattern across the region. Single-position selectors automatically use site mode. In strict region mode, samples with heterozygous or missing calls are excluded unless `--impute` is used.

</details>

### Gene Annotation and Haplotype Table Figure

```bash
# With gene model annotation
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out

# Compact table theme
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --table-theme compact --output-file out
```

<details>
<summary><b>Details</b></summary>

When a GFF3/GTF file is supplied, the figure draws a pyGenomeTracks-style gene model above the table (backbone/introns, CDS, UTR, and a terminal arrow for strand), with SNP ticks and guide lines connecting each variant to its allele column, plus a CDS/UTR/intron legend. Without `--gff` only the table is drawn. Output also includes `gff_ann_summary.tsv`.

`--table-theme` selects the table look: `detailed` (default; square cells with white gridlines) or `compact` (flat wide cells, flush, shorter header). The gene model is shared by both themes.

</details>

<img src="data/figure/haplotype_table.png" alt="Haplotype summary table" width="800">

### Population Grouping

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

<details>
<summary><b>Details</b></summary>

`popgroup.txt` is a two-column tab-separated file:

```text
sample  population
C1      wild
C2      wild
C13     landrace
```

Population groups are shown as per-haplotype count columns in the output table and as grouped counts in the figure.

</details>

### Geographic Distribution

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
```

<details>
<summary><b>Details</b></summary>

Coordinate input is tab-separated:

```text
ID    longitude  latitude
C1    116.40     39.90
C2    116.40     39.90
```

Use `--show-counts` to draw sample-count labels at map pie centers, or `--hide-counts` to keep them hidden.

World map example resources are bundled under `data/`:

- `sample_world_geo.txt`
- `world_countries.shp`, `world_countries.shx`, `world_countries.dbf`
- `data/figure/haplotype_map_world.png`

</details>

<img src="data/figure/haplotype_map_china.png" alt="Haplotype geographic distribution" width="600">

<img src="data/figure/haplotype_map_world.png" alt="World haplotype geographic distribution" width="600">

### Haplotype Network

```bash
# TCS network (default)
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out

# Median-joining network
haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
```

| Method | Name | Description |
| --- | --- | --- |
| `tcs` | Statistical parsimony network | Templeton, Crandall & Sing (1992) |
| `msn` | Minimum spanning network | Based on Hamming distances |
| `mjn` | Median-joining network | Bandelt et al. (1999) |

<details>
<summary><b>Details</b></summary>

Network figures follow PopART conventions: node area reflects haplotype count, pie slices show population composition, edge ticks show mutation steps, and small black vertices indicate inferred intermediates.

</details>

![Network algorithms comparison - MSN / TCS / MJN](data/figure/haplotype_network_algorithms.png)

## Phenotype Statistics

### Basic Usage

```bash
haplokit phenotype \
  --hapresult out/hapresult.tsv \
  --phenotypes phenotype.csv \
  --population popgroup.txt \
  --trait yield \
  --min-hap-size 5 \
  --method welch \
  --output yield_stats.tsv \
  --summary-output yield_summary.tsv
```

### With Boxplot

```bash
haplokit phenotype \
  -H data/example_phenotype_haplotypes.tsv \
  -P data/example_phenotype.csv \
  -p data/popgroup.txt \
  -t yield \
  -m 4 \
  --method welch \
  --plot-box \
  -F png \
  -T "Yield by haplotype and population" \
  -b data/figure/phenotype_population_boxplot.png
```

<img src="data/figure/phenotype_population_boxplot.png" alt="Population-stratified phenotype boxplot" width="900">

<details>
<summary><b>Statistical scenarios</b></summary>

| Scenario | Grouping used for tests | Pairwise comparison reported | Boxplot annotation |
| --- | --- | --- | --- |
| No population file, multiple haplotypes | `trait x haplotype` | All retained haplotype pairs for each trait | Haplotype-pair significance labels |
| Population file, multiple haplotypes | `trait x population x haplotype` | Haplotype pairs inside each population | Within-population haplotype comparisons and between-population comparisons for the same haplotype |
| Population file, one retained haplotype | `trait x haplotype x population` | Population pairs within that haplotype | Between-population labels only |
| Multiple traits | Each selected trait is analyzed independently | One result block per trait | Plotting requires exactly one `--trait` |
| Missing phenotype values | Non-numeric or missing values are ignored per trait | Counts reflect only numeric observations | `effective_n` records the usable sample count |
| IQR outlier preprocessing | Optional Tukey IQR k=1.5 within each `trait x population x haplotype` group | Tests use records remaining after outlier removal | Plot uses the same filtered records; summary records removed counts |

</details>

<details>
<summary><b>Pairwise test methods</b></summary>

Hypothesis tests use `scipy.stats`.

| `--method` | Test | Typical use | P-value adjustment |
| --- | --- | --- | --- |
| `welch` | Welch two-sample t-test | Default when variances may differ | Bonferroni by default |
| `student` | Student two-sample t-test | Similar variance assumption | Bonferroni by default |
| `mannwhitney` | Mann-Whitney U test | Non-parametric rank comparison | Bonferroni by default |
| `tukey` | Tukey HSD | Multi-group post-hoc comparison | Uses Tukey HSD p-values directly |

</details>

<details>
<summary><b>Outlier preprocessing</b></summary>

Use `--remove-outliers` to remove extreme phenotype values before statistics and boxplot rendering:

```bash
haplokit phenotype \
  -H out/hapresult.tsv \
  -P phenotype.csv \
  -t yield \
  --remove-outliers \
  -o yield_stats.tsv \
  -s yield_summary.tsv
```

The rule is Tukey IQR with `k=1.5`: values outside `[Q1 - 1.5 x IQR, Q3 + 1.5 x IQR]` are removed. Filtering is performed separately within each `trait x population x haplotype` group. Groups with fewer than four numeric values are left unchanged.

Summary output records preprocessing with:

| Column | Meaning |
| --- | --- |
| `raw_count` | Numeric values before outlier removal |
| `raw_min`, `raw_max` | Raw group range before removal |
| `outlier_removed` | Number of removed values in the summary group |
| `outlier_method` | `none` or `iqr` |
| `outlier_iqr_k` | IQR multiplier, currently `1.5` when enabled |

</details>

<details>
<summary><b>Phenotype output files</b></summary>

| File | Content |
| --- | --- |
| `phenotype_stats.tsv` | Pairwise comparison rows with group counts, means, standard deviations, ANOVA result, pairwise statistic, raw p-value, adjusted p-value, significance label, and `effective_n` |
| summary TSV (`--summary-output`) | Per-trait/per-population/per-haplotype summary statistics, including outlier accounting when preprocessing is enabled |
| boxplot (`--plot-box`) | One selected trait visualized with the same filtering, grouping, and comparison logic used by the statistics |

</details>

## Other Workflows

### BED Batch Processing

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

<details>
<summary><b>Details</b></summary>

`regions.bed` requires at least three tab-separated columns:

```text
chr1  1000  2000
chr2  5000  6000
```

Each BED row is processed independently. Output files are named with a region suffix such as `_chr1_1000_2000`.

</details>

### Approximate Grouping

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff` clusters haplotypes that differ at no more than the given fraction of variant positions. For example, `--max-diff 0.2` merges haplotypes differing at ≤20% of variant sites.

### Sample Subset and Missing-Call Imputation

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list` contains one sample ID per line. `--impute` treats missing genotypes as reference (`0/0`) to increase sample retention.

### Gene ID and Gene List Selectors

```bash
# Single gene ID (requires --gff)
haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 --output-file out

# Gene list file (one gene ID per line)
haplokit view in.vcf.gz --gene-list gene_list.txt --gff genes.gff3 --output-file out

# With upstream/downstream extensions
haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 \
  --upstream 2000 --downstream 1000 --strand-aware --output-file out
```

## Command Reference

### `haplokit view`

```text
haplokit view <input.vcf.gz|input.bcf> (-r <region> | -R <regions.bed> | -t <targets> | -T <targets.txt> | --gene-id <id> | --gene-list <file>) [options]
```

**Selector options** (exactly one required):

| Option | Description |
| --- | --- |
| `-r, --region` | Single region: `chr:start-end` or `chr:pos` |
| `-R, --regions-file` | BED file with multiple regions |
| `-t, --targets` | Comma-separated target regions on one chromosome |
| `-T, --targets-file` | File with one target region per line (same chromosome) |
| `-G, --gene-id` | Single gene ID (requires `--gff/--gff3`) |
| `-l, --gene-list` | File with one gene ID per line (requires `--gff/--gff3`) |

**Calling options**:

| Option | Default | Description |
| --- | --- | --- |
| `-b, --by` | `auto` | Grouping mode: `auto`, `region`, or `site` |
| `-i, --impute` | off | Treat missing genotypes as reference |
| `-m, --max-diff` | off | Merge haplotypes with difference ratio ≤ threshold `[0,1]` |

**Annotation options**:

| Option | Default | Description |
| --- | --- | --- |
| `-g, --gff3, --gff` | off | GFF3/GTF file for gene selectors and annotation |
| `-u, --upstream` | `0` | Upstream bases added to gene selectors |
| `-d, --downstream` | `0` | Downstream bases added to gene selectors |
| `-a, --strand-aware` | off | Apply upstream/downstream relative to gene strand |

**Sample/population options**:

| Option | Default | Description |
| --- | --- | --- |
| `-S, --samples-file` | off | Restrict to sample IDs in a file |
| `-p, --population` | off | Sample-to-population table (2 columns) |

**Output options**:

| Option | Default | Description |
| --- | --- | --- |
| `-o, --output` | `summary` | JSONL payload mode: `summary` or `detail` |
| `-f, --output-format` | `tsv` | Output format: `tsv` or `jsonl` |
| `-O, --output-file` | current directory | Output directory or prefix |

**Visualization options**:

| Option | Default | Description |
| --- | --- | --- |
| `-P, --plot` | off | Render haplotype table figure |
| `-F, --plot-format` | `png` | `png`, `pdf`, `svg`, or `tiff` |
| `--table-theme` | `detailed` | Table theme: `detailed` or `compact` |
| `-z, --figsize` | auto | Figure size as `WIDTH,HEIGHT` (inches) |
| `-e, --geo` | off | Sample coordinates for map plotting |
| `--show-counts` | off | Show sample-count labels at map pie centers |
| `--hide-counts` | off | Hide sample-count labels (default) |
| `-n, --network` | off | Render haplotype network |
| `-N, --network-method` | `tcs` | Network algorithm: `tcs`, `msn`, or `mjn` |

**Labeling options**:

| Option | Default | Description |
| --- | --- | --- |
| `-H, --hap-prefix` | `Hap` | Haplotype label prefix |
| `-D, --hap-pad` | `2` | Zero-padding width for labels |

**Notes**:
- Exactly one selector is required: `-r`, `-R`, `-t`, `-T`, `--gene-id`, or `--gene-list`
- Targets supplied with `-t` or `-T` must all be on the same chromosome
- `--gene-id` and `--gene-list` require `--gff/--gff3`

### `haplokit phenotype`

```text
haplokit phenotype -H <hapresult.tsv|sample_hap.tsv> -P <phenotype.tsv|phenotype.csv> [options]
```

**Input options**:

| Option | Default | Description |
| --- | --- | --- |
| `-H, --hapresult, --haplotypes` | required | `hapresult.tsv` or two-column sample-haplotype table |
| `-P, --phenotypes, --phenotype, --pheno-file` | required | Phenotype table; first column is sample ID |
| `-p, --population, --pop-group` | off | Sample-to-population table |

**Test options**:

| Option | Default | Description |
| --- | --- | --- |
| `-t, --trait` | all numeric traits | Trait to analyze; repeatable for multiple traits |
| `-m, --min-hap-size` | `5` | Minimum numeric observations per test group |
| `-M, --method` | `welch` | Test method: `welch`, `student`, `mannwhitney`, or `tukey` |
| `-a, --adjust` | `bonferroni` | P-value adjustment: `bonferroni` or `none` (non-Tukey) |
| `--remove-outliers` | off | Remove Tukey IQR k=1.5 outliers before statistics |

**Output options**:

| Option | Default | Description |
| --- | --- | --- |
| `-o, --output` | `phenotype_stats.tsv` | Pairwise statistics TSV |
| `-s, --summary-output` | off | Per-haplotype summary TSV |

**Boxplot options**:

| Option | Default | Description |
| --- | --- | --- |
| `-B, --plot-box` | off | Render phenotype boxplot |
| `-b, --box-output` | `phenotype_box.png` | Boxplot output path |
| `-F, --plot-format` | output suffix | `png`, `pdf`, `svg`, or `tiff` |
| `-z, --figsize` | auto | Boxplot size as `WIDTH,HEIGHT` (inches) |
| `-T, --title` | auto | Boxplot title |
| `-c, --comparison` | all available pairs | Haplotype pair to annotate (e.g., `Hap01,Hap02`); repeatable |

**Delimiter options**:

| Option | Default | Description |
| --- | --- | --- |
| `-d, --delimiter` | `auto` | Haplotype input delimiter: `auto`, `tab`, or `comma` |
| `-D, --phenotype-delimiter` | `auto` | Phenotype input delimiter |
| `-G, --population-delimiter` | `auto` | Population input delimiter |

**Notes**:
- `--plot-box` requires exactly one selected trait (`-t`)

<details>
<summary><b>Backend</b></summary>

Backend binary: `haplokit_cpp`.

Backend discovery order:

1. `HAPLOKIT_CPP_BIN` environment variable
2. Packaged binary: `haplokit/_bin/haplokit_cpp`
3. Local builds: `build-wsl/haplokit_cpp`, `build/haplokit_cpp`, `build-haplokit-python/haplokit_cpp`
4. Source-tree CMake build fallback

Native components:

- [htslib](https://github.com/samtools/htslib) for indexed VCF/BCF reading
- [gffsub](https://github.com/WWz33/gffsub) for GFF3/GTF parsing and interval queries

</details>

<details>
<summary><b>Development</b></summary>

```bash
# Build C++ backend
cmake -S . -B build-wsl && cmake --build build-wsl -j12

# Run Python tests
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python

# Run C++ tests
ctest --test-dir build-wsl --output-on-failure
```

</details>

## Input/Output Formats

### Input Formats

| Format | Extension | Description | Notes |
| --- | --- | --- | --- |
| VCF | `.vcf.gz` | Variant Call Format | Must be indexed (`.tbi` or `.csi`) |
| BCF | `.bcf` | Binary VCF | Must be indexed |
| BED | `.bed` | Genomic intervals | 3+ tab-separated columns: chr, start, end |
| GFF3 | `.gff3` | Gene annotations | Standard GFF3 format |
| GTF | `.gtf` | Gene transfer format | Standard GTF format |
| Population map | `.txt` | Sample-to-population mapping | 2 tab-separated columns: sample ID, population name |
| Phenotype table | `.csv`, `.tsv` | Trait values | First column: sample ID; remaining columns: numeric traits |
| Gene list | `.txt` | Gene IDs | One gene ID per line |

### Output Formats

#### `hapresult.tsv`

Main output from `haplokit view`. Contains haplotype definitions and sample assignments.

| Column | Type | Description |
| --- | --- | --- |
| `Hap` | string | Haplotype label (e.g., Hap01, Hap02) |
| `Count` | integer | Number of samples with this haplotype |
| `Freq` | float | Frequency (0-1) |
| `State` | string | Allele pattern, pipe-delimited (e.g., `0\|0\|1\|0`) |
| `Samples` | string | Comma-separated sample IDs |
| `{Pop1}`, `{Pop2}`, ... | integer | Per-population counts (only if `-p` provided) |

#### `hap_summary.tsv`

Simplified summary of haplotype counts.

| Column | Type | Description |
| --- | --- | --- |
| `Hap` | string | Haplotype label |
| `Count` | integer | Total count |
| `Freq` | float | Frequency |
| `Samples` | string | Sample list |

#### `gff_ann_summary.tsv`

Gene annotation summary (only if `--gff` provided).

| Column | Type | Description |
| --- | --- | --- |
| `Position` | integer | Variant position |
| `Feature` | string | Feature type (CDS, 5'UTR, 3'UTR, intron, upstream, downstream, intergenic) |
| `Gene` | string | Gene ID |
| `Distance` | integer | Distance to feature (0 if within) |

#### `phenotype_stats.tsv`

Pairwise statistical comparisons from `haplokit phenotype`.

| Column | Type | Description |
| --- | --- | --- |
| `trait` | string | Trait name |
| `group1`, `group2` | string | Compared groups (haplotype or population labels) |
| `n1`, `n2` | integer | Sample counts per group |
| `mean1`, `mean2` | float | Group means |
| `std1`, `std2` | float | Group standard deviations |
| `statistic` | float | Test statistic |
| `pvalue` | float | Raw p-value |
| `pvalue_adjusted` | float | Adjusted p-value (Bonferroni or Tukey) |
| `significance` | string | Significance label: `***`, `**`, `*`, `ns` |
| `effective_n` | integer | Usable sample count (excludes missing values) |

## Citation

If you use haplokit in your research, please cite:

> WWz33. (2024). haplokit: Command-line haplotype analysis for indexed VCF/BCF data. https://github.com/WWz33/haplokit

`haplokit` is inspired by geneHapR:

> Zhang, R., Jia, G. & Diao, X. (2023). geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199. https://doi.org/10.1186/s12859-023-05318-9

Network visualization follows PopART conventions:

> Leigh, J. W. & Bryant, D. (2015). popart: full-feature software for haplotype network construction. Methods in Ecology and Evolution 6, 1110-1116. https://doi.org/10.1111/2041-210X.12410

## Support

- **Issues**: [GitHub Issues](https://github.com/WWz33/haplokit/issues)
- **Discussions**: [GitHub Discussions](https://github.com/WWz33/haplokit/discussions)

## License

GPL-3.0-or-later
