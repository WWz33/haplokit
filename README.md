# haplokit

Command-line haplotype analysis for indexed VCF/BCF data, with a C++ data plane and Python statistics/plotting layer.

<!-- README-I18N:START -->

**English** | [Chinese](./README.zh-CN.md)

<!-- README-I18N:END -->

`haplokit` is designed for gene- or interval-level haplotype analysis in population genomic studies. It reads indexed VCF/BCF files, builds haplotype tables, annotates variants with GFF3/GTF gene models, summarizes population composition, renders maps and haplotype networks, and tests phenotype differences between haplotype groups.

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

Source builds require Linux/WSL, Python 3.10+, a C++17 compiler, CMake 3.22+, `make`, and native libraries used by the vendored htslib build.

Conda/mamba example:

```bash
mamba install -c conda-forge compilers make cmake libcurl zlib bzip2 xz
python -m pip install --no-cache-dir haplokit
```

Ubuntu/Debian example:

```bash
sudo apt-get update
sudo apt-get install -y build-essential make cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
python -m pip install --no-cache-dir haplokit
```

From a source checkout:

```bash
pip install .
```

For development:

```bash
pip install -e .
```

If the C++ backend is built outside the package, point the CLI to it:

```bash
export HAPLOKIT_CPP_BIN=/path/to/haplokit_cpp
```

Common linker errors map to missing native libraries:

| Error | Install |
| --- | --- |
| `cannot find -lcurl` | `libcurl` / `libcurl4-openssl-dev` |
| `cannot find -lbz2` | `bzip2` / `libbz2-dev` |
| `cannot find -llzma` | `xz` / `liblzma-dev` |
| `cannot find -lz` | `zlib` / `zlib1g-dev` |

## Quick Start

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

Main outputs:

| File | Meaning |
| --- | --- |
| `out/hapresult.tsv` | Haplotype allele pattern and sample accessions |
| `out/hap_summary.tsv` | Haplotype counts and frequencies |

## Haplotype Workflows

### Interval or Single-Site Haplotype Calling

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

Interval selectors group by the full allele pattern across the region. Single-position selectors automatically use site mode. In strict region mode, samples with heterozygous or missing calls are excluded unless `--impute` is used.

### Gene Annotation and Haplotype Table Figure

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
```

The GFF3/GTF file is used for gene selectors and for the functional category strip in the figure. Output includes the table figure and `gff_ann_summary.tsv`.

<img src="data/figure/haplotype_table.png" alt="Haplotype summary table" width="800">

### Population Grouping

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

`popgroup.txt` is a two-column tab-separated file:

```text
sample  population
C1      wild
C2      wild
C13     landrace
```

Population groups are shown as per-haplotype count columns in the output table and as grouped counts in the figure.

### Geographic Distribution

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
```

Coordinate input is tab-separated:

```text
ID    longitude  latitude
C1    116.40     39.90
C2    116.40     39.90
```

Use `--show-counts` to draw sample-count labels at map pie centers, or `--hide-counts` to keep them hidden.

<img src="data/figure/haplotype_map_china.png" alt="Haplotype geographic distribution" width="600">

World map example resources are bundled under `data/`:

- `sample_world_geo.txt`
- `world_countries.shp`, `world_countries.shx`, `world_countries.dbf`
- `data/figure/haplotype_map_world.png`

<img src="data/figure/haplotype_map_world.png" alt="World haplotype geographic distribution" width="600">

### Haplotype Network

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out
haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
```

Supported network methods:

| Method | Meaning |
| --- | --- |
| `msn` | Minimum spanning network |
| `tcs` | Statistical parsimony network |
| `mjn` | Median-joining network |

Network figures follow PopART conventions: node area reflects haplotype count, pie slices show population composition, edge ticks show mutation steps, and small black vertices indicate inferred intermediates.

![Network algorithms comparison - MSN / TCS / MJN](data/figure/haplotype_network_algorithms.png)

## Phenotype Statistics

`haplokit phenotype` joins haplotype assignments with numeric phenotype traits. Input haplotypes can be the `hapresult.tsv` produced by `haplokit view` or a simple two-column sample-to-haplotype table. The phenotype table uses the first column as sample ID and all remaining selected columns as numeric traits.

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

Boxplot example:

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

### Statistical Scenarios

| Scenario | Grouping used for tests | Pairwise comparison reported | Boxplot annotation |
| --- | --- | --- | --- |
| No population file, multiple haplotypes | `trait x haplotype` | All retained haplotype pairs for each trait | Haplotype-pair significance labels |
| Population file, multiple haplotypes | `trait x population x haplotype` | Haplotype pairs inside each population | Within-population haplotype comparisons and between-population comparisons for the same haplotype |
| Population file, one retained haplotype | `trait x haplotype x population` | Population pairs within that haplotype | Between-population labels only |
| Multiple traits | Each selected trait is analyzed independently | One result block per trait | Plotting requires exactly one `--trait` |
| Missing phenotype values | Non-numeric or missing values are ignored per trait | Counts reflect only numeric observations | `effective_n` records the usable sample count |
| IQR outlier preprocessing | Optional Tukey IQR k=1.5 within each `trait x population x haplotype` group | Tests use records remaining after outlier removal | Plot uses the same filtered records; summary records removed counts |

### Pairwise Test Methods

Hypothesis tests use `scipy.stats`.

| `--method` | Test | Typical use | P-value adjustment |
| --- | --- | --- | --- |
| `welch` | Welch two-sample t-test | Default when variances may differ | Bonferroni by default |
| `student` | Student two-sample t-test | Similar variance assumption | Bonferroni by default |
| `mannwhitney` | Mann-Whitney U test | Non-parametric rank comparison | Bonferroni by default |
| `tukey` | Tukey HSD | Multi-group post-hoc comparison | Uses Tukey HSD p-values directly |

### Outlier Preprocessing

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

### Phenotype Output Files

| File | Content |
| --- | --- |
| `phenotype_stats.tsv` | Pairwise comparison rows with group counts, means, standard deviations, ANOVA result, pairwise statistic, raw p-value, adjusted p-value, significance label, and `effective_n` |
| summary TSV (`--summary-output`) | Per-trait/per-population/per-haplotype summary statistics, including outlier accounting when preprocessing is enabled |
| boxplot (`--plot-box`) | One selected trait visualized with the same filtering, grouping, and comparison logic used by the statistics |

## Other Workflows

### BED Batch Processing

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

`regions.bed` requires at least three tab-separated columns:

```text
chr1  1000  2000
chr2  5000  6000
```

Each BED row is processed independently. Output files are named with a region suffix such as `_chr1_1000_2000`.

### Approximate Grouping

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff` clusters haplotypes that differ at no more than the given fraction of variant positions.

### Sample Subset and Missing-Call Imputation

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list` contains one sample ID per line. `--impute` treats missing genotypes as reference (`0/0`) to increase sample retention.

## Command Reference

### `haplokit view`

```text
haplokit view <input.vcf.gz|input.bcf> (-r <region> | -R <regions.bed> | -t <targets> | -T <targets.txt> | --gene-id <id> | --gene-list <file>) [options]
```

| Option | Default | Description |
| --- | --- | --- |
| `-r, --region` | required selector | `chr:start-end` or `chr:pos` |
| `-R, --regions-file` | required selector | BED file |
| `-t, --targets` | required selector | Comma-separated target regions on one chromosome (`chr:pos` or `chr:start-end`) |
| `-T, --targets-file` | required selector | File containing one target region per line on one chromosome; `-` is not accepted |
| `-G, --gene-id` | required selector | Resolve one gene through `--gff/--gff3` |
| `-l, --gene-list` | required selector | One gene ID per line; requires `--gff/--gff3` |
| `-S, --samples-file` | off | Restrict to sample IDs in a file |
| `-b, --by` | `auto` | `auto`, `region`, or `site` |
| `-i, --impute` | off | Treat missing genotypes as reference |
| `-m, --max-diff` | off | Approximate grouping threshold in `[0,1]` |
| `-g, --gff3, --gff` | off | GFF3/GTF file for gene selectors and annotation |
| `-u, --upstream` | `0` | Upstream bases for gene selectors |
| `-d, --downstream` | `0` | Downstream bases for gene selectors |
| `-a, --strand-aware` | off | Apply upstream/downstream relative to gene strand |
| `-o, --output` | `summary` | JSONL payload mode: `summary` or `detail` |
| `-f, --output-format` | `tsv` | Output format: `tsv` or `jsonl` |
| `-O, --output-file` | current directory | Output directory, prefix, or JSONL file |
| `-P, --plot` | off | Render haplotype table figure |
| `-F, --plot-format` | `png` | `png`, `pdf`, `svg`, or `tiff` |
| `-z, --figsize` | auto | Figure size as `WIDTH,HEIGHT` |
| `-p, --population` | off | Sample-to-population table |
| `-e, --geo` | off | Sample coordinates for map plotting |
| `--show-counts`, `--hide-counts` | hidden | Control map count labels |
| `-n, --network` | off | Render haplotype network |
| `-N, --network-method` | `tcs` | `tcs`, `msn`, or `mjn` |
| `-H, --hap-prefix` | `Hap` | Haplotype label prefix |
| `-D, --hap-pad` | `2` | Zero-padding width for labels |

Exactly one selector is required: `-r`, `-R`, `-t`, `-T`, `--gene-id`, or `--gene-list`.
Targets supplied with `-t` or `-T` must all be on the same chromosome.

### `haplokit phenotype`

```text
haplokit phenotype -H <hapresult.tsv|sample_hap.tsv> -P <phenotype.tsv|phenotype.csv> [options]
```

| Option | Default | Description |
| --- | --- | --- |
| `-H, --hapresult, --haplotypes` | required | `hapresult.tsv` or two-column sample-haplotype table |
| `-P, --phenotypes, --phenotype, --pheno-file` | required | Phenotype table; first column is sample ID |
| `-p, --population, --pop-group` | off | Sample-to-population table |
| `-t, --trait` | all numeric traits | Trait to analyze; repeatable |
| `-m, --min-hap-size` | `5` | Minimum numeric observations per test group |
| `-M, --method` | `welch` | `welch`, `student`, `mannwhitney`, or `tukey` |
| `-a, --adjust` | `bonferroni` | `bonferroni` or `none` for non-Tukey tests |
| `--remove-outliers` | off | Remove Tukey IQR k=1.5 outliers before statistics and plot |
| `-o, --output` | `phenotype_stats.tsv` | Pairwise statistics TSV |
| `-s, --summary-output` | off | Per-haplotype summary TSV |
| `-B, --plot-box` | off | Render phenotype boxplot |
| `-b, --box-output` | `phenotype_box.png` | Boxplot output path |
| `-F, --plot-format` | output suffix | `png`, `pdf`, `svg`, or `tiff` |
| `-z, --figsize` | auto | Boxplot size as `WIDTH,HEIGHT` |
| `-T, --title` | auto | Boxplot title |
| `-c, --comparison` | all available pairs | Explicit haplotype pair such as `Hap01,Hap02`; repeatable |
| `-d, --delimiter` | `auto` | Haplotype input delimiter: `auto`, `tab`, or `comma` |
| `-D, --phenotype-delimiter` | `auto` | Phenotype input delimiter |
| `-G, --population-delimiter` | `auto` | Population input delimiter |

`--plot-box` requires exactly one selected trait.

## Backend

Backend binary: `haplokit_cpp`.

Backend discovery order:

1. `HAPLOKIT_CPP_BIN`
2. Packaged binary: `haplokit/_bin/haplokit_cpp`
3. Local builds: `build-wsl/haplokit_cpp`, `build/haplokit_cpp`, `build-haplokit-python/haplokit_cpp`
4. Source-tree CMake build fallback

Native components:

- [htslib](https://github.com/samtools/htslib) for indexed VCF/BCF reading
- [gffsub](https://github.com/WWz33/gffsub) for GFF3/GTF parsing and interval queries

## Development

```bash
cmake -S . -B build-wsl && cmake --build build-wsl -j12
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python
ctest --test-dir build-wsl --output-on-failure
```

## References

`haplokit` is inspired by geneHapR:

> Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9

Network visualization follows PopART conventions:

> Leigh, J. W. & Bryant, D. popart: full-feature software for haplotype network construction. Methods in Ecology and Evolution 6, 1110-1116 (2015). https://doi.org/10.1111/2041-210X.12410

## License

GPL-3.0-or-later
