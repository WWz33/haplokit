1|# haplokit
2|
3|Command-line haplotype analysis for indexed VCF/BCF data, with a C++ data plane and Python statistics/plotting layer.
4|
5|<!-- README-I18N:START -->
6|
7|**English** | [Chinese](./README.zh-CN.md)
8|
9|<!-- README-I18N:END -->
10|
11|## Table of Contents
12|
13|- [Capabilities](#capabilities)
14|- [Installation](#installation)
15|- [Quick Start](#quick-start)
16|- [Haplotype Workflows](#haplotype-workflows)
17|  - [Interval or Single-Site Haplotype Calling](#interval-or-single-site-haplotype-calling)
18|  - [Gene Annotation and Haplotype Table Figure](#gene-annotation-and-haplotype-table-figure)
19|  - [Population Grouping](#population-grouping)
20|  - [Geographic Distribution](#geographic-distribution)
21|  - [Haplotype Network](#haplotype-network)
22|- [Phenotype Statistics](#phenotype-statistics)
23|- [Other Workflows](#other-workflows)
24|  - [BED Batch Processing](#bed-batch-processing)
25|  - [Approximate Grouping](#approximate-grouping)
26|  - [Sample Subset and Missing-Call Imputation](#sample-subset-and-missing-call-imputation)
27|- [Command Reference](#command-reference)
28|- [Input/Output Formats](#inputoutput-formats)
29|- [Citation](#citation)
30|- [Support](#support)
31|- [License](#license)
32|
33|## Capabilities
34|
35|| Module | Purpose | Typical output |
36|| --- | --- | --- |
37|| `view` | Extract haplotypes from a genomic interval, single site, gene ID, gene list, or BED file | `hapresult.tsv`, `hap_summary.tsv` |
38|| Annotation | Resolve gene selectors and annotate variant positions with gene model context | `gff_ann_summary.tsv`, annotated haplotype table figure |
39|| Population summary | Count haplotypes by population group | population columns in tables and figures |
40|| Geographic map | Draw haplotype composition at sampling locations | map figure with pie charts and count scale |
41|| Network | Build PopART-style haplotype networks with MSN, TCS, or MJN | network figure with population pies and mutation ticks |
42|| `phenotype` | Join haplotypes with numeric traits, run pairwise tests, and draw boxplots | `phenotype_stats.tsv`, phenotype summary TSV, boxplot |
43|
44|## Installation
45|
46|```bash
47|pip install haplokit
48|```
49|
50|<details>
51|<summary><b>From Source (Advanced)</b></summary>
52|
53|#### Requirements
54|
55|- Linux/WSL
56|- Python 3.10+
57|- C++17 compiler
58|- CMake 3.22+
59|- `make`
60|- Native libraries for htslib
61|
62|#### Conda/mamba
63|
64|```bash
65|mamba install -c conda-forge compilers make cmake libcurl zlib bzip2 xz
66|python -m pip install --no-cache-dir haplokit
67|```
68|
69|#### Ubuntu/Debian
70|
71|```bash
72|sudo apt-get update
73|sudo apt-get install -y build-essential make cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
74|python -m pip install --no-cache-dir haplokit
75|```
76|
77|#### From source checkout
78|
79|```bash
80|pip install .
81|```
82|
83|#### Development install
84|
85|```bash
86|pip install -e .
87|```
88|
89|#### Custom backend path
90|
91|If the C++ backend is built outside the package, point the CLI to it:
92|
93|```bash
94|export HAPLOKIT_CPP_BIN=/path/to/haplokit_cpp
95|```
96|
97|#### Troubleshooting linker errors
98|
99|| Error | Install |
100|| --- | --- |
101|| `cannot find -lcurl` | `libcurl` / `libcurl4-openssl-dev` |
102|| `cannot find -lbz2` | `bzip2` / `libbz2-dev` |
103|| `cannot find -llzma` | `xz` / `liblzma-dev` |
104|| `cannot find -lz` | `zlib` / `zlib1g-dev` |
105|
106|</details>
107|
108|## Quick Start
109|
110|```bash
111|haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
112|```
113|
114|This extracts haplotypes from the region `scaffold_1:4300-5000` and writes results to the `out/` directory:
115|
116|| File | Meaning |
117|| --- | --- |
118|| `out/hapresult.tsv` | Haplotype allele pattern and sample accessions |
119|| `out/hap_summary.tsv` | Haplotype counts and frequencies |
120|
121|## Haplotype Workflows
122|
123|### Interval or Single-Site Haplotype Calling
124|
125|```bash
126|# Region mode: group by full allele pattern across the interval
127|haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
128|
129|# Site mode: group by alleles at a single position
130|haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
131|```
132|
133|<details>
134|<summary><b>Details</b></summary>
135|
136|Interval selectors group by the full allele pattern across the region. Single-position selectors automatically use site mode. In strict region mode, samples with heterozygous or missing calls are excluded unless `--impute` is used.
137|
138|</details>
139|
140|### Gene Annotation and Haplotype Table Figure
141|
142|```bash
143|# With gene model annotation
144|haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
145|
146|# Compact table theme
147|haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --table-theme compact --output-file out
148|```
149|
150|<details>
151|<summary><b>Details</b></summary>
152|
153|When a GFF3/GTF file is supplied, the figure draws a pyGenomeTracks-style gene model above the table (backbone/introns, CDS, UTR, and a terminal arrow for strand), with SNP ticks and guide lines connecting each variant to its allele column, plus a CDS/UTR/intron legend. Without `--gff` only the table is drawn. Output also includes `gff_ann_summary.tsv`.
154|
155|`--table-theme` selects the table look: `detailed` (default; square cells with white gridlines) or `compact` (flat wide cells, flush, shorter header). The gene model is shared by both themes.
156|
157|</details>
158|
159|<img src="data/figure/haplotype_table.png" alt="Haplotype summary table" width="800">
160|
161|### Population Grouping
162|
163|```bash
164|haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
165|```
166|
167|<details>
168|<summary><b>Details</b></summary>
169|
170|`popgroup.txt` is a two-column tab-separated file:
171|
172|```text
173|sample  population
174|C1      wild
175|C2      wild
176|C13     landrace
177|```
178|
179|Population groups are shown as per-haplotype count columns in the output table and as grouped counts in the figure.
180|
181|</details>
182|
183|### Geographic Distribution
184|
185|```bash
186|haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
187|```
188|
189|<details>
190|<summary><b>Details</b></summary>
191|
192|Coordinate input is tab-separated:
193|
194|```text
195|ID    longitude  latitude
196|C1    116.40     39.90
197|C2    116.40     39.90
198|```
199|
200|Use `--show-counts` to draw sample-count labels at map pie centers, or `--hide-counts` to keep them hidden.
201|
202|World map example resources are bundled under `data/`:
203|
204|- `sample_world_geo.txt`
205|- `world_countries.shp`, `world_countries.shx`, `world_countries.dbf`
206|- `data/figure/haplotype_map_world.png`
207|
208|</details>
209|
210|<img src="data/figure/haplotype_map_china.png" alt="Haplotype geographic distribution" width="600">
211|
212|<img src="data/figure/haplotype_map_world.png" alt="World haplotype geographic distribution" width="600">
213|
214|### Haplotype Network
215|
216|```bash
217|# TCS network (default)
218|haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out
219|
220|# Median-joining network
221|haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
222|```
223|
224|| Method | Name | Description |
225|| --- | --- | --- |
226|| `tcs` | Statistical parsimony network | Templeton, Crandall & Sing (1992) |
227|| `msn` | Minimum spanning network | Based on Hamming distances |
228|| `mjn` | Median-joining network | Bandelt et al. (1999) |
229|
230|<details>
231|<summary><b>Details</b></summary>
232|
233|Network figures follow PopART conventions: node area reflects haplotype count, pie slices show population composition, edge ticks show mutation steps, and small black vertices indicate inferred intermediates.
234|
235|</details>
236|
237|![Network algorithms comparison - MSN / TCS / MJN](data/figure/haplotype_network_algorithms.png)
238|
239|## Phenotype Statistics
240|
241|### Basic Usage
242|
243|```bash
244|haplokit phenotype \
245|  --hapresult out/hapresult.tsv \
246|  --phenotypes phenotype.csv \
247|  --population popgroup.txt \
248|  --trait yield \
249|  --min-hap-size 5 \
250|  --method welch \
251|  --output yield_stats.tsv \
252|  --summary-output yield_summary.tsv
253|```
254|
255|### With Boxplot
256|
257|```bash
258|haplokit phenotype \
259|  -H data/example_phenotype_haplotypes.tsv \
260|  -P data/example_phenotype.csv \
261|  -p data/popgroup.txt \
262|  -t yield \
263|  -m 4 \
264|  --method welch \
265|  --plot-box \
266|  -F png \
267|  -T "Yield by haplotype and population" \
268|  -b data/figure/phenotype_population_boxplot.png
269|```
270|
271|<img src="data/figure/phenotype_population_boxplot.png" alt="Population-stratified phenotype boxplot" width="900">
272|
273|<details>
274|<summary><b>Statistical scenarios</b></summary>
275|
276|| Scenario | Grouping used for tests | Pairwise comparison reported | Boxplot annotation |
277|| --- | --- | --- | --- |
278|| No population file, multiple haplotypes | `trait x haplotype` | All retained haplotype pairs for each trait | Haplotype-pair significance labels |
279|| Population file, multiple haplotypes | `trait x population x haplotype` | Haplotype pairs inside each population | Within-population haplotype comparisons and between-population comparisons for the same haplotype |
280|| Population file, one retained haplotype | `trait x haplotype x population` | Population pairs within that haplotype | Between-population labels only |
281|| Multiple traits | Each selected trait is analyzed independently | One result block per trait | Plotting requires exactly one `--trait` |
282|| Missing phenotype values | Non-numeric or missing values are ignored per trait | Counts reflect only numeric observations | `effective_n` records the usable sample count |
283|| IQR outlier preprocessing | Optional Tukey IQR k=1.5 within each `trait x population x haplotype` group | Tests use records remaining after outlier removal | Plot uses the same filtered records; summary records removed counts |
284|
285|</details>
286|
287|<details>
288|<summary><b>Pairwise test methods</b></summary>
289|
290|Hypothesis tests use `scipy.stats`.
291|
292|| `--method` | Test | Typical use | P-value adjustment |
293|| --- | --- | --- | --- |
294|| `welch` | Welch two-sample t-test | Default when variances may differ | Bonferroni by default |
295|| `student` | Student two-sample t-test | Similar variance assumption | Bonferroni by default |
296|| `mannwhitney` | Mann-Whitney U test | Non-parametric rank comparison | Bonferroni by default |
297|| `tukey` | Tukey HSD | Multi-group post-hoc comparison | Uses Tukey HSD p-values directly |
298|
299|</details>
300|
301|<details>
302|<summary><b>Outlier preprocessing</b></summary>
303|
304|Use `--remove-outliers` to remove extreme phenotype values before statistics and boxplot rendering:
305|
306|```bash
307|haplokit phenotype \
308|  -H out/hapresult.tsv \
309|  -P phenotype.csv \
310|  -t yield \
311|  --remove-outliers \
312|  -o yield_stats.tsv \
313|  -s yield_summary.tsv
314|```
315|
316|The rule is Tukey IQR with `k=1.5`: values outside `[Q1 - 1.5 x IQR, Q3 + 1.5 x IQR]` are removed. Filtering is performed separately within each `trait x population x haplotype` group. Groups with fewer than four numeric values are left unchanged.
317|
318|Summary output records preprocessing with:
319|
320|| Column | Meaning |
321|| --- | --- |
322|| `raw_count` | Numeric values before outlier removal |
323|| `raw_min`, `raw_max` | Raw group range before removal |
324|| `outlier_removed` | Number of removed values in the summary group |
325|| `outlier_method` | `none` or `iqr` |
326|| `outlier_iqr_k` | IQR multiplier, currently `1.5` when enabled |
327|
328|</details>
329|
330|<details>
331|<summary><b>Phenotype output files</b></summary>
332|
333|| File | Content |
334|| --- | --- |
335|| `phenotype_stats.tsv` | Pairwise comparison rows with group counts, means, standard deviations, ANOVA result, pairwise statistic, raw p-value, adjusted p-value, significance label, and `effective_n` |
336|| summary TSV (`--summary-output`) | Per-trait/per-population/per-haplotype summary statistics, including outlier accounting when preprocessing is enabled |
337|| boxplot (`--plot-box`) | One selected trait visualized with the same filtering, grouping, and comparison logic used by the statistics |
338|
339|</details>
340|
341|## Other Workflows
342|
343|### BED Batch Processing
344|
345|```bash
346|haplokit view in.vcf.gz -R regions.bed --output-file out_batch
347|```
348|
349|<details>
350|<summary><b>Details</b></summary>
351|
352|`regions.bed` requires at least three tab-separated columns:
353|
354|```text
355|chr1  1000  2000
356|chr2  5000  6000
357|```
358|
359|Each BED row is processed independently. Output files are named with a region suffix such as `_chr1_1000_2000`.
360|
361|</details>
362|
363|### Approximate Grouping
364|
365|```bash
366|haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
367|```
368|
369|`--max-diff` clusters haplotypes that differ at no more than the given fraction of variant positions. For example, `--max-diff 0.2` merges haplotypes differing at ≤20% of variant sites.
370|
371|### Sample Subset and Missing-Call Imputation
372|
373|```bash
374|haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
375|```
376|
377|`samples.list` contains one sample ID per line. `--impute` treats missing genotypes as reference (`0/0`) to increase sample retention.
378|
379|### Gene ID and Gene List Selectors
380|
381|```bash
382|# Single gene ID (requires --gff)
383|haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 --output-file out
384|
385|# Gene list file (one gene ID per line)
386|haplokit view in.vcf.gz --gene-list gene_list.txt --gff genes.gff3 --output-file out
387|
388|# With upstream/downstream extensions
389|haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 \
390|  --upstream 2000 --downstream 1000 --strand-aware --output-file out
391|```
392|
393|## Command Reference
394|
395|### `haplokit view`
396|
397|```text
398|haplokit view <input.vcf.gz|input.bcf> (-r <region> | -R <regions.bed> | -t <targets> | -T <targets.txt> | --gene-id <id> | --gene-list <file>) [options]
399|```
400|
401|**Selector options** (exactly one required):
402|
403|| Option | Description |
404|| --- | --- |
405|| `-r, --region` | Single region: `chr:start-end` or `chr:pos` |
406|| `-R, --regions-file` | BED file with multiple regions |
407|| `-t, --targets` | Comma-separated target regions on one chromosome |
408|| `-T, --targets-file` | File with one target region per line (same chromosome) |
409|| `-G, --gene-id` | Single gene ID (requires `--gff/--gff3`) |
410|| `-l, --gene-list` | File with one gene ID per line (requires `--gff/--gff3`) |
411|
412|**Calling options**:
413|
414|| Option | Default | Description |
415|| --- | --- | --- |
416|| `-b, --by` | `auto` | Grouping mode: `auto`, `region`, or `site` |
417|| `-i, --impute` | off | Treat missing genotypes as reference |
418|| `-m, --max-diff` | off | Merge haplotypes with difference ratio ≤ threshold `[0,1]` |
419|
420|**Annotation options**:
421|
422|| Option | Default | Description |
423|| --- | --- | --- |
424|| `-g, --gff3, --gff` | off | GFF3/GTF file for gene selectors and annotation |
425|| `-u, --upstream` | `0` | Upstream bases added to gene selectors |
426|| `-d, --downstream` | `0` | Downstream bases added to gene selectors |
427|| `-a, --strand-aware` | off | Apply upstream/downstream relative to gene strand |
428|
429|**Sample/population options**:
430|
431|| Option | Default | Description |
432|| --- | --- | --- |
433|| `-S, --samples-file` | off | Restrict to sample IDs in a file |
434|| `-p, --population` | off | Sample-to-population table (2 columns) |
435|
436|**Output options**:
437|
438|| Option | Default | Description |
439|| --- | --- | --- |
440|| `-o, --output` | `summary` | JSONL payload mode: `summary` or `detail` |
441|| `-f, --output-format` | `tsv` | Output format: `tsv` or `jsonl` |
442|| `-O, --output-file` | current directory | Output directory or prefix |
443|
444|**Visualization options**:
445|
446|| Option | Default | Description |
447|| --- | --- | --- |
448|| `-P, --plot` | off | Render haplotype table figure |
449|| `-F, --plot-format` | `png` | `png`, `pdf`, `svg`, or `tiff` |
450|| `--table-theme` | `detailed` | Table theme: `detailed` or `compact` |
451|| `-z, --figsize` | auto | Figure size as `WIDTH,HEIGHT` (inches) |
452|| `-e, --geo` | off | Sample coordinates for map plotting |
453|| `--show-counts` | off | Show sample-count labels at map pie centers |
454|| `--hide-counts` | off | Hide sample-count labels (default) |
455|| `-n, --network` | off | Render haplotype network |
456|| `-N, --network-method` | `tcs` | Network algorithm: `tcs`, `msn`, or `mjn` |
457|
458|**Labeling options**:
459|
460|| Option | Default | Description |
461|| --- | --- | --- |
462|| `-H, --hap-prefix` | `Hap` | Haplotype label prefix |
463|| `-D, --hap-pad` | `2` | Zero-padding width for labels |
464|
465|**Notes**:
466|- Exactly one selector is required: `-r`, `-R`, `-t`, `-T`, `--gene-id`, or `--gene-list`
467|- Targets supplied with `-t` or `-T` must all be on the same chromosome
468|- `--gene-id` and `--gene-list` require `--gff/--gff3`
469|
470|### `haplokit phenotype`
471|
472|```text
473|haplokit phenotype -H <hapresult.tsv|sample_hap.tsv> -P <phenotype.tsv|phenotype.csv> [options]
474|```
475|
476|**Input options**:
477|
478|| Option | Default | Description |
479|| --- | --- | --- |
480|| `-H, --hapresult, --haplotypes` | required | `hapresult.tsv` or two-column sample-haplotype table |
481|| `-P, --phenotypes, --phenotype, --pheno-file` | required | Phenotype table; first column is sample ID |
482|| `-p, --population, --pop-group` | off | Sample-to-population table |
483|
484|**Test options**:
485|
486|| Option | Default | Description |
487|| --- | --- | --- |
488|| `-t, --trait` | all numeric traits | Trait to analyze; repeatable for multiple traits |
489|| `-m, --min-hap-size` | `5` | Minimum numeric observations per test group |
490|| `-M, --method` | `welch` | Test method: `welch`, `student`, `mannwhitney`, or `tukey` |
491|| `-a, --adjust` | `bonferroni` | P-value adjustment: `bonferroni` or `none` (non-Tukey) |
492|| `--remove-outliers` | off | Remove Tukey IQR k=1.5 outliers before statistics |
493|
494|**Output options**:
495|
496|| Option | Default | Description |
497|| --- | --- | --- |
498|| `-o, --output` | `phenotype_stats.tsv` | Pairwise statistics TSV |
499|| `-s, --summary-output` | off | Per-haplotype summary TSV |
500|
501|