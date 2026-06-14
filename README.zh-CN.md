1|# haplokit
2|
3|面向 indexed VCF/BCF 的命令行单倍型分析工具，使用 C++ 处理数据平面，使用 Python 完成统计分析与绘图。
4|
5|<!-- README-I18N:START -->
6|
7|[English](./README.md) | **中文**
8|
9|<!-- README-I18N:END -->
10|
11|## 目录
12|
13|- [功能概览](#功能概览)
14|- [安装](#安装)
15|- [快速开始](#快速开始)
16|- [单倍型分析流程](#单倍型分析流程)
17|  - [区间或单个位点单倍型识别](#区间或单个位点单倍型识别)
18|  - [基因注释与单倍型表图](#基因注释与单倍型表图)
19|  - [群体分组](#群体分组)
20|  - [地理分布](#地理分布)
21|  - [单倍型网络](#单倍型网络)
22|- [表型统计模块](#表型统计模块)
23|- [其他流程](#其他流程)
24|  - [BED 批量处理](#bed-批量处理)
25|  - [近似分组](#近似分组)
26|  - [样本子集与缺失填补](#样本子集与缺失填补)
27|- [参数速查](#参数速查)
28|- [输入/输出格式](#输入输出格式)
29|- [引用](#引用)
30|- [支持](#支持)
31|- [许可证](#许可证)
32|
33|## 功能概览
34|
35|| 模块 | 用途 | 典型输出 |
36|| --- | --- | --- |
37|| `view` | 从区间、单个位点、基因 ID、基因列表或 BED 文件提取单倍型 | `hapresult.tsv`, `hap_summary.tsv` |
38|| 基因注释 | 解析 gene selector，并为变异位置添加基因结构上下文 | `gff_ann_summary.tsv`, 注释版单倍型表图 |
39|| 群体统计 | 按群体统计每个单倍型的样本数 | 表格和图中的 population columns |
40|| 地理分布图 | 在采样地点绘制单倍型组成 | 带饼图和样本数比例尺的地图 |
41|| 单倍型网络 | 使用 MSN、TCS 或 MJN 构建 PopART 风格网络 | 带群体饼图和突变刻度的 network figure |
42|| `phenotype` | 连接单倍型与数值表型，输出检验结果和箱线图 | `phenotype_stats.tsv`, 表型 summary TSV, boxplot |
43|
44|## 安装
45|
46|```bash
47|pip install haplokit
48|```
49|
50|<details>
51|<summary><b>源码构建（高级）</b></summary>
52|
53|#### 系统要求
54|
55|- Linux/WSL
56|- Python 3.10+
57|- C++17 编译器
58|- CMake 3.22+
59|- `make`
60|- htslib 所需的本地库
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
77|#### 从源码目录安装
78|
79|```bash
80|pip install .
81|```
82|
83|#### 开发模式安装
84|
85|```bash
86|pip install -e .
87|```
88|
89|#### 自定义后端路径
90|
91|如果 C++ 后端构建在其他位置，可以显式指定：
92|
93|```bash
94|export HAPLOKIT_CPP_BIN=/path/to/haplokit_cpp
95|```
96|
97|#### 常见链接错误排查
98|
99|| 错误 | 需要安装 |
100|| --- | --- |
101|| `cannot find -lcurl` | `libcurl` / `libcurl4-openssl-dev` |
102|| `cannot find -lbz2` | `bzip2` / `libbz2-dev` |
103|| `cannot find -llzma` | `xz` / `liblzma-dev` |
104|| `cannot find -lz` | `zlib` / `zlib1g-dev` |
105|
106|</details>
107|
108|## 快速开始
109|
110|```bash
111|haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
112|```
113|
114|从区间 `scaffold_1:4300-5000` 提取单倍型，结果写入 `out/` 目录：
115|
116|| 文件 | 含义 |
117|| --- | --- |
118|| `out/hapresult.tsv` | 单倍型等位基因模式与样本列表 |
119|| `out/hap_summary.tsv` | 单倍型计数和频率 |
120|
121|## 单倍型分析流程
122|
123|### 区间或单个位点单倍型识别
124|
125|```bash
126|# 区间模式：按整个区域内的等位基因模式分组
127|haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
128|
129|# 单点模式：按单个位点的等位基因分组
130|haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
131|```
132|
133|<details>
134|<summary><b>详细说明</b></summary>
135|
136|区间 selector 会按整个区域内的等位基因模式分组。单点 selector 会自动进入 site mode。严格区间模式下，带杂合或缺失调用的样本会被排除；如需保留缺失样本，可以使用 `--impute`。
137|
138|</details>
139|
140|### 基因注释与单倍型表图
141|
142|```bash
143|# 带基因模型注释
144|haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
145|
146|# 紧凑表格样式
147|haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --table-theme compact --output-file out
148|```
149|
150|<details>
151|<summary><b>详细说明</b></summary>
152|
153|提供 GFF3/GTF 文件时，图中会在表格上方绘制 pyGenomeTracks 风格的基因模型（backbone/内含子、CDS、UTR，末端箭头表示链方向），并用 SNP 刻度和引导线把每个变异连到对应的等位基因列，同时附带 CDS/UTR/intron 图例。不提供 `--gff` 时只绘制表格。输出还包括 `gff_ann_summary.tsv`。
154|
155|`--table-theme` 选择表格样式：`detailed`（默认，方形单元格、白色网格线）或 `compact`（扁宽单元格、无边贴合、更矮的表头）。基因模型在两种样式下通用。
156|
157|</details>
158|
159|<img src="data/figure/haplotype_table.png" alt="单倍型汇总表" width="800">
160|
161|### 群体分组
162|
163|```bash
164|haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
165|```
166|
167|<details>
168|<summary><b>详细说明</b></summary>
169|
170|`popgroup.txt` 是两列 tab 分隔文件：
171|
172|```text
173|sample  population
174|C1      wild
175|C2      wild
176|C13     landrace
177|```
178|
179|群体信息会进入输出表格的计数列，也会进入图中的群体统计。
180|
181|</details>
182|
183|### 地理分布
184|
185|```bash
186|haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
187|```
188|
189|<details>
190|<summary><b>详细说明</b></summary>
191|
192|坐标文件为 tab 分隔：
193|
194|```text
195|ID    longitude  latitude
196|C1    116.40     39.90
197|C2    116.40     39.90
198|```
199|
200|使用 `--show-counts` 在地图饼图中心显示样本数，或使用 `--hide-counts` 显式隐藏。
201|
202|`data/` 中附带世界地图示例资源：
203|
204|- `sample_world_geo.txt`
205|- `world_countries.shp`, `world_countries.shx`, `world_countries.dbf`
206|- `data/figure/haplotype_map_world.png`
207|
208|</details>
209|
210|<img src="data/figure/haplotype_map_china.png" alt="单倍型地理分布图" width="600">
211|
212|<img src="data/figure/haplotype_map_world.png" alt="世界单倍型地理分布图" width="600">
213|
214|### 单倍型网络
215|
216|```bash
217|# TCS 网络（默认）
218|haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out
219|
220|# 中值连接网络
221|haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
222|```
223|
224|| 方法 | 名称 | 说明 |
225|| --- | --- | --- |
226|| `tcs` | 统计简约性网络 | Templeton, Crandall & Sing (1992) |
227|| `msn` | 最小生成树网络 | 基于 Hamming 距离 |
228|| `mjn` | 中值连接网络 | Bandelt et al. (1999) |
229|
230|<details>
231|<summary><b>详细说明</b></summary>
232|
233|网络图遵循 PopART 风格：节点面积表示单倍型样本数，饼图扇区表示群体组成，边上的刻度表示突变步数，小黑点表示推断出的中间节点。
234|
235|</details>
236|
237|![网络算法对比 - MSN / TCS / MJN](data/figure/haplotype_network_algorithms.png)
238|
239|## 表型统计模块
240|
241|### 基本用法
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
255|### 带箱线图
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
271|<img src="data/figure/phenotype_population_boxplot.png" alt="群体分层表型箱线图" width="900">
272|
273|<details>
274|<summary><b>统计场景</b></summary>
275|
276|| 场景 | 检验分组 | result 中的两两比较 | 箱线图显著性标注 |
277|| --- | --- | --- | --- |
278|| 不提供群体文件，保留多个单倍型 | `trait x haplotype` | 每个性状内所有保留单倍型两两比较 | 单倍型之间的显著性 |
279|| 提供群体文件，保留多个单倍型 | `trait x population x haplotype` | 每个群体内部的单倍型两两比较 | 群体内单倍型比较；同一单倍型的群体间比较 |
280|| 提供群体文件，只保留一个单倍型 | `trait x haplotype x population` | 该单倍型内部的群体两两比较 | 仅显示群体间比较 |
281|| 多个性状 | 每个性状独立分析 | 每个性状输出独立结果块 | 绘图必须用 `--trait` 选择一个性状 |
282|| 表型缺失值 | 按性状忽略非数值和缺失值 | 计数只包含数值样本 | `effective_n` 记录进入该分层的有效样本数 |
283|| IQR 极值预处理 | 可选；在每个 `trait x population x haplotype` 内执行 Tukey IQR k=1.5 | 检验使用删除极值后的样本 | 图使用同一批过滤后的样本；summary 记录删除数量 |
284|
285|</details>
286|
287|<details>
288|<summary><b>两两检验方法</b></summary>
289|
290|假设检验使用 `scipy.stats`。
291|
292|| `--method` | 检验 | 适用场景 | P 值校正 |
293|| --- | --- | --- | --- |
294|| `welch` | Welch two-sample t-test | 默认；不假设方差相等 | 非 Tukey 检验默认 Bonferroni |
295|| `student` | Student two-sample t-test | 可接受方差相等假设时 | 非 Tukey 检验默认 Bonferroni |
296|| `mannwhitney` | Mann-Whitney U test | 非参数秩检验 | 非 Tukey 检验默认 Bonferroni |
297|| `tukey` | Tukey HSD | 多组 post-hoc 比较 | 直接使用 Tukey HSD p 值 |
298|
299|</details>
300|
301|<details>
302|<summary><b>极值预处理</b></summary>
303|
304|使用 `--remove-outliers` 在统计和绘图前删除极端表型值：
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
316|规则为 Tukey IQR，`k=1.5`：删除超出 `[Q1 - 1.5 x IQR, Q3 + 1.5 x IQR]` 的值。过滤在每个 `trait x population x haplotype` 分组内独立执行。少于 4 个数值样本的分组不执行删除。
317|
318|summary 输出中会记录预处理信息：
319|
320|| 列名 | 含义 |
321|| --- | --- |
322|| `raw_count` | 删除极值前的数值样本数 |
323|| `raw_min`, `raw_max` | 删除前的原始范围 |
324|| `outlier_removed` | 该 summary 分组中删除的值数量 |
325|| `outlier_method` | `none` 或 `iqr` |
326|| `outlier_iqr_k` | IQR 倍数；启用时为 `1.5` |
327|
328|</details>
329|
330|<details>
331|<summary><b>表型输出文件</b></summary>
332|
333|| 文件 | 内容 |
334|| --- | --- |
335|| `phenotype_stats.tsv` | 两两比较结果，包括分组样本数、均值、标准差、ANOVA、两两统计量、原始 P 值、校正后 P 值、显著性标签和 `effective_n` |
336|| summary TSV (`--summary-output`) | 每个性状、群体、单倍型的 summary 统计；启用极值预处理时记录删除数量 |
337|| boxplot (`--plot-box`) | 对一个选定性状绘制箱线图，使用与统计结果一致的过滤、分组和比较逻辑 |
338|
339|</details>
340|
341|## 其他流程
342|
343|### BED 批量处理
344|
345|```bash
346|haplokit view in.vcf.gz -R regions.bed --output-file out_batch
347|```
348|
349|<details>
350|<summary><b>详细说明</b></summary>
351|
352|`regions.bed` 至少包含三列 tab 分隔字段：
353|
354|```text
355|chr1  1000  2000
356|chr2  5000  6000
357|```
358|
359|每个 BED 行独立处理。输出文件按区间 suffix 命名，例如 `_chr1_1000_2000`。
360|
361|</details>
362|
363|### 近似分组
364|
365|```bash
366|haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
367|```
368|
369|`--max-diff` 会把差异比例不高于阈值的单倍型聚为同一组。例如，`--max-diff 0.2` 会合并差异位点 ≤20% 的单倍型。
370|
371|### 样本子集与缺失填补
372|
373|```bash
374|haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
375|```
376|
377|`samples.list` 每行一个样本 ID。`--impute` 将缺失基因型视为参考型 `0/0`，以提高样本保留率。
378|
379|### 基因 ID 和基因列表选择器
380|
381|```bash
382|# 单个基因 ID（需要 --gff）
383|haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 --output-file out
384|
385|# 基因列表文件（每行一个基因 ID）
386|haplokit view in.vcf.gz --gene-list gene_list.txt --gff genes.gff3 --output-file out
387|
388|# 带上下游扩展
389|haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 \
390|  --upstream 2000 --downstream 1000 --strand-aware --output-file out
391|```
392|
393|## 参数速查
394|
395|### `haplokit view`
396|
397|```text
398|haplokit view <input.vcf.gz|input.bcf> (-r <region> | -R <regions.bed> | -t <targets> | -T <targets.txt> | --gene-id <id> | --gene-list <file>) [options]
399|```
400|
401|**Selector 选项**（必须且只能提供一个）：
402|
403|| 参数 | 说明 |
404|| --- | --- |
405|| `-r, --region` | 单个区间：`chr:start-end` 或 `chr:pos` |
406|| `-R, --regions-file` | BED 文件，包含多个区间 |
407|| `-t, --targets` | 同一染色体上逗号分隔的 target 区间 |
408|| `-T, --targets-file` | 同一染色体上每行一个 target 区间的文件 |
409|| `-G, --gene-id` | 单个基因 ID（需要 `--gff/--gff3`） |
410|| `-l, --gene-list` | 每行一个基因 ID 的文件（需要 `--gff/--gff3`） |
411|
412|**识别选项**：
413|
414|| 参数 | 默认 | 说明 |
415|| --- | --- | --- |
416|| `-b, --by` | `auto` | 分组模式：`auto`, `region`, `site` |
417|| `-i, --impute` | 关闭 | 将缺失基因型视为参考型 |
418|| `-m, --max-diff` | 关闭 | 合并差异比例 ≤ 阈值的单倍型，范围 `[0,1]` |
419|
420|**注释选项**：
421|
422|| 参数 | 默认 | 说明 |
423|| --- | --- | --- |
424|| `-g, --gff3, --gff` | 关闭 | 用于 gene selector 和注释的 GFF3/GTF |
425|| `-u, --upstream` | `0` | gene selector 上游扩展长度 |
426|| `-d, --downstream` | `0` | gene selector 下游扩展长度 |
427|| `-a, --strand-aware` | 关闭 | 按基因链方向解释 upstream/downstream |
428|
429|**样本/群体选项**：
430|
431|| 参数 | 默认 | 说明 |
432|| --- | --- | --- |
433|| `-S, --samples-file` | 关闭 | 限制到指定样本 |
434|| `-p, --population` | 关闭 | 样本到群体的映射表（2 列） |
435|
436|**输出选项**：
437|
438|| 参数 | 默认 | 说明 |
439|| --- | --- | --- |
440|| `-o, --output` | `summary` | JSONL 输出内容：`summary` 或 `detail` |
441|| `-f, --output-format` | `tsv` | 输出格式：`tsv` 或 `jsonl` |
442|| `-O, --output-file` | 当前目录 | 输出目录或前缀 |
443|
444|**可视化选项**：
445|
446|| 参数 | 默认 | 说明 |
447|| --- | --- | --- |
448|| `-P, --plot` | 关闭 | 绘制单倍型表图 |
449|| `-F, --plot-format` | `png` | `png`, `pdf`, `svg`, `tiff` |
450|| `--table-theme` | `detailed` | 表格样式：`detailed` 或 `compact` |
451|| `-z, --figsize` | 自动 | 图尺寸，格式 `WIDTH,HEIGHT`（英寸） |
452|| `-e, --geo` | 关闭 | 地图绘图坐标文件 |
453|| `--show-counts` | 关闭 | 在地图饼图中心显示样本数 |
454|| `--hide-counts` | 关闭 | 隐藏样本数标签（默认） |
455|| `-n, --network` | 关闭 | 绘制单倍型网络 |
456|| `-N, --network-method` | `tcs` | 网络算法：`tcs`, `msn`, `mjn` |
457|
458|**标签选项**：
459|
460|| 参数 | 默认 | 说明 |
461|| --- | --- | --- |
462|| `-H, --hap-prefix` | `Hap` | 单倍型标签前缀 |
463|| `-D, --hap-pad` | `2` | 标签数字补零宽度 |
464|
465|**注意事项**：
466|- 必须且只能提供一个 selector：`-r`、`-R`、`-t`、`-T`、`--gene-id` 或 `--gene-list`
467|- 通过 `-t` 或 `-T` 提供的 targets 必须位于同一条染色体
468|- `--gene-id` 和 `--gene-list` 需要 `--gff/--gff3`
469|
470|### `haplokit phenotype`
471|
472|```text
473|haplokit phenotype -H <hapresult.tsv|sample_hap.tsv> -P <phenotype.tsv|phenotype.csv> [options]
474|```
475|
476|**输入选项**：
477|
478|| 参数 | 默认 | 说明 |
479|| --- | --- | --- |
480|| `-H, --hapresult, --haplotypes` | 必选 | `hapresult.tsv` 或两列 sample-haplotype 表 |
481|| `-P, --phenotypes, --phenotype, --pheno-file` | 必选 | 表型表；第一列是样本 ID |
482|| `-p, --population, --pop-group` | 关闭 | 样本到群体的映射表 |
483|
484|**检验选项**：
485|
486|| 参数 | 默认 | 说明 |
487|| --- | --- | --- |
488|| `-t, --trait` | 所有数值性状 | 选择分析性状；可重复 |
489|| `-m, --min-hap-size` | `5` | 每个检验分组所需的最小数值样本数 |
490|| `-M, --method` | `welch` | 检验方法：`welch`, `student`, `mannwhitney`, `tukey` |
491|| `-a, --adjust` | `bonferroni` | 非 Tukey 检验的 P 值校正：`bonferroni` 或 `none` |
492|| `--remove-outliers` | 关闭 | 统计和绘图前删除 Tukey IQR k=1.5 极值 |
493|
494|**输出选项**：
495|
496|| 参数 | 默认 | 说明 |
497|| --- | --- | --- |
498|| `-o, --output` | `phenotype_stats.tsv` | 两两统计结果 TSV |
499|| `-s, --summary-output` | 关闭 | 每个单倍型的 summary TSV |
500|
501|