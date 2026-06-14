# haplokit

面向 indexed VCF/BCF 的命令行单倍型分析工具，使用 C++ 处理数据平面，使用 Python 完成统计分析与绘图。

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

## 目录

- [功能概览](#功能概览)
- [安装](#安装)
- [快速开始](#快速开始)
- [单倍型分析流程](#单倍型分析流程)
  - [区间或单个位点单倍型识别](#区间或单个位点单倍型识别)
  - [基因注释与单倍型表图](#基因注释与单倍型表图)
  - [群体分组](#群体分组)
  - [地理分布](#地理分布)
  - [单倍型网络](#单倍型网络)
- [表型统计模块](#表型统计模块)
- [其他流程](#其他流程)
  - [BED 批量处理](#bed-批量处理)
  - [近似分组](#近似分组)
  - [样本子集与缺失填补](#样本子集与缺失填补)
- [参数速查](#参数速查)
- [输入/输出格式](#输入输出格式)
- [引用](#引用)
- [支持](#支持)
- [许可证](#许可证)

## 功能概览

| 模块 | 用途 | 典型输出 |
| --- | --- | --- |
| `view` | 从区间、单个位点、基因 ID、基因列表或 BED 文件提取单倍型 | `hapresult.tsv`, `hap_summary.tsv` |
| 基因注释 | 解析 gene selector，并为变异位置添加基因结构上下文 | `gff_ann_summary.tsv`, 注释版单倍型表图 |
| 群体统计 | 按群体统计每个单倍型的样本数 | 表格和图中的 population columns |
| 地理分布图 | 在采样地点绘制单倍型组成 | 带饼图和样本数比例尺的地图 |
| 单倍型网络 | 使用 MSN、TCS 或 MJN 构建 PopART 风格网络 | 带群体饼图和突变刻度的 network figure |
| `phenotype` | 连接单倍型与数值表型，输出检验结果和箱线图 | `phenotype_stats.tsv`, 表型 summary TSV, boxplot |

## 安装

```bash
pip install haplokit
```

<details>
<summary><b>源码构建（高级）</b></summary>

#### 系统要求

- Linux/WSL
- Python 3.10+
- C++17 编译器
- CMake 3.22+
- `make`
- htslib 所需的本地库

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

#### 从源码目录安装

```bash
pip install .
```

#### 开发模式安装

```bash
pip install -e .
```

#### 自定义后端路径

如果 C++ 后端构建在其他位置，可以显式指定：

```bash
export HAPLOKIT_CPP_BIN=/path/to/haplokit_cpp
```

#### 常见链接错误排查

| 错误 | 需要安装 |
| --- | --- |
| `cannot find -lcurl` | `libcurl` / `libcurl4-openssl-dev` |
| `cannot find -lbz2` | `bzip2` / `libbz2-dev` |
| `cannot find -llzma` | `xz` / `liblzma-dev` |
| `cannot find -lz` | `zlib` / `zlib1g-dev` |

</details>

## 快速开始

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

从区间 `scaffold_1:4300-5000` 提取单倍型，结果写入 `out/` 目录：

| 文件 | 含义 |
| --- | --- |
| `out/hapresult.tsv` | 单倍型等位基因模式与样本列表 |
| `out/hap_summary.tsv` | 单倍型计数和频率 |

## 单倍型分析流程

### 区间或单个位点单倍型识别

```bash
# 区间模式：按整个区域内的等位基因模式分组
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out

# 单点模式：按单个位点的等位基因分组
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

<details>
<summary><b>详细说明</b></summary>

区间 selector 会按整个区域内的等位基因模式分组。单点 selector 会自动进入 site mode。严格区间模式下，带杂合或缺失调用的样本会被排除；如需保留缺失样本，可以使用 `--impute`。

</details>

### 基因注释与单倍型表图

```bash
# 带基因模型注释
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out

# 紧凑表格样式
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --table-theme compact --output-file out
```

<details>
<summary><b>详细说明</b></summary>

提供 GFF3/GTF 文件时，图中会在表格上方绘制 pyGenomeTracks 风格的基因模型（backbone/内含子、CDS、UTR，末端箭头表示链方向），并用 SNP 刻度和引导线把每个变异连到对应的等位基因列，同时附带 CDS/UTR/intron 图例。不提供 `--gff` 时只绘制表格。输出还包括 `gff_ann_summary.tsv`。

`--table-theme` 选择表格样式：`detailed`（默认，方形单元格、白色网格线）或 `compact`（扁宽单元格、无边贴合、更矮的表头）。基因模型在两种样式下通用。

</details>

<img src="data/figure/haplotype_table.png" alt="单倍型汇总表" width="800">

### 群体分组

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

<details>
<summary><b>详细说明</b></summary>

`popgroup.txt` 是两列 tab 分隔文件：

```text
sample  population
C1      wild
C2      wild
C13     landrace
```

群体信息会进入输出表格的计数列，也会进入图中的群体统计。

</details>

### 地理分布

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
```

<details>
<summary><b>详细说明</b></summary>

坐标文件为 tab 分隔：

```text
ID    longitude  latitude
C1    116.40     39.90
C2    116.40     39.90
```

使用 `--show-counts` 在地图饼图中心显示样本数，或使用 `--hide-counts` 显式隐藏。

`data/` 中附带世界地图示例资源：

- `sample_world_geo.txt`
- `world_countries.shp`, `world_countries.shx`, `world_countries.dbf`
- `data/figure/haplotype_map_world.png`

</details>

<img src="data/figure/haplotype_map_china.png" alt="单倍型地理分布图" width="600">

<img src="data/figure/haplotype_map_world.png" alt="世界单倍型地理分布图" width="600">

### 单倍型网络

```bash
# TCS 网络（默认）
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out

# 中值连接网络
haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
```

| 方法 | 名称 | 说明 |
| --- | --- | --- |
| `tcs` | 统计简约性网络 | Templeton, Crandall & Sing (1992) |
| `msn` | 最小生成树网络 | 基于 Hamming 距离 |
| `mjn` | 中值连接网络 | Bandelt et al. (1999) |

<details>
<summary><b>详细说明</b></summary>

网络图遵循 PopART 风格：节点面积表示单倍型样本数，饼图扇区表示群体组成，边上的刻度表示突变步数，小黑点表示推断出的中间节点。

</details>

![网络算法对比 - MSN / TCS / MJN](data/figure/haplotype_network_algorithms.png)

## 表型统计模块

### 基本用法

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

### 带箱线图

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

<img src="data/figure/phenotype_population_boxplot.png" alt="群体分层表型箱线图" width="900">

<details>
<summary><b>统计场景</b></summary>

| 场景 | 检验分组 | result 中的两两比较 | 箱线图显著性标注 |
| --- | --- | --- | --- |
| 不提供群体文件，保留多个单倍型 | `trait x haplotype` | 每个性状内所有保留单倍型两两比较 | 单倍型之间的显著性 |
| 提供群体文件，保留多个单倍型 | `trait x population x haplotype` | 每个群体内部的单倍型两两比较 | 群体内单倍型比较；同一单倍型的群体间比较 |
| 提供群体文件，只保留一个单倍型 | `trait x haplotype x population` | 该单倍型内部的群体两两比较 | 仅显示群体间比较 |
| 多个性状 | 每个性状独立分析 | 每个性状输出独立结果块 | 绘图必须用 `--trait` 选择一个性状 |
| 表型缺失值 | 按性状忽略非数值和缺失值 | 计数只包含数值样本 | `effective_n` 记录进入该分层的有效样本数 |
| IQR 极值预处理 | 可选；在每个 `trait x population x haplotype` 内执行 Tukey IQR k=1.5 | 检验使用删除极值后的样本 | 图使用同一批过滤后的样本；summary 记录删除数量 |

</details>

<details>
<summary><b>两两检验方法</b></summary>

假设检验使用 `scipy.stats`。

| `--method` | 检验 | 适用场景 | P 值校正 |
| --- | --- | --- | --- |
| `welch` | Welch two-sample t-test | 默认；不假设方差相等 | 非 Tukey 检验默认 Bonferroni |
| `student` | Student two-sample t-test | 可接受方差相等假设时 | 非 Tukey 检验默认 Bonferroni |
| `mannwhitney` | Mann-Whitney U test | 非参数秩检验 | 非 Tukey 检验默认 Bonferroni |
| `tukey` | Tukey HSD | 多组 post-hoc 比较 | 直接使用 Tukey HSD p 值 |

</details>

<details>
<summary><b>极值预处理</b></summary>

使用 `--remove-outliers` 在统计和绘图前删除极端表型值：

```bash
haplokit phenotype \
  -H out/hapresult.tsv \
  -P phenotype.csv \
  -t yield \
  --remove-outliers \
  -o yield_stats.tsv \
  -s yield_summary.tsv
```

规则为 Tukey IQR，`k=1.5`：删除超出 `[Q1 - 1.5 x IQR, Q3 + 1.5 x IQR]` 的值。过滤在每个 `trait x population x haplotype` 分组内独立执行。少于 4 个数值样本的分组不执行删除。

summary 输出中会记录预处理信息：

| 列名 | 含义 |
| --- | --- |
| `raw_count` | 删除极值前的数值样本数 |
| `raw_min`, `raw_max` | 删除前的原始范围 |
| `outlier_removed` | 该 summary 分组中删除的值数量 |
| `outlier_method` | `none` 或 `iqr` |
| `outlier_iqr_k` | IQR 倍数；启用时为 `1.5` |

</details>

<details>
<summary><b>表型输出文件</b></summary>

| 文件 | 内容 |
| --- | --- |
| `phenotype_stats.tsv` | 两两比较结果，包括分组样本数、均值、标准差、ANOVA、两两统计量、原始 P 值、校正后 P 值、显著性标签和 `effective_n` |
| summary TSV (`--summary-output`) | 每个性状、群体、单倍型的 summary 统计；启用极值预处理时记录删除数量 |
| boxplot (`--plot-box`) | 对一个选定性状绘制箱线图，使用与统计结果一致的过滤、分组和比较逻辑 |

</details>

## 其他流程

### BED 批量处理

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

<details>
<summary><b>详细说明</b></summary>

`regions.bed` 至少包含三列 tab 分隔字段：

```text
chr1  1000  2000
chr2  5000  6000
```

每个 BED 行独立处理。输出文件按区间 suffix 命名，例如 `_chr1_1000_2000`。

</details>

### 近似分组

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff` 会把差异比例不高于阈值的单倍型聚为同一组。例如，`--max-diff 0.2` 会合并差异位点 ≤20% 的单倍型。

### 样本子集与缺失填补

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list` 每行一个样本 ID。`--impute` 将缺失基因型视为参考型 `0/0`，以提高样本保留率。

### 基因 ID 和基因列表选择器

```bash
# 单个基因 ID（需要 --gff）
haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 --output-file out

# 基因列表文件（每行一个基因 ID）
haplokit view in.vcf.gz --gene-list gene_list.txt --gff genes.gff3 --output-file out

# 带上下游扩展
haplokit view in.vcf.gz --gene-id Glyma.01G001000 --gff genes.gff3 \
  --upstream 2000 --downstream 1000 --strand-aware --output-file out
```

## 参数速查

### `haplokit view`

```text
haplokit view <input.vcf.gz|input.bcf> (-r <region> | -R <regions.bed> | -t <targets> | -T <targets.txt> | --gene-id <id> | --gene-list <file>) [options]
```

**Selector 选项**（必须且只能提供一个）：

| 参数 | 说明 |
| --- | --- |
| `-r, --region` | 单个区间：`chr:start-end` 或 `chr:pos` |
| `-R, --regions-file` | BED 文件，包含多个区间 |
| `-t, --targets` | 同一染色体上逗号分隔的 target 区间 |
| `-T, --targets-file` | 同一染色体上每行一个 target 区间的文件 |
| `-G, --gene-id` | 单个基因 ID（需要 `--gff/--gff3`） |
| `-l, --gene-list` | 每行一个基因 ID 的文件（需要 `--gff/--gff3`） |

**识别选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-b, --by` | `auto` | 分组模式：`auto`, `region`, `site` |
| `-i, --impute` | 关闭 | 将缺失基因型视为参考型 |
| `-m, --max-diff` | 关闭 | 合并差异比例 ≤ 阈值的单倍型，范围 `[0,1]` |

**注释选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-g, --gff3, --gff` | 关闭 | 用于 gene selector 和注释的 GFF3/GTF |
| `-u, --upstream` | `0` | gene selector 上游扩展长度 |
| `-d, --downstream` | `0` | gene selector 下游扩展长度 |
| `-a, --strand-aware` | 关闭 | 按基因链方向解释 upstream/downstream |

**样本/群体选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-S, --samples-file` | 关闭 | 限制到指定样本 |
| `-p, --population` | 关闭 | 样本到群体的映射表（2 列） |

**输出选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-o, --output` | `summary` | JSONL 输出内容：`summary` 或 `detail` |
| `-f, --output-format` | `tsv` | 输出格式：`tsv` 或 `jsonl` |
| `-O, --output-file` | 当前目录 | 输出目录或前缀 |

**可视化选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-P, --plot` | 关闭 | 绘制单倍型表图 |
| `-F, --plot-format` | `png` | `png`, `pdf`, `svg`, `tiff` |
| `--table-theme` | `detailed` | 表格样式：`detailed` 或 `compact` |
| `-z, --figsize` | 自动 | 图尺寸，格式 `WIDTH,HEIGHT`（英寸） |
| `-e, --geo` | 关闭 | 地图绘图坐标文件 |
| `--show-counts` | 关闭 | 在地图饼图中心显示样本数 |
| `--hide-counts` | 关闭 | 隐藏样本数标签（默认） |
| `-n, --network` | 关闭 | 绘制单倍型网络 |
| `-N, --network-method` | `tcs` | 网络算法：`tcs`, `msn`, `mjn` |

**标签选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-H, --hap-prefix` | `Hap` | 单倍型标签前缀 |
| `-D, --hap-pad` | `2` | 标签数字补零宽度 |

**注意事项**：
- 必须且只能提供一个 selector：`-r`、`-R`、`-t`、`-T`、`--gene-id` 或 `--gene-list`
- 通过 `-t` 或 `-T` 提供的 targets 必须位于同一条染色体
- `--gene-id` 和 `--gene-list` 需要 `--gff/--gff3`

### `haplokit phenotype`

```text
haplokit phenotype -H <hapresult.tsv|sample_hap.tsv> -P <phenotype.tsv|phenotype.csv> [options]
```

**输入选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-H, --hapresult, --haplotypes` | 必选 | `hapresult.tsv` 或两列 sample-haplotype 表 |
| `-P, --phenotypes, --phenotype, --pheno-file` | 必选 | 表型表；第一列是样本 ID |
| `-p, --population, --pop-group` | 关闭 | 样本到群体的映射表 |

**检验选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-t, --trait` | 所有数值性状 | 选择分析性状；可重复 |
| `-m, --min-hap-size` | `5` | 每个检验分组所需的最小数值样本数 |
| `-M, --method` | `welch` | 检验方法：`welch`, `student`, `mannwhitney`, `tukey` |
| `-a, --adjust` | `bonferroni` | 非 Tukey 检验的 P 值校正：`bonferroni` 或 `none` |
| `--remove-outliers` | 关闭 | 统计和绘图前删除 Tukey IQR k=1.5 极值 |

**输出选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-o, --output` | `phenotype_stats.tsv` | 两两统计结果 TSV |
| `-s, --summary-output` | 关闭 | 每个单倍型的 summary TSV |

**箱线图选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-B, --plot-box` | 关闭 | 绘制表型箱线图 |
| `-b, --box-output` | `phenotype_box.png` | 箱线图输出路径 |
| `-F, --plot-format` | 输出后缀 | `png`, `pdf`, `svg`, `tiff` |
| `-z, --figsize` | 自动 | 箱线图尺寸，格式 `WIDTH,HEIGHT`（英寸） |
| `-T, --title` | 自动 | 箱线图标题 |
| `-c, --comparison` | 所有可用 pair | 指定标注的单倍型 pair（例如 `Hap01,Hap02`）；可重复 |

**分隔符选项**：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-d, --delimiter` | `auto` | 单倍型输入分隔符：`auto`, `tab`, `comma` |
| `-D, --phenotype-delimiter` | `auto` | 表型输入分隔符 |
| `-G, --population-delimiter` | `auto` | 群体输入分隔符 |

**注意事项**：
- `--plot-box` 需要恰好选择一个性状（`-t`）

<details>
<summary><b>后端</b></summary>

后端二进制：`haplokit_cpp`。

后端发现顺序：

1. `HAPLOKIT_CPP_BIN` 环境变量
2. 包内二进制：`haplokit/_bin/haplokit_cpp`
3. 本地构建：`build-wsl/haplokit_cpp`, `build/haplokit_cpp`, `build-haplokit-python/haplokit_cpp`
4. 源码树 CMake 构建兜底

本地组件：

- [htslib](https://github.com/samtools/htslib)：indexed VCF/BCF 读取
- [gffsub](https://github.com/WWz33/gffsub)：GFF3/GTF 解析和区间查询

</details>

<details>
<summary><b>开发</b></summary>

```bash
# 构建 C++ 后端
cmake -S . -B build-wsl && cmake --build build-wsl -j12

# 运行 Python 测试
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python

# 运行 C++ 测试
ctest --test-dir build-wsl --output-on-failure
```

</details>

## 输入/输出格式

### 输入格式

| 格式 | 扩展名 | 说明 | 备注 |
| --- | --- | --- | --- |
| VCF | `.vcf.gz` | Variant Call Format | 必须建立索引（`.tbi` 或 `.csi`） |
| BCF | `.bcf` | 二进制 VCF | 必须建立索引 |
| BED | `.bed` | 基因组区间 | 3+ 列 tab 分隔：chr, start, end |
| GFF3 | `.gff3` | 基因注释 | 标准 GFF3 格式 |
| GTF | `.gtf` | 基因转移格式 | 标准 GTF 格式 |
| 群体映射表 | `.txt` | 样本到群体的映射 | 2 列 tab 分隔：样本 ID, 群体名称 |
| 表型表 | `.csv`, `.tsv` | 性状值 | 第一列：样本 ID；其余列：数值性状 |
| 基因列表 | `.txt` | 基因 ID | 每行一个基因 ID |

### 输出格式

#### `hapresult.tsv`

`haplokit view` 的主要输出，包含单倍型定义和样本分配。

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `Hap` | string | 单倍型标签（如 Hap01, Hap02） |
| `Count` | integer | 该单倍型的样本数 |
| `Freq` | float | 频率（0-1） |
| `State` | string | 等位基因模式，竖线分隔（如 `0\|0\|1\|0`） |
| `Samples` | string | 逗号分隔的样本 ID |
| `{Pop1}`, `{Pop2}`, ... | integer | 每个群体的计数（仅在提供 `-p` 时出现） |

#### `hap_summary.tsv`

单倍型计数的简化汇总。

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `Hap` | string | 单倍型标签 |
| `Count` | integer | 总计数 |
| `Freq` | float | 频率 |
| `Samples` | string | 样本列表 |

#### `gff_ann_summary.tsv`

基因注释汇总（仅在提供 `--gff` 时生成）。

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `Position` | integer | 变异位置 |
| `Feature` | string | 特征类型（CDS, 5'UTR, 3'UTR, intron, upstream, downstream, intergenic） |
| `Gene` | string | 基因 ID |
| `Distance` | integer | 到特征的距离（0 表示在特征内） |

#### `phenotype_stats.tsv`

`haplokit phenotype` 的两两统计比较结果。

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `trait` | string | 性状名称 |
| `group1`, `group2` | string | 比较的两组（单倍型或群体标签） |
| `n1`, `n2` | integer | 每组的样本数 |
| `mean1`, `mean2` | float | 组均值 |
| `std1`, `std2` | float | 组标准差 |
| `statistic` | float | 检验统计量 |
| `pvalue` | float | 原始 P 值 |
| `pvalue_adjusted` | float | 校正后 P 值（Bonferroni 或 Tukey） |
| `significance` | string | 显著性标签：`***`, `**`, `*`, `ns` |
| `effective_n` | integer | 有效样本数（排除缺失值后） |

## 引用

如果您在研究中使用 haplokit，请引用：

> WWz33. (2024). haplokit: Command-line haplotype analysis for indexed VCF/BCF data. https://github.com/WWz33/haplokit

`haplokit` 受 geneHapR 启发：

> Zhang, R., Jia, G. & Diao, X. (2023). geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199. https://doi.org/10.1186/s12859-023-05318-9

网络图遵循 PopART 的可视化约定：

> Leigh, J. W. & Bryant, D. (2015). popart: full-feature software for haplotype network construction. Methods in Ecology and Evolution 6, 1110-1116. https://doi.org/10.1111/2041-210X.12410

## 支持

- **问题反馈**: [GitHub Issues](https://github.com/WWz33/haplokit/issues)
- **讨论交流**: [GitHub Discussions](https://github.com/WWz33/haplokit/discussions)

## 许可证

GPL-3.0-or-later
