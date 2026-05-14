# haplokit

面向 CLI 的单倍型查看工具，提供 bcftools 风格选择器、C++ 后端加速和 Python 绘图。

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

## 安装

```bash
pip install haplokit
```

> 源码构建需要 Linux/WSL、Python 3.10+、C++17 工具链、CMake 3.22+ — 见[贡献开发](#贡献开发)。

## 快速开始

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

输出：

- `out/hapresult.tsv` — 逐样本单倍型详情
- `out/hap_summary.tsv` — 单倍型计数汇总

## 使用场景

### 1. 区域查询 — 严格精确分组

识别基因组区域内的所有不同单倍型。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
```

在 `out/` 中生成 `hapresult.tsv` + `hap_summary.tsv`。每行单倍型展示精确等位基因模式；含杂合或缺失呼叫的样本被排除。

### 2. 单位点查询

分析单个变异位点的单倍型。

```bash
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

`--by` 对 `chr:pos` 选择器自动推断为 `site`。

### 3. 基因注释 + 图表

在单倍型表上叠加基因结构。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
```

添加 SnpEff 风格功能分类色带（CDS、UTR、exon、intron、intergenic）到变异位点上。输出图片（`out/*.png`）+ `gff_ann_summary.tsv`。

<img src="plottable.png" alt="单倍型汇总表" width="800">

图表组件：

- **标题**：区域 + 重叠基因名（提供 `--gff` 时）
- **功能色带**（仅 `--gff`）：彩色条按功能类别标注每个变异
- **POS / ALLELE 行**：变异位置和替代等位基因
- **单倍型行**（H001、H002、...）：每位点等位基因；空 = 参考
- **群体列**（`--population`）：各群体各单倍型样本数
- **n/N**：单倍型频率
- **图例**（仅 `--gff`）：功能类别颜色
- **Indel 脚注**：多等位基因 Indel 用上标标记标注

### 4. 群体分组

比较不同群体间的单倍型分布。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

`popgroup.txt` 为 Tab 分隔：`sample<TAB>population`。在表和图中添加群体列。

### 5. 地理分布图

在采样地点上映射单倍型组成。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo sample_geo.txt --plot --output-file out
```

<img src="plotmap.png" alt="单倍型地理分布图" width="600">

图表组件：

- **饼图**：每位置单倍型组成；大小 ∝ √(样本数)
- **计数标签**：饼图中心显示样本总数
- **图例**：单倍型颜色键
- **底图**：GeoJSON 省界多边形（中国）

### 6. BED 批处理

一次运行处理多个区域。

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

每行 BED 独立处理。输出文件按区域 slug 加后缀（`_chr1_1000_2000`）。

### 7. 近似分组

在容差范围内聚类相似单倍型。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff`（0–1）：差异 ≤ 20% 位点的单倍型归为一组。分组模式从 `strict-region` 变为 `approx-region`。

### 8. 样本子集 + 填补

限制分析到特定样本；缺失呼叫按参考等位基因填补。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list`：每行一个样本 ID。`--impute` 将缺失 GT 视为 `0/0`，提高样本保留率。

## 完整参数

```
haplokit view <input_vcf> (-r <region> | -R <regions.bed>) [options]
```

`<input_vcf>` 须为已索引 VCF/BCF（`.vcf.gz` + `.tbi`，或 BCF 索引）。

| 参数 | 类型 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `-r, --region` | 字符串 | — | `chr:start-end` 或 `chr:pos` |
| `-R, --regions-file` | 路径 | — | BED 文件（≥3 列，Tab 分隔） |
| `-S, --samples-file` | 路径 | — | 每行一个样本 ID |
| `--by` | `auto\|region\|site` | `auto` | 分组模式；auto 按选择器形态推断 |
| `--impute` | 开关 | 关 | 缺失 GT 按参考等位基因填补 |
| `-g, --gff` | 路径 | — | GFF3/GTF 基因注释 |
| `-p, --population` | 路径 | — | Tab 分隔的 样本→群体 映射 |
| `--output` | `summary\|detail` | `summary` | 仅 JSONL 模式生效；TSV 总是双表 |
| `--output-format` | `tsv\|jsonl` | `tsv` | 输出格式 |
| `--output-file` | 路径 | — | 输出目录、前缀或 JSONL 文件 |
| `--plot` | 开关 | 关 | 生成单倍型表图 |
| `--plot-format` | `png\|pdf\|svg\|tiff` | `png` | 图片格式 |
| `--max-diff` | 浮点 [0,1] | — | 近似分组阈值 |
| `--geo` | 路径 | — | 样本地理坐标（用于地图） |

选择器规则：`-r` 与 `-R` 互斥且必须提供其一。`--by site` 仅可用于 `-r chr:pos`。

## 后端

C++ 后端（`haplokit_cpp`）处理 VCF 读取和单倍型分组。发现顺序：

1. `HAPLOKIT_CPP_BIN` 环境变量
2. 打包内二进制：`haplokit/_bin/haplokit_cpp`
3. 仓库构建产物：`build-wsl/haplokit_cpp` → `build/haplokit_cpp`
4. 回退：自动运行 `cmake` 构建

供应商依赖库：

- **[htslib](https://github.com/samtools/htslib)** — VCF/BCF 读取，索引随机访问
- **[gffsub](https://github.com/WWz33/gffsub)** — GFF3/GTF 解析，overlap/nearest-gene 查询

## 贡献开发

```bash
cmake -S . -B build-wsl && cmake --build build-wsl -j12
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python
ctest --test-dir build-wsl --output-on-failure
```

## 致谢

设计灵感来自 geneHapR：

> Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9

## 许可证

GPL-3.0-or-later
