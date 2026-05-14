# haplokit

面向 CLI 的单倍型查看工具，提供 bcftools 风格选择器、C++ 后端加速和 Python 绘图。

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

## 安装

### 运行要求

- Linux 或 WSL（发布链路当前为 Linux-first）
- Python 3.10+
- 源码构建需要 C++17 工具链和 CMake 3.22+

### 从 PyPI 安装

```bash
python -m pip install haplokit
```

### 从源码安装

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .[test]
```

## 快速开始

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

该命令会生成：

- `out/hapresult.tsv`
- `out/hap_summary.tsv`

## 图表

### 单倍型汇总表

<img src="plottable.png" alt="单倍型汇总表" width="800">

单倍型汇总表展示选定基因组区域内所有变异位点处识别到的单倍型模式。每个视觉组件传递特定信息：

- **标题**：显示基因组区域（`CHR:start-end`），当提供 `--gff` 时还显示重叠的基因名称。
- **功能分类色带**（仅 `--gff`）：位于 POS 行上方的细长彩色条。根据 GFF3 注释将每个变异位点归类为 SnpEff 风格的功能类别（CDS、UTR、exon、intron、upstream/downstream、intergenic），一目了然地展示变异落在哪些功能区域。
- **POS 行**：所选区域内每个变异位点的物理位置。
- **ALLELE 行**：每个位点处的替代等位基因，按等位基因着色。
- **单倍型行**（H001、H002、...）：每行代表一种独特的单倍型模式。单元格显示该位点携带的等位基因，空单元格表示参考等位基因。
- **群体列**（仅 `--population`）：每个群体（如野生种、地方品种、栽培种）对应一列，显示属于各单倍型的样本数。
- **n/N 列**：各单倍型频率，格式为计数/总数。
- **SnpEff 风格图例**（仅 `--gff`）：底部图例展示功能类别颜色（CDS、UTR、exon、intron、intergenic）。
- **Indel 脚注**：多等位基因 Indel（如 `T/A,GG`）用上标标记标注，完整序列在脚注区域解释。

### 单倍型地理分布图

<img src="plotmap.png" alt="单倍型地理分布图" width="600">

地理分布图将单倍型组成饼图叠加在底图上，揭示采样地点间的单倍型空间变异模式。

- **底图**：从 GeoJSON 渲染的省界多边形（中国底图来自阿里云 DataV API）。
- **饼图**：在每个采样位置，饼图展示单倍型组成，每个扇形代表一种单倍型。
- **√频率缩放**：符号大小与 √(样本总数) 成正比，与 R 的 `symbol.lim` 逻辑一致。
- **计数标签**：每个位置的样本总数显示在饼图中心。
- **坐标轴**：经度（x 轴）和纬度（y 轴），刻度标记为淡灰色。
- **图例**：左上角的单倍型颜色图例标识每种单倍型。
- **标题**：可选的图表标题。

## 与 bcftools 对齐的选择器语义

`haplokit view` 采用与 `bcftools` 工作流同形态的选择器词汇。

- `-r/--region chr:start-end`：区间选择
- `-r/--region chr:pos`：单个位点选择
- `-R/--regions-file regions.bed`：按 BED 每行独立处理
- `-S/--samples-file samples.list`：仅保留样本文件中的样本

校验规则：

- `-r` 与 `-R` 二选一，且必须提供其一
- `-r` 与 `-R` 互斥
- `--by site` 仅可用于 `-r chr:pos`
- `--by region` 与 `-r chr:pos` 冲突
- `--by site` 与 `-r chr:start-end` 冲突

## C++ 后端加速

Python CLI 会把核心单倍型分组计算委托给 `haplokit_cpp`。

### 供应商依赖库

- **[htslib](https://github.com/samtools/htslib)** — 高通量测序数据读写 C 库。原生支持 VCF 和 BCF 格式，提供索引随机访问和高效基因型解码。构建时链接为静态库。

- **[gffsub](https://github.com/WWz33/gffsub)** — 轻量级 GFF3/GTF 解析与过滤库。解析基因注释文件，支持 feature 类型过滤（最长转录本选择），提供 overlap/nearest-gene 查询用于单倍型注释。

后端发现顺序：

1. `HAPLOKIT_CPP_BIN`
2. 打包内二进制：`haplokit/_bin/haplokit_cpp`
3. 仓库构建产物：`build-wsl/haplokit_cpp`，然后 `build/haplokit_cpp`
4. 回退本地构建：`cmake -S . -B build-wsl` 与 `cmake --build build-wsl --clean-first -j1`

若发现与回退构建后仍找不到后端，CLI 会报错退出。

## 命令

```text
haplokit view <input_vcf> (-r <region> | -R <regions.bed>) [options]
```

`<input_vcf>` 应为已建立索引的 VCF/BCF（`.vcf.gz` + `.tbi`，或 BCF 索引）。

## 参数与默认值

| 参数 | 类型 | 默认值 | 行为 |
| --- | --- | --- | --- |
| `input_vcf` | 位置参数路径 | 解析器默认 `None`（实际必须提供） | 输入索引化 VCF/BCF |
| `-r, --region` | 选择器字符串 | `None` | `chr:start-end` 或 `chr:pos` |
| `-R, --regions-file` | BED 路径 | `None` | BED 每行至少 3 列（tab 分隔） |
| `-S, --samples-file` | 样本列表路径 | `None` | 每行一个样本 ID |
| `--by` | `auto \| region \| site` | `auto` | 按选择器形态自动推断；解析器会做一致性校验 |
| `--impute` | 开关 | `False` | 分组前将缺失基因型按参考态填补 |
| `-g, --gff3, --gff` | GFF/GFF3 路径 | `None` | 启用基因 overlap/nearest 注释 |
| `-p, --population` | 群体分组文件路径 | `None` | Tab 分隔文件，映射样本 → 群体分组 |
| `--output` | `summary \| detail` | `summary` | 在 JSONL 模式下生效；TSV 模式总是写两类表 |
| `--output-format` | `tsv \| jsonl` | `tsv` | 默认契约是 `tsv` |
| `--output-file` | 路径 | `None` | 可作为输出目录、前缀文件或 JSONL 目标文件 |
| `--plot` | 开关 | `False` | 每个 selector 输出一张单倍型表图 |
| `--plot-format` | `png \| pdf \| svg \| tiff` | `png` | 输出图片格式 |
| `--max-diff` | `[0,1]` 浮点数 | `None` | 启用阈值近似分组 |

## 输出行为（按当前实现）

### `--output-format tsv`（默认）

- 始终同时写出：
  - `hapresult*.tsv`
  - `hap_summary*.tsv`
- `--output` 不改变 TSV 双文件行为。
- 若启用 `--plot`，额外写每个 selector 的图片文件（格式由 `--plot-format` 设定）。
- 若启用 `--gff/--gff3`，额外写 `gff_ann_summary.tsv`。

命名规则：

- 单 selector 且不传 `--output-file`：当前目录写 `hapresult.tsv` 与 `hap_summary.tsv`
- `--output-file <dir>`（无后缀）：写入该目录
- `--output-file <path/custom.tsv>`：走前缀命名：
  - `custom.hapresult.tsv`
  - `custom.hap_summary.tsv`
- BED 多 selector：按 selector slug 命名：
  - `hapresult_<chrom>_<start>_<end>.tsv`
  - `hap_summary_<chrom>_<start>_<end>.tsv`
  - 若重名，追加 `_<NNN>`

### `--output-format jsonl`（兼容模式）

- 不传 `--output-file`：输出到 stdout。
- `--output-file` 为目录路径：写 `<dir>/result.jsonl`。
- `--output-file` 为文件路径：写到该精确路径。
- JSONL 模式下会遵循 `--output summary|detail`。

## 示例

### 区间模式（严格精确分组）

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
```

### 单位点模式

```bash
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

### BED 批处理 + 样本子集 + 绘图 + 注释

```bash
haplokit view in.vcf.gz -R regions.bed -S samples.list --plot --gff genes.gff3 --output-file out_bed
```

### 带群体分组和 PNG 图片

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --plot-format png --output-file out
```

### JSONL detail + 近似分组

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-format jsonl --output detail --output-file result.jsonl
```

## 升级说明

- 默认输出已是 TSV（`hapresult/hap_summary` 双文件）。
- `--output-format jsonl` 仅建议用于兼容旧流程。
- 默认图片格式已改为 PNG（原为 SVG+PDF+TIFF 三格式同时输出）。使用 `--plot-format` 选择格式。

## 贡献开发

Linux/WSL 验证流程：

```bash
cmake -S . -B build-wsl
cmake --build build-wsl -j12
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python
ctest --test-dir build-wsl --output-on-failure
```

参考文档：

- `docs/specs/haplokit-view-cli.md`
- `docs/specs/haplokit-result-schema.md`
- `docs/development/haplokit-linux-workflow.md`
- `docs/release/pypi-release.md`

## 致谢

haplokit 的设计灵感来自 geneHapR：

> Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9

## 许可证

GPL-3.0-or-later。
