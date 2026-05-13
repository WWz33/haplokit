# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Structure

This repository contains two independent projects that share test data:

- **haptools** — 活跃开发的 Python+C++ 命令行单倍型查看器
- **geneHapR** — 遗留 R 包，仅保留作为科学参考，不再活跃开发

新功能必须在 haptools (Python/C++) 中实现，不要修改 R 代码。

---

## haptools (Python + C++)

### 概述

CLI 单倍型查看器。读取索引的 VCF/BCF 文件，按指定基因组区域的单倍型模式对样本分组，生成汇总表、详情表和可选图表。支持 bcftools 风格的区域选择器、BED 批处理、样本子集、GFF 基因注释和近似单倍型分组。

### 命令

**构建 C++ 后端：**
```bash
cmake -S . -B build-wsl
cmake --build build-wsl --parallel
```

**安装 Python 包（自动编译 C++ 后端）：**
```bash
pip install -e ".[test]"
```
`setup.py` 中的自定义 `BuildPyWithCpp` 命令会在安装时运行 CMake，将 `haptools_cpp` 打包到 `haptools/_bin/`。

**运行：**
```bash
haptools view <input.vcf.gz> -r chr:start-end --output-file out
python -m haptools view <input.vcf.gz> -r chr:start-end --output-file out
```

**测试：**
```bash
# Python 测试
HAPTOOLS_CPP_BIN=$PWD/build-wsl/haptools_cpp python -m pytest -q tests/python

# C++ 测试
ctest --test-dir build-wsl --output-on-failure
```

**构建 sdist：**
```bash
python -m build --sdist
```

### 架构

#### 两层设计：Python CLI + C++ 子进程后端

- **Python CLI** (`haptools/cli.py`)：参数解析、TSV/JSONL 输出、绘图 (matplotlib)
- **C++ 后端** (`src/cpp/`)：通过 htslib 读取 VCF、按单倍型对样本分组、通过 gffsub 解析 GFF3 基因注释

通信方式：Python 通过 `subprocess.run()` 调用 C++ 后端，结果以 JSON 通过 stdout 传递。两个子命令：`view-json`（单区域）和 `view-bed-jsonl`（BED 文件，JSONL 输出）。

后端发现链：`HAPTOOLS_CPP_BIN` 环境变量 → `haptools/_bin/haptools_cpp` 打包二进制 → 仓库构建目录 → 自动尝试 CMake 构建。

#### 目录结构

| 路径 | 用途 |
|---|---|
| `haptools/` | Python 包（cli.py, plot.py, __main__.py） |
| `src/cpp/` | C++17 后端（vcf_reader, selector, view_backend, gff_annotator） |
| `tests/python/` | pytest 测试套件 |
| `tests/cpp/` | C++ 测试可执行文件（CTest） |
| `deps/htslib/` | 供应商 htslib，CMake 构建时编译为静态库 |
| `deps/gffsub/` | 供应商 gffsub，GFF3 解析/过滤库，编译进 haptools_core |
| `docs/` | 架构文档、规范、开发流程 |
| `inst/extdata/` | 共享测试数据（VCF、GFF、BED、FASTA 等） |

#### 选择器语义（bcftools 风格）

- `-r chr:start-end` — 区域模式
- `-r chr:pos` — 位点模式
- `-R regions.bed` — BED 批处理
- `-S samples.list` — 样本子集
- `--by auto|region|site` — 从选择器形状自动推断

#### 打包机制

`setup.py` 的 `BuildPyWithCpp` 命令编译 C++ 后端并打包进 Python 包，使终端用户从 sdist 安装时无需 C++ 工具链。

#### 供应商 htslib

C++ 构建默认链接 `deps/htslib/` 作为静态库。可通过 `HAPTOOLS_USE_VENDORED_HTSLIB=OFF` 使用系统 htslib。

### 系统依赖（C++ 构建）

`build-essential`, `cmake`, `pkg-config`, `libbz2-dev`, `liblzma-dev`, `libcurl4-openssl-dev`, `zlib1g-dev`

---

## geneHapR (R 包) — 遗留参考

### 概述

geneHapR 是一个 R 包（v1.2.6），用于基因单倍型分析、可视化和表型关联。已不再活跃开发，保留作为科学和历史参考。

### 目录结构

整个 R 包位于 `geneHapR/` 子目录下，是一个标准 R 包结构：

| 路径 | 用途 |
|---|---|
| `geneHapR/R/` | 27 个 R 源文件（导入、过滤、单倍型计算、可视化、GUI） |
| `geneHapR/tests/testthat/` | R 测试套件（testthat） |
| `geneHapR/vignettes/` | R vignette（Rmd + 图表） |
| `geneHapR/man/` | R 文档（.Rd 文件） |
| `geneHapR/data/` | R 数据文件（.rda, .RData） |
| `geneHapR/DESCRIPTION` | R 包元数据 |
| `geneHapR/NAMESPACE` | R 包导出/导入 |
| `geneHapR/inst/extdata` | symlink → `../../inst/extdata`（共享测试数据） |

### 共享测试数据

根目录 `inst/extdata/` 包含 VCF 文件、GFF 注释、BED、FASTA、PLINK ped/map 等测试数据，供 haptools 和 geneHapR 共用。geneHapR 通过 `geneHapR/inst/extdata` symlink 访问。

### R 包命令

```bash
# 在 R 中安装开发版本（需在 geneHapR/ 目录下）
cd geneHapR && Rscript -e "devtools::install()"

# 运行测试
cd geneHapR && Rscript -e "devtools::test()"

# 构建文档
cd geneHapR && Rscript -e "devtools::document()"
```
