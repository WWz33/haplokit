# gffsub

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

一款使用 C++ 编写的快速 GFF3/GTF 提取和过滤工具。

## 功能特性

- **区域提取**: 按基因组区域提取特征
- **BED 文件支持**: 使用 BED 文件过滤区域
- **最长转录本**: 保留每个基因的最长转录本亚型
- **多格式输出**: 支持 GFF3、GTF2、GTF3 和 BED 格式

## 安装

### 环境要求

- C++17 编译器（g++ 或 clang++）

### 编译

```bash
cd gffsub_dev
make
```

## 使用方法

```
./gffsub <input.gff3> [options]
```

### 输入/区域选项

| 短选项 | 长选项 | 说明 |
|--------|--------|------|
| `-r` | `--region CHR:START-END` | 提取区域（1-based，闭合区间） |
| `-b` | `--bed FILE` | BED 文件格式的区域（0-based → 1-based） |

### 特征过滤选项

| 短选项 | 长选项 | 说明 |
|--------|--------|------|
| `-f` | `--feature TYPE` | 按特征类型过滤（gene、mRNA、exon、CDS...） |
| `-L` | `--longest` | 保留每个基因的最长转录本 |

### 输出选项

| 短选项 | 长选项 | 说明 |
|--------|--------|------|
| `-t` | `--output-format FMT` | 格式：gff3、gtf2、gtf3、bed |
| `-o` | `--output FILE` | 输出文件（默认：stdout） |
| `-h` | `--help` | 显示帮助信息 |

## 使用示例

```bash
# 提取区域内的基因
./gffsub annotation.gff3 -r chr1:1-100000 -f gene

# 使用 BED 文件提取
./gffsub annotation.gff3 -b regions.bed -f exon

# 保留最长转录本
./gffsub annotation.gff3 --longest

# 输出为 GTF3 格式
./gffsub annotation.gff3 -r chr1:1-100000 -t gtf3 -o out.gtf
```

## 输出格式说明

| 格式 | 头部 | 坐标系统 |
|------|------|----------|
| gff3 | `##gff-version 3` | 1-based 闭合区间 |
| gtf2 | `##gtf-version 2` | 1-based 闭合区间 |
| gtf3 | `##gtf-version 2.2.1` | 1-based 闭合区间 |
| bed | （无头部） | 0-based 半开区间 |

## 分发说明

单一静态二进制文件 — 复制到任意 Linux 系统即可运行：

```bash
./gffsub annotation.gff3 -r chr1:1-100000 -f gene
```

**运行环境要求:** glibc 系 Linux（x86_64）

## 许可证

MIT License
