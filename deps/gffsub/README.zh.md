# gffsub

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

`gffsub` 是一个面向日常基因组注释处理的命令行工具，适用于 GFF3/GTF 风格文件。它可以帮助你按区域提取注释、从标识符恢复完整 gene model、构建上下游窗口、保留代表性转录本，并在下游 pipeline 前快速运行 QC。

当普通区间过滤不够用，而你需要理解 GFF3 语义的行为时，可以使用 `gffsub`：例如第 9 列属性查询、`Parent`/child 遍历、gene model 提取，以及注释文件专用的质量检查。

## 从任务开始

| 我想要... | 使用 |
|----------|------|
| 提取某个基因组区间内的 gene | `gffsub annotation.gff3 -r chr1:1-100000 -f gene` |
| 提取某条染色体或 contig 上的记录 | `gffsub annotation.gff3 --seqid chr1` |
| 提取某个注释来源的记录 | `gffsub annotation.gff3 --source Gnomon` |
| 按 score 列过滤记录 | `gffsub annotation.gff3 --score 42.5` |
| 按 strand 列过滤记录 | `gffsub annotation.gff3 --strand -` |
| 按 phase 列过滤 CDS 记录 | `gffsub annotation.gff3 --phase 0 -f CDS` |
| 使用 BED 区间作为输入 | `gffsub annotation.gff3 -b regions.bed -f exon` |
| 按精确 `ID` 查找一个 feature | `gffsub annotation.gff3 --id GeneA` |
| 批量提取精确 ID | `gffsub annotation.gff3 --ids genes.txt` |
| 按常见命名键查找基因 | `gffsub annotation.gff3 --name GeneA` |
| 按任意精确属性值查找记录 | `gffsub annotation.gff3 --where biotype=protein_coding` |
| 用模式文件搜索字段或属性 | `gffsub annotation.gff3 --grep-file genes.txt --grep-field ID` |
| 用正则搜索字段或属性 | `gffsub annotation.gff3 --grep-regex 'ID:^Glyma\.01G'` |
| 用表达式组合语义过滤 | `gffsub annotation.gff3 -I 'type=="gene" && attr.biotype=="protein_coding"'` |
| 用表达式排除记录 | `gffsub annotation.gff3 -E 'attr.Note~"transposon|retroelement"'` |
| 查找离某个区间最近的基因 | `gffsub annotation.gff3 --nearest chr1:1000-2000` |
| 包含匹配记录的后代 | `gffsub annotation.gff3 --id GeneA -C` |
| 包含匹配记录的祖先 | `gffsub annotation.gff3 --id ExonA --parents` |
| 从任意 feature 恢复完整 gene model | `gffsub annotation.gff3 --id ExonA --model` |
| 输出适合 pipeline 读取的 summary | `gffsub annotation.gff3 --id GeneA --summary tsv` |
| 提取指定属性值 | `gffsub annotation.gff3 --id GeneA --out-attrs ID,Name,Parent` |
| 提取上下游背景区域 | `gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware` |
| 每个基因保留最长转录本 | `gffsub annotation.gff3 --longest` |
| 检查注释语法和 graph 问题 | `gffsub annotation.gff3 --qc` |

## 支持的输入和输出

| 类型 | 支持内容 |
|------|----------|
| 注释输入 | GFF3/GTF 风格 feature records |
| 区间输入 | `CHR:START-END` 字符串和 BED 文件 |
| 标识符输入 | 可重复的 `--id` 值，或 `--ids` 指定的一行一个 ID 文件 |
| 模式输入 | `--grep-file` 与 `--grep-field` 指定的一行一个 pattern 文件 |
| 注释输出 | `gff3`, `gtf` (= `gtf2`), `gtf2`, `gtf3`, `bed` |
| 表格输出 | TSV/JSON summary 和 TSV QC report |

`gffsub` 处理的是注释记录。它不会把序列 FASTA 或变异 VCF 内容作为主输入来处理。

## 安装

### 环境要求

- C++17 编译器（`g++` 或 `clang++`）
- CMake 可选；仓库也提供 `Makefile`

### 使用 Make 编译

```bash
cd gffsub
make
```

### 使用 CMake 编译并测试

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## 场景：按基因组上下文提取注释

当问题从坐标区间、染色体/contig、source 列或 BED 文件开始时，使用这类过滤方式。

```bash
./gffsub annotation.gff3 -r chr1:1-100000 -f gene
./gffsub annotation.gff3 --seqid chr1
./gffsub annotation.gff3 --source Gnomon
./gffsub annotation.gff3 --score 42.5
./gffsub annotation.gff3 --strand -
./gffsub annotation.gff3 --phase 0 -f CDS
./gffsub annotation.gff3 -b regions.bed -f exon
./gffsub annotation.gff3 -r chr1:1-100000 -t bed
./gffsub annotation.gff3 -r chr1:1-100000 -o subset.gff3
```

坐标规则是显式的：

| 输入或输出 | 坐标系统 |
|-----------|----------|
| GFF3/GTF records | 1-based 闭合区间 |
| `CHR:START-END` regions | 1-based 闭合区间 |
| BED input/output | 0-based 半开区间 |

选项会在当前记录集合上继续过滤。例如，`-r chr1:1-100000 -f gene -t bed` 会先保留与区间重叠的记录，再限制为 `gene`，最后输出 BED 坐标。

## 场景：查找基因并恢复 Gene Model

当问题从 feature ID、基因名、属性或附近 locus 开始时，使用 selector 选项。

```bash
./gffsub annotation.gff3 --id Glyma.01G000100
./gffsub annotation.gff3 --ids genes.txt
./gffsub annotation.gff3 --name ABC1
./gffsub annotation.gff3 --where biotype=protein_coding
./gffsub annotation.gff3 --where Dbxref=GeneID:123
./gffsub annotation.gff3 --nearest chr1:1000-2000
./gffsub annotation.gff3 --id Glyma.01G000100 -C
./gffsub annotation.gff3 --id ExonA --parents
./gffsub annotation.gff3 --id ExonA --model
```

`--nearest` 会在同一 seqid 上查找距离 1-based 闭合区间最近的 gene。与区间重叠的 gene 距离为 `0`；并列时按输入文件顺序。

在批处理 pipeline 中，可以输出 summary，而不是原始 GFF3：

```bash
./gffsub annotation.gff3 --ids genes.txt --summary tsv
./gffsub annotation.gff3 --id GeneA --summary json
./gffsub annotation.gff3 --id gene0001 --out-attrs ID,Name,Alias,Dbxref
```

Summary 字段包括 `query_id`, `matched_id`, `matched_by`, `seqid`, `start`, `end`, `strand`, `type`, `parent_id`, `child_count`, `transcript_count`, `exon_count`, `cds_length` 和 `status`。如果使用 `--out-attrs`，选中的第 9 列键会追加为 TSV 列，或在 JSON 中输出到 `attrs` 下。

### 搜索和输出中的属性键

GFF3 在第 9 列用分号分隔的 `KEY=VALUE` 对保存记录属性：

```gff3
chr1	src	gene	100	400	.	+	.	ID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123
```

使用 `--id` 做精确 `ID` 查询，使用 `--name` 在常见基因命名键中查找，使用 `--where KEY=VALUE` 做任意精确属性值过滤。

| 任务 | 命令 | 使用的键 |
|------|------|----------|
| 精确 feature 查询 | `--id gene0001` | `ID` |
| 批量精确 feature 查询 | `--ids genes.txt` | 每行一个 `ID` 值 |
| 基因查询 | `--name ABC1` | gene 记录中的 `ID`, `gene_id`, `Name`, `locus_tag`, `Alias` 或完整 `Dbxref` 值 |
| 任意精确属性过滤 | `--where Parent=gene0001` | 任意第 9 列 `KEY=VALUE`，包括 `ID`, `Name`, `Alias`, `Parent`, `Dbxref`, `Accession` 或 `Parent_Accession` |
| 最近基因查询 | `--nearest chr1:1000-2000` | 同一 seqid 上距离 1-based 闭合区间最近的 gene |
| 包含匹配记录后代 | `-C`, `--children` | 通过 `Parent` 连接的 child 记录；`--include-children` 是较长别名 |
| 包含匹配记录祖先 | `--parents` | 沿 `Parent` 链向上找到的 parent 记录；`--include-parents` 是较长别名 |
| 提取完整 gene model | `--model`, `--gene-model` | 所属 gene 以及 transcript/exon/CDS/UTR 后代 |
| 打印指定属性 | `--out-attrs ID,Name,Parent` | 记录匹配后选择的第 9 列键 |

`--attr KEY=VALUE` 是 `--where KEY=VALUE` 的兼容别名。`--output-attrs` 是 `--out-attrs` 的较长别名。`--attrs` 保留为已弃用的兼容别名。

### Grep 和表达式过滤

`gffsub` 按 GFF 语义 subtract 记录：列、属性、ID、`Parent`/child 关系、gene model、转录本结构和 QC 状态。字段级过滤可以用 grep 风格模式快速完成；逻辑更复杂时，用表达式把筛选条件写清楚。

```bash
./gffsub annotation.gff3 --grep ID:Glyma.01G
./gffsub annotation.gff3 --grep-file genes.txt --grep-field ID
./gffsub annotation.gff3 --grep-regex 'ID:^Glyma\.01G'
./gffsub annotation.gff3 --grep-regex 'seqid:^chr[0-9]+$' -f gene
./gffsub annotation.gff3 -I 'type=="gene" && attr.biotype=="protein_coding"'
./gffsub annotation.gff3 -I '(type=="gene" && length>=1000) || attr.ID~"^Glyma\.01G"'
./gffsub annotation.gff3 -E 'attr.Note~"transposon|retroelement"'
```

Grep 字段可以是 GFF 核心列（`seqid`, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`, `length`, `attrs`），也可以是属性（`ID`, `Name`, `Parent`, `Alias`, `Dbxref`, `Note`, `biotype`, `gene_id`, `transcript_id`, `locus_tag` 或 `attr.KEY`）。`--grep` 做子串匹配，`--grep-regex` 使用 ECMAScript 正则，`--grep-file` 读取每个非空行作为一个 pattern，`-v` 反选 grep 结果，`--ignore-case` 对 grep 和表达式中的字符串/正则匹配生效。

表达式过滤使用同一组字段名，支持 `==`, `!=`, `~`, `!~`, `<`, `<=`, `>`, `>=`, `&&`, `||`, `!` 和括号。缺失值按 `.` 比较。

## 场景：提取上游或下游窗口

当你需要某个基因或 feature 周围的局部注释上下文时使用 window 选项，例如启动子区域查看或邻近 feature 审查。

```bash
./gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware
```

不加 `--strand-aware` 时，upstream 表示更小的基因组坐标，downstream 表示更大的基因组坐标。加上 `--strand-aware` 后，上下游按 feature 的链方向解释。

## 场景：每个基因保留一个转录本

当下游工具要求每个基因只有一个代表性转录本时，使用 `--longest`。

```bash
./gffsub annotation.gff3 --longest
./gffsub annotation.gff3 --longest -@ 6
```

最长 isoform 逻辑遵循本项目现有的 AGAT 风格规则：如果基因存在 CDS isoform，则比较 CDS 长度；否则比较 exon 长度。

## 场景：在 Pipeline 前检查注释质量

在把注释送入 graph-aware workflow、基于 ID 的提取或格式转换前，可以先运行 `--qc`。

```bash
./gffsub annotation.gff3 --qc
```

QC 输出 TSV 表，字段为 `severity`, `code`, `line_idx`, `id` 和 `message`。

| 检查类别 | Codes |
|----------|-------|
| Header 和记录形状 | `invalid_gff_version`, `invalid_column_count` |
| Attributes | `invalid_attribute_syntax`, `invalid_attribute_value`, `invalid_attribute_escape`, `duplicate_attribute_tag`, `invalid_attribute_multivalue`, `invalid_percent_encoding` |
| GFF3 核心列 | `invalid_seqid`, `invalid_source`, `invalid_feature_type`, `invalid_coordinate`, `invalid_range`, `invalid_score`, `invalid_strand`, `invalid_phase`, `invalid_cds_phase` |
| Sequence-region directives | `invalid_sequence_region`, `duplicate_sequence_region`, `outside_sequence_region` |
| Feature graph | `duplicate_id`, `duplicate_parent`, `parent_cycle`, `missing_derives_from`, `missing_parent`, `child_outside_parent` |
| Structured attributes | `invalid_dbxref`, `invalid_gap`, `invalid_is_circular`, `invalid_ontology_term`, `invalid_target` |

在严格 GFF3 QC 中，attribute 列必须是 `.` 或以分号分隔的 `tag=value` 字段；每个 `tag=value` 属性都必须有非空值；作为属性内容出现的 ampersand 和双引号必须分别 URL-escape 为 `%26` 和 `%22`；逗号分隔多值只接受 `Parent`, `Alias`, `Note`, `Dbxref` 和 `Ontology_term`；未知 source 列应写成 `.`，不要留空。带 `Is_circular=true` 的 `region` feature 可以让同一 seqid 上的 feature 越过其 `##sequence-region` 终点。

## CLI 参考

`gffsub` 以顶层入口为主：常用 GFF3 工作流从 `gffsub <input.gff3> [options]` 开始。`query`, `window` 和 `qc` 子命令保留为兼容的进阶入口，并共享同一套输出语义。

### 顶层模式

```bash
gffsub <input.gff3> [options]
```

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入 GFF3/GTF 风格注释文件。 |
| `--id` | ID | 保留精确 feature `ID`。该选项可以重复使用。 |
| `--ids`, `--id-list` | 文件 | 每个非空行读取一个精确 feature ID。`--id-list` 是较长别名。 |
| `--name` | key | 保留一个按 `ID`, `Name`, `gene_id`, `locus_tag`, `Alias` 或完整 `Dbxref` 值找到的基因。 |
| `--where`, `--attr` | `KEY=VALUE` | 保留精确 GFF3 属性值匹配的 feature。该选项可以重复使用。 |
| `--grep` | `FIELD:PATTERN` | 保留字段或属性中包含 `PATTERN` 的记录。该选项可以重复使用。 |
| `--grep-regex` | `FIELD:REGEX` | 保留字段或属性匹配 ECMAScript 正则的记录。该选项可以重复使用。 |
| `--grep-file` | 文件 | 每个非空行读取一个 grep pattern。与 `--grep-field` 组合使用。 |
| `--grep-field` | 字段 | `--grep-file` 使用的字段，例如 `ID`, `Name`, `seqid`, `type` 或 `attr.KEY`。 |
| `--grep-file-regex` | 标志 | 将 `--grep-file` 中的每行按正则表达式处理，而不是子串 pattern。 |
| `-I`, `--include-expr` | 表达式 | 保留匹配 GFF 语义表达式的记录。该选项可以重复使用。 |
| `-E`, `--exclude-expr` | 表达式 | 删除匹配 GFF 语义表达式的记录。该选项可以重复使用。 |
| `-v`, `--invert-match` | 标志 | 反转 `--grep`, `--grep-regex` 或 `--grep-file` 的匹配结果。 |
| `--ignore-case` | 标志 | 对 grep 和表达式中的字符串/正则匹配启用大小写不敏感。 |
| `-C`, `--children`, `--include-children` | 标志 | 包含由 `--id`, `--ids`, `--name`, `--where` 或 `--nearest` 匹配记录的后代。 |
| `--parents`, `--include-parents` | 标志 | 包含由 `--id`, `--ids`, `--name`, `--where` 或 `--nearest` 匹配记录的祖先。 |
| `--model`, `--gene-model` | 标志 | 包含匹配记录所属的完整 gene model。 |
| `--nearest`, `--nearest-gene` | `CHR:START-END` | 保留同一 seqid 上离 1-based 闭合区间最近的 gene。 |
| `--out-attrs`, `--output-attrs` | `KEY1,KEY2,...` | 将指定第 9 列属性作为额外 TSV/JSON 字段输出。只与 query-style selector 组合。 |
| `--attrs` | `KEY1,KEY2,...` | `--out-attrs` 的已弃用兼容别名。 |
| `--summary`, `--summary-format` | `tsv`, `json` | 输出 summary 行，而不是 GFF3 记录。只与 query-style selector 组合。 |
| `--up`, `--upstream` | 整数 | 与 `--id` 配合，提取与目标上游扩展窗口重叠的记录。 |
| `--down`, `--downstream` | 整数 | 与 `--id` 配合，提取与目标下游扩展窗口重叠的记录。 |
| `--strand-aware` | 标志 | 窗口提取时，按 feature 链方向解释 biological upstream/downstream。 |
| `--qc` | 标志 | 运行注释 QC。 |
| `--seqid` | seqid | 保留第 1 列 GFF3 seqid 精确等于该值的记录。 |
| `--source` | source | 保留第 2 列 GFF3 source 精确等于该值的记录。 |
| `--score` | number, `.` | 保留第 6 列 GFF3 score 匹配该数值，或 score 列为 `.` 的记录。 |
| `--strand` | `+`, `-`, `.`, `?` | 保留第 7 列 GFF3 strand 精确等于该值的记录。 |
| `--phase` | `0`, `1`, `2`, `.` | 保留第 8 列 GFF3 phase 精确等于该值的记录。 |
| `-r`, `--region` | `CHR:START-END` | 保留与 1-based 闭合区间重叠的 feature。 |
| `-b`, `--bed` | 文件 | 保留与 BED 区间重叠的 feature；BED 按 0-based 半开区间读取。 |
| `-f`, `--feature`, `--type` | 类型 | 只保留第三列匹配该 feature type 的记录，例如 `gene`, `mRNA`, `transcript`, `exon` 或 `CDS`。 |
| `-L`, `--longest` | 标志 | 每个基因保留最长 transcript isoform。 |
| `-@`, `--threads` | 整数 | 设置 `--longest` 使用的线程数；超过 256 会被限制为 256。 |
| `-t`, `--format`, `--output-format` | `gff3`, `gtf`, `gtf2`, `gtf3`, `bed` | 选择输出格式。`gtf` 按 `gtf2` 处理；默认是 `gff3`。 |
| `-o`, `--output` | 文件 | 写入文件，而不是 stdout。 |
| `-h`, `--help` | 标志 | 显示顶层模式帮助。 |

### 兼容子命令

```bash
gffsub query <input.gff3> [options]
gffsub window <input.gff3> --id ID [options]
gffsub qc <input.gff3>
```

使用 `query` 表达显式 query-style 命令行，使用 `window` 兼容旧的上下游窗口脚本，使用 `qc` 作为 `--qc` 的子命令形式。多数工作流都可以用上面的顶层选项完成。

兼容说明：

- `query` 接受上面的 selector 选项，也支持用 `--type` 按 feature type 限制查询输出。
- `window` 必须提供 `--id`；命令会先尝试精确 `ID`，再做基因查找。`--up` 和 `--down` 默认是 `0`，且必须为非负整数。
- `qc` 等价于顶层 `--qc`，并输出 QC 场景中描述的同一组 TSV 字段。

## 输出格式

| 格式 | Header | 坐标 |
|------|--------|------|
| gff3 | `##gff-version 3` | 1-based 闭合区间 |
| gtf2 | `##gtf-version 2` | 1-based 闭合区间 |
| gtf3 | `##gtf-version 2.2.1` | 1-based 闭合区间 |
| bed | 无 header | 0-based 半开区间 |

## C++ API

公开 C++ API 以 `AnnotationIndex` 为中心：

```cpp
auto index = gffsub::AnnotationIndex::from_gff3("annotation.gff3");
auto gene = index.find_gene("GeneA");
auto model = index.gene_model("GeneA");
auto hits = index.overlap("chr1", 1000, 2000);
```

可用查询包括 `find_by_id`、`find_gene`、`parents_of`、`children_of`、`descendants_of`、`gene_model`、`overlap`、`nearest_gene` 和 `with_attribute`。

## 分发

`gffsub` 会构建为单个命令行二进制文件。将二进制复制到 glibc-based Linux x86_64 机器后，即可配合注释文件运行。

## 许可证

MIT License
