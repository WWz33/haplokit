# gffsub

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

`gffsub` is a command-line tool for day-to-day genome annotation work with GFF3/GTF-style files. It helps you subset annotations, recover complete gene models from identifiers, build upstream/downstream windows, keep representative transcripts, and run quick QC before downstream pipelines.

Use it when a plain interval filter is not enough and you need GFF3-aware behavior such as column-9 attribute lookup, `Parent`/child traversal, gene model extraction, and annotation-specific quality checks.

## Start With Your Task

| I want to... | Use this |
|-------------|----------|
| Extract genes in a genomic interval | `gffsub annotation.gff3 -r chr1:1-100000 -f gene` |
| Extract records from one chromosome or contig | `gffsub annotation.gff3 --seqid chr1` |
| Extract records from one annotation source | `gffsub annotation.gff3 --source Gnomon` |
| Filter records by score | `gffsub annotation.gff3 --score 42.5` |
| Filter records by strand | `gffsub annotation.gff3 --strand -` |
| Filter CDS records by phase | `gffsub annotation.gff3 --phase 0 -f CDS` |
| Use BED intervals as input | `gffsub annotation.gff3 -b regions.bed -f exon` |
| Find one feature by exact `ID` | `gffsub annotation.gff3 --id GeneA` |
| Extract many exact IDs | `gffsub annotation.gff3 --ids genes.txt` |
| Find a gene by common naming keys | `gffsub annotation.gff3 --name GeneA` |
| Find records by any exact attribute value | `gffsub annotation.gff3 --where biotype=protein_coding` |
| Grep a field or attribute with a pattern file | `gffsub annotation.gff3 --grep-file genes.txt --grep-field ID` |
| Grep a field or attribute with regex | `gffsub annotation.gff3 --grep-regex 'ID:^Glyma\.01G'` |
| Combine semantic filters in an expression | `gffsub annotation.gff3 -I 'type=="gene" && attr.biotype=="protein_coding"'` |
| Exclude records by a semantic expression | `gffsub annotation.gff3 -E 'attr.Note~"transposon|retroelement"'` |
| Find the nearest gene to a region | `gffsub annotation.gff3 --nearest chr1:1000-2000` |
| Include descendants of matched records | `gffsub annotation.gff3 --id GeneA -C` |
| Include ancestors of matched records | `gffsub annotation.gff3 --id ExonA --parents` |
| Recover the full gene model from any feature | `gffsub annotation.gff3 --id ExonA --model` |
| Produce a pipeline-friendly summary | `gffsub annotation.gff3 --id GeneA --summary tsv` |
| Extract selected attribute values | `gffsub annotation.gff3 --id GeneA --out-attrs ID,Name,Parent` |
| Extract upstream/downstream context | `gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware` |
| Keep the longest transcript per gene | `gffsub annotation.gff3 --longest` |
| Check annotation syntax and graph problems | `gffsub annotation.gff3 --qc` |

## Supported Inputs And Outputs

| Kind | Supported |
|------|-----------|
| Annotation input | GFF3/GTF-style feature records |
| Region input | `CHR:START-END` strings and BED files |
| Identifier input | repeated `--id` values or one-ID-per-line files with `--ids` |
| Pattern input | one-pattern-per-line files with `--grep-file` and `--grep-field` |
| Annotation output | `gff3`, `gtf` (= `gtf2`), `gtf2`, `gtf3`, `bed` |
| Tabular output | TSV/JSON summaries and TSV QC reports |

`gffsub` works on annotation records. It does not process sequence FASTA or variant VCF content as primary input.

## Install

### Requirements

- C++17 compiler (`g++` or `clang++`)
- CMake is optional; the repository also includes a `Makefile`

### Build With Make

```bash
cd gffsub
make
```

### Build And Test With CMake

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Scenario: Subset Annotation By Genomic Context

Use this mode when your question starts from a coordinate interval, chromosome/contig, source column, or BED file.

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

Coordinate rules are explicit:

| Input or output | Coordinate system |
|-----------------|-------------------|
| GFF3/GTF records | 1-based inclusive |
| `CHR:START-END` regions | 1-based inclusive |
| BED input/output | 0-based half-open |

Options compose left to right by filtering the current record set. For example, `-r chr1:1-100000 -f gene -t bed` keeps region-overlapping records, limits them to `gene`, then prints BED coordinates.

## Scenario: Find Genes And Recover Gene Models

Use selector options when your question starts from a feature ID, gene name, attribute, or nearby locus.

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

`--nearest` searches genes on the same seqid as a 1-based inclusive region. Overlapping genes have distance `0`; ties use input file order.

For batch pipelines, ask for summaries instead of raw GFF3:

```bash
./gffsub annotation.gff3 --ids genes.txt --summary tsv
./gffsub annotation.gff3 --id GeneA --summary json
./gffsub annotation.gff3 --id gene0001 --out-attrs ID,Name,Alias,Dbxref
```

Summary fields include `query_id`, `matched_id`, `matched_by`, `seqid`, `start`, `end`, `strand`, `type`, `parent_id`, `child_count`, `transcript_count`, `exon_count`, `cds_length`, and `status`. If `--out-attrs` is present, selected column-9 keys are appended as TSV columns or emitted under `attrs` in JSON.

### Attribute Keys In Search And Output

GFF3 stores record attributes in column 9 as semicolon-separated `KEY=VALUE` pairs:

```gff3
chr1	src	gene	100	400	.	+	.	ID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123
```

Use `--id` for exact `ID` lookup, `--name` for gene lookup across common naming keys, and `--where KEY=VALUE` for any exact attribute-value filter.

| Task | Command | Keys used |
|------|---------|-----------|
| Exact feature lookup | `--id gene0001` | `ID` |
| Batch exact feature lookup | `--ids genes.txt` | `ID` values, one per line |
| Gene lookup | `--name ABC1` | gene records by `ID`, `gene_id`, `Name`, `locus_tag`, `Alias`, or full `Dbxref` value |
| Any exact attribute filter | `--where Parent=gene0001` | any column-9 `KEY=VALUE`, including `ID`, `Name`, `Alias`, `Parent`, `Dbxref`, `Accession`, or `Parent_Accession` |
| Nearest gene lookup | `--nearest chr1:1000-2000` | same-seqid gene with the shortest distance to a 1-based inclusive region |
| Include matched descendants | `-C`, `--children` | child records linked by `Parent`; `--include-children` is a verbose alias |
| Include matched ancestors | `--parents` | parent records reached by walking `Parent` links upward; `--include-parents` is a verbose alias |
| Extract full gene model | `--model`, `--gene-model` | containing gene plus transcript/exon/CDS/UTR descendants |
| Print selected attributes | `--out-attrs ID,Name,Parent` | selected column-9 keys after records are matched |

`--attr KEY=VALUE` is a compatibility alias for `--where KEY=VALUE`. `--output-attrs` is a verbose alias for `--out-attrs`. `--attrs` remains as a deprecated compatibility alias.

### Grep And Expression Filters

`gffsub` subtracts records by GFF semantics: columns, attributes, IDs, `Parent`/child links, gene models, transcript structure, and QC status. For field-level filtering, use grep-style patterns for quick tasks and expression filters when the logic needs to be explicit.

```bash
./gffsub annotation.gff3 --grep ID:Glyma.01G
./gffsub annotation.gff3 --grep-file genes.txt --grep-field ID
./gffsub annotation.gff3 --grep-regex 'ID:^Glyma\.01G'
./gffsub annotation.gff3 --grep-regex 'seqid:^chr[0-9]+$' -f gene
./gffsub annotation.gff3 -I 'type=="gene" && attr.biotype=="protein_coding"'
./gffsub annotation.gff3 -I '(type=="gene" && length>=1000) || attr.ID~"^Glyma\.01G"'
./gffsub annotation.gff3 -E 'attr.Note~"transposon|retroelement"'
```

Grep fields can be core GFF columns (`seqid`, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`, `length`, `attrs`) or attributes (`ID`, `Name`, `Parent`, `Alias`, `Dbxref`, `Note`, `biotype`, `gene_id`, `transcript_id`, `locus_tag`, or `attr.KEY`). `--grep` does substring matching, `--grep-regex` uses ECMAScript regular expressions, `--grep-file` reads one pattern per non-empty line, `-v` inverts grep matches, and `--ignore-case` applies to grep and expression string matches.

Expression filters use the same field names and support `==`, `!=`, `~`, `!~`, `<`, `<=`, `>`, `>=`, `&&`, `||`, `!`, and parentheses. Missing values compare as `.`.

## Scenario: Extract Upstream Or Downstream Windows

Use window options when you need local annotation context around a gene or feature, such as promoter inspection or neighboring-feature review.

```bash
./gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware
```

Without `--strand-aware`, upstream means lower genomic coordinates and downstream means higher genomic coordinates. With `--strand-aware`, upstream/downstream follows the feature strand.

## Scenario: Keep One Transcript Per Gene

Use `--longest` when a downstream tool expects one representative transcript per gene.

```bash
./gffsub annotation.gff3 --longest
./gffsub annotation.gff3 --longest -@ 6
```

The longest isoform logic follows the existing AGAT-style rule in this project: if a gene has CDS isoforms, compare CDS length; otherwise compare exon length.

## Scenario: Check Annotation Quality Before A Pipeline

Run `--qc` before feeding annotations into graph-aware workflows, ID-based extraction, or format conversion.

```bash
./gffsub annotation.gff3 --qc
```

QC writes a TSV table with `severity`, `code`, `line_idx`, `id`, and `message`.

| Check family | Codes |
|--------------|-------|
| Header and record shape | `invalid_gff_version`, `invalid_column_count` |
| Attributes | `invalid_attribute_syntax`, `invalid_attribute_value`, `invalid_attribute_escape`, `duplicate_attribute_tag`, `invalid_attribute_multivalue`, `invalid_percent_encoding` |
| Core GFF3 columns | `invalid_seqid`, `invalid_source`, `invalid_feature_type`, `invalid_coordinate`, `invalid_range`, `invalid_score`, `invalid_strand`, `invalid_phase`, `invalid_cds_phase` |
| Sequence-region directives | `invalid_sequence_region`, `duplicate_sequence_region`, `outside_sequence_region` |
| Feature graph | `duplicate_id`, `duplicate_parent`, `parent_cycle`, `missing_derives_from`, `missing_parent`, `child_outside_parent` |
| Structured attributes | `invalid_dbxref`, `invalid_gap`, `invalid_is_circular`, `invalid_ontology_term`, `invalid_target` |

In strict GFF3 QC, the attribute column must be `.` or semicolon-separated `tag=value` fields; each `tag=value` attribute must have a non-empty value; ampersands and double quotes used as attribute content must be URL-escaped as `%26` and `%22`; comma-separated values are accepted only for `Parent`, `Alias`, `Note`, `Dbxref`, and `Ontology_term`; and an unknown source column should be written as `.` rather than left empty. A `region` feature marked `Is_circular=true` may make features on that seqid wrap past the end of their `##sequence-region`.

## CLI Reference

`gffsub` is top-level first: common GFF3 work starts as `gffsub <input.gff3> [options]`. The `query`, `window`, and `qc` subcommands remain as compatible advanced entry points and share the same output semantics.

### Top-Level Mode

```bash
gffsub <input.gff3> [options]
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input GFF3/GTF-style annotation file. |
| `--id` | ID | Keep the exact feature `ID`. This option can be repeated. |
| `--ids`, `--id-list` | file | Read one exact feature ID per non-empty line. `--id-list` is a verbose alias. |
| `--name` | key | Keep one gene found by `ID`, `Name`, `gene_id`, `locus_tag`, `Alias`, or full `Dbxref` value. |
| `--where`, `--attr` | `KEY=VALUE` | Keep features with an exact GFF3 attribute value. This option can be repeated. |
| `--grep` | `FIELD:PATTERN` | Keep records whose field or attribute contains `PATTERN`. This option can be repeated. |
| `--grep-regex` | `FIELD:REGEX` | Keep records whose field or attribute matches an ECMAScript regular expression. This option can be repeated. |
| `--grep-file` | file | Read one grep pattern per non-empty line. Combine with `--grep-field`. |
| `--grep-field` | field | Field used by `--grep-file`, such as `ID`, `Name`, `seqid`, `type`, or `attr.KEY`. |
| `--grep-file-regex` | flag | Treat `--grep-file` lines as regular expressions instead of substring patterns. |
| `-I`, `--include-expr` | expression | Keep records matching a GFF semantic expression. This option can be repeated. |
| `-E`, `--exclude-expr` | expression | Drop records matching a GFF semantic expression. This option can be repeated. |
| `-v`, `--invert-match` | flag | Invert `--grep`, `--grep-regex`, or `--grep-file` matches. |
| `--ignore-case` | flag | Apply case-insensitive matching to grep and expression string/regex matches. |
| `-C`, `--children`, `--include-children` | flag | Include descendants of records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`. |
| `--parents`, `--include-parents` | flag | Include ancestors of records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`. |
| `--model`, `--gene-model` | flag | Include the full gene model containing matched records. |
| `--nearest`, `--nearest-gene` | `CHR:START-END` | Keep the nearest gene on the same seqid as a 1-based inclusive region. |
| `--out-attrs`, `--output-attrs` | `KEY1,KEY2,...` | Print selected column-9 attributes as extra TSV/JSON fields. Combine only with query-style selectors. |
| `--attrs` | `KEY1,KEY2,...` | Deprecated compatibility alias for `--out-attrs`. |
| `--summary`, `--summary-format` | `tsv`, `json` | Print summary rows instead of GFF3 records. Combine only with query-style selectors. |
| `--up`, `--upstream` | integer | With `--id`, extract records overlapping the upstream-expanded target window. |
| `--down`, `--downstream` | integer | With `--id`, extract records overlapping the downstream-expanded target window. |
| `--strand-aware` | flag | With window extraction, interpret upstream/downstream biologically by feature strand. |
| `--qc` | flag | Run annotation QC. |
| `--seqid` | seqid | Keep records whose first GFF3 column exactly matches the value. |
| `--source` | source | Keep records whose second GFF3 column exactly matches the value. |
| `--score` | number, `.` | Keep records whose sixth GFF3 column matches the numeric score, or whose score column is `.`. |
| `--strand` | `+`, `-`, `.`, `?` | Keep records whose seventh GFF3 column exactly matches the value. |
| `--phase` | `0`, `1`, `2`, `.` | Keep records whose eighth GFF3 column exactly matches the value. |
| `-r`, `--region` | `CHR:START-END` | Keep features overlapping a 1-based inclusive region. |
| `-b`, `--bed` | file | Keep features overlapping BED intervals; BED is read as 0-based half-open. |
| `-f`, `--feature`, `--type` | type | Keep only records whose third column matches the feature type, such as `gene`, `mRNA`, `transcript`, `exon`, or `CDS`. |
| `-L`, `--longest` | flag | Keep the longest transcript isoform per gene. |
| `-@`, `--threads` | integer | Set worker threads for `--longest`; values above 256 are capped. |
| `-t`, `--format`, `--output-format` | `gff3`, `gtf`, `gtf2`, `gtf3`, `bed` | Select output format. `gtf` is accepted as `gtf2`; default is `gff3`. |
| `-o`, `--output` | file | Write output to a file instead of stdout. |
| `-h`, `--help` | flag | Show help for top-level mode. |

### Compatibility Subcommands

```bash
gffsub query <input.gff3> [options]
gffsub window <input.gff3> --id ID [options]
gffsub qc <input.gff3>
```

Use `query` for explicit query-style command lines, `window` for older upstream/downstream scripts, and `qc` for the subcommand form of `--qc`. Most workflows can be written with the top-level options shown above.

Compatibility notes:

- `query` accepts the selector options above and also supports `--type` for restricting query output by feature type.
- `window` requires `--id`; it first tries an exact `ID`, then gene lookup. `--up` and `--down` default to `0` and must be non-negative.
- `qc` is equivalent to top-level `--qc` and writes the same TSV fields described in the QC scenario.

## Output Formats

| Format | Header | Coordinate |
|--------|--------|------------|
| gff3 | `##gff-version 3` | 1-based inclusive |
| gtf2 | `##gtf-version 2` | 1-based inclusive |
| gtf3 | `##gtf-version 2.2.1` | 1-based inclusive |
| bed | no header | 0-based half-open |

## C++ API

The public C++ API centers on `AnnotationIndex`:

```cpp
auto index = gffsub::AnnotationIndex::from_gff3("annotation.gff3");
auto gene = index.find_gene("GeneA");
auto model = index.gene_model("GeneA");
auto hits = index.overlap("chr1", 1000, 2000);
```

Available queries include `find_by_id`, `find_gene`, `parents_of`, `children_of`, `descendants_of`, `gene_model`, `overlap`, `nearest_gene`, and `with_attribute`.

## Distribution

`gffsub` builds as a single command-line binary. Copy the binary to a glibc-based Linux x86_64 machine and run it with your annotation files.

## License

MIT License
