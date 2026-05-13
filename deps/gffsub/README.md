# gffsub

<!-- README-I18N:START -->

**English** | [汉语](./README.zh.md)

<!-- README-I18N:END -->

A fast GFF3/GTF extraction and filtering tool written in C++.

## Features

- **Region extraction**: Extract features by genomic region
- **BED support**: Filter using BED files
- **Longest isoform**: Keep only the longest transcript per gene
- **Multi-format output**: GFF3, GTF2, GTF3, and BED formats

## Installation

### Requirements

- C++17 compiler (g++ or clang++)

### Build

```bash
cd gffsub
make
```

## Usage

```
./gffsub <input.gff3> [options]
```

### Input/Region Options

| Short | Long | Description |
|-------|------|-------------|
| `-r` | `--region CHR:START-END` | Extract region (1-based inclusive) |
| `-b` | `--bed FILE` | BED file for regions (0-based → 1-based) |

### Feature Filter Options

| Short | Long | Description |
|-------|------|-------------|
| `-f` | `--feature TYPE` | Filter by feature type (gene, mRNA, exon, CDS...) |
| `-L` | `--longest` | Keep longest transcript per gene |

### Output Options

| Short | Long | Description |
|-------|------|-------------|
| `-t` | `--output-format FMT` | Format: gff3, gtf2, gtf3, bed |
| `-o` | `--output FILE` | Output file (default: stdout) |
| `-h` | `--help` | Show help |

## Examples

```bash
# Extract genes in region
./gffsub annotation.gff3 -r chr1:1-100000 -f gene

# Extract using BED file
./gffsub annotation.gff3 -b regions.bed -f exon

# Keep only longest isoform
./gffsub annotation.gff3 --longest

# Output as GTF3
./gffsub annotation.gff3 -r chr1:1-100000 -t gtf3 -o out.gtf
```

## Distribution

Single static binary - copy to any Linux machine and run:

```bash
./gffsub annotation.gff3 -r chr1:1-100000 -f gene
```

**Requirements:** glibc-based Linux (x86_64)

## Output Formats

| Format | Header | Coordinate |
|--------|--------|------------|
| gff3 | `##gff-version 3` | 1-based inclusive |
| gtf2 | `##gtf-version 2` | 1-based inclusive |
| gtf3 | `##gtf-version 2.2.1` | 1-based inclusive |
| bed | (no header) | 0-based half-open |

## License

MIT License
