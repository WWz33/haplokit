# Requirements

Implement a Python-only haplotype phenotype module for haplokit, using `repo/hastat` and `repo/geneHapR` as functional references while matching haplokit's current file contracts.

Scope:
- Add a `haplokit phenotype` command, with `pheno` alias.
- Support `stat` for haplotype-vs-phenotype statistics.
- Support `box` for phenotype distribution plots by haplotype.
- Read haplotype assignments from existing `hapresult.tsv` or simple two-column `samples,haplotypes` CSV/TSV files.
- Read phenotype tables whose first column is sample ID and remaining columns are numeric traits.
- Filter haplotypes per trait with `--min-hap-size`, default `5`.
- Output TSV statistics with per-pair counts, means, standard deviations, ANOVA p-value, pairwise p-value, Bonferroni-adjusted p-value, and significance labels.
- Export plotting through `haplokit.plot` without R runtime calls.

Out of scope for this pass:
- R-compatible plotting style clones.
- Population composition bar/pie commands.
- Integrating phenotype logic into `view`.

