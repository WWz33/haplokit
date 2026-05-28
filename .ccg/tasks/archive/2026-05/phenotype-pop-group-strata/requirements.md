# Requirements

Current phenotype statistics compare haplotypes globally. Add population-group support so users can pass a sample-to-population table and run the same haplotype-vs-phenotype tests within each population stratum.

Scope:
- Add `--population/--pop-group` to `haplokit phenotype stat`.
- Parse two-column sample group files compatible with existing `view -p/--population` input.
- When population is provided, run ANOVA and pairwise haplotype tests separately per population.
- Add a `population` column to statistics and summary output.
- Keep global behavior when no population file is provided.
- Update tests and README.

