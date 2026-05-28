# Requirements

- Keep the project split where backend phenotype loading/statistics prepares structured records and plotting renders from those records.
- Add a population-stratified phenotype example based on the existing `data/popgroup.txt` samples.
- Generate a representative boxplot image through the official haplokit plotting API.
- Keep existing phenotype boxplot behavior compatible when no population file is provided.
- Treat missing phenotype values like GCTA GWAS workflows: ignore them per trait and report effective non-missing sample counts.
