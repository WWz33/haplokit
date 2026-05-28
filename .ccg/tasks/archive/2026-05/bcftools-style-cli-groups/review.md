# Review

Critical: none found.

Implementation notes:
- Kept haplokit command names: `view`, `phenotype/pheno`, `stat`, `box`.
- Used domain-specific help groups instead of copied external module names.
- `view --help` now groups target, calling, gene annotation, output, visualization, network, and labeling options.
- `phenotype stat/box --help` now groups haplotype/phenotype input, phenotype tests/plots, output, comparison annotation, and parsing options.

Verification:
- Windows parser/help targeted tests: 8 passed, 1 skipped.
- WSL full `tests/python`: 58 passed.

