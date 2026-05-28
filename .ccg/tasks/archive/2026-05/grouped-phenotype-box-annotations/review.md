# Review

## Automated checks

- `python -m pytest -q tests\python\test_phenotype.py --basetemp=.pytest-tmp-phenotype-grouped-box -p no:cacheprovider`: 17 passed.
- `python -m pytest -q tests\python --basetemp=.pytest-tmp-python-full -p no:cacheprovider`: 44 passed, 2 skipped. Existing Windows GBK subprocess decode warnings remain in network tests.
- `git diff --check`: passed.

## Review findings addressed

- CLI help still described population boxplots as facets; updated to grouped population boxes.
- Tests only checked for any star annotation; added focused annotation-count tests for grouped multiple-haplotype and single-haplotype population paths.
