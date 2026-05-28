# Review

## Automated/Agent Review

Test-review agent reported no Critical findings.

Resolved warnings:
- Made `tests/python/test_phenotype.py` self-contained by inserting repo root into `sys.path`.
- Added coverage for summary-style `Accession/freq` parsing.
- Added CLI rejection coverage for malformed `--comparison`.

Noted but not applicable:
- The agent's WSL warning used system `pytest`; validation here used `/home/ww/miniforge3/bin/python -m pytest` and passed.

## Manual Review

Critical: none found.

Warnings addressed:
- Tukey HSD p-values are already adjusted; `pairwise_statistics()` now avoids applying Bonferroni again when `method="tukey"`.
- `--phenotype` was ambiguous with `--trait`; canonical CLI now uses `--phenotypes` / `--pheno-file`, with `--phenotype` retained as a compatibility alias.

Residual risk:
- Bar/pie composition plots from the reference projects are intentionally out of scope for this pass.
- `ruff` was not available in either Windows PATH or the configured WSL Python environment, so lint was not run.

Verification:
- Windows: `python -m pytest -q tests\python\test_phenotype.py tests\python\test_haplokit_cli_contract.py tests\python\test_python_only_plotting_contract.py tests\python\test_packaging_contract.py --basetemp=.pytest-tmp-phenotype -p no:cacheprovider` -> 22 passed, 1 skipped.
- WSL: `/home/ww/miniforge3/bin/python -m pytest -q tests/python --basetemp=.pytest-tmp-python-wsl -p no:cacheprovider` -> 54 passed.

