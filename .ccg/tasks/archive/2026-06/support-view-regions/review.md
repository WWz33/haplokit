# Review

## External Review

- Critical: none found.
- Warning: `-T/--targets-file` was initially expanded into a comma-joined argv string, risking command-line length limits for large files.
  - Resolution: changed Python CLI to pass targets-file path to the backend via `view-targets-file-json`; backend now reads target files directly.
- Warning: docs did not state that target selectors must be on one chromosome.
  - Resolution: README and README.zh-CN now document the same-chromosome requirement.
- Info: direct C++ backend behavior should match Python same-chromosome semantics.
  - Resolution: `view-targets-json` computes a covering region for all target runs, enforcing same-chromosome targets.

Second short review agent timed out twice and was closed without output.

## Verification

- `wsl bash -lc "cd /mnt/f/codex/genehapr-master && cmake --build build-wsl-targets --parallel 4 && ctest --test-dir build-wsl-targets --output-on-failure"`
- `python -m py_compile haplokit\cli.py tests\python\test_haplokit_cli_contract.py`
- `pytest tests\python\test_summary_contract.py -q`
- WSL Python `haplokit.cli.main` checks for `-t scaffold_1:4300,scaffold_1:4950` and `-T /tmp/haplokit_targets.txt`, both returning sites `[4300, 4950]`.
