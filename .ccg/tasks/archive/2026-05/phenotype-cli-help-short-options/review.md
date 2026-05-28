# Review

Critical: none found.

Changes verified:
- `view`, `phenotype`, `phenotype stat`, and `phenotype box` now expose argparse descriptions/help.
- Existing `view` short options are preserved.
- Added `-C` for `view --map-facecolor`.
- Added short options for phenotype file/stat/box arguments, including `-H`, `-P`, `-s`, `-m`, `-M`, `-a`, `-d`, `-D`, `-G`, `-c`, and `-T` where applicable.
- `--min-hap-size` remains a phenotype-module shared option, not a top-level global option.

Verification:
- Windows parser/help targeted tests: 8 passed, 1 skipped.
- WSL full `tests/python`: 58 passed.

