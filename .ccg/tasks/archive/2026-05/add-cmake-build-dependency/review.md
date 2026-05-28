# Review

## Result

- Added `cmake>=3.22` to `build-system.requires`.
- Bumped `haplokit` to `0.1.3`.
- Updated backend build helper to prefer the real CMake binary from the Python `cmake` package when available, avoiding broken PATH stubs.

## Validation

- Local `pytest tests/python/test_packaging_contract.py -q` passed.
- WSL packaging contract tests passed.
- `python -m build --sdist --wheel` passed and showed isolated install of `cmake>=3.22`.
- `twine check` passed.
- Clean WSL venv install from `haplokit-0.1.3.tar.gz` passed and `haplokit view` smoke test produced expected output.
