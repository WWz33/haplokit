# Review

## Result

- Fixed regular and editable pip installs to build the C++ backend automatically.
- Fixed runtime backend discovery to search packaged binaries, repo build dirs, `build-python-package`, and `build-haplokit-python`.
- Fixed CMake cache path mismatch by falling back to a Python-specific build dir.
- Marked wheels as platform-specific when native binaries are bundled.

## Validation

- `pytest tests/python/test_packaging_contract.py -q` passed locally.
- WSL `/home/ww/miniforge3/bin/python -m pytest tests/python/test_packaging_contract.py -q` passed.
- WSL `pip install -e . --no-deps` passed.
- WSL sdist install from `haplokit-0.1.2.tar.gz` passed and `haplokit view` produced expected outputs.
