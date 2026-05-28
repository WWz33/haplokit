# Review

## Critical

- None.

## Warning

- README now documents both conda/mamba and distro-native dependency installs, but the project is still source-build only on Linux/WSL. Users on macOS or Windows without WSL should not expect `pip install haplokit` to work.

## Info

- `haplokit/_backend.py` now injects `CONDA_PREFIX`-derived include and library paths into the CMake subprocesses used by pip builds.
- `cmake/htslib.cmake` now resolves `z`, `m`, `bz2`, `lzma`, and `curl` with `find_library()` when vendored htslib is used, which avoids bare `-lcurl` failures in conda/mamba environments.
- Verified with:
  - `pytest tests/python/test_packaging_contract.py -q`
  - `wsl bash -lc "cd /mnt/f/codex/genehapr-master && /home/ww/miniforge3/bin/python -m pytest -q tests/python"`
  - `wsl bash -lc "cd /mnt/f/codex/genehapr-master && cmake -S . -B build-conda-deps-test -DCMAKE_BUILD_TYPE=Release && cmake --build build-conda-deps-test --parallel 2 && ctest --test-dir build-conda-deps-test --output-on-failure"`
