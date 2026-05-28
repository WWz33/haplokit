# Review

## Critical

- None.

## Warning

- Runtime path injection depends on `CONDA_PREFIX`; users must run haplokit from an activated conda/mamba environment for conda-installed shared libraries to be discovered automatically.

## Info

- `native_runtime_environment()` now prepends conda runtime library directories for packaged backend subprocesses.
- `haplokit view` backend calls and network backend calls now pass the runtime environment.
- README documents the old-install workaround for `libbz2.so.1.0` runtime loader failures.
- Verified with:
  - `pytest tests/python/test_packaging_contract.py -q`
  - `wsl bash -lc "cd /mnt/f/codex/genehapr-master && /home/ww/miniforge3/bin/python -m pytest -q tests/python"`
  - `wsl bash -lc "cd /mnt/f/codex/genehapr-master && cmake -S . -B build-runtime-path-test -DCMAKE_BUILD_TYPE=Release && cmake --build build-runtime-path-test --parallel 2 && ctest --test-dir build-runtime-path-test --output-on-failure"`
  - `wsl bash -lc "cd /mnt/f/codex/genehapr-master && /home/ww/miniforge3/bin/python -m twine check dist/haplokit-0.1.5.tar.gz"`
