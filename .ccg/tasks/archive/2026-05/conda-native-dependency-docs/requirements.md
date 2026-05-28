# Requirements

- Document Linux native build dependencies needed by haplokit source installs.
- Explain conda/mamba packages and linker-path exports for `-lcurl`, `-lbz2`, `-llzma`, and `-lz` errors.
- Make pip auto-builds prefer active conda/mamba native include and library paths when `CONDA_PREFIX` is set.
- Preserve Linux-only package scope and existing C++ backend discovery behavior.
