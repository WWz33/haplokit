# pip auto-build backend

## User Request

After `git clone`, `pip` installation should automatically build the `haplokit_cpp`
C++ backend so `haplokit view ...` works without a manual `cmake` command.

## Acceptance Criteria

- Regular `pip install .` continues to build and package `haplokit_cpp`.
- Editable/source-tree usage can build `haplokit_cpp` automatically or provide a clear actionable error.
- Runtime fallback does not silently swallow CMake failures.
- Installation documentation explains `pip install .`, `pip install -e .`, and required system toolchain.
- Packaging contract tests cover the new behavior.

## Constraints

- Source builds target POSIX/Linux; Windows users should use WSL/Linux.
- Do not introduce new build backends unless needed.
- Keep changes surgical and compatible with existing setuptools configuration.
