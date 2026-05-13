# haptools Repository Layout

## Purpose

This document defines which parts of this repository belong to the active `haptools` delivery lane and which parts are retained as legacy reference material.

## Active haptools Surface

These paths are part of the default build, test, CI, and release contract:

- `haptools/`
- `src/cpp/`
- `tests/python/`
- `tests/cpp/`
- `CMakeLists.txt`
- `docs/specs/`
- `docs/plans/`
- `.github/workflows/` (when added)
- Python packaging files at repo root (`pyproject.toml`, `MANIFEST.in`, etc., when added)

## Legacy Reference Surface

These paths are retained for scientific and historical context, but are not part of default `haptools` runtime/CI lanes:

- `R/`
- `tests/testthat/`
- `vignettes/`
- legacy R-package metadata files (`DESCRIPTION`, `NAMESPACE`, `man/`, etc.)

## Default Quality Gates

For `haptools` changes, required verification should come from:

- Python tests in `tests/python/`
- C++ tests in `tests/cpp/`
- packaging/build checks for the Python/C++ install path

R-based tests and R-package tooling are optional reference validation, not default gates for `haptools` delivery.

## Boundary Rules

- New runtime features for `haptools` must land in Python/C++ surfaces, not in `R/`.
- New CI checks for merge and release should target `haptools` active surfaces only.
- Legacy files may be cited for behavior parity, but should not be treated as required runtime dependencies.
