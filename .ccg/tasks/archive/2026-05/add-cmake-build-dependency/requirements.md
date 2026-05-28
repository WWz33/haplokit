# Add CMake build dependency

## User Report

`pip install haplokit` from PyPI fails while building the sdist because the target
Linux environment does not have a `cmake` executable on `PATH`.

## Acceptance Criteria

- PEP 517 build isolation installs a suitable CMake executable automatically.
- Release version is bumped because PyPI versions are immutable.
- New sdist builds and passes `twine check`.
- GitHub and PyPI are updated.
