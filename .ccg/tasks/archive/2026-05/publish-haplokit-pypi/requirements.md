# Publish haplokit to PyPI

## User Request

Publish the `haplokit` Python package to PyPI.

## Acceptance Criteria

- Verify package metadata and version before upload.
- Build source distribution and wheel from the repository state intended for release.
- Run package validation (`twine check`) and relevant project tests where feasible.
- Confirm PyPI credentials without exposing secrets.
- Upload only if the target version is not already published and credentials are available.
- Record release result and archive this CCG task.

## Constraints

- Publishing to PyPI is irreversible for a version; do not upload a duplicate or unverified artifact.
- Project metadata targets POSIX/Linux and the build invokes CMake for a C++ backend.
- Prefer WSL/Linux for build validation and release packaging.
- Do not print or persist PyPI tokens.
