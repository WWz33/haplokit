# Review

## Result

- Synced `deps/gffsub` source, tests, build files, and README/license from `F:/codex/gffsub`.
- Updated haplokit top-level CMake integration from old `gff3_*` source names to the refactored `annotation_*`, `attributes`, `region`, and `gtf_parser` units.
- Verified haplokit call sites still align with `gffsub::AnnotationIndex`, `GffRecord`, and `window_region`.

## Validation

- WSL CMake build in `build-gffsub-align` passed.
- `ctest --test-dir build-gffsub-align --output-on-failure` passed, 3/3.
- WSL Python tests `test_haplokit_real_vcf.py`, `test_haplokit_cli_contract.py`, and `test_packaging_contract.py` passed, 33/33.
