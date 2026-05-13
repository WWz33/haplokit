# haptools PyPI Release (Linux sdist lane)

## Scope

- Target package name: `haptools` (fallback naming policy must be decided before first public publish if name is unavailable).
- Release platform scope: Linux-first.
- Artifact type: source distribution (`sdist`).

## One-Time Repository Setup

1. Enable GitHub Actions for the repository.
2. Configure PyPI trusted publishing for this repository/workflow:
   - project: `haptools`
   - owner/repo: your GitHub repo
   - workflow file: `.github/workflows/release.yml`
   - environment (optional): `pypi`
3. Ensure `deps/htslib/htscodecs` is available in the repository state used by CI (`submodules: recursive` is already enabled in workflow files).

## Standard Pre-Release Checklist

1. Version in `pyproject.toml` is updated.
2. Local Linux validation passes:
   - `cmake -S . -B build-wsl`
   - `cmake --build build-wsl --parallel`
   - `HAPTOOLS_CPP_BIN=$PWD/build-wsl/haptools_cpp python -m pytest -q tests/python`
   - `ctest --test-dir build-wsl --output-on-failure`
3. Source distribution builds and passes install smoke:
   - `python -m build --sdist`
   - `python -m venv .venv-sdist-smoke && source .venv-sdist-smoke/bin/activate`
   - `python -m pip install -U pip`
   - `python -m pip install dist/*.tar.gz`
   - `REPO_ROOT=$PWD && TMP_DIR=$(mktemp -d) && cd "$TMP_DIR"`
   - `haptools view "$REPO_ROOT/data/var.sorted.vcf.gz" -r scaffold_1:4300-5000 --output-file smoke_out`
   - `test -f smoke_out/hapresult.tsv && test -f smoke_out/hap_summary.tsv`
4. No default runtime/CI lane depends on R.

## Release Steps

1. Create and push a semver-style tag:

```bash
git tag v0.1.0
git push origin v0.1.0
```

2. GitHub `release` workflow runs automatically:
   - builds backend
   - runs Python and C++ tests
   - builds `sdist`
   - smoke-installs from built `sdist` and runs a real `haptools view` command
   - runs `twine check`
   - publishes to PyPI via trusted publishing

3. Verify package on PyPI and perform smoke install on Linux:

```bash
python -m pip install haptools
```

## Release Close-out Checklist

- [ ] Target commit `ci` workflow is green (includes `sdist` install smoke).
- [ ] Local Linux e2e regression completed (`pytest`, `ctest`, `build --sdist`, `twine check`, `sdist` install smoke).
- [ ] Tag pushed and `release` workflow finished successfully.
- [ ] PyPI package page and install smoke verification completed.

## Failure Handling

- Build or test failure: fix on branch, retag with next version.
- `sdist` smoke failure: inspect packaged `_bin/haptools_cpp` presence and rebuild from clean tree before retag.
- Publish failure due to PyPI setup: fix trusted publishing configuration and re-run workflow.
- Package-name conflict: apply the agreed fallback package name policy and update metadata/workflows before retry.
