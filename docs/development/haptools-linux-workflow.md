# haptools Linux Development Workflow

## Goal

Provide one reproducible Linux/WSL workflow for contributors to build and validate `haptools` without R dependencies.

## Prerequisites

- Linux or WSL
- Python 3.12+
- C++17 compiler toolchain
- CMake 3.22+

Example dependency install (Ubuntu-like systems):

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake python3 python3-venv
```

Vendored dependency note:

`haptools` builds `htslib` from `deps/htslib`. Ensure the nested `htscodecs` sources are present:

```bash
git -C deps/htslib submodule update --init --recursive
```

## Environment Setup

```bash
cd /path/to/genehapr-master
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install gffutils pysam pytest
```

## Build C++ Backend

```bash
cmake -S . -B build-wsl
cmake --build build-wsl -j12
```

## Run Tests

Python lane:

```bash
HAPTOOLS_CPP_BIN=$PWD/build-wsl/haptools_cpp python -m pytest -q tests/python
```

C++ lane:

```bash
ctest --test-dir build-wsl --output-on-failure
```

## Run CLI Manually

```bash
HAPTOOLS_CPP_BIN=$PWD/build-wsl/haptools_cpp \
python -m haptools.cli view inst/extdata/var.sorted.vcf.gz \
  -r scaffold_1:4300-5000 \
  --plot \
  --output-file out
```

## Notes

- `haptools view --plot` is Python-side plotting and should run without R.
- `R/` and `tests/testthat/` are legacy/reference surfaces and are not required for default `haptools` validation.
