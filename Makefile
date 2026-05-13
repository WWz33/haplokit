.PHONY: all build test clean install

CMAKE_BUILD_DIR := build-wsl
CMAKE_JOBS := $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

all: build

build:
	cmake -S . -B $(CMAKE_BUILD_DIR)
	cmake --build $(CMAKE_BUILD_DIR) --parallel $(CMAKE_JOBS)

test: build
	ctest --test-dir $(CMAKE_BUILD_DIR) --output-on-failure
	HAPTOOLS_CPP_BIN=$$PWD/$(CMAKE_BUILD_DIR)/haptools_cpp python -m pytest -q tests/python

clean:
	cmake --build $(CMAKE_BUILD_DIR) --target clean 2>/dev/null || true
	rm -rf $(CMAKE_BUILD_DIR)

install: build
	cmake --install $(CMAKE_BUILD_DIR)
