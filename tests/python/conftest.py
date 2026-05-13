from __future__ import annotations

import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
BUILD_DIR = ROOT / "build-wsl"
CPP_BIN = BUILD_DIR / "haptools_cpp"


@pytest.fixture(scope="session", autouse=True)
def ensure_cpp_backend() -> None:
    os.environ["HAPTOOLS_CPP_BIN"] = str(CPP_BIN)
