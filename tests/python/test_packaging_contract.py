from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import haptools


def test_pyproject_declares_haptools_console_entrypoint_and_linux_scope() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'name = "haptools"' in pyproject
    assert 'haptools = "haptools.cli:main"' in pyproject
    assert "Operating System :: POSIX :: Linux" in pyproject


def test_setup_py_defines_cpp_build_hook() -> None:
    setup_py = (ROOT / "setup.py").read_text(encoding="utf-8")
    assert "class BuildPyWithCpp" in setup_py
    assert "cmdclass={\"build_py\": BuildPyWithCpp}" in setup_py


def test_project_version_matches_package_version() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, flags=re.MULTILINE)
    assert match is not None
    assert match.group(1) == haptools.__version__


def test_module_entrypoint_invokes_cli() -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "haptools"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode != 0
    stderr = completed.stderr
    assert "a subcommand is required" in stderr or "usage:" in stderr


def test_cli_checks_packaged_backend_path() -> None:
    cli_py = (ROOT / "haptools" / "cli.py").read_text(encoding="utf-8")
    assert 'Path(__file__).resolve().parent / "_bin" / "haptools_cpp"' in cli_py


def test_vendored_htscodecs_version_header_is_present_for_sdist_builds() -> None:
    configure_ac = (ROOT / "deps" / "htslib" / "htscodecs" / "configure.ac").read_text(encoding="utf-8")
    match = re.search(r"AC_INIT\(htscodecs,\s*([0-9]+\.[0-9]+\.[0-9]+)\)", configure_ac)
    assert match is not None
    expected = match.group(1)

    version_h = (ROOT / "deps" / "htslib" / "htscodecs" / "htscodecs" / "version.h").read_text(
        encoding="utf-8"
    )
    assert f'#define HTSCODECS_VERSION_TEXT "{expected}"' in version_h
