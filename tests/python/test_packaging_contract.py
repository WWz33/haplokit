from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import haplokit
import haplokit.cli as cli_module
from haplokit._backend import (
    CppBackendBuildError,
    _native_build_environment,
    compatible_build_dir,
    find_haplokit_cpp,
    native_runtime_environment,
)


def test_pyproject_declares_haplokit_console_entrypoint_and_linux_scope() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'name = "haplokit"' in pyproject
    assert 'haplokit = "haplokit.cli:main"' in pyproject
    assert "Operating System :: POSIX :: Linux" in pyproject
    assert '"cmake>=3.22"' in pyproject


def test_setup_py_defines_cpp_build_hook() -> None:
    setup_py = (ROOT / "setup.py").read_text(encoding="utf-8")
    assert "class BinaryDistribution" in setup_py
    assert "def has_ext_modules" in setup_py
    assert "class BuildPyWithCpp" in setup_py
    assert 'cmdclass = {"build_py": BuildPyWithCpp}' in setup_py
    assert "class EditableWheelWithCpp" in setup_py
    assert 'cmdclass["editable_wheel"] = EditableWheelWithCpp' in setup_py
    assert "distclass=BinaryDistribution" in setup_py


def test_project_version_matches_package_version() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, flags=re.MULTILINE)
    assert match is not None
    assert match.group(1) == haplokit.__version__


def test_module_entrypoint_invokes_cli() -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "haplokit"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode != 0
    stderr = completed.stderr
    assert "a subcommand is required" in stderr or "usage:" in stderr


def test_cli_checks_packaged_backend_path() -> None:
    cli_py = (ROOT / "haplokit" / "cli.py").read_text(encoding="utf-8")
    assert "find_haplokit_cpp" in cli_py


def test_network_backend_checks_python_build_dir() -> None:
    network_py = (ROOT / "haplokit" / "network.py").read_text(encoding="utf-8")
    assert "build-python-package" in network_py
    assert "PYTHON_BUILD_DIR" in network_py
    assert "native_runtime_environment" in network_py


def test_find_haplokit_cpp_auto_builds_source_tree(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    package_dir = repo_root / "haplokit"
    build_dir = repo_root / "build-wsl"
    package_dir.mkdir(parents=True)
    (repo_root / "CMakeLists.txt").write_text("cmake_minimum_required(VERSION 3.22)\n", encoding="utf-8")
    monkeypatch.delenv("HAPLOKIT_CPP_BIN", raising=False)

    calls: list[list[str]] = []

    def fake_run(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(cmd)
        if cmd[:2] == ["cmake", "--build"]:
            build_dir.mkdir()
            (build_dir / "haplokit_cpp").write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr("haplokit._backend.subprocess.run", fake_run)

    assert find_haplokit_cpp(repo_root, package_dir) == build_dir / "haplokit_cpp"
    assert "-S" in calls[0]
    assert str(repo_root.resolve()) in calls[0]
    assert "--build" in calls[1]


def test_find_haplokit_cpp_reports_cmake_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    package_dir = repo_root / "haplokit"
    package_dir.mkdir(parents=True)
    (repo_root / "CMakeLists.txt").write_text("cmake_minimum_required(VERSION 3.22)\n", encoding="utf-8")
    monkeypatch.delenv("HAPLOKIT_CPP_BIN", raising=False)

    def fake_run(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(cmd, 1, "", "missing zlib")

    monkeypatch.setattr("haplokit._backend.subprocess.run", fake_run)

    with pytest.raises(CppBackendBuildError, match="missing zlib"):
        find_haplokit_cpp(repo_root, package_dir)


def test_native_build_environment_exposes_conda_native_paths(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    conda_prefix = tmp_path / "conda-env"
    (conda_prefix / "include").mkdir(parents=True)
    (conda_prefix / "lib" / "pkgconfig").mkdir(parents=True)
    (conda_prefix / "share" / "pkgconfig").mkdir(parents=True)
    monkeypatch.setenv("CONDA_PREFIX", str(conda_prefix))
    monkeypatch.setenv("LIBRARY_PATH", "/existing/lib")

    env = _native_build_environment()

    assert env["CPATH"].split(os.pathsep)[0] == str(conda_prefix / "include")
    assert env["LIBRARY_PATH"].split(os.pathsep)[:2] == [str(conda_prefix / "lib"), "/existing/lib"]
    assert env["CMAKE_PREFIX_PATH"].split(os.pathsep)[0] == str(conda_prefix)
    assert env["PKG_CONFIG_PATH"].split(os.pathsep)[:2] == [
        str(conda_prefix / "lib" / "pkgconfig"),
        str(conda_prefix / "share" / "pkgconfig"),
    ]


def test_native_runtime_environment_exposes_conda_library_path(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    conda_prefix = tmp_path / "conda-env"
    (conda_prefix / "lib").mkdir(parents=True)
    (conda_prefix / "bin").mkdir()
    (conda_prefix / "Library" / "bin").mkdir(parents=True)
    monkeypatch.setenv("CONDA_PREFIX", str(conda_prefix))
    monkeypatch.setenv("LD_LIBRARY_PATH", "/existing/runtime")

    env = native_runtime_environment()

    if os.name == "nt":
        assert str(conda_prefix / "bin") in env["PATH"].split(os.pathsep)
    else:
        assert env["LD_LIBRARY_PATH"].split(os.pathsep)[:2] == [
            str(conda_prefix / "lib"),
            "/existing/runtime",
        ]


def test_cli_backend_subprocess_uses_native_runtime_environment(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    captured: dict[str, dict[str, str]] = {}

    def fake_run(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        captured["env"] = kwargs["env"]  # type: ignore[assignment]
        return subprocess.CompletedProcess(cmd, 0, "scaffold_1:1-2", "")

    conda_prefix = tmp_path / "conda-env"
    (conda_prefix / "lib").mkdir(parents=True)
    (conda_prefix / "bin").mkdir()
    (conda_prefix / "Library" / "bin").mkdir(parents=True)
    monkeypatch.setenv("CONDA_PREFIX", str(conda_prefix))
    monkeypatch.setattr(cli_module, "_cpp_backend_path", lambda: Path("/tmp/haplokit_cpp"))
    monkeypatch.setattr(cli_module.subprocess, "run", fake_run)

    assert cli_module._resolve_gene_region("anno.gff", "gene1") == "scaffold_1:1-2"
    env = captured["env"]
    if os.name == "nt":
        assert str(conda_prefix / "bin") in env["PATH"].split(os.pathsep)
    else:
        assert env["LD_LIBRARY_PATH"].split(os.pathsep)[0] == str(conda_prefix / "lib")


def test_incompatible_cmake_cache_uses_python_build_dir(tmp_path: Path) -> None:
    repo_root = (tmp_path / "repo").resolve()
    build_dir = repo_root / "build-wsl"
    build_dir.mkdir(parents=True)
    (build_dir / "CMakeCache.txt").write_text(
        "CMAKE_HOME_DIRECTORY:INTERNAL=/old/source/tree\n",
        encoding="utf-8",
    )

    assert compatible_build_dir(repo_root, build_dir) == repo_root / "build-haplokit-python"


def test_vendored_htscodecs_version_header_is_present_for_sdist_builds() -> None:
    configure_ac = (ROOT / "deps" / "htslib" / "htscodecs" / "configure.ac").read_text(encoding="utf-8")
    match = re.search(r"AC_INIT\(htscodecs,\s*([0-9]+\.[0-9]+\.[0-9]+)\)", configure_ac)
    assert match is not None
    expected = match.group(1)

    version_h = (ROOT / "deps" / "htslib" / "htscodecs" / "htscodecs" / "version.h").read_text(
        encoding="utf-8"
    )
    assert f'#define HTSCODECS_VERSION_TEXT "{expected}"' in version_h


def test_sdist_includes_network_json_header() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "recursive-include src/cpp/third_party *.hpp" in manifest
    assert (ROOT / "src" / "cpp" / "third_party" / "json.hpp").exists()
