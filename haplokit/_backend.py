from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


CPP_BACKEND = "haplokit_cpp"
NETWORK_BACKEND = "haplokit_network_backend"
BACKEND_BINARIES = (CPP_BACKEND, NETWORK_BACKEND)
PYTHON_BUILD_DIR = "build-haplokit-python"


class CppBackendBuildError(RuntimeError):
    """Raised when the C++ backend cannot be built from the source tree."""


def executable_candidates(directory: Path, name: str) -> list[Path]:
    candidates = [directory / name]
    if not name.endswith(".exe"):
        candidates.append(directory / f"{name}.exe")
    return candidates


def existing_executable(directory: Path, name: str) -> Path | None:
    for candidate in executable_candidates(directory, name):
        if candidate.exists():
            return candidate
    return None


def backend_search_paths(repo_root: Path, package_dir: Path, env_path: str | None) -> list[Path]:
    candidates: list[Path] = []
    if env_path:
        candidates.append(Path(env_path))
    candidates.extend(
        [
            package_dir / "_bin" / CPP_BACKEND,
            package_dir / "_bin" / f"{CPP_BACKEND}.exe",
            repo_root / "build-wsl" / CPP_BACKEND,
            repo_root / "build-wsl" / f"{CPP_BACKEND}.exe",
            repo_root / "build" / CPP_BACKEND,
            repo_root / "build" / f"{CPP_BACKEND}.exe",
            repo_root / "build-python-package" / CPP_BACKEND,
            repo_root / "build-python-package" / f"{CPP_BACKEND}.exe",
            repo_root / PYTHON_BUILD_DIR / CPP_BACKEND,
            repo_root / PYTHON_BUILD_DIR / f"{CPP_BACKEND}.exe",
        ]
    )
    return candidates


def _cmake_cache_home_directory(build_dir: Path) -> str | None:
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return None
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("CMAKE_HOME_DIRECTORY:INTERNAL="):
            return line.split("=", 1)[1]
    return None


def compatible_build_dir(repo_root: Path, preferred_build_dir: Path) -> Path:
    cache_home = _cmake_cache_home_directory(preferred_build_dir)
    if cache_home is None:
        return preferred_build_dir

    current_home = str(repo_root)
    if cache_home.replace("\\", "/") == current_home.replace("\\", "/"):
        return preferred_build_dir

    return repo_root / PYTHON_BUILD_DIR


def _format_failed_command(cmd: list[str], completed: subprocess.CompletedProcess[str]) -> str:
    output = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)
    if len(output) > 4000:
        output = output[-4000:]
    detail = output or f"command exited with code {completed.returncode}"
    return f"{' '.join(cmd)} failed:\n{detail}"


def _run_checked(cmd: list[str], repo_root: Path) -> None:
    try:
        completed = subprocess.run(
            cmd,
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise CppBackendBuildError(f"{' '.join(cmd)} failed: {exc}") from exc
    if completed.returncode != 0:
        raise CppBackendBuildError(_format_failed_command(cmd, completed))


def build_cpp_backends(repo_root: Path, build_dir: Path | None = None) -> dict[str, Path]:
    repo_root = repo_root.resolve()
    build_dir = compatible_build_dir(repo_root, (build_dir or repo_root / "build-wsl").resolve())
    if not (repo_root / "CMakeLists.txt").exists():
        raise CppBackendBuildError(
            f"cannot auto-build {CPP_BACKEND}: CMakeLists.txt not found at {repo_root}"
        )

    _run_checked(
        ["cmake", "-S", str(repo_root), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"],
        repo_root,
    )
    _run_checked(
        ["cmake", "--build", str(build_dir), "--parallel", str(os.cpu_count() or 1)],
        repo_root,
    )

    built = {
        name: path
        for name in BACKEND_BINARIES
        if (path := existing_executable(build_dir, name)) is not None
    }
    if CPP_BACKEND not in built:
        raise CppBackendBuildError(f"expected built backend at {build_dir / CPP_BACKEND}")
    return built


def copy_backend_binaries(repo_root: Path, build_dir: Path, out_dir: Path) -> None:
    built = build_cpp_backends(repo_root, build_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, path in built.items():
        shutil.copy2(path, out_dir / path.name)


def find_haplokit_cpp(repo_root: Path, package_dir: Path, auto_build: bool = True) -> Path:
    env_path = os.environ.get("HAPLOKIT_CPP_BIN")
    candidates = backend_search_paths(repo_root, package_dir, env_path)
    for candidate in candidates:
        if candidate.exists():
            return candidate

    if auto_build and (repo_root / "CMakeLists.txt").exists():
        build_cpp_backends(repo_root, repo_root / "build-wsl")
        for candidate in backend_search_paths(repo_root, package_dir, env_path):
            if candidate.exists():
                return candidate

    checked = "\n".join(f"  - {candidate}" for candidate in candidates)
    raise FileNotFoundError(
        f"{CPP_BACKEND} backend not found.\n"
        "Run `pip install .` from the repository root, or set HAPLOKIT_CPP_BIN to the compiled backend.\n"
        f"Checked:\n{checked}"
    )
