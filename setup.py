from __future__ import annotations

import importlib.util
from pathlib import Path

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py

try:
    from setuptools.command.editable_wheel import editable_wheel
except Exception:  # pragma: no cover - setuptools>=69 provides this for PEP 660
    editable_wheel = None


def _load_backend_module():
    backend_path = Path(__file__).resolve().parent / "haplokit" / "_backend.py"
    spec = importlib.util.spec_from_file_location("haplokit_build_backend", backend_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load backend build helper from {backend_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_backend = _load_backend_module()


class BinaryDistribution(Distribution):
    def has_ext_modules(self) -> bool:
        return True


class BuildPyWithCpp(build_py):
    def run(self) -> None:
        super().run()

        repo_root = Path(__file__).resolve().parent
        build_dir = repo_root / "build-python-package"
        out_dir = Path(self.build_lib) / "haplokit" / "_bin"

        _backend.copy_backend_binaries(repo_root, build_dir, out_dir)


cmdclass = {"build_py": BuildPyWithCpp}


if editable_wheel is not None:

    class EditableWheelWithCpp(editable_wheel):
        def run(self) -> None:
            repo_root = Path(__file__).resolve().parent
            _backend.build_cpp_backends(repo_root, repo_root / "build-wsl")
            super().run()


    cmdclass["editable_wheel"] = EditableWheelWithCpp


setup(cmdclass=cmdclass, distclass=BinaryDistribution)
