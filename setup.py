from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildPyWithCpp(build_py):
    def run(self) -> None:
        super().run()

        repo_root = Path(__file__).resolve().parent
        build_dir = repo_root / "build-python-package"

        env = os.environ.copy()
        subprocess.check_call(["cmake", "-S", str(repo_root), "-B", str(build_dir)], cwd=repo_root, env=env)
        subprocess.check_call(["cmake", "--build", str(build_dir), "--parallel"], cwd=repo_root, env=env)

        built_bin = build_dir / "haplokit_cpp"
        if not built_bin.exists():
            raise RuntimeError(f"expected built backend at {built_bin}")

        out_dir = Path(self.build_lib) / "haplokit" / "_bin"
        out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built_bin, out_dir / "haplokit_cpp")


setup(
    cmdclass={"build_py": BuildPyWithCpp},
)

