"""Build the native extension and install py_kraken in editable mode."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def run(command: list[str], cwd: Path) -> None:
    print("+", subprocess.list2cmdline(command))
    subprocess.run(command, cwd=cwd, check=True)


def main() -> int:
    wrapper_root = Path(__file__).resolve().parent
    repository = wrapper_root.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=repository / "build_python")
    parser.add_argument("--compiler", type=Path)
    parser.add_argument("--config", default="Release")
    args = parser.parse_args()

    configure = [
        "cmake",
        "-S",
        str(repository),
        "-B",
        str(args.build_dir),
        "-DBUILD_AS_PYTHON=ON",
        "-DBUILD_AS_EXE=ON",
        f"-DPython3_EXECUTABLE={sys.executable}",
    ]
    if args.compiler:
        configure.append(f"-DCMAKE_CXX_COMPILER={args.compiler}")
    run(configure, repository)
    run(["cmake", "--build", str(args.build_dir), "--config", args.config, "--parallel"], repository)
    run([sys.executable, "-m", "pip", "install", "-e", str(wrapper_root), "--no-deps"], repository)
    run(
        [
            sys.executable,
            "-c",
            "import py_kraken; from py_kraken import OpenOceanKraken_interface; "
            "assert py_kraken.native is not None; print('py_kraken import OK')",
        ],
        repository,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
