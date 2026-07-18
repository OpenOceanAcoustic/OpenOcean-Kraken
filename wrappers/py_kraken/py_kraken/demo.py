"""Headless end-to-end py_kraken example."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from .ook_interface import OpenOceanKraken_interface
from .ook_plot import plot_field, plot_modes


def main() -> int:
    repository = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture",
        type=Path,
        default=(repository / "../test/MunkK.env").resolve(),
    )
    parser.add_argument("--output-dir", type=Path, default=repository / "tmp/py_kraken_demo")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    interface = OpenOceanKraken_interface(thread_num=1)
    interface.ook_load_env(args.fixture)
    interface.ook_run()
    pressure = interface.ook_get_pressure()
    modes = interface.ook_get_modes()

    field_axes = plot_field(pressure)
    field_axes.figure.savefig(args.output_dir / "pressure.png", dpi=150, bbox_inches="tight")
    mode_figure = plot_modes(modes)
    mode_figure.savefig(args.output_dir / "modes.png", dpi=150, bbox_inches="tight")
    print(f"pressure shape: {pressure.values.shape}")
    print(f"mode profiles: {len(modes.profiles)}")
    print(f"plots written to: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
