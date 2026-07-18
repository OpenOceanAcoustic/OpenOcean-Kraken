from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
import numpy as np

from for_test.shd_reader import ShadeFile


matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


PRESSURE_FLOOR = 1.0e-12
VALID_RELATIVE_AMPLITUDE = 1.0e-8


@dataclass(frozen=True)
class PressureCube:
    title: str
    frequency: float
    source_depths: np.ndarray
    receiver_depths: np.ndarray
    ranges_metres: np.ndarray
    pressure: np.ndarray


@dataclass(frozen=True)
class ComparisonMetrics:
    complex_relative_l2: float
    tl_mae_db: float
    tl_rmse_db: float
    tl_p95_db: float
    tl_max_db: float
    valid_fraction: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class SpatialDiagnostics:
    positive_range_relative_l2: float
    interior_positive_range_relative_l2: float
    zero_range_reference_energy_fraction: float
    boundary_error_fraction: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


def shade_to_cube(shade: ShadeFile) -> PressureCube:
    shape = (
        len(shade.source_depths),
        len(shade.receiver_depths),
        len(shade.ranges_metres),
    )
    pressure = np.asarray(shade.pressure, dtype=np.complex128)
    if pressure.size != int(np.prod(shape)):
        raise ValueError(
            f"pressure size {pressure.size} does not match grid {shape}"
        )
    return PressureCube(
        title=shade.title,
        frequency=float(shade.frequency),
        source_depths=np.asarray(shade.source_depths, dtype=float),
        receiver_depths=np.asarray(shade.receiver_depths, dtype=float),
        ranges_metres=np.asarray(shade.ranges_metres, dtype=float),
        pressure=pressure.reshape(shape),
    )


def compute_tl(pressure: np.ndarray) -> np.ndarray:
    return -20.0 * np.log10(
        np.maximum(np.abs(np.asarray(pressure)), PRESSURE_FLOOR)
    )


def valid_tl_mask(ookc: np.ndarray, krakenc: np.ndarray) -> np.ndarray:
    joint_amplitude = np.maximum(np.abs(ookc), np.abs(krakenc))
    peak = float(np.max(joint_amplitude, initial=0.0))
    if peak == 0.0:
        return np.ones(joint_amplitude.shape, dtype=bool)
    return joint_amplitude >= peak * VALID_RELATIVE_AMPLITUDE


def compute_metrics(
    ookc: np.ndarray, krakenc: np.ndarray
) -> ComparisonMetrics:
    ookc = np.asarray(ookc, dtype=np.complex128)
    krakenc = np.asarray(krakenc, dtype=np.complex128)
    if ookc.shape != krakenc.shape:
        raise ValueError("pressure shape mismatch")
    difference = ookc - krakenc
    denominator = float(np.linalg.norm(krakenc.ravel()))
    relative_l2 = (
        float(np.linalg.norm(difference.ravel()) / denominator)
        if denominator > 0.0
        else float("inf")
    )
    mask = valid_tl_mask(ookc, krakenc)
    db_error = np.abs(compute_tl(ookc)[mask] - compute_tl(krakenc)[mask])
    return ComparisonMetrics(
        complex_relative_l2=relative_l2,
        tl_mae_db=float(np.mean(db_error)),
        tl_rmse_db=float(np.sqrt(np.mean(np.square(db_error)))),
        tl_p95_db=float(np.percentile(db_error, 95.0)),
        tl_max_db=float(np.max(db_error)),
        valid_fraction=float(np.mean(mask)),
    )


def _relative_l2(ookc: np.ndarray, krakenc: np.ndarray) -> float:
    denominator = float(np.linalg.norm(krakenc.ravel()))
    if denominator == 0.0:
        return float("inf")
    return float(np.linalg.norm((ookc - krakenc).ravel()) / denominator)


def spatial_diagnostics(
    ookc: np.ndarray,
    krakenc: np.ndarray,
    ranges_metres: np.ndarray,
) -> SpatialDiagnostics:
    ookc = np.asarray(ookc, dtype=np.complex128)
    krakenc = np.asarray(krakenc, dtype=np.complex128)
    ranges_metres = np.asarray(ranges_metres, dtype=float)
    if ookc.shape != krakenc.shape or ookc.ndim != 2:
        raise ValueError("spatial diagnostics require matching 2D fields")
    if ookc.shape[1] != ranges_metres.size:
        raise ValueError("range grid size does not match pressure field")
    positive_range = ranges_metres > 0.0
    if not np.any(positive_range):
        raise ValueError("spatial diagnostics require a positive range")
    positive_l2 = _relative_l2(
        ookc[:, positive_range], krakenc[:, positive_range]
    )
    if ookc.shape[0] > 2:
        interior_ookc = ookc[1:-1, positive_range]
        interior_krakenc = krakenc[1:-1, positive_range]
    else:
        interior_ookc = ookc[:, positive_range]
        interior_krakenc = krakenc[:, positive_range]
    interior_l2 = _relative_l2(interior_ookc, interior_krakenc)
    zero_range = np.isclose(ranges_metres, 0.0, rtol=0.0, atol=1.0e-9)
    reference_norm = float(np.linalg.norm(krakenc.ravel()))
    zero_range_fraction = (
        float(
            np.square(np.linalg.norm(krakenc[:, zero_range].ravel()))
            / np.square(reference_norm)
        )
        if reference_norm > 0.0 and np.any(zero_range)
        else 0.0
    )
    error = ookc - krakenc
    error_norm = float(np.linalg.norm(error.ravel()))
    if error_norm == 0.0:
        boundary_fraction = 0.0
    elif error.shape[0] == 1:
        boundary_fraction = 1.0
    else:
        boundary_error = np.concatenate([error[0].ravel(), error[-1].ravel()])
        boundary_fraction = float(
            np.square(np.linalg.norm(boundary_error))
            / np.square(error_norm)
        )
    return SpatialDiagnostics(
        positive_range_relative_l2=positive_l2,
        interior_positive_range_relative_l2=interior_l2,
        zero_range_reference_energy_fraction=zero_range_fraction,
        boundary_error_fraction=boundary_fraction,
    )


def assert_matching_grids(left: PressureCube, right: PressureCube) -> None:
    if not np.isclose(left.frequency, right.frequency, rtol=0.0, atol=1.0e-8):
        raise ValueError("frequency mismatch")
    for label, left_grid, right_grid in (
        ("source-depth", left.source_depths, right.source_depths),
        ("receiver-depth", left.receiver_depths, right.receiver_depths),
        ("range", left.ranges_metres, right.ranges_metres),
    ):
        if left_grid.shape != right_grid.shape or not np.allclose(
            left_grid, right_grid, rtol=0.0, atol=1.0e-8
        ):
            raise ValueError(f"{label} grid mismatch")
    if left.pressure.shape != right.pressure.shape:
        raise ValueError("pressure grid mismatch")


def comparison_color_limits(
    ookc: np.ndarray, krakenc: np.ndarray
) -> tuple[float, float, float]:
    ookc_tl = compute_tl(ookc)
    krakenc_tl = compute_tl(krakenc)
    combined = np.concatenate([ookc_tl.ravel(), krakenc_tl.ravel()])
    finite = combined[np.isfinite(combined)]
    if finite.size == 0:
        raise ValueError("transmission-loss arrays contain no finite values")
    low, high = np.percentile(finite, [2.0, 98.0])
    if not high > low:
        high = low + 1.0
    absolute_difference = np.abs(ookc_tl - krakenc_tl)
    difference_half_width = max(
        0.1, float(np.percentile(absolute_difference, 98.0))
    )
    return float(low), float(high), difference_half_width


def nearest_receiver_index(
    receiver_depths: np.ndarray, source_depth: float
) -> int:
    receiver_depths = np.asarray(receiver_depths, dtype=float)
    if receiver_depths.size == 0:
        raise ValueError("receiver-depth grid is empty")
    return int(np.argmin(np.abs(receiver_depths - source_depth)))


def plot_comparison(
    case_name: str,
    ookc: PressureCube,
    krakenc: PressureCube,
    source_index: int,
    output_path: Path,
) -> ComparisonMetrics:
    assert_matching_grids(ookc, krakenc)
    if not 0 <= source_index < len(ookc.source_depths):
        raise IndexError(f"source index out of range: {source_index}")
    ookc_pressure = ookc.pressure[source_index]
    krakenc_pressure = krakenc.pressure[source_index]
    ookc_tl = compute_tl(ookc_pressure)
    krakenc_tl = compute_tl(krakenc_pressure)
    difference = ookc_tl - krakenc_tl
    metrics = compute_metrics(ookc_pressure, krakenc_pressure)
    diagnostics = spatial_diagnostics(
        ookc_pressure, krakenc_pressure, ookc.ranges_metres
    )
    tl_min, tl_max, difference_half_width = comparison_color_limits(
        ookc_pressure, krakenc_pressure
    )
    ranges_km = ookc.ranges_metres / 1000.0
    depths = ookc.receiver_depths
    source_depth = float(ookc.source_depths[source_index])
    figure, axes = plt.subplots(
        2, 2, figsize=(15, 10), constrained_layout=True
    )
    for axis, values, title in (
        (axes[0, 0], ookc_tl, "OOKc transmission loss"),
        (axes[0, 1], krakenc_tl, "KrakenC transmission loss"),
    ):
        mesh = axis.pcolormesh(
            ranges_km,
            depths,
            values,
            shading="auto",
            cmap="viridis",
            vmin=tl_min,
            vmax=tl_max,
        )
        axis.invert_yaxis()
        axis.set(
            xlabel="Range (km)",
            ylabel="Receiver depth (m)",
            title=title,
        )
        figure.colorbar(mesh, ax=axis, label="TL (dB)")
    difference_norm = TwoSlopeNorm(
        vmin=-difference_half_width,
        vcenter=0.0,
        vmax=difference_half_width,
    )
    difference_mesh = axes[1, 0].pcolormesh(
        ranges_km,
        depths,
        difference,
        shading="auto",
        cmap="RdBu_r",
        norm=difference_norm,
    )
    axes[1, 0].invert_yaxis()
    axes[1, 0].set(
        xlabel="Range (km)",
        ylabel="Receiver depth (m)",
        title="OOKc - KrakenC TL difference",
    )
    figure.colorbar(difference_mesh, ax=axes[1, 0], label="Difference (dB)")
    receiver_index = nearest_receiver_index(depths, source_depth)
    axes[1, 1].plot(
        ranges_km,
        ookc_tl[receiver_index],
        label="OOKc",
        linewidth=1.2,
    )
    axes[1, 1].plot(
        ranges_km,
        krakenc_tl[receiver_index],
        label="KrakenC",
        linewidth=1.0,
        linestyle="--",
    )
    axes[1, 1].invert_yaxis()
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()
    axes[1, 1].set(
        xlabel="Range (km)",
        ylabel="TL (dB)",
        title=(
            "Range cut at receiver depth "
            f"{depths[receiver_index]:.1f} m"
        ),
    )
    axes[1, 1].text(
        0.02,
        0.02,
        f"relative L2 = {metrics.complex_relative_l2:.3e}\n"
        f"interior r>0 L2 = "
        f"{diagnostics.interior_positive_range_relative_l2:.3e}\n"
        f"TL RMSE = {metrics.tl_rmse_db:.3f} dB\n"
        f"TL P95 = {metrics.tl_p95_db:.3f} dB",
        transform=axes[1, 1].transAxes,
        va="bottom",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"},
    )
    figure.suptitle(
        f"{case_name} | {ookc.frequency:g} Hz | "
        f"source depth {source_depth:g} m"
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return metrics


def create_overview(
    figure_paths: list[Path], output_path: Path
) -> None:
    if not figure_paths:
        raise ValueError("at least one figure is required")
    columns = min(3, len(figure_paths))
    rows = (len(figure_paths) + columns - 1) // columns
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(6 * columns, 4 * rows),
        squeeze=False,
        constrained_layout=True,
    )
    for axis, path in zip(axes.ravel(), figure_paths):
        axis.imshow(plt.imread(path))
        axis.set_title(Path(path).stem, fontsize=9)
        axis.axis("off")
    for axis in axes.ravel()[len(figure_paths) :]:
        axis.axis("off")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)


def failure_figures(
    case_name: str,
    source_count: int,
    message: str,
    figures_directory: Path,
) -> list[Path]:
    if source_count < 1:
        raise ValueError("source count must be positive")
    figures_directory = Path(figures_directory)
    figures_directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for source_index in range(source_count):
        suffix = f"source{source_index + 1}"
        path = figures_directory / f"{case_name}_{suffix}_failed.png"
        figure, axis = plt.subplots(
            figsize=(12, 8), constrained_layout=True
        )
        axis.axis("off")
        axis.text(
            0.5,
            0.55,
            f"{case_name}\ncomparison unavailable",
            ha="center",
            va="center",
            fontsize=22,
        )
        axis.text(
            0.5,
            0.42,
            message,
            ha="center",
            va="center",
            fontsize=13,
            wrap=True,
        )
        figure.savefig(
            path, dpi=200, bbox_inches="tight", facecolor="white"
        )
        plt.close(figure)
        paths.append(path)
    return paths
