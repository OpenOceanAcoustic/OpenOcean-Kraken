"""Metrics and array helpers for OOKc pressure/particle-velocity fields."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from for_test.shd_reader import ShadeFile


@dataclass(frozen=True)
class DeltaStats:
    valid_points: int
    correlation: float
    mean_abs_delta_db: float
    p50_abs_delta_db: float
    p90_abs_delta_db: float
    p95_abs_delta_db: float
    p98_abs_delta_db: float
    max_abs_delta_db: float


@dataclass(frozen=True)
class SurfaceMetric:
    has_surface_row: bool
    surface_depth_m: float | None
    surface_max_abs: float | None
    global_max_abs: float
    surface_to_global_max_ratio: float | None
    first_interior_to_global_max_ratio: float | None
    is_near_zero: bool
    zero_ratio_limit: float


@dataclass(frozen=True)
class VelocitySliceMetrics:
    interior_point_count: int
    line_receiver_depth_m: float
    pressure_horizontal_interior: DeltaStats
    pressure_horizontal_source_line: DeltaStats
    pressure_vertical_interior: DeltaStats
    pressure_surface: SurfaceMetric
    horizontal_surface: SurfaceMetric
    vertical_surface: SurfaceMetric
    rescaled_horizontal_p50_abs_delta_db: float


def reshape_shade(shade: ShadeFile) -> np.ndarray:
    """Return SHD samples ordered as source, receiver depth, and range."""

    shape = (
        len(shade.source_depths),
        len(shade.receiver_depths),
        len(shade.ranges_metres),
    )
    expected = int(np.prod(shape, dtype=np.int64))
    values = np.asarray(shade.pressure, dtype=np.complex128)
    if values.size != expected:
        raise ValueError(
            f"SHD sample count {values.size} does not match grid size {expected}"
        )
    return values.reshape(shape)


def _matching_grid(label: str, left: list[float], right: list[float]) -> None:
    left_array = np.asarray(left, dtype=float)
    right_array = np.asarray(right, dtype=float)
    if left_array.shape != right_array.shape or not np.allclose(
        left_array, right_array, rtol=0.0, atol=1.0e-8
    ):
        raise ValueError(f"{label} grid mismatch")


def assert_matching_shade_grids(
    pressure: ShadeFile, horizontal: ShadeFile, vertical: ShadeFile
) -> None:
    """Require P/H/V SHD files to describe the same physical grid."""

    for candidate in (horizontal, vertical):
        if not np.isclose(
            pressure.frequency, candidate.frequency, rtol=0.0, atol=1.0e-8
        ):
            raise ValueError("frequency mismatch")
        _matching_grid(
            "source-depth", pressure.source_depths, candidate.source_depths
        )
        _matching_grid(
            "receiver-depth", pressure.receiver_depths, candidate.receiver_depths
        )
        _matching_grid("range", pressure.ranges_metres, candidate.ranges_metres)
    reshape_shade(pressure)
    reshape_shade(horizontal)
    reshape_shade(vertical)


def tl_db(values: np.ndarray, scale: float = 1.0) -> np.ndarray:
    amplitude = np.maximum(np.abs(np.asarray(values)) * scale, 1.0e-30)
    return -20.0 * np.log10(amplitude)


def _delta_stats(left: np.ndarray, right: np.ndarray) -> DeltaStats:
    left = np.asarray(left, dtype=float).ravel()
    right = np.asarray(right, dtype=float).ravel()
    mask = np.isfinite(left) & np.isfinite(right)
    if not np.any(mask):
        raise ValueError("velocity comparison contains no finite overlap")
    left = left[mask]
    right = right[mask]
    delta = np.abs(left - right)
    correlation = float("nan")
    if left.size > 1 and np.std(left) > 0.0 and np.std(right) > 0.0:
        correlation = float(np.corrcoef(left, right)[0, 1])
    return DeltaStats(
        valid_points=int(left.size),
        correlation=correlation,
        mean_abs_delta_db=float(np.mean(delta)),
        p50_abs_delta_db=float(np.percentile(delta, 50.0)),
        p90_abs_delta_db=float(np.percentile(delta, 90.0)),
        p95_abs_delta_db=float(np.percentile(delta, 95.0)),
        p98_abs_delta_db=float(np.percentile(delta, 98.0)),
        max_abs_delta_db=float(np.max(delta)),
    )


def surface_boundary_metrics(
    values: np.ndarray,
    receiver_depths: np.ndarray,
    *,
    zero_ratio_limit: float = 1.0e-8,
) -> SurfaceMetric:
    values = np.asarray(values, dtype=np.complex128)
    receiver_depths = np.asarray(receiver_depths, dtype=float)
    if values.ndim != 2 or values.shape[0] != receiver_depths.size:
        raise ValueError("surface diagnostic grid mismatch")
    has_surface = bool(
        receiver_depths.size
        and np.isclose(receiver_depths[0], 0.0, rtol=0.0, atol=1.0e-8)
    )
    global_max = float(np.max(np.abs(values))) if values.size else 0.0
    surface_max = float(np.max(np.abs(values[0]))) if has_surface else None
    first_interior = (
        float(np.max(np.abs(values[1])))
        if has_surface and values.shape[0] > 1
        else None
    )
    surface_ratio = (
        surface_max / global_max
        if has_surface and global_max > 0.0
        else (0.0 if has_surface and surface_max == 0.0 else None)
    )
    interior_ratio = (
        first_interior / global_max
        if first_interior is not None and global_max > 0.0
        else None
    )
    return SurfaceMetric(
        has_surface_row=has_surface,
        surface_depth_m=float(receiver_depths[0]) if has_surface else None,
        surface_max_abs=surface_max,
        global_max_abs=global_max,
        surface_to_global_max_ratio=surface_ratio,
        first_interior_to_global_max_ratio=interior_ratio,
        is_near_zero=bool(
            has_surface
            and surface_ratio is not None
            and surface_ratio <= zero_ratio_limit
        ),
        zero_ratio_limit=zero_ratio_limit,
    )


def compute_velocity_slice_metrics(
    pressure: np.ndarray,
    horizontal: np.ndarray,
    vertical: np.ndarray,
    receiver_depths: np.ndarray,
    ranges_metres: np.ndarray,
    source_depth: float,
    *,
    rho_c0: float = 1500.0,
) -> VelocitySliceMetrics:
    """Compare one source slice using the validated equivalent-pressure scale."""

    pressure = np.asarray(pressure, dtype=np.complex128)
    horizontal = np.asarray(horizontal, dtype=np.complex128)
    vertical = np.asarray(vertical, dtype=np.complex128)
    receiver_depths = np.asarray(receiver_depths, dtype=float)
    ranges_metres = np.asarray(ranges_metres, dtype=float)
    if pressure.shape != horizontal.shape or pressure.shape != vertical.shape:
        raise ValueError("P/H/V field shape mismatch")
    if pressure.ndim != 2 or pressure.shape != (
        receiver_depths.size,
        ranges_metres.size,
    ):
        raise ValueError("P/H/V spatial grid mismatch")
    if not np.isfinite(rho_c0) or rho_c0 <= 0.0:
        raise ValueError("rho*c0 must be positive and finite")

    positive_range = ranges_metres > 0.0
    depth_mask = np.ones(receiver_depths.size, dtype=bool)
    if receiver_depths.size and np.isclose(
        receiver_depths[0], 0.0, rtol=0.0, atol=1.0e-8
    ):
        depth_mask[0] = False
    spatial_mask = depth_mask[:, np.newaxis] & positive_range[np.newaxis, :]
    if not np.any(spatial_mask):
        raise ValueError("velocity comparison has no interior positive-range points")

    pressure_tl = tl_db(pressure)
    horizontal_tl = tl_db(horizontal)
    vertical_tl = tl_db(vertical)
    line_index = int(np.argmin(np.abs(receiver_depths - source_depth)))
    line_mask = positive_range
    if not np.any(line_mask):
        raise ValueError("velocity source-depth line has no positive ranges")

    rescaled_horizontal_tl = tl_db(horizontal, scale=rho_c0)
    return VelocitySliceMetrics(
        interior_point_count=int(np.count_nonzero(spatial_mask)),
        line_receiver_depth_m=float(receiver_depths[line_index]),
        pressure_horizontal_interior=_delta_stats(
            pressure_tl[spatial_mask], horizontal_tl[spatial_mask]
        ),
        pressure_horizontal_source_line=_delta_stats(
            pressure_tl[line_index, line_mask],
            horizontal_tl[line_index, line_mask],
        ),
        pressure_vertical_interior=_delta_stats(
            pressure_tl[spatial_mask], vertical_tl[spatial_mask]
        ),
        pressure_surface=surface_boundary_metrics(pressure, receiver_depths),
        horizontal_surface=surface_boundary_metrics(horizontal, receiver_depths),
        vertical_surface=surface_boundary_metrics(vertical, receiver_depths),
        rescaled_horizontal_p50_abs_delta_db=_delta_stats(
            pressure_tl[spatial_mask], rescaled_horizontal_tl[spatial_mask]
        ).p50_abs_delta_db,
    )
