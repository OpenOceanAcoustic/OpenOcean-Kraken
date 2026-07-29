"""Acceptance grading and anomaly attribution for OOKc/KrakenC metrics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ConsistencyThresholds:
    pass_interior_relative_l2: float = 1.0e-3
    pass_tl_p95_db: float = 0.05
    watch_interior_relative_l2: float = 5.0e-3
    watch_tl_p95_db: float = 0.5
    boundary_error_fraction: float = 0.95
    zero_range_energy_fraction: float = 0.5
    zero_range_l2_ratio: float = 5.0
    local_spike_min_db: float = 1.0
    local_spike_to_p95_ratio: float = 10.0


DEFAULT_THRESHOLDS = ConsistencyThresholds()


@dataclass(frozen=True)
class ConsistencyResult:
    grade: str
    anomaly_regions: tuple[str, ...]
    rationale: str


def _number(metric: dict, key: str) -> float:
    try:
        return float(metric[key])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"missing or invalid consistency metric: {key}") from error


def classify_slice(
    metric: dict,
    thresholds: ConsistencyThresholds = DEFAULT_THRESHOLDS,
) -> ConsistencyResult:
    """Classify one source slice using approved acceptance thresholds."""

    interior_l2 = _number(metric, "interior_positive_range_relative_l2")
    tl_p95 = _number(metric, "tl_p95_db")
    if interior_l2 <= thresholds.pass_interior_relative_l2 and tl_p95 <= thresholds.pass_tl_p95_db:
        grade = "pass"
    elif interior_l2 <= thresholds.watch_interior_relative_l2 and tl_p95 <= thresholds.watch_tl_p95_db:
        grade = "watch"
    else:
        grade = "fail"

    full_l2 = _number(metric, "complex_relative_l2")
    positive_l2 = _number(metric, "positive_range_relative_l2")
    zero_energy = _number(metric, "zero_range_reference_energy_fraction")
    boundary_fraction = _number(metric, "boundary_error_fraction")
    tl_max = _number(metric, "tl_max_db")
    anomalies: list[str] = []
    if (
        full_l2 > thresholds.zero_range_l2_ratio * max(positive_l2, 1.0e-30)
        and zero_energy >= thresholds.zero_range_energy_fraction
    ):
        anomalies.append("zero_range")
    if boundary_fraction >= thresholds.boundary_error_fraction:
        anomalies.append("boundary")
    if (
        interior_l2 > thresholds.pass_interior_relative_l2
        or tl_p95 > thresholds.pass_tl_p95_db
    ):
        anomalies.append("interior_positive_range")
    if (
        tl_p95 <= thresholds.pass_tl_p95_db
        and tl_max >= thresholds.local_spike_min_db
        and tl_max >= thresholds.local_spike_to_p95_ratio * max(tl_p95, 1.0e-12)
    ):
        anomalies.append("local_interference_spike")
    if not anomalies:
        anomalies.append("none")

    rationale = (
        f"interior L2={interior_l2:.3e}; TL P95={tl_p95:.6f} dB; "
        f"grade={grade}"
    )
    return ConsistencyResult(grade, tuple(anomalies), rationale)
