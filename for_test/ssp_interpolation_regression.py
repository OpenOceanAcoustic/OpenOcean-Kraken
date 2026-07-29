"""Check modal-loss and transmission-loss non-regression for ALG-022."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import statistics
import struct
import tempfile
from typing import Any, Sequence

import numpy as np

from .mod_reader import ModeFileError, read_first_wavenumber_set
from .shd_reader import read_shade_file


GRID_ATOL_METRES = 1.0e-8
FREQUENCY_ATOL_HZ = 1.0e-8
E_IM_REGRESSION_ATOL = 1.0e-12
TL_REGRESSION_ATOL_DB = 1.0e-9
KNOWN_BEFORE_METRICS = {
    "e_im": {
        "median": 5.289488948137276e-4,
        "p95": 7.564144223952983e-4,
    },
    "tl_db": {
        "median": 8.269118904209449e-4,
        "p95": 3.7592181663796964e-3,
    },
}


class RegressionError(RuntimeError):
    """A deterministic validation failure with a machine-readable status."""

    def __init__(self, status: str, detail: str) -> None:
        self.status = status
        self.detail = detail
        super().__init__(f"{status}: {detail}")


@dataclass(frozen=True)
class AcousticFields:
    wavenumbers: np.ndarray
    frequency_hz: float
    source_depths_metres: np.ndarray
    receiver_depths_metres: np.ndarray
    ranges_metres: np.ndarray
    pressure: np.ndarray


@dataclass(frozen=True)
class RunPair:
    reference: AcousticFields
    actual: AcousticFields


@dataclass(frozen=True)
class ArtifactPaths:
    case_stem: str
    reference_mod: Path
    reference_shd: Path
    actual_mod: Path
    actual_shd: Path


def relative_imaginary_error(actual: complex, reference: complex) -> float:
    return abs(actual.imag - reference.imag) / max(
        abs(reference.imag), 1.0e-12
    )


def nearest_rank(values: list[float], quantile: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("EMPTY_METRIC_DOMAIN")
    index = math.ceil(quantile * len(ordered)) - 1
    return ordered[index]


def transmission_loss(pressure: np.ndarray) -> np.ndarray:
    return -20.0 * np.log10(np.maximum(np.abs(pressure), 1.0e-30))


def metric_summary(values: Sequence[float]) -> dict[str, float]:
    finite_values = [float(value) for value in values]
    if not finite_values:
        raise ValueError("EMPTY_METRIC_DOMAIN")
    if not all(math.isfinite(value) for value in finite_values):
        raise ValueError("NONFINITE_VALUE")
    return {
        "median": float(statistics.median(finite_values)),
        "p95": nearest_rank(finite_values, 0.95),
    }


def _all_finite(fields: AcousticFields) -> bool:
    return math.isfinite(fields.frequency_hz) and all(
        bool(np.all(np.isfinite(values)))
        for values in (
            fields.wavenumbers,
            fields.source_depths_metres,
            fields.receiver_depths_metres,
            fields.ranges_metres,
            fields.pressure,
        )
    )


def _internal_shape_is_valid(fields: AcousticFields) -> bool:
    if (
        fields.wavenumbers.ndim != 1
        or fields.source_depths_metres.ndim != 1
        or fields.receiver_depths_metres.ndim != 1
        or fields.ranges_metres.ndim != 1
        or fields.pressure.ndim != 3
    ):
        return False
    expected_pressure_shape = (
        fields.source_depths_metres.size,
        fields.receiver_depths_metres.size,
        fields.ranges_metres.size,
    )
    return fields.pressure.shape == expected_pressure_shape


def _validation_status(
    reference: AcousticFields, actual: AcousticFields
) -> str:
    if reference.wavenumbers.size != actual.wavenumbers.size:
        return "MODE_COUNT_MISMATCH"
    if not _all_finite(reference) or not _all_finite(actual):
        return "NONFINITE_VALUE"
    if not _internal_shape_is_valid(reference) or not _internal_shape_is_valid(
        actual
    ):
        return "SHAPE_MISMATCH"
    shape_pairs = (
        (reference.wavenumbers, actual.wavenumbers),
        (reference.source_depths_metres, actual.source_depths_metres),
        (reference.receiver_depths_metres, actual.receiver_depths_metres),
        (reference.ranges_metres, actual.ranges_metres),
        (reference.pressure, actual.pressure),
    )
    if any(left.shape != right.shape for left, right in shape_pairs):
        return "SHAPE_MISMATCH"
    coordinate_pairs = (
        (reference.source_depths_metres, actual.source_depths_metres),
        (reference.receiver_depths_metres, actual.receiver_depths_metres),
        (reference.ranges_metres, actual.ranges_metres),
    )
    if any(
        not np.allclose(
            left, right, rtol=0.0, atol=GRID_ATOL_METRES
        )
        for left, right in coordinate_pairs
    ):
        return "GRID_MISMATCH"
    if not np.isclose(
        reference.frequency_hz,
        actual.frequency_hz,
        rtol=0.0,
        atol=FREQUENCY_ATOL_HZ,
    ):
        return "GRID_MISMATCH"
    positive_ranges = reference.ranges_metres > 0.0
    receiver_slice = (
        slice(1, -1)
        if reference.receiver_depths_metres.size > 2
        else slice(None)
    )
    tl_pressure = reference.pressure[:, receiver_slice, :][
        :, :, positive_ranges
    ]
    if tl_pressure.size == 0:
        return "EMPTY_TL_DOMAIN"
    if reference.wavenumbers.size == 0:
        return "EMPTY_METRIC_DOMAIN"
    return "OK"


def _candidate_metrics(
    reference: AcousticFields, actual: AcousticFields
) -> dict[str, dict[str, float]]:
    modal_errors = [
        relative_imaginary_error(actual_value, reference_value)
        for actual_value, reference_value in zip(
            actual.wavenumbers.tolist(),
            reference.wavenumbers.tolist(),
            strict=True,
        )
    ]
    positive_ranges = reference.ranges_metres > 0.0
    receiver_slice = (
        slice(1, -1)
        if reference.receiver_depths_metres.size > 2
        else slice(None)
    )
    reference_pressure = reference.pressure[:, receiver_slice, :][
        :, :, positive_ranges
    ]
    actual_pressure = actual.pressure[:, receiver_slice, :][
        :, :, positive_ranges
    ]
    tl_errors = np.abs(
        transmission_loss(actual_pressure)
        - transmission_loss(reference_pressure)
    ).ravel()
    return {
        "e_im": metric_summary(modal_errors),
        "tl_db": metric_summary(tl_errors.tolist()),
    }


def evaluate_candidate(
    reference: AcousticFields, actual: AcousticFields
) -> dict[str, Any]:
    status = _validation_status(reference, actual)
    if status != "OK":
        return {"status": status}
    return {
        "status": "OK",
        "metrics": _candidate_metrics(reference, actual),
    }


def _cross_run_status(before: RunPair, after: RunPair) -> str:
    canonical = before.reference
    fields = (
        before.actual,
        after.reference,
        after.actual,
    )
    for current in fields:
        arrays = (
            (canonical.wavenumbers, current.wavenumbers),
            (canonical.source_depths_metres, current.source_depths_metres),
            (
                canonical.receiver_depths_metres,
                current.receiver_depths_metres,
            ),
            (canonical.ranges_metres, current.ranges_metres),
            (canonical.pressure, current.pressure),
        )
        if any(first.shape != second.shape for first, second in arrays):
            return "SHAPE_MISMATCH"
        coordinates = (
            (
                canonical.source_depths_metres,
                current.source_depths_metres,
            ),
            (
                canonical.receiver_depths_metres,
                current.receiver_depths_metres,
            ),
            (canonical.ranges_metres, current.ranges_metres),
        )
        if any(
            not np.allclose(
                first, second, rtol=0.0, atol=GRID_ATOL_METRES
            )
            for first, second in coordinates
        ):
            return "GRID_MISMATCH"
        if not np.isclose(
            canonical.frequency_hz,
            current.frequency_hz,
            rtol=0.0,
            atol=FREQUENCY_ATOL_HZ,
        ):
            return "GRID_MISMATCH"
    return "OK"


def _baseline_self_check(
    actual: dict[str, dict[str, float]],
    expected: dict[str, dict[str, float]],
) -> dict[str, Any]:
    checks = {
        "e_im.median": (
            actual["e_im"]["median"],
            expected["e_im"]["median"],
            E_IM_REGRESSION_ATOL,
        ),
        "e_im.p95": (
            actual["e_im"]["p95"],
            expected["e_im"]["p95"],
            E_IM_REGRESSION_ATOL,
        ),
        "tl_db.median": (
            actual["tl_db"]["median"],
            expected["tl_db"]["median"],
            TL_REGRESSION_ATOL_DB,
        ),
        "tl_db.p95": (
            actual["tl_db"]["p95"],
            expected["tl_db"]["p95"],
            TL_REGRESSION_ATOL_DB,
        ),
    }
    passed = all(
        abs(actual_value - expected_value) <= tolerance
        for actual_value, expected_value, tolerance in checks.values()
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "expected": expected,
    }


def metrics_regressed(
    before: dict[str, dict[str, float]],
    after: dict[str, dict[str, float]],
) -> bool:
    return any(
        (
            after["e_im"][statistic]
            > before["e_im"][statistic] + E_IM_REGRESSION_ATOL
        )
        or (
            after["tl_db"][statistic]
            > before["tl_db"][statistic] + TL_REGRESSION_ATOL_DB
        )
        for statistic in ("median", "p95")
    )


def compare_runs(
    before: RunPair,
    after: RunPair,
    *,
    expected_before_metrics: dict[str, dict[str, float]] | None = None,
) -> dict[str, Any]:
    for phase, pair in (("before", before), ("after", after)):
        status = _validation_status(pair.reference, pair.actual)
        if status != "OK":
            return {"status": status, "verdict": "FAIL", "phase": phase}
    cross_status = _cross_run_status(before, after)
    if cross_status != "OK":
        return {
            "status": cross_status,
            "verdict": "FAIL",
            "phase": "before_after",
        }

    before_metrics = _candidate_metrics(before.reference, before.actual)
    after_metrics = _candidate_metrics(after.reference, after.actual)
    result: dict[str, Any] = {
        "status": "OK",
        "verdict": "PASS",
        "before": before_metrics,
        "after": after_metrics,
        "thresholds": {
            "e_im_absolute": E_IM_REGRESSION_ATOL,
            "tl_db_absolute": TL_REGRESSION_ATOL_DB,
        },
    }
    if expected_before_metrics is not None:
        self_check = _baseline_self_check(
            before_metrics, expected_before_metrics
        )
        result["baseline_self_check"] = self_check
        if self_check["status"] != "PASS":
            result["status"] = "BASELINE_SELF_CHECK_FAILED"
            result["verdict"] = "FAIL"
            return result

    if metrics_regressed(before_metrics, after_metrics):
        result["status"] = "REGRESSION"
        result["verdict"] = "FAIL"
    return result


def _unique_artifact(directory: Path, suffix: str) -> Path:
    candidates = sorted(directory.glob(f"*{suffix}"))
    if len(candidates) != 1 or not candidates[0].is_file():
        raise RegressionError(
            "INVALID_RUN_LAYOUT",
            f"expected exactly one {suffix} file in {directory}",
        )
    return candidates[0]


def discover_artifacts(run_dir: Path) -> ArtifactPaths:
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise RegressionError(
            "INVALID_RUN_LAYOUT", f"run directory does not exist: {run_dir}"
        )
    reference_directory = run_dir / "krakenc"
    actual_directory = run_dir / "ookc"
    if not reference_directory.is_dir() or not actual_directory.is_dir():
        raise RegressionError(
            "INVALID_RUN_LAYOUT",
            "run directory must contain krakenc and ookc directories",
        )
    reference_mod = _unique_artifact(reference_directory, ".mod")
    reference_shd = _unique_artifact(reference_directory, ".shd")
    actual_mod = _unique_artifact(actual_directory, ".mod")
    actual_shd = _unique_artifact(actual_directory, ".shd")
    stems = {
        path.stem
        for path in (
            reference_mod,
            reference_shd,
            actual_mod,
            actual_shd,
        )
    }
    if len(stems) != 1:
        raise RegressionError(
            "INVALID_RUN_LAYOUT",
            "krakenc and ookc MOD/SHD stems must match",
        )
    return ArtifactPaths(
        case_stem=stems.pop(),
        reference_mod=reference_mod,
        reference_shd=reference_shd,
        actual_mod=actual_mod,
        actual_shd=actual_shd,
    )


def _read_fields(mod_path: Path, shd_path: Path) -> AcousticFields:
    wavenumbers = read_first_wavenumber_set(mod_path)
    shade = read_shade_file(shd_path)
    expected_pressure_size = (
        len(shade.source_depths)
        * len(shade.receiver_depths)
        * len(shade.ranges_metres)
    )
    if len(shade.pressure) != expected_pressure_size:
        raise RegressionError(
            "INVALID_SHD_LAYOUT",
            f"unsupported pressure dimensions in {shd_path}",
        )
    pressure = np.asarray(shade.pressure, dtype=np.complex128).reshape(
        (
            len(shade.source_depths),
            len(shade.receiver_depths),
            len(shade.ranges_metres),
        )
    )
    return AcousticFields(
        wavenumbers=np.asarray(
            wavenumbers, dtype=np.complex128
        ),
        frequency_hz=float(shade.frequency),
        source_depths_metres=np.asarray(
            shade.source_depths, dtype=np.float64
        ),
        receiver_depths_metres=np.asarray(
            shade.receiver_depths, dtype=np.float64
        ),
        ranges_metres=np.asarray(shade.ranges_metres, dtype=np.float64),
        pressure=pressure,
    )


def load_run_pair(run_dir: Path) -> tuple[str, RunPair]:
    artifacts = discover_artifacts(run_dir)
    return (
        artifacts.case_stem,
        RunPair(
            reference=_read_fields(
                artifacts.reference_mod, artifacts.reference_shd
            ),
            actual=_read_fields(artifacts.actual_mod, artifacts.actual_shd),
        ),
    )


def _atomic_write_json(output: Path, payload: dict[str, Any]) -> None:
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=output.name + ".",
            suffix=".tmp",
            dir=output.parent,
            delete=False,
        ) as stream:
            temporary_path = Path(stream.name)
            json.dump(
                payload,
                stream,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            stream.write("\n")
        temporary_path.replace(output)
    except BaseException:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-run-dir", required=True, type=Path)
    parser.add_argument("--after-run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    before_path = arguments.before_run_dir.resolve()
    after_path = arguments.after_run_dir.resolve()
    base_payload: dict[str, Any] = {
        "schema": "OpenOcean-Krakenc.ssp-interpolation-regression",
        "schema_version": 1,
        "before_run_dir": str(before_path),
        "after_run_dir": str(after_path),
    }
    try:
        before_case, before = load_run_pair(before_path)
        after_case, after = load_run_pair(after_path)
        if before_case != after_case:
            raise RegressionError(
                "CASE_MISMATCH",
                f"before case {before_case!r} != after case {after_case!r}",
            )
        comparison = compare_runs(
            before,
            after,
            expected_before_metrics=KNOWN_BEFORE_METRICS,
        )
        payload = {**base_payload, "case": before_case, **comparison}
    except RegressionError as error:
        payload = {
            **base_payload,
            "status": error.status,
            "verdict": "FAIL",
            "detail": error.detail,
        }
    except (OSError, ValueError, ModeFileError, struct.error) as error:
        payload = {
            **base_payload,
            "status": "ARTIFACT_READ_ERROR",
            "verdict": "FAIL",
            "detail": str(error),
        }
    _atomic_write_json(arguments.output, payload)
    return 0 if payload["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
