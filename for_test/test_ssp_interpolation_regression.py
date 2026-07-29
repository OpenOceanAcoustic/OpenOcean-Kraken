from __future__ import annotations

import json
import math
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from for_test import ssp_interpolation_regression as regression
from for_test.shd_reader import read_shade_file


def _fields(
    *,
    wavenumbers: list[complex] | None = None,
    frequency_hz: float = 100.0,
    source_depths: list[float] | None = None,
    receiver_depths: list[float] | None = None,
    ranges: list[float] | None = None,
    pressure: np.ndarray | None = None,
) -> regression.AcousticFields:
    source_values = [25.0] if source_depths is None else source_depths
    receiver_values = (
        [0.0, 50.0, 100.0, 150.0]
        if receiver_depths is None
        else receiver_depths
    )
    range_values = [-10.0, 0.0, 1000.0] if ranges is None else ranges
    if pressure is None:
        pressure = np.ones(
            (
                len(source_values),
                len(receiver_values),
                len(range_values),
            ),
            dtype=np.complex128,
        )
    return regression.AcousticFields(
        wavenumbers=np.asarray(
            [1.0 + 1.0j, 2.0 + 2.0j]
            if wavenumbers is None
            else wavenumbers,
            dtype=np.complex128,
        ),
        frequency_hz=frequency_hz,
        source_depths_metres=np.asarray(source_values, dtype=np.float64),
        receiver_depths_metres=np.asarray(
            receiver_values, dtype=np.float64
        ),
        ranges_metres=np.asarray(range_values, dtype=np.float64),
        pressure=np.asarray(pressure, dtype=np.complex128),
    )


def _fields_with_frequency(
    frequency_hz: float, **kwargs: object
) -> regression.AcousticFields:
    return _fields(frequency_hz=frequency_hz, **kwargs)


def _write_test_mod(path: Path) -> None:
    record_words = 25
    record_bytes = 4 * record_words
    data = bytearray(29 * record_bytes)
    struct.pack_into("<i", data, 0, record_words)
    struct.pack_into("<iiii", data, 84, 0, 0, 1, 1)
    struct.pack_into("<i", data, 5 * record_bytes, 2)
    struct.pack_into("<ff", data, 9 * record_bytes, 1.25, -0.5)
    struct.pack_into("<ff", data, 9 * record_bytes + 8, 2.5, -1.0)
    path.write_bytes(data)


class FormulaTests(unittest.TestCase):
    def test_relative_imaginary_error_uses_reference_floor(self) -> None:
        self.assertEqual(
            regression.relative_imaginary_error(3.0 + 2.0e-13j, 0.0j),
            0.2,
        )

    def test_metric_summary_uses_standard_median_for_odd_and_even(self) -> None:
        self.assertEqual(regression.metric_summary([9.0, 1.0, 5.0])["median"], 5.0)
        self.assertEqual(
            regression.metric_summary([8.0, 2.0, 4.0, 6.0])["median"],
            5.0,
        )

    def test_nearest_rank_p95_uses_ceil_index(self) -> None:
        values = [float(value) for value in range(20, 0, -1)]
        self.assertEqual(regression.nearest_rank(values, 0.95), 19.0)

    def test_nearest_rank_rejects_empty_metric_domain(self) -> None:
        with self.assertRaisesRegex(ValueError, "^EMPTY_METRIC_DOMAIN$"):
            regression.nearest_rank([], 0.95)

    def test_transmission_loss_uses_pressure_floor(self) -> None:
        actual = regression.transmission_loss(
            np.asarray([1.0, 0.1, 0.0], dtype=np.complex128)
        )
        np.testing.assert_allclose(actual, [0.0, 20.0, 600.0], atol=1.0e-12)


class CandidateEvaluationTests(unittest.TestCase):
    def test_same_run_frequency_mismatch_has_no_statistics(self) -> None:
        result = regression.evaluate_candidate(
            _fields_with_frequency(100.0),
            _fields_with_frequency(100.0 + 2.0e-8),
        )
        self.assertEqual(result["status"], "GRID_MISMATCH")
        self.assertNotIn("metrics", result)

    def test_same_run_frequency_within_tolerance_is_accepted(self) -> None:
        result = regression.evaluate_candidate(
            _fields_with_frequency(100.0),
            _fields_with_frequency(100.0 + 1.0e-9),
        )
        self.assertEqual(result["status"], "OK")

    def test_nonfinite_frequency_has_no_statistics(self) -> None:
        result = regression.evaluate_candidate(
            _fields_with_frequency(100.0),
            _fields_with_frequency(math.nan),
        )
        self.assertEqual(result["status"], "NONFINITE_VALUE")
        self.assertNotIn("metrics", result)

    def test_mode_count_mismatch_has_no_error_statistics(self) -> None:
        result = regression.evaluate_candidate(
            _fields(wavenumbers=[1.0 + 1.0j]),
            _fields(wavenumbers=[1.0 + 1.0j, 2.0 + 2.0j]),
        )
        self.assertEqual(result["status"], "MODE_COUNT_MISMATCH")
        self.assertNotIn("metrics", result)

    def test_empty_strictly_positive_tl_domain_has_no_statistics(self) -> None:
        reference = _fields(ranges=[-1.0, 0.0])
        result = regression.evaluate_candidate(reference, reference)
        self.assertEqual(result["status"], "EMPTY_TL_DOMAIN")
        self.assertNotIn("metrics", result)

    def test_nonfinite_arrays_fail_without_statistics(self) -> None:
        nonfinite_pressure = np.ones((1, 4, 3), dtype=np.complex128)
        nonfinite_pressure[0, 2, 2] = complex(math.nan, 0.0)
        cases = {
            "wavenumber": _fields(
                wavenumbers=[1.0 + math.inf * 1j, 2.0 + 2.0j]
            ),
            "source coordinate": _fields(source_depths=[math.nan]),
            "receiver coordinate": _fields(receiver_depths=[0.0, math.inf]),
            "range coordinate": _fields(ranges=[0.0, math.nan]),
            "pressure": _fields(pressure=nonfinite_pressure),
        }
        reference = _fields()
        for name, actual in cases.items():
            with self.subTest(name=name):
                result = regression.evaluate_candidate(reference, actual)
                self.assertEqual(result["status"], "NONFINITE_VALUE")
                self.assertNotIn("metrics", result)

    def test_grid_difference_of_one_nanometre_is_accepted(self) -> None:
        reference = _fields()
        actual = _fields(ranges=[-10.0, 0.0, 1000.0 + 1.0e-9])
        result = regression.evaluate_candidate(reference, actual)
        self.assertEqual(result["status"], "OK")

    def test_grid_difference_of_twenty_nanometres_is_rejected(self) -> None:
        reference = _fields()
        actual = _fields(ranges=[-10.0, 0.0, 1000.0 + 2.0e-8])
        result = regression.evaluate_candidate(reference, actual)
        self.assertEqual(result["status"], "GRID_MISMATCH")
        self.assertNotIn("metrics", result)

    def test_different_pressure_grid_is_rejected_without_resampling(self) -> None:
        reference = _fields()
        actual = _fields(
            receiver_depths=[0.0, 50.0, 150.0],
            pressure=np.ones((1, 3, 3), dtype=np.complex128),
        )
        result = regression.evaluate_candidate(reference, actual)
        self.assertEqual(result["status"], "SHAPE_MISMATCH")
        self.assertNotIn("metrics", result)

    def test_tl_domain_uses_positive_ranges_and_interior_receiver_depths(self) -> None:
        reference = _fields()
        pressure = np.ones((1, 4, 3), dtype=np.complex128)
        pressure[:, 1:3, 2] = 0.1
        actual = _fields(pressure=pressure)
        result = regression.evaluate_candidate(reference, actual)
        self.assertEqual(
            result["metrics"]["tl_db"],
            {"median": 20.0, "p95": 20.0},
        )


class AcceptanceTests(unittest.TestCase):
    def test_cross_run_frequency_mismatch_uses_canonical_anchor(self) -> None:
        before = regression.RunPair(
            reference=_fields_with_frequency(100.0),
            actual=_fields_with_frequency(100.0),
        )
        after = regression.RunPair(
            reference=_fields_with_frequency(100.0 + 9.0e-9),
            actual=_fields_with_frequency(100.0 + 18.0e-9),
        )

        result = regression.compare_runs(before, after)

        self.assertEqual(result["status"], "GRID_MISMATCH")
        self.assertNotIn("before", result)
        self.assertNotIn("after", result)

    def test_all_frequency_roles_within_canonical_tolerance_are_accepted(
        self,
    ) -> None:
        before = regression.RunPair(
            reference=_fields_with_frequency(100.0),
            actual=_fields_with_frequency(100.0 - 9.0e-9),
        )
        after = regression.RunPair(
            reference=_fields_with_frequency(100.0 + 1.0e-9),
            actual=_fields_with_frequency(100.0 + 9.0e-9),
        )

        result = regression.compare_runs(before, after)

        self.assertEqual(result["status"], "OK")
        self.assertEqual(result["verdict"], "PASS")

    def test_known_before_metrics_are_computed_by_controlled_data(self) -> None:
        known = {
            "e_im": {
                "median": 5.289488948137276e-4,
                "p95": 7.564144223952983e-4,
            },
            "tl_db": {
                "median": 8.269118904209449e-4,
                "p95": 3.7592181663796964e-3,
            },
        }
        self.assertEqual(regression.KNOWN_BEFORE_METRICS, known)
        modal_median = known["e_im"]["median"]
        modal_p95 = known["e_im"]["p95"]
        tl_median = known["tl_db"]["median"]
        tl_p95 = known["tl_db"]["p95"]

        modal_errors = [0.0] * 9 + [modal_median] * 9 + [modal_p95] * 2
        reference_wavenumbers = [complex(index + 1.0, 1.0) for index in range(20)]
        actual_wavenumbers = [
            complex(value.real, value.imag + error)
            for value, error in zip(
                reference_wavenumbers, modal_errors, strict=True
            )
        ]
        tl_errors = [0.0] * 9 + [tl_median] * 9 + [tl_p95] * 2
        actual_pressure = np.asarray(
            [10.0 ** (-error / 20.0) for error in tl_errors],
            dtype=np.complex128,
        ).reshape(1, 2, 10)
        reference = _fields(
            wavenumbers=reference_wavenumbers,
            receiver_depths=[25.0, 75.0],
            ranges=[float(index + 1) for index in range(10)],
            pressure=np.ones((1, 2, 10), dtype=np.complex128),
        )
        before = _fields(
            wavenumbers=actual_wavenumbers,
            receiver_depths=[25.0, 75.0],
            ranges=[float(index + 1) for index in range(10)],
            pressure=actual_pressure,
        )

        result = regression.compare_runs(
            regression.RunPair(reference=reference, actual=before),
            regression.RunPair(reference=reference, actual=before),
            expected_before_metrics=known,
        )

        self.assertEqual(result["status"], "OK")
        self.assertEqual(result["verdict"], "PASS")
        self.assertEqual(result["baseline_self_check"]["status"], "PASS")
        self.assertAlmostEqual(
            result["before"]["e_im"]["median"], modal_median, delta=1.0e-12
        )
        self.assertAlmostEqual(
            result["before"]["e_im"]["p95"], modal_p95, delta=1.0e-12
        )
        self.assertAlmostEqual(
            result["before"]["tl_db"]["median"], tl_median, delta=1.0e-9
        )
        self.assertAlmostEqual(
            result["before"]["tl_db"]["p95"], tl_p95, delta=1.0e-9
        )

    def test_after_thresholds_are_inclusive_and_regressions_fail(self) -> None:
        reference = _fields()
        before = _fields(
            wavenumbers=[1.0 + 1.1j, 2.0 + 2.2j],
            pressure=np.full((1, 4, 3), 0.9, dtype=np.complex128),
        )
        same = regression.compare_runs(
            regression.RunPair(reference, before),
            regression.RunPair(reference, before),
        )
        self.assertEqual(same["verdict"], "PASS")

        worse = _fields(
            wavenumbers=[1.0 + 1.2j, 2.0 + 2.4j],
            pressure=np.full((1, 4, 3), 0.8, dtype=np.complex128),
        )
        regressed = regression.compare_runs(
            regression.RunPair(reference, before),
            regression.RunPair(reference, worse),
        )
        self.assertEqual(regressed["status"], "REGRESSION")
        self.assertEqual(regressed["verdict"], "FAIL")

    def test_metric_regression_thresholds_are_exact_and_inclusive(self) -> None:
        regression_check = getattr(regression, "metrics_regressed", None)
        self.assertIsNotNone(
            regression_check,
            "metrics_regressed must expose the exact threshold decision",
        )
        before = {
            "e_im": {"median": 1.0, "p95": 1.0},
            "tl_db": {"median": 1.0, "p95": 1.0},
        }
        tolerances = {
            "e_im": 1.0e-12,
            "tl_db": 1.0e-9,
        }
        for metric, tolerance in tolerances.items():
            for statistic in ("median", "p95"):
                with self.subTest(
                    metric=metric, statistic=statistic, position="boundary"
                ):
                    after = {
                        name: dict(values)
                        for name, values in before.items()
                    }
                    boundary = before[metric][statistic] + tolerance
                    after[metric][statistic] = boundary
                    self.assertFalse(
                        regression_check(before, after)
                    )
                with self.subTest(
                    metric=metric, statistic=statistic, position="above"
                ):
                    after[metric][statistic] = math.nextafter(
                        boundary, math.inf
                    )
                    self.assertTrue(
                        regression_check(before, after)
                    )

    def test_all_four_coordinate_grids_use_one_canonical_anchor(self) -> None:
        canonical_range = 1000.0 - 9.0e-9
        centre_range = 1000.0
        far_range = 1000.0 + 9.0e-9
        before = regression.RunPair(
            reference=_fields(ranges=[canonical_range]),
            actual=_fields(ranges=[centre_range]),
        )
        after = regression.RunPair(
            reference=_fields(ranges=[centre_range]),
            actual=_fields(ranges=[far_range]),
        )

        result = regression.compare_runs(before, after)

        self.assertEqual(result["status"], "GRID_MISMATCH")
        self.assertNotIn("before", result)
        self.assertNotIn("after", result)

    def test_validation_failure_emits_no_before_or_after_statistics(self) -> None:
        reference = _fields()
        mismatched = _fields(wavenumbers=[1.0 + 1.0j])
        result = regression.compare_runs(
            regression.RunPair(reference, mismatched),
            regression.RunPair(reference, reference),
        )
        self.assertEqual(result["status"], "MODE_COUNT_MISMATCH")
        self.assertNotIn("before", result)
        self.assertNotIn("after", result)


class ArtifactLayoutTests(unittest.TestCase):
    def test_field_reader_uses_first_mode_set_when_mod_has_trailing_records(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            mod_path = Path(temporary) / "case.mod"
            shd_path = Path(temporary) / "case.shd"
            _write_test_mod(mod_path)
            shd_data = bytearray(
                Path("for_test/fixtures/two_profile_small.shd").read_bytes()
            )
            (record_words,) = struct.unpack_from("<i", shd_data, 0)
            struct.pack_into(
                "<d", shd_data, 3 * 4 * record_words, 73.25
            )
            shd_path.write_bytes(shd_data)

            fields = regression._read_fields(
                mod_path, shd_path
            )
            shade = read_shade_file(shd_path)

            np.testing.assert_array_equal(
                fields.wavenumbers,
                np.asarray([1.25 - 0.5j, 2.5 - 1.0j]),
            )
            self.assertEqual(shade.frequency, 73.25)
            self.assertEqual(fields.frequency_hz, shade.frequency)

    def test_run_layout_requires_unique_matching_case_stem(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            for engine in ("krakenc", "ookc"):
                directory = run_dir / engine
                directory.mkdir()
                (directory / "case.mod").touch()
                (directory / "case.shd").touch()

            paths = regression.discover_artifacts(run_dir)

            self.assertEqual(paths.case_stem, "case")
            self.assertEqual(paths.reference_mod, run_dir / "krakenc" / "case.mod")
            self.assertEqual(paths.actual_shd, run_dir / "ookc" / "case.shd")

    def test_run_layout_rejects_multiple_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            for engine in ("krakenc", "ookc"):
                directory = run_dir / engine
                directory.mkdir()
                (directory / "case.mod").touch()
                (directory / "case.shd").touch()
            (run_dir / "ookc" / "extra.mod").touch()

            with self.assertRaisesRegex(
                regression.RegressionError, "^INVALID_RUN_LAYOUT:"
            ):
                regression.discover_artifacts(run_dir)


class CliTests(unittest.TestCase):
    def test_truncated_shd_writes_failure_json_without_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid_shd = Path(
                "for_test/fixtures/two_profile_small.shd"
            ).read_bytes()
            (record_words,) = struct.unpack_from("<i", valid_shd, 0)
            truncated_shd = valid_shd[: 10 * 4 * record_words]
            for phase in ("before", "after"):
                run_dir = root / phase
                for engine in ("krakenc", "ookc"):
                    directory = run_dir / engine
                    directory.mkdir(parents=True)
                    _write_test_mod(directory / "case.mod")
                    (directory / "case.shd").write_bytes(truncated_shd)
            output = root / "result.json"

            return_code = regression.main(
                [
                    "--before-run-dir",
                    str(root / "before"),
                    "--after-run-dir",
                    str(root / "after"),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(return_code, 1)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "ARTIFACT_READ_ERROR")
            self.assertEqual(payload["verdict"], "FAIL")
            self.assertNotIn("before", payload)
            self.assertNotIn("after", payload)

    def test_corrupt_artifacts_write_failure_json_without_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for phase in ("before", "after"):
                run_dir = root / phase
                for engine in ("krakenc", "ookc"):
                    directory = run_dir / engine
                    directory.mkdir(parents=True)
                    (directory / "case.mod").touch()
                    (directory / "case.shd").touch()
            output = root / "result.json"

            return_code = regression.main(
                [
                    "--before-run-dir",
                    str(root / "before"),
                    "--after-run-dir",
                    str(root / "after"),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(return_code, 1)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "ARTIFACT_READ_ERROR")
            self.assertEqual(payload["verdict"], "FAIL")
            self.assertNotIn("before", payload)
            self.assertNotIn("after", payload)


if __name__ == "__main__":
    unittest.main()
