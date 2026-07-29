import unittest
import json
from pathlib import Path
import tempfile
import warnings

from for_test.generate_numerical_consistency import (
    build_records,
    generate_artifacts,
    summary_payload,
)
from for_test.numerical_consistency import classify_slice


def metric(**overrides):
    values = {
        "case": "synthetic",
        "source_depth_m": 50.0,
        "complex_relative_l2": 1.0e-5,
        "positive_range_relative_l2": 1.0e-5,
        "interior_positive_range_relative_l2": 1.0e-5,
        "zero_range_reference_energy_fraction": 0.01,
        "boundary_error_fraction": 0.01,
        "tl_p95_db": 0.001,
        "tl_max_db": 0.01,
        "tl_rmse_db": 0.001,
        "max_difference_range_km": 1.0,
        "max_difference_depth_m": 20.0,
    }
    values.update(overrides)
    return values


class NumericalConsistencyTests(unittest.TestCase):
    def test_classification_uses_approved_pass_watch_fail_boundaries(self):
        passing = classify_slice(
            metric(
                interior_positive_range_relative_l2=1.0e-3,
                tl_p95_db=0.05,
            )
        )
        watching = classify_slice(
            metric(
                interior_positive_range_relative_l2=5.0e-3,
                tl_p95_db=0.5,
            )
        )
        failing = classify_slice(
            metric(interior_positive_range_relative_l2=5.01e-3)
        )

        self.assertEqual(passing.grade, "pass")
        self.assertEqual(watching.grade, "watch")
        self.assertEqual(failing.grade, "fail")

    def test_anomaly_attribution_supports_multiple_regions(self):
        result = classify_slice(
            metric(
                complex_relative_l2=0.4,
                positive_range_relative_l2=2.0e-3,
                interior_positive_range_relative_l2=2.0e-3,
                zero_range_reference_energy_fraction=0.9,
                boundary_error_fraction=0.97,
                tl_p95_db=0.1,
                tl_max_db=5.0,
            )
        )

        self.assertEqual(result.grade, "watch")
        self.assertEqual(
            set(result.anomaly_regions),
            {"zero_range", "boundary", "interior_positive_range"},
        )

    def test_local_spike_is_distinguished_from_broad_tl_error(self):
        spike = classify_slice(metric(tl_p95_db=0.01, tl_max_db=2.0))
        clean = classify_slice(metric(tl_p95_db=0.01, tl_max_db=0.05))

        self.assertEqual(spike.anomaly_regions, ("local_interference_spike",))
        self.assertEqual(clean.anomaly_regions, ("none",))

    def test_generator_writes_all_acceptance_artifacts(self):
        metrics = [
            metric(case="pass_case"),
            metric(
                case="watch_case",
                interior_positive_range_relative_l2=2.0e-3,
                tl_p95_db=0.1,
            ),
        ]
        records = build_records(metrics)
        summary = summary_payload(records)

        self.assertEqual(summary["slice_count"], 2)
        self.assertEqual(summary["grade_counts"], {"pass": 1, "watch": 1, "fail": 0})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps({"case_count": 2, "metrics": metrics}),
                encoding="utf-8",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                generated = generate_artifacts(manifest, root)

            self.assertEqual(len(generated), 4)
            self.assertFalse(
                [warning for warning in caught if "Glyph" in str(warning.message)]
            )
            for path in generated:
                self.assertTrue(path.is_file(), path)
                self.assertGreater(path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
