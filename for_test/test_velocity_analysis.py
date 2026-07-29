import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

from for_test.shd_reader import ShadeFile
from for_test.shd_reader import read_shade_file
from for_test.velocity_analysis import (
    assert_matching_shade_grids,
    compute_velocity_slice_metrics,
    reshape_shade,
    surface_boundary_metrics,
)
from for_test.plot_velocity_comparison import (
    analyze_velocity_shades,
    discover_velocity_cases,
    plot_velocity_slice,
    safe_reset_velocity_output,
    write_velocity_analysis,
    write_velocity_metrics,
    velocity_manifest_payload,
    velocity_figure_name,
)


class VelocityAnalysisTests(unittest.TestCase):
    def test_reshape_shade_expands_sources_receivers_and_ranges(self):
        shade = ShadeFile(
            title="synthetic",
            frequency=50.0,
            source_depths=[10.0, 20.0],
            receiver_depths=[0.0, 5.0],
            ranges_metres=[0.0, 1000.0],
            pressure=[complex(value) for value in range(8)],
        )

        cube = reshape_shade(shade)

        self.assertEqual(cube.shape, (2, 2, 2))
        self.assertEqual(cube[1, 1, 1], 7.0 + 0.0j)

    def test_velocity_metrics_exclude_surface_and_zero_range(self):
        pressure = np.ones((3, 3), dtype=np.complex128)
        horizontal = pressure * (10.0 ** (-0.1 / 20.0))
        vertical = pressure * 0.5

        metrics = compute_velocity_slice_metrics(
            pressure,
            horizontal,
            vertical,
            receiver_depths=np.array([0.0, 1.0, 2.0]),
            ranges_metres=np.array([0.0, 100.0, 200.0]),
            source_depth=1.0,
        )

        self.assertEqual(metrics.interior_point_count, 4)
        self.assertAlmostEqual(
            metrics.pressure_horizontal_interior.p50_abs_delta_db,
            0.1,
            places=10,
        )
        self.assertAlmostEqual(
            metrics.pressure_vertical_interior.p50_abs_delta_db,
            20.0 * np.log10(2.0),
            places=10,
        )
        self.assertEqual(metrics.line_receiver_depth_m, 1.0)
        self.assertEqual(metrics.pressure_horizontal_source_line.valid_points, 2)

    def test_surface_boundary_and_repeated_scaling_are_diagnostics(self):
        pressure = np.ones((3, 2), dtype=np.complex128)
        horizontal = pressure.copy()
        vertical = pressure.copy()
        pressure[0, :] = 0.0
        horizontal[0, :] = 0.0

        metrics = compute_velocity_slice_metrics(
            pressure,
            horizontal,
            vertical,
            receiver_depths=np.array([0.0, 1.0, 2.0]),
            ranges_metres=np.array([100.0, 200.0]),
            source_depth=1.0,
            rho_c0=1500.0,
        )

        self.assertTrue(metrics.pressure_surface.is_near_zero)
        self.assertTrue(metrics.horizontal_surface.is_near_zero)
        self.assertFalse(metrics.vertical_surface.is_near_zero)
        self.assertAlmostEqual(
            metrics.rescaled_horizontal_p50_abs_delta_db,
            20.0 * np.log10(1500.0),
            places=10,
        )

    def test_surface_metrics_without_a_surface_row_are_unavailable(self):
        metrics = compute_velocity_slice_metrics(
            np.ones((2, 2), dtype=np.complex128),
            np.ones((2, 2), dtype=np.complex128),
            np.ones((2, 2), dtype=np.complex128),
            receiver_depths=np.array([5.0, 10.0]),
            ranges_metres=np.array([100.0, 200.0]),
            source_depth=5.0,
        )

        for surface in (
            metrics.pressure_surface,
            metrics.horizontal_surface,
            metrics.vertical_surface,
        ):
            self.assertFalse(surface.has_surface_row)
            self.assertIsNone(surface.surface_depth_m)
            self.assertIsNone(surface.surface_max_abs)
            self.assertIsNone(surface.surface_to_global_max_ratio)
            self.assertIsNone(surface.first_interior_to_global_max_ratio)
            self.assertEqual(surface.global_max_abs, 1.0)
            self.assertEqual(surface.zero_ratio_limit, 1.0e-8)
            self.assertFalse(surface.is_near_zero)

    def test_surface_metrics_without_an_interior_row_leave_interior_ratio_unavailable(self):
        surface = surface_boundary_metrics(
            np.ones((1, 2), dtype=np.complex128), np.array([0.0])
        )

        self.assertTrue(surface.has_surface_row)
        self.assertEqual(surface.surface_depth_m, 0.0)
        self.assertEqual(surface.surface_max_abs, 1.0)
        self.assertEqual(surface.surface_to_global_max_ratio, 1.0)
        self.assertIsNone(surface.first_interior_to_global_max_ratio)

    def test_unavailable_surface_metrics_are_blank_in_csv_and_safe_for_artifacts(self):
        pressure = ShadeFile(
            "p",
            50.0,
            [10.0],
            [5.0, 10.0],
            [0.0, 1000.0, 2000.0],
            [
                1.0 + 0.0j,
                2.0 + 0.0j,
                3.0 + 0.0j,
                4.0 + 0.0j,
                5.0 + 0.0j,
                6.0 + 0.0j,
            ],
        )
        horizontal = ShadeFile(
            "h",
            50.0,
            [10.0],
            [5.0, 10.0],
            [0.0, 1000.0, 2000.0],
            [
                0.9 + 0.0j,
                1.8 + 0.0j,
                2.7 + 0.0j,
                3.6 + 0.0j,
                4.5 + 0.0j,
                5.4 + 0.0j,
            ],
        )
        vertical = ShadeFile(
            "v",
            50.0,
            [10.0],
            [5.0, 10.0],
            [0.0, 1000.0, 2000.0],
            [
                0.5 + 0.0j,
                1.0 + 0.0j,
                1.5 + 0.0j,
                2.0 + 0.0j,
                2.5 + 0.0j,
                3.0 + 0.0j,
            ],
        )
        rows = analyze_velocity_shades("no_surface", pressure, horizontal, vertical)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            csv_path = root / "metrics.csv"
            figure_path = root / "slice.png"
            analysis_path = root / "analysis.md"
            write_velocity_metrics(rows, csv_path)
            plot_velocity_slice(
                "no_surface", pressure, horizontal, vertical, 0, rows[0], figure_path
            )
            write_velocity_analysis(
                [{"case": "no_surface", "status": "success"}], rows, analysis_path
            )

            with csv_path.open(encoding="utf-8-sig", newline="") as stream:
                csv_row = next(csv.DictReader(stream))
            self.assertEqual(csv_row["pressure_surface_depth_m"], "")
            self.assertEqual(csv_row["pressure_surface_max_abs"], "")
            self.assertEqual(csv_row["pressure_surface_to_global_max_ratio"], "")
            self.assertEqual(
                csv_row["pressure_first_interior_to_global_max_ratio"], ""
            )
            self.assertTrue(figure_path.is_file())
            self.assertIn("n/a/n/a/n/a", analysis_path.read_text(encoding="utf-8"))

    def test_no_surface_manifest_payload_is_strict_json_serializable(self):
        shade = ShadeFile(
            "p",
            50.0,
            [10.0],
            [5.0, 10.0],
            [0.0, 1000.0, 2000.0],
            [
                1.0 + 0.0j,
                2.0 + 0.0j,
                3.0 + 0.0j,
                4.0 + 0.0j,
                5.0 + 0.0j,
                6.0 + 0.0j,
            ],
        )
        rows = analyze_velocity_shades("no_surface", shade, shade, shade)

        encoded = json.dumps(
            velocity_manifest_payload(
                [{"case": "no_surface", "status": "success", "metrics": rows}], rows
            ),
            allow_nan=False,
        )

        self.assertIn('"metrics"', encoded)

    def test_matching_shade_grids_rejects_range_mismatch(self):
        base = ShadeFile(
            "x", 50.0, [10.0], [0.0, 5.0], [0.0, 1000.0], [0j] * 4
        )
        changed = ShadeFile(
            "x", 50.0, [10.0], [0.0, 5.0], [0.0, 1100.0], [0j] * 4
        )

        with self.assertRaisesRegex(ValueError, "range grid mismatch"):
            assert_matching_shade_grids(base, base, changed)

    def test_velocity_case_discovery_is_recursive_and_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for folder, name in (("z_group", "beta"), ("a_group", "alpha")):
                directory = root / folder
                directory.mkdir()
                (directory / f"{name}.env").write_text("env", encoding="utf-8")
                (directory / f"{name}.flp").write_text("flp", encoding="utf-8")

            cases = discover_velocity_cases(root)

            self.assertEqual([case.name for case in cases], ["alpha", "beta"])

    def test_velocity_case_discovery_rejects_duplicate_stems(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for folder in ("one", "two"):
                directory = root / folder
                directory.mkdir()
                (directory / "same.env").write_text("env", encoding="utf-8")
                (directory / "same.flp").write_text("flp", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "duplicate ENV stem"):
                discover_velocity_cases(root)

    def test_velocity_analysis_expands_every_source_slice(self):
        pressure = ShadeFile(
            "p",
            50.0,
            [10.0, 20.0],
            [0.0, 10.0],
            [0.0, 1000.0],
            [1.0 + 0.0j] * 8,
        )
        horizontal = ShadeFile(
            "h",
            50.0,
            [10.0, 20.0],
            [0.0, 10.0],
            [0.0, 1000.0],
            [0.9 + 0.0j] * 8,
        )
        vertical = ShadeFile(
            "v",
            50.0,
            [10.0, 20.0],
            [0.0, 10.0],
            [0.0, 1000.0],
            [0.5 + 0.0j] * 8,
        )

        rows = analyze_velocity_shades("multi", pressure, horizontal, vertical)

        self.assertEqual([row["source_index"] for row in rows], [0, 1])
        self.assertEqual([row["source_depth_m"] for row in rows], [10.0, 20.0])
        self.assertEqual(
            velocity_figure_name("multi", 1, 20.0),
            "multi_src1_sd20m.png",
        )

    def test_velocity_output_reset_requires_ownership_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            (workspace / "test").mkdir()
            output = workspace / "velocity"
            output.mkdir()
            (output / "user.txt").write_text("keep", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "ownership marker"):
                safe_reset_velocity_output(output, workspace)

    def test_velocity_output_reset_replaces_owned_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            (workspace / "test").mkdir()
            output = workspace / "velocity"
            output.mkdir()
            (output / ".ookc_velocity_comparison_output").write_text(
                "owned", encoding="utf-8"
            )
            (output / "old.txt").write_text("old", encoding="utf-8")

            safe_reset_velocity_output(output, workspace)

            self.assertFalse((output / "old.txt").exists())
            self.assertTrue(
                (output / ".ookc_velocity_comparison_output").is_file()
            )

    def test_velocity_manifest_counts_actual_cases_and_slices(self):
        records = [
            {"case": "a", "status": "success", "source_count": 2},
            {"case": "b", "status": "velocity_failed", "source_count": 0},
        ]
        rows = [{"case": "a"}, {"case": "a"}]

        manifest = velocity_manifest_payload(records, rows)

        self.assertEqual(manifest["case_count"], 2)
        self.assertEqual(manifest["successful_case_count"], 1)
        self.assertEqual(manifest["source_slice_count"], 2)
        self.assertEqual(manifest["successful_source_slice_count"], 2)


@unittest.skipUnless(
    os.environ.get("OPENOCEAN_KRAKENC_VELOCITY_RUNNER"),
    "OPENOCEAN_KRAKENC_VELOCITY_RUNNER is not configured",
)
class VelocityRunnerTests(unittest.TestCase):
    def test_runner_exports_pressure_horizontal_and_vertical_shade_files(self):
        runner = Path(os.environ["OPENOCEAN_KRAKENC_VELOCITY_RUNNER"])
        self.assertTrue(runner.is_file(), f"missing velocity runner: {runner}")
        workspace = Path(__file__).resolve().parents[2]
        source_dir = workspace / "test" / "krakenc_fortran_env"
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            for suffix in (".env", ".flp"):
                shutil.copy2(
                    source_dir / f"pekeris_isovelocity{suffix}",
                    run_dir / f"pekeris_isovelocity{suffix}",
                )
            output_root = run_dir / "velocity"
            completed = subprocess.run(
                [
                    str(runner),
                    str(run_dir / "pekeris_isovelocity.env"),
                    str(output_root),
                    "2",
                ],
                cwd=run_dir,
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            record = json.loads(completed.stdout)
            self.assertTrue(record["passes"])
            self.assertEqual(record["source_count"], 1)
            pressure = read_shade_file(Path(f"{output_root}_P.shd"))
            horizontal = read_shade_file(Path(f"{output_root}_H.shd"))
            vertical = read_shade_file(Path(f"{output_root}_V.shd"))
            assert_matching_shade_grids(pressure, horizontal, vertical)
            self.assertGreater(len(horizontal.pressure), 0)
            self.assertGreater(len(vertical.pressure), 0)


if __name__ == "__main__":
    unittest.main()
