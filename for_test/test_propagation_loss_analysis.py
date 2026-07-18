import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from for_test.propagation_loss_analysis import (
    PressureCube,
    assert_matching_grids,
    comparison_color_limits,
    compute_metrics,
    compute_tl,
    create_overview,
    failure_figures,
    nearest_receiver_index,
    plot_comparison,
    shade_to_cube,
    spatial_diagnostics,
    valid_tl_mask,
)
from for_test.plot_propagation_loss_comparison import (
    CaseInput,
    OUTPUT_MARKER,
    discover_cases,
    expected_source_depths,
    process_cases,
    run_case,
    run_logged_command,
    safe_reset_output,
    stage_case_inputs,
    stage_fortran_case_inputs,
    write_analysis_report,
    write_manifest,
    write_metrics_csv,
)
from for_test.shd_reader import ShadeFile


def synthetic_cube() -> PressureCube:
    pressure = np.array(
        [
            [1.0 + 0.0j, 0.8 + 0.1j, 0.5 + 0.2j, 0.3 + 0.1j],
            [0.9 + 0.0j, 0.6 + 0.1j, 0.4 + 0.2j, 0.2 + 0.1j],
            [0.7 + 0.0j, 0.5 + 0.1j, 0.3 + 0.2j, 0.1 + 0.1j],
        ],
        dtype=np.complex128,
    )
    return PressureCube(
        title="synthetic",
        frequency=50.0,
        source_depths=np.array([25.0]),
        receiver_depths=np.array([0.0, 25.0, 50.0]),
        ranges_metres=np.array([100.0, 1000.0, 2000.0, 3000.0]),
        pressure=pressure[np.newaxis, :, :],
    )


def workspace_test_directory() -> Path:
    for ancestor in Path(__file__).resolve().parents:
        candidate = ancestor / "test"
        if any(candidate.rglob("MunkKleaky.env")):
            return candidate
    raise RuntimeError("unable to locate workspace test directory")


class PropagationLossAnalysisTests(unittest.TestCase):
    def test_shade_to_cube_preserves_source_receiver_range_order(self):
        shade = ShadeFile(
            "synthetic",
            50.0,
            [25.0, 250.0],
            [10.0, 20.0],
            [0.0, 1000.0, 2000.0],
            [complex(value, 0.0) for value in range(1, 13)],
        )

        cube = shade_to_cube(shade)

        self.assertEqual(cube.pressure.shape, (2, 2, 3))
        np.testing.assert_array_equal(cube.pressure[1, 0], [7, 8, 9])

    def test_compute_tl_applies_pressure_floor(self):
        actual = compute_tl(np.array([1.0, 0.1, 0.0]))

        np.testing.assert_allclose(actual, [0.0, 20.0, 240.0])

    def test_valid_tl_mask_uses_joint_peak(self):
        ookc = np.array([1.0, 1.0e-9, 0.0])
        krakenc = np.array([0.5, 0.0, 1.0e-7])

        np.testing.assert_array_equal(
            valid_tl_mask(ookc, krakenc), [True, False, True]
        )

    def test_compute_metrics_reports_complex_l2_and_db_statistics(self):
        reference = np.array([1.0 + 0.0j, 0.1 + 0.0j])
        actual = np.array([1.0 + 0.0j, 0.01 + 0.0j])

        metrics = compute_metrics(actual, reference)

        self.assertAlmostEqual(
            metrics.complex_relative_l2, 0.09 / np.sqrt(1.01)
        )
        self.assertAlmostEqual(metrics.tl_mae_db, 10.0)
        self.assertAlmostEqual(metrics.tl_rmse_db, np.sqrt(200.0))
        self.assertAlmostEqual(metrics.tl_p95_db, 19.0)
        self.assertAlmostEqual(metrics.tl_max_db, 20.0)
        self.assertEqual(metrics.valid_fraction, 1.0)

    def test_spatial_diagnostics_separate_zero_range_and_boundaries(self):
        reference = np.ones((3, 3), dtype=np.complex128)
        actual = reference.copy()
        actual[:, 0] = 2.0

        diagnostics = spatial_diagnostics(
            actual, reference, np.array([0.0, 1000.0, 2000.0])
        )

        self.assertEqual(diagnostics.positive_range_relative_l2, 0.0)
        self.assertEqual(
            diagnostics.interior_positive_range_relative_l2, 0.0
        )
        self.assertAlmostEqual(
            diagnostics.zero_range_reference_energy_fraction,
            1.0 / 3.0,
        )
        self.assertAlmostEqual(
            diagnostics.boundary_error_fraction, 2.0 / 3.0
        )

    def test_assert_matching_grids_rejects_receiver_depth_difference(self):
        left = ShadeFile("a", 50.0, [25.0], [10.0], [1000.0], [1 + 0j])
        right = ShadeFile("b", 50.0, [25.0], [11.0], [1000.0], [1 + 0j])

        with self.assertRaisesRegex(ValueError, "receiver-depth grid mismatch"):
            assert_matching_grids(shade_to_cube(left), shade_to_cube(right))

    def test_comparison_color_limits_keep_zero_difference_visible(self):
        pressure = np.array([1.0, 0.1])

        tl_min, tl_max, difference_half_width = comparison_color_limits(
            pressure, pressure
        )

        self.assertLess(tl_min, tl_max)
        self.assertEqual(difference_half_width, 0.1)

    def test_nearest_receiver_index_selects_source_depth(self):
        actual = nearest_receiver_index(np.array([0.0, 20.0, 40.0]), 31.0)

        self.assertEqual(actual, 2)

    def test_plot_comparison_and_overview_write_nonempty_png_files(self):
        cube = synthetic_cube()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            comparison = directory / "comparison.png"
            overview = directory / "overview.png"

            metrics = plot_comparison("synthetic", cube, cube, 0, comparison)
            create_overview([comparison], overview)

            self.assertEqual(metrics.complex_relative_l2, 0.0)
            self.assertGreater(comparison.stat().st_size, 10_000)
            self.assertGreater(overview.stat().st_size, 10_000)

    def test_failure_figures_write_one_png_per_source(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = failure_figures(
                "failed_case", 2, "solver failed", Path(directory)
            )

            self.assertEqual(len(paths), 2)
            self.assertTrue(all(path.stat().st_size > 10_000 for path in paths))


class BatchOrchestrationTests(unittest.TestCase):
    def test_cli_can_be_invoked_by_script_path(self):
        project = Path(__file__).resolve().parents[1]

        completed = subprocess.run(
            [
                sys.executable,
                str(project / "for_test" / "plot_propagation_loss_comparison.py"),
                "--help",
            ],
            cwd=project,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--test-dir", completed.stdout)

    def test_discover_cases_finds_all_workspace_env_files(self):
        cases = discover_cases(workspace_test_directory())

        self.assertEqual(len(cases), 11)
        self.assertEqual(sum(case.name == "MunkKleaky" for case in cases), 1)
        self.assertTrue(all(case.env.parent == case.root.parent for case in cases))
        self.assertEqual(
            [case.env.relative_to(workspace_test_directory()).as_posix()
             for case in cases],
            sorted(
                [case.env.relative_to(workspace_test_directory()).as_posix()
                 for case in cases],
                key=str.casefold,
            ),
        )

    def test_discover_cases_rejects_duplicate_stems_in_nested_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            test_dir = Path(directory)
            for group in ("a", "b"):
                nested = test_dir / group
                nested.mkdir()
                (nested / "duplicate.env").write_text("env", encoding="utf-8")
                (nested / "duplicate.flp").write_text("flp", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "duplicate ENV stem"):
                discover_cases(test_dir)

    def test_stage_inputs_maps_missing_neggradk_reflection_file(self):
        case = next(
            case
            for case in discover_cases(workspace_test_directory())
            if case.name == "neggradK_brc"
        )
        with tempfile.TemporaryDirectory() as directory:
            record = stage_case_inputs(case, Path(directory))

            self.assertTrue(
                (Path(directory) / "neggradK_brc.brc").is_file()
            )
            self.assertEqual(
                record["neggradK_brc.brc"]["source"],
                "neggradC_brc.brc",
            )

    def test_fortran_staging_uses_equivalent_irc_for_brc_case(self):
        case = next(
            case
            for case in discover_cases(workspace_test_directory())
            if case.name == "neggradC_brc"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)

            record = stage_fortran_case_inputs(case, directory)

            staged_env = (directory / "neggradC_brc.env").read_text(
                encoding="utf-8"
            )
            original_env = case.env.read_text(encoding="utf-8")
            self.assertIn("'P'", staged_env)
            self.assertIn("'F'", original_env)
            self.assertTrue((directory / "neggradC_brc.irc").is_file())
            self.assertEqual(record["_adaptation"]["kind"], "brc_to_irc")

    def test_safe_reset_output_rejects_test_directory(self):
        test_directory = workspace_test_directory()

        with self.assertRaisesRegex(ValueError, "unsafe output directory"):
            safe_reset_output(test_directory, test_directory.parent)

    def test_safe_reset_output_rejects_unmarked_existing_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory)
            unrelated = workspace / "unrelated"
            unrelated.mkdir()
            protected = unrelated / "keep.txt"
            protected.write_text("keep", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "ownership marker"):
                safe_reset_output(unrelated, workspace)

            self.assertEqual(protected.read_text(encoding="utf-8"), "keep")

    def test_safe_reset_output_replaces_only_marked_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory)
            output = workspace / "generated"
            output.mkdir()
            (output / OUTPUT_MARKER).write_text("owned", encoding="utf-8")
            old = output / "old.txt"
            old.write_text("old", encoding="utf-8")

            safe_reset_output(output, workspace)

            self.assertFalse(old.exists())
            self.assertTrue((output / OUTPUT_MARKER).is_file())

    def test_process_cases_records_failure_and_continues(self):
        def runner(name):
            if name == "first":
                raise RuntimeError("intentional")
            return {"case": name, "status": "success"}

        records = process_cases(["first", "second"], runner)

        self.assertEqual(
            [record["status"] for record in records],
            ["analysis_failed", "success"],
        )

    def test_run_logged_command_retains_stdout_and_stderr(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            record = run_logged_command(
                [
                    "python",
                    "-c",
                    "import sys; print('out'); print('err', file=sys.stderr)",
                ],
                directory,
                "sample",
            )

            self.assertEqual(record["returncode"], 0)
            self.assertIn(
                b"out", Path(record["stdout"]["path"]).read_bytes()
            )
            self.assertIn(
                b"err", Path(record["stderr"]["path"]).read_bytes()
            )

    def test_expected_source_depths_include_both_munk_sources(self):
        self.assertEqual(expected_source_depths("MunkKleaky"), [25.0, 250.0])
        self.assertEqual(expected_source_depths("wedge"), [None])

    def test_report_writers_emit_required_artifacts(self):
        rows = [
            {
                "case": "synthetic",
                "frequency_hz": 50.0,
                "source_depth_m": 25.0,
                "complex_relative_l2": 0.0,
                "positive_range_relative_l2": 0.0,
                "interior_positive_range_relative_l2": 0.0,
                "zero_range_reference_energy_fraction": 0.0,
                "boundary_error_fraction": 0.0,
                "tl_mae_db": 0.0,
                "tl_rmse_db": 0.0,
                "tl_p95_db": 0.0,
                "tl_max_db": 0.0,
                "valid_fraction": 1.0,
                "max_difference_range_km": 1.0,
                "max_difference_depth_m": 20.0,
                "figure": "figures/synthetic_sd25m.png",
            }
        ]
        cases = [
            {
                "case": "synthetic",
                "status": "success",
                "source_count": 1,
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)

            write_metrics_csv(rows, output / "metrics.csv")
            write_manifest(cases, rows, output / "manifest.json")
            write_analysis_report(cases, rows, output / "analysis.md")

            self.assertIn(
                "complex_relative_l2",
                (output / "metrics.csv").read_text(encoding="utf-8"),
            )
            manifest = json.loads(
                (output / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["case_count"], 1)
            self.assertEqual(manifest["source_task_count"], 1)
            report = (output / "analysis.md").read_text(encoding="utf-8")
            self.assertIn("数值一致性", report)
            self.assertIn("声传播物理特征", report)
            self.assertIn("总体结论", report)
            self.assertIn("overview.png", report)
            self.assertIn("synthetic", report)
            self.assertIn("等价刚性 IRC", report)

    def test_run_case_generates_metrics_for_a_successful_source_field(self):
        shade = ShadeFile(
            "synthetic",
            50.0,
            [25.0],
            [0.0, 25.0, 50.0],
            [100.0, 1000.0, 2000.0, 3000.0],
            list(synthetic_cube().pressure.ravel()),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            input_directory = directory / "inputs"
            input_directory.mkdir()
            env = input_directory / "synthetic.env"
            flp = input_directory / "synthetic.flp"
            env.write_text("synthetic", encoding="utf-8")
            flp.write_text("synthetic", encoding="utf-8")
            case = CaseInput("synthetic", env.with_suffix(""), env, flp)

            def fake_run(command, cwd, label, timeout_seconds=600):
                suffix = ".mod" if label in {"ookc_mod", "krakenc"} else ".shd"
                (Path(cwd) / f"synthetic{suffix}").write_bytes(b"artifact")
                return {
                    "command": [str(value) for value in command],
                    "cwd": str(cwd),
                    "returncode": 0,
                    "timed_out": False,
                    "seconds": 0.01,
                }

            with mock.patch(
                "for_test.plot_propagation_loss_comparison.run_logged_command",
                side_effect=fake_run,
            ), mock.patch(
                "for_test.plot_propagation_loss_comparison.read_shade_file",
                side_effect=[shade, shade],
            ):
                record = run_case(
                    case,
                    Path("ookc.exe"),
                    Path("krakenc.exe"),
                    Path("field.exe"),
                    directory / "output",
                    1,
                )

            self.assertEqual(record["status"], "success")
            self.assertEqual(record["source_count"], 1)
            self.assertEqual(len(record["metrics"]), 1)
            self.assertEqual(record["metrics"][0]["complex_relative_l2"], 0.0)
            self.assertTrue(Path(record["figures"][0]["path"]).is_file())

    def test_run_case_creates_placeholder_after_ookc_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            env = directory / "failed.env"
            flp = directory / "failed.flp"
            env.write_text("failed", encoding="utf-8")
            flp.write_text("failed", encoding="utf-8")
            case = CaseInput("failed", env.with_suffix(""), env, flp)
            failed_run = {
                "command": ["ookc"],
                "cwd": str(directory),
                "returncode": 1,
                "timed_out": False,
                "seconds": 0.01,
            }

            with mock.patch(
                "for_test.plot_propagation_loss_comparison.run_logged_command",
                return_value=failed_run,
            ):
                record = run_case(
                    case,
                    Path("ookc.exe"),
                    Path("krakenc.exe"),
                    Path("field.exe"),
                    directory / "output",
                    1,
                )

            self.assertEqual(record["status"], "ookc_failed")
            self.assertEqual(len(record["figures"]), 1)
            self.assertTrue(Path(record["figures"][0]["path"]).is_file())


if __name__ == "__main__":
    unittest.main()
