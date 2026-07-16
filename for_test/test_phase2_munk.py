import json
import math
import os
import subprocess
import unittest
from pathlib import Path

from for_test.mod_reader import read_first_mode_set, read_first_wavenumber_set
from for_test.reference_pipeline import run_reference_pipeline


class Phase2MunkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.reference_dir = cls.project / "build" / "reference" / "MunkKleaky"
        run_reference_pipeline(
            cls.workspace / "test" / "MunkKleaky",
            cls.workspace / "krakenFortran" / "build_mingw" / "out" / "krakenc.exe",
            cls.workspace / "krakenFortran" / "build_mingw" / "out" / "field.exe",
            cls.reference_dir,
        )

    def test_mod_reader_returns_complete_first_wavenumber_set(self):
        wavenumbers = read_first_wavenumber_set(
            self.reference_dir / "MunkKleaky.mod"
        )
        self.assertEqual(len(wavenumbers), 329)
        self.assertAlmostEqual(wavenumbers[0].real, 0.20936152338981628)
        self.assertAlmostEqual(wavenumbers[-1].real, 0.016605129465460777)

    def test_mod_reader_returns_depths_and_complex_modes(self):
        mode_set = read_first_mode_set(self.reference_dir / "MunkKleaky.mod")
        self.assertEqual(len(mode_set.depths), 1001)
        self.assertEqual(len(mode_set.modes), 329)
        self.assertEqual(len(mode_set.modes[0]), 1001)
        self.assertAlmostEqual(mode_set.depths[0], 0.0)
        self.assertAlmostEqual(mode_set.depths[-1], 5000.0)
        self.assertTrue(all(math.isfinite(value.real) and math.isfinite(value.imag)
                            for value in mode_set.modes[0]))

    def _run_cpp(self):
        executable = self.binary_dir / "OpenOceanKrakenc_acoustic_case_runner.exe"
        completed = subprocess.run(
            [str(executable), str(self.workspace / "test" / "MunkKleaky.env")],
            cwd=self.project,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(completed.stdout)

    def _run_cli(self):
        executable = self.binary_dir / "OpenOcean-Krakenc.exe"
        completed = subprocess.run(
            [
                str(executable),
                "--eigen",
                str(self.workspace / "test" / "MunkKleaky.env"),
            ],
            cwd=self.project,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(completed.stdout)

    def test_cpp_munk_matches_fortran_mod_and_is_deterministic(self):
        reference = read_first_wavenumber_set(
            self.reference_dir / "MunkKleaky.mod"
        )
        first = self._run_cpp()
        second = self._run_cpp()
        self.assertTrue(first["passes"], first)
        self.assertEqual(first, second)
        self.assertEqual(first["mode_count"], len(reference))
        actual = [complex(real, imag) for real, imag in first["wavenumbers"]]
        self.assertTrue(all(math.isfinite(value.real) and math.isfinite(value.imag)
                            for value in actual))
        absolute_errors = [abs(left - right) for left, right in zip(actual, reference)]
        relative_errors = [
            error / max(abs(expected), 1.0e-30)
            for error, expected in zip(absolute_errors, reference)
        ]
        self.assertLessEqual(max(absolute_errors), 1.0e-7)
        self.assertLessEqual(max(relative_errors), 1.0e-6)
        report = {
            "passes": True,
            "case": "MunkKleaky",
            "mode_count": len(reference),
            "mesh_sets_used": first["mesh_sets_used"],
            "max_absolute_wavenumber_error_1_per_m": max(absolute_errors),
            "max_relative_wavenumber_error": max(relative_errors),
            "deterministic_repeat": first == second,
        }
        (self.reference_dir / "phase2_comparison.json").write_text(
            json.dumps(report, indent=2), encoding="utf-8"
        )

    def test_public_cli_exposes_verified_acoustic_solver(self):
        self.assertEqual(self._run_cli(), self._run_cpp())


if __name__ == "__main__":
    unittest.main()
