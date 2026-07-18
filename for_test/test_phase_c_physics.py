import json
import math
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from for_test.mod_reader import read_wavenumber_sets
from for_test.oracle_config import oracle_identity, resolve_oracle_binary
from for_test.fixture_paths import case_root


class PhaseCPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.oracle = resolve_oracle_binary("krakenc.exe")

    def run_fortran(self, env: Path, directory: Path) -> Path:
        target = directory / env.name
        shutil.copy2(env, target)
        root = target.stem
        completed = subprocess.run(
            [str(self.oracle), root], cwd=directory,
            capture_output=True, text=True, timeout=240, check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = directory / f"{root}.mod"
        self.assertTrue(result.is_file() and result.stat().st_size > 0)
        return result

    @staticmethod
    def relative_errors(actual, expected):
        return [abs(left - right) / max(abs(right), 1.0e-30)
                for left, right in zip(actual, expected)]

    def test_scholte_floor_matches_frozen_oracle(self):
        env = self.project / "for_test" / "fixtures" / "g3_scholte.env"
        with tempfile.TemporaryDirectory(prefix="ookc_phase_c_scholte_") as raw:
            directory = Path(raw)
            expected_sets = read_wavenumber_sets(
                self.run_fortran(env, directory))
            completed = subprocess.run(
                [str(self.binary_dir / "OpenOceanKrakenc_acoustic_case_runner.exe"),
                 str(env)], cwd=self.project, capture_output=True, text=True,
                timeout=240, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(completed.stdout)
            actual = [complex(*value) for value in payload["wavenumbers"]]
            expected = expected_sets[0]
            self.assertEqual(len(actual), len(expected))
            self.assertEqual(len(expected), 7)
            errors = self.relative_errors(actual, expected)
            self.assertLessEqual(max(errors, default=0.0), 2.0e-5)
            minimum_speed = min(
                2.0 * math.pi * 50.0 / value.real for value in actual)
            self.assertAlmostEqual(minimum_speed, 446.7375608441422, delta=0.05)

    def test_solve3_profiles_match_frozen_oracle(self):
        env = case_root("solve3_mode_gain").with_suffix(".env")
        with tempfile.TemporaryDirectory(prefix="ookc_phase_c_solve3_") as raw:
            directory = Path(raw)
            expected = read_wavenumber_sets(self.run_fortran(env, directory))
            output = directory / "target.mod"
            completed = subprocess.run(
                [str(self.binary_dir / "OpenOcean-Krakenc.exe"), "--mod",
                 str(env), str(output), "--threads", "2"],
                cwd=self.project, capture_output=True, text=True,
                timeout=600, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(completed.stdout)
            actual = read_wavenumber_sets(output)
            self.assertEqual(payload["profile_count"], len(expected))
            self.assertEqual(len(actual), len(expected))
            maximum_error = 0.0
            for actual_set, expected_set in zip(actual, expected):
                self.assertGreaterEqual(len(actual_set), len(expected_set))
                maximum_error = max(
                    maximum_error,
                    max(self.relative_errors(
                        actual_set[:len(expected_set)], expected_set), default=0.0))
            self.assertLessEqual(maximum_error, 2.0e-5)

            evidence = {
                "passes": True,
                "oracle": oracle_identity(self.oracle),
                "profile_count": len(actual),
                "target_mode_counts": [len(values) for values in actual],
                "oracle_mode_counts": [len(values) for values in expected],
                "maximum_relative_wavenumber_error": maximum_error,
            }
            evidence_path = (self.binary_dir / "reference" /
                             "phase_c_solve3_differential.json")
            evidence_path.parent.mkdir(parents=True, exist_ok=True)
            evidence_path.write_text(
                json.dumps(evidence, indent=2, ensure_ascii=False),
                encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
