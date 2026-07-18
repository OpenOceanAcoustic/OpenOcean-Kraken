import json
import math
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

from for_test.mod_reader import read_first_wavenumber_set
from for_test.oracle_config import oracle_identity, resolve_oracle_binary
from for_test.fixture_paths import case_root


class PhaseBPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.runner = cls.binary_dir / "OpenOceanKrakenc_acoustic_case_runner.exe"
        cls.krakenc = resolve_oracle_binary("krakenc.exe")
        cls.base_text = case_root("MunkKleaky").with_suffix(".env").read_text(
            encoding="utf-8"
        )
        cls.smooth = cls._run_cpp(case_root("MunkKleaky").with_suffix(".env"))
        cls.results = []

    @classmethod
    def tearDownClass(cls):
        report_dir = cls.project / "build" / "reference"
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "phase_b_physics.json").write_text(
            json.dumps({
                "passes": True,
                "oracle": oracle_identity(cls.krakenc),
                "cases": cls.results,
            }, indent=2),
            encoding="utf-8",
        )

    @classmethod
    def _run_cpp(cls, env_path):
        completed = subprocess.run(
            [str(cls.runner), str(env_path)],
            cwd=cls.project,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
        if completed.returncode != 0:
            raise AssertionError(completed.stderr)
        return [complex(real, imaginary)
                for real, imaginary in json.loads(completed.stdout)["wavenumbers"]]

    def _check_variant(self, name, text, minimum_effect):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            env_path = directory / f"{name}.env"
            env_path.write_text(text, encoding="utf-8")
            completed = subprocess.run(
                [str(self.krakenc), name],
                cwd=directory,
                capture_output=True,
                text=True,
                timeout=300,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            expected = read_first_wavenumber_set(directory / f"{name}.mod")
            actual = self._run_cpp(env_path)

        self.assertEqual(len(actual), len(expected))
        self.assertEqual(len(actual), len(self.smooth))
        self.assertTrue(all(math.isfinite(value.real) and math.isfinite(value.imag)
                            for value in actual))
        relative_errors = [
            abs(left - right) / max(abs(right), 1.0e-30)
            for left, right in zip(actual, expected)
        ]
        effect = max(abs(left - right)
                     for left, right in zip(actual, self.smooth))
        self.assertGreater(effect, minimum_effect)
        self.assertLessEqual(max(relative_errors, default=0.0), 2.0e-5)
        self.results.append({
            "case": name,
            "mode_count": len(actual),
            "max_relative_wavenumber_error": max(relative_errors, default=0.0),
            "max_effect_from_smooth_1_per_m": effect,
        })

    def test_thorp_volume_absorption_matches_fortran(self):
        self._check_variant(
            "MunkThorp",
            self.base_text.replace("'NVW .'", "'NVWT.'", 1),
            1.0e-12,
        )

    def test_francois_garrison_volume_absorption_matches_fortran(self):
        self._check_variant(
            "MunkFrancGarr",
            self.base_text.replace(
                "'NVW .'", "'NVWF.'\n10.0 35.0 8.0 1000.0", 1
            ),
            1.0e-12,
        )

    def test_biological_volume_absorption_matches_fortran(self):
        self._check_variant(
            "MunkBiological",
            self.base_text.replace(
                "'NVW .'",
                "'NVWB.'\n1\n0.0 5000.0 100.0 5.0 2.0",
                1,
            ),
            1.0e-12,
        )

    def test_kuperman_ingenito_roughness_matches_fortran(self):
        self._check_variant(
            "MunkRough",
            self.base_text.replace("'A' 0.0", "'A' 0.05", 1),
            1.0e-12,
        )


if __name__ == "__main__":
    unittest.main()
