import json
import math
import os
import subprocess
import unittest
from pathlib import Path

from for_test.mod_reader import read_first_wavenumber_set
from for_test.oracle_config import resolve_oracle_binary
from for_test.reference_pipeline import run_reference_pipeline
from for_test.fixture_paths import case_root


class Phase3ReflectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.krakenc = resolve_oracle_binary("krakenc.exe")
        cls.field = resolve_oracle_binary("field.exe")
        cls.runner = cls.binary_dir / "OpenOceanKrakenc_acoustic_case_runner.exe"

    def _solve(self, case):
        reference_dir = self.project / "build" / "reference" / case
        run_reference_pipeline(
            case_root(case),
            self.krakenc,
            self.field,
            reference_dir,
            run_field=False,
        )
        expected = read_first_wavenumber_set(reference_dir / f"{case}.mod")
        completed = subprocess.run(
            [str(self.runner), str(case_root(case).with_suffix(".env"))],
            cwd=self.project,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        actual = [complex(real, imag)
                  for real, imag in json.loads(completed.stdout)["wavenumbers"]]
        self.assertEqual(len(actual), len(expected))
        self.assertTrue(all(math.isfinite(value.real) and math.isfinite(value.imag)
                            for value in actual))
        errors = [abs(left - right) / max(abs(right), 1.0e-30)
                  for left, right in zip(actual, expected)]
        self.assertLessEqual(max(errors, default=0.0), 2.0e-5)
        return actual

    def _solve_cpp(self, case):
        completed = subprocess.run(
            [str(self.runner), str(case_root(case).with_suffix(".env"))],
            cwd=self.project,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return [complex(real, imag)
                for real, imag in json.loads(completed.stdout)["wavenumbers"]]

    def test_brc_and_irc_match_fortran_and_each_other(self):
        irc = self._solve("neggradC_irc")
        brc = self._solve_cpp("neggradC_brc")
        self.assertEqual(len(brc), len(irc))
        self.assertLessEqual(max((abs(left - right) for left, right in zip(brc, irc)),
                                 default=0.0), 1.0e-8)


if __name__ == "__main__":
    unittest.main()
