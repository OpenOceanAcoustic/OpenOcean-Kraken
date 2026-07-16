import json
import math
import os
import subprocess
import unittest
from pathlib import Path

from for_test.mod_reader import read_first_wavenumber_set
from for_test.reference_pipeline import run_reference_pipeline


class Phase3ElasticTests(unittest.TestCase):
    cases = (
        "elastic_fd_two_layer",
        "multilayer_elastic_stack",
        "multilayer_mud_sand",
    )

    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.krakenc = cls.workspace / "krakenFortran" / "build_mingw" / "out" / "krakenc.exe"
        cls.field = cls.workspace / "krakenFortran" / "build_mingw" / "out" / "field.exe"

    def test_elastic_wavenumbers_match_fortran(self):
        executable = self.binary_dir / "OpenOceanKrakenc_acoustic_case_runner.exe"
        for case in self.cases:
            with self.subTest(case=case):
                reference_dir = self.project / "build" / "reference" / case
                run_reference_pipeline(
                    self.workspace / "test" / case,
                    self.krakenc,
                    self.field,
                    reference_dir,
                )
                expected = read_first_wavenumber_set(reference_dir / f"{case}.mod")
                completed = subprocess.run(
                    [str(executable), str(self.workspace / "test" / f"{case}.env")],
                    cwd=self.project,
                    capture_output=True,
                    text=True,
                    timeout=240,
                    check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                result = json.loads(completed.stdout)
                actual = [complex(real, imag) for real, imag in result["wavenumbers"]]
                self.assertTrue(all(math.isfinite(value.real) and math.isfinite(value.imag)
                                    for value in actual))
                if case == "multilayer_elastic_stack":
                    # The Fortran coarse mesh terminates at 21 roots while its fine mesh
                    # continues, so only the first 19 Richardson pairings are complete.
                    # C++ bidirectional continuation intentionally returns the finite tail.
                    comparison_count = 19
                    self.assertGreaterEqual(len(actual), comparison_count)
                else:
                    comparison_count = len(expected)
                    self.assertEqual(len(actual), len(expected))
                errors = [abs(actual[index] - expected[index]) /
                          max(abs(expected[index]), 1.0e-30)
                          for index in range(comparison_count)]
                self.assertLessEqual(max(errors, default=0.0), 2.0e-5)


if __name__ == "__main__":
    unittest.main()
