import json
import math
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from for_test.mod_reader import read_first_mode_set, read_first_wavenumber_set
from for_test.oracle_config import resolve_oracle_binary
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
        cls.krakenc = resolve_oracle_binary("krakenc.exe")
        cls.field = resolve_oracle_binary("field.exe")

    def case_root(self, case):
        test_root = Path(os.environ["OPENOCEANKRAKENC_TEST_ROOT"])
        matches = list(test_root.rglob(f"{case}.env"))
        self.assertEqual(len(matches), 1, f"expected one nested ENV for {case}")
        return matches[0].with_suffix("")

    def test_elastic_wavenumbers_match_fortran(self):
        executable = self.binary_dir / "OpenOceanKrakenc_acoustic_case_runner.exe"
        for case in self.cases:
            with self.subTest(case=case):
                reference_dir = self.project / "build" / "reference" / case
                case_root = self.case_root(case)
                run_reference_pipeline(
                    case_root,
                    self.krakenc,
                    self.field,
                    reference_dir,
                )
                expected = read_first_wavenumber_set(reference_dir / f"{case}.mod")
                completed = subprocess.run(
                    [str(executable), str(case_root.with_suffix(".env")),
                     "--diagnostics"],
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
                diagnostics = result["mode_diagnostics"]
                self.assertEqual(len(diagnostics), len(actual))
                self.assertTrue(all(item["matched_mesh_sets"] == result["mesh_sets_used"]
                                    and math.isfinite(item["mesh_relative_spread"])
                                    and item["relative_correction"] <= 1.0e-9
                                    and "Fortran ordinal" in item["decision"]
                                    for item in diagnostics))
                expected_counts = {
                    "elastic_fd_two_layer": 20,
                    "multilayer_elastic_stack": 25,
                    "multilayer_mud_sand": 4,
                }
                self.assertEqual(len(expected), expected_counts[case])
                self.assertEqual(len(actual), expected_counts[case])
                comparison_indices = range(len(expected))
                errors = [abs(actual[index] - expected[index]) /
                          max(abs(expected[index]), 1.0e-30)
                          for index in comparison_indices]
                self.assertLessEqual(max(errors, default=0.0), 2.0e-5)

                with tempfile.TemporaryDirectory(prefix=f"ookc_{case}_mod_") as temp:
                    target_mod = Path(temp) / f"{case}.mod"
                    cli = self.binary_dir / "OpenOcean-Krakenc.exe"
                    generated = subprocess.run(
                        [str(cli), "--mod", str(case_root.with_suffix(".env")),
                         str(target_mod), "--threads", "2"],
                        cwd=self.project,
                        capture_output=True,
                        text=True,
                        timeout=300,
                        check=False,
                    )
                    self.assertEqual(generated.returncode, 0, generated.stderr)
                    actual_modes = read_first_mode_set(target_mod)
                    expected_modes = read_first_mode_set(
                        reference_dir / f"{case}.mod")
                    correlations = []
                    for left, right in zip(actual_modes.modes,
                                           expected_modes.modes):
                        left = np.asarray(left)
                        right = np.asarray(right)
                        correlations.append(
                            abs(np.vdot(left, right)) /
                            (np.linalg.norm(left) * np.linalg.norm(right))
                        )
                    self.assertGreaterEqual(min(correlations), 0.999)


if __name__ == "__main__":
    unittest.main()
