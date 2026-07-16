import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from for_test.reference_pipeline import run_reference_pipeline
from for_test.shd_reader import read_shade_file


class Phase4RangeDependentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.executable = cls.binary_dir / "OpenOcean-Krakenc.exe"
        cls.fortran_field = (
            cls.workspace / "krakenFortran" / "build_mingw" / "out" / "field.exe"
        )
        cls.reference_dir = cls.project / "build" / "reference" / "stepK_rd"
        run_reference_pipeline(
            cls.workspace / "test" / "stepK_rd",
            cls.workspace / "krakenFortran" / "build_mingw" / "out" / "krakenc.exe",
            cls.fortran_field,
            cls.reference_dir,
        )
        cls.wedge_reference_dir = cls.project / "build" / "reference" / "wedge"
        run_reference_pipeline(
            cls.workspace / "test" / "wedge",
            cls.workspace / "krakenFortran" / "build_mingw" / "out" / "krakenc.exe",
            cls.fortran_field,
            cls.wedge_reference_dir,
        )

    def test_stepk_adiabatic_field_matches_fortran(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cpp_dir = root / "cpp"
            fortran_dir = root / "fortran"
            cpp_dir.mkdir()
            fortran_dir.mkdir()
            flp = (self.workspace / "test" / "stepK_rd.flp").read_text(
                encoding="utf-8"
            ).replace("'RC'", "'RA'", 1)
            for directory in (cpp_dir, fortran_dir):
                shutil.copy2(self.reference_dir / "stepK_rd.mod",
                             directory / "stepK_rd.mod")
                (directory / "stepK_rd.flp").write_text(flp, encoding="utf-8")

            fortran_run = subprocess.run(
                [str(self.fortran_field), "stepK_rd"], cwd=fortran_dir,
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(fortran_run.returncode, 0, fortran_run.stderr)
            cpp_run = subprocess.run(
                [str(self.executable), "--field", str(cpp_dir / "stepK_rd")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(cpp_run.returncode, 0, cpp_run.stderr)

            actual = read_shade_file(cpp_dir / "stepK_rd.shd")
            expected = read_shade_file(fortran_dir / "stepK_rd.shd")
            self.assertEqual(actual.source_depths, expected.source_depths)
            self.assertEqual(actual.receiver_depths, expected.receiver_depths)
            self.assertLessEqual(max((abs(left - right)
                                      for left, right in zip(actual.ranges_metres,
                                                             expected.ranges_metres)),
                                     default=0.0), 1.0e-8)
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 2.0e-4)

    def test_stepk_adiabatic_beam_pattern_matches_fortran(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cpp_dir = root / "cpp"
            fortran_dir = root / "fortran"
            cpp_dir.mkdir()
            fortran_dir.mkdir()
            flp = (self.workspace / "test" / "stepK_rd.flp").read_text(
                encoding="utf-8").replace("'RC'", "'RA*'", 1)
            sbp = "2\n-90.0 -6.020599913\n90.0 -6.020599913\n"
            for directory in (cpp_dir, fortran_dir):
                shutil.copy2(self.reference_dir / "stepK_rd.mod",
                             directory / "stepK_rd.mod")
                (directory / "stepK_rd.flp").write_text(flp, encoding="utf-8")
                (directory / "stepK_rd.sbp").write_text(sbp, encoding="utf-8")
            fortran_run = subprocess.run(
                [str(self.fortran_field), "stepK_rd"], cwd=fortran_dir,
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(fortran_run.returncode, 0, fortran_run.stderr)
            cpp_run = subprocess.run(
                [str(self.executable), "--field", str(cpp_dir / "stepK_rd")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(cpp_run.returncode, 0, cpp_run.stderr)
            actual = read_shade_file(cpp_dir / "stepK_rd.shd")
            expected = read_shade_file(fortran_dir / "stepK_rd.shd")
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 2.0e-4)

    def test_stepk_scaled_coupled_phase_matches_fortran(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cpp_dir = root / "cpp"
            fortran_dir = root / "fortran"
            cpp_dir.mkdir()
            fortran_dir.mkdir()
            flp = (self.workspace / "test" / "stepK_rd.flp").read_text(
                encoding="utf-8").replace("'RC'", "'SC'", 1)
            for directory in (cpp_dir, fortran_dir):
                shutil.copy2(self.reference_dir / "stepK_rd.mod",
                             directory / "stepK_rd.mod")
                (directory / "stepK_rd.flp").write_text(flp, encoding="utf-8")
            fortran_run = subprocess.run(
                [str(self.fortran_field), "stepK_rd"], cwd=fortran_dir,
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(fortran_run.returncode, 0, fortran_run.stderr)
            cpp_run = subprocess.run(
                [str(self.executable), "--field", str(cpp_dir / "stepK_rd")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(cpp_run.returncode, 0, cpp_run.stderr)
            actual = read_shade_file(cpp_dir / "stepK_rd.shd")
            expected = read_shade_file(fortran_dir / "stepK_rd.shd")
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 2.0e-4)

    def test_stepk_coupled_field_matches_fortran(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cpp_dir = root / "cpp"
            cpp_dir.mkdir()
            shutil.copy2(self.reference_dir / "stepK_rd.mod",
                         cpp_dir / "stepK_rd.mod")
            shutil.copy2(self.workspace / "test" / "stepK_rd.flp",
                         cpp_dir / "stepK_rd.flp")
            cpp_run = subprocess.run(
                [str(self.executable), "--field", str(cpp_dir / "stepK_rd")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(cpp_run.returncode, 0, cpp_run.stderr)

            actual = read_shade_file(cpp_dir / "stepK_rd.shd")
            expected = read_shade_file(self.reference_dir / "stepK_rd.shd")
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 2.0e-4)

    def test_stepk_cpp_krakenc_and_field_end_to_end(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            case_root = root / "stepK_rd"
            shutil.copy2(self.workspace / "test" / "stepK_rd.flp",
                         case_root.with_suffix(".flp"))
            mod_run = subprocess.run(
                [str(self.executable), "--mod",
                 str(self.workspace / "test" / "stepK_rd.env"),
                 str(case_root.with_suffix(".mod"))],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(mod_run.returncode, 0, mod_run.stderr)
            field_run = subprocess.run(
                [str(self.executable), "--field", str(case_root)],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(field_run.returncode, 0, field_run.stderr)
            actual = read_shade_file(case_root.with_suffix(".shd"))
            expected = read_shade_file(self.reference_dir / "stepK_rd.shd")
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 5.0e-3)

    def _compare_wedge(self, propagation_type):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cpp_dir = root / "cpp"
            fortran_dir = root / "fortran"
            cpp_dir.mkdir()
            fortran_dir.mkdir()
            original = (self.workspace / "test" / "wedge.flp").read_text(
                encoding="utf-8")
            flp = original.replace("'RC'", f"'R{propagation_type}'", 1)
            for directory in (cpp_dir, fortran_dir):
                shutil.copy2(self.wedge_reference_dir / "wedge.mod",
                             directory / "wedge.mod")
                (directory / "wedge.flp").write_text(flp, encoding="utf-8")

            fortran_run = subprocess.run(
                [str(self.fortran_field), "wedge"], cwd=fortran_dir,
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(fortran_run.returncode, 0, fortran_run.stderr)
            cpp_run = subprocess.run(
                [str(self.executable), "--field", str(cpp_dir / "wedge")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(cpp_run.returncode, 0, cpp_run.stderr)
            actual = read_shade_file(cpp_dir / "wedge.shd")
            expected = read_shade_file(fortran_dir / "wedge.shd")
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 2.0e-4)

    def test_wedge_coupled_field_matches_fortran(self):
        self._compare_wedge("C")

    def test_wedge_cpp_krakenc_and_field_end_to_end(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            case_root = root / "wedge"
            shutil.copy2(self.workspace / "test" / "wedge.flp",
                         case_root.with_suffix(".flp"))
            mod_run = subprocess.run(
                [str(self.executable), "--mod",
                 str(self.workspace / "test" / "wedge.env"),
                 str(case_root.with_suffix(".mod"))],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(mod_run.returncode, 0, mod_run.stderr)
            field_run = subprocess.run(
                [str(self.executable), "--field", str(case_root)],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(field_run.returncode, 0, field_run.stderr)
            actual = read_shade_file(case_root.with_suffix(".shd"))
            expected = read_shade_file(self.wedge_reference_dir / "wedge.shd")
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure,
                                                      expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 5.0e-3)



if __name__ == "__main__":
    unittest.main()
