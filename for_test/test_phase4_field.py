import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from for_test.oracle_config import resolve_oracle_binary
from for_test.shd_reader import read_shade_file
from for_test.fixture_paths import case_root


class Phase4FieldTests(unittest.TestCase):
    def test_munk_cpp_field_matches_fortran_field(self):
        project = Path(__file__).resolve().parents[1]
        workspace = project.parent
        binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", project / "build"))
        executable = binary_dir / "OpenOcean-Krakenc.exe"
        fortran_field = resolve_oracle_binary("field.exe")
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cpp_dir = root / "cpp"
            fortran_dir = root / "fortran"
            cpp_dir.mkdir()
            fortran_dir.mkdir()
            for directory in (cpp_dir, fortran_dir):
                for source in case_root("MunkKleaky").parent.glob("MunkKleaky.*"):
                    shutil.copy2(source, directory / source.name)

            mod_run = subprocess.run(
                [str(executable), "--mod", str(case_root("MunkKleaky").with_suffix(".env")),
                 str(cpp_dir / "MunkKleaky.mod")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(mod_run.returncode, 0, mod_run.stderr)
            shutil.copy2(cpp_dir / "MunkKleaky.mod", fortran_dir / "MunkKleaky.mod")

            cpp_run = subprocess.run(
                [str(executable), "--field", str(cpp_dir / "MunkKleaky")],
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(cpp_run.returncode, 0, cpp_run.stderr)
            fortran_run = subprocess.run(
                [str(fortran_field), "MunkKleaky"], cwd=fortran_dir,
                capture_output=True, text=True, timeout=300, check=False)
            self.assertEqual(fortran_run.returncode, 0, fortran_run.stderr)

            actual = read_shade_file(cpp_dir / "MunkKleaky.shd")
            expected = read_shade_file(fortran_dir / "MunkKleaky.shd")
            self.assertEqual(actual.source_depths, expected.source_depths)
            self.assertEqual(actual.receiver_depths, expected.receiver_depths)
            self.assertLessEqual(max((abs(left - right)
                                      for left, right in zip(actual.ranges_metres,
                                                             expected.ranges_metres)),
                                     default=0.0), 1.0e-8)
            self.assertEqual(len(actual.pressure), len(expected.pressure))
            error_energy = sum(abs(left - right) ** 2
                               for left, right in zip(actual.pressure, expected.pressure))
            reference_energy = sum(abs(value) ** 2 for value in expected.pressure)
            self.assertLessEqual((error_energy / reference_energy) ** 0.5, 2.0e-5)


if __name__ == "__main__":
    unittest.main()
