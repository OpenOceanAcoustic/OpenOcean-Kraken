import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from for_test.mod_reader import read_first_mode_set
from for_test.oracle_config import resolve_oracle_binary
from for_test.reference_pipeline import run_reference_pipeline
from for_test.fixture_paths import case_root


class Phase3MunkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[1]
        cls.workspace = cls.project.parent
        cls.binary_dir = Path(os.environ.get(
            "OPENOCEAN_KRAKENC_BINARY_DIR", cls.project / "build"))
        cls.reference_dir = (
            cls.project / "build" / "reference" / "phase3_munk" / "MunkKleaky"
        )
        run_reference_pipeline(
            case_root("MunkKleaky"),
            resolve_oracle_binary("krakenc.exe"),
            resolve_oracle_binary("field.exe"),
            cls.reference_dir,
        )
        cls.reference = read_first_mode_set(cls.reference_dir / "MunkKleaky.mod")

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as stream:
            cls.output_path = Path(stream.name)
        executable = cls.binary_dir / "OpenOceanKrakenc_phase3_munk_runner.exe"
        completed = subprocess.run(
            [str(executable), str(case_root("MunkKleaky").with_suffix(".env")),
             str(cls.output_path)],
            cwd=cls.project,
            capture_output=True,
            text=True,
            timeout=240,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr)
        cls.actual = json.loads(cls.output_path.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "output_path"):
            cls.output_path.unlink(missing_ok=True)

    def test_all_modes_are_finite_and_match_mod_shapes(self):
        self.assertEqual(self.actual["mode_count"], 329)
        self.assertEqual(self.actual["depths"], self.reference.depths)
        self.assertEqual(len(self.actual["modes"]), 329)
        worst_relative_l2 = 0.0
        for actual_pairs, expected in zip(self.actual["modes"], self.reference.modes):
            actual = [complex(real, imag) for real, imag in actual_pairs]
            self.assertTrue(all(math.isfinite(value.real) and math.isfinite(value.imag)
                                for value in actual))
            inner = sum(value.conjugate() * target
                        for value, target in zip(actual, expected))
            phase = inner / abs(inner)
            error = math.sqrt(sum(abs(phase * value - target) ** 2
                                  for value, target in zip(actual, expected)))
            scale = max(math.sqrt(sum(abs(value) ** 2 for value in expected)), 1.0e-30)
            worst_relative_l2 = max(worst_relative_l2, error / scale)
        self.assertLessEqual(worst_relative_l2, 2.0e-3)

    def test_group_velocities_match_fortran_report(self):
        report = (self.reference_dir / "MunkKleaky.prt").read_text(
            encoding="utf-8", errors="replace")
        expected = {}
        pattern = re.compile(
            r"^\s*(\d+)\s+[-+0-9.E]+\s+[-+0-9.E]+\s+[-+0-9.E]+\s+([-+0-9.E]+)\s*$"
        )
        for line in report.splitlines():
            match = pattern.match(line)
            if match:
                expected[int(match.group(1)) - 1] = float(match.group(2))
        self.assertGreaterEqual(len(expected), 30)
        velocities = self.actual["group_velocities"]
        for index, target in expected.items():
            self.assertTrue(math.isfinite(velocities[index]))
            self.assertLessEqual(abs(velocities[index] - target) / abs(target), 1.0e-4)

    def test_cpp_mod_is_fortran_field_compatible(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for source in case_root("MunkKleaky").parent.glob("MunkKleaky.*"):
                shutil.copy2(source, directory / source.name)
            output_mod = directory / "MunkKleaky.mod"
            executable = self.binary_dir / "OpenOcean-Krakenc.exe"
            completed = subprocess.run(
                [str(executable), "--mod",
                 str(case_root("MunkKleaky").with_suffix(".env")),
                 str(output_mod)],
                cwd=self.project,
                capture_output=True,
                text=True,
                timeout=300,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            generated = read_first_mode_set(output_mod)
            self.assertEqual(generated.depths, self.reference.depths)
            self.assertEqual(len(generated.modes), len(self.reference.modes))
            self.assertLessEqual(
                max((abs(left - right) / max(abs(right), 1.0e-30)
                     for left, right in zip(generated.wavenumbers,
                                            self.reference.wavenumbers)),
                    default=0.0),
                1.0e-6,
            )
            worst_mode_error = 0.0
            for actual, expected in zip(generated.modes, self.reference.modes):
                inner = sum(value.conjugate() * target
                            for value, target in zip(actual, expected))
                phase = inner / abs(inner)
                error = math.sqrt(sum(abs(phase * value - target) ** 2
                                      for value, target in zip(actual, expected)))
                scale = max(math.sqrt(sum(abs(value) ** 2 for value in expected)),
                            1.0e-30)
                worst_mode_error = max(worst_mode_error, error / scale)
            self.assertLessEqual(worst_mode_error, 2.0e-3)

            field = resolve_oracle_binary("field.exe")
            field_run = subprocess.run(
                [str(field), "MunkKleaky"],
                cwd=directory,
                capture_output=True,
                text=True,
                timeout=180,
                check=False,
            )
            self.assertEqual(field_run.returncode, 0, field_run.stderr)
            self.assertTrue((directory / "MunkKleaky.shd").stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
