import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from for_test.shd_reader import read_shade_file


CASES = (
    "MunkKleaky",
    "elastic_fd_two_layer",
    "multilayer_elastic_stack",
    "multilayer_mud_sand",
    "neggradC_brc",
    "neggradC_irc",
    "stepK_rd",
    "wedge",
)

INPUT_SUFFIXES = {".env", ".flp", ".sbp", ".brc", ".irc", ".trc"}


class EnvJsonEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        binary_dir = Path(os.environ["OPENOCEAN_KRAKENC_BINARY_DIR"])
        cls.binary = binary_dir / "OpenOcean-Krakenc.exe"
        if not cls.binary.exists():
            cls.binary = binary_dir / "OpenOcean-Krakenc"
        cls.fixtures = Path(
            os.environ["OPENOCEANKRAKENC_TEST_ROOT"]
        ).resolve()

    def run_cli(self, *arguments):
        completed = subprocess.run(
            [str(self.binary), *map(str, arguments)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(completed.stdout)

    def copy_case(self, case, target):
        for source in self.fixtures.rglob(case + ".*"):
            if source.suffix.lower() in INPUT_SUFFIXES:
                shutil.copy2(source, target / source.name)

    def test_representative_matrix_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory(prefix="ookc_env_json_matrix_") as temp:
            root = Path(temp)
            for case in CASES:
                with self.subTest(case=case):
                    case_dir = root / case
                    case_dir.mkdir()
                    self.copy_case(case, case_dir)
                    env_path = case_dir / f"{case}.env"
                    json_path = case_dir / f"{case}_json.json"

                    self.run_cli("--convert-env", env_path, json_path)

                    env_eigen = self.run_cli("--eigen", env_path)
                    json_eigen = self.run_cli("--eigen", json_path)
                    self.assertEqual(
                        env_eigen["mode_count"], json_eigen["mode_count"])
                    self.assertEqual(
                        env_eigen["wavenumbers"], json_eigen["wavenumbers"])

                    mod_outputs = []
                    reference_mod = None
                    for source_name, source_path in (
                            ("env", env_path), ("json", json_path)):
                        for repetition in range(2):
                            mod_path = case_dir / (
                                f"{case}_{source_name}_{repetition}.mod")
                            self.run_cli(
                                "--mod", source_path, mod_path,
                                "--threads", "4")
                            if reference_mod is None:
                                reference_mod = mod_path
                            mod_outputs.append(mod_path.read_bytes())
                    self.assertTrue(mod_outputs)
                    self.assertTrue(all(
                        data == mod_outputs[0] for data in mod_outputs[1:]))

                    shutil.copy2(reference_mod, case_dir / f"{case}.mod")
                    legacy_field = self.run_cli("--field", case_dir / case)
                    legacy_shade = read_shade_file(Path(legacy_field["path"]))
                    env_field = self.run_cli("--field", env_path)
                    json_field = self.run_cli("--field", json_path)
                    env_shd = Path(env_field["path"])
                    json_shd = Path(json_field["path"])
                    self.assertEqual(env_field["source_depth_count"],
                                     json_field["source_depth_count"])
                    self.assertEqual(env_field["receiver_depth_count"],
                                     json_field["receiver_depth_count"])
                    self.assertEqual(env_field["range_count"],
                                     json_field["range_count"])
                    self.assertEqual(env_shd.read_bytes(), json_shd.read_bytes())
                    memory_shade = read_shade_file(env_shd)
                    self.assertEqual(len(memory_shade.pressure),
                                     len(legacy_shade.pressure))
                    error_energy = sum(
                        abs(left - right) ** 2
                        for left, right in zip(memory_shade.pressure,
                                               legacy_shade.pressure))
                    reference_energy = sum(
                        abs(value) ** 2 for value in legacy_shade.pressure)
                    relative_error = (
                        (error_energy / reference_energy) ** 0.5
                        if reference_energy else error_energy ** 0.5)
                    self.assertLessEqual(relative_error, 2.0e-5)


if __name__ == "__main__":
    unittest.main()
