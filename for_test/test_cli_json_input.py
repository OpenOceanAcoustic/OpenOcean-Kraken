import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class CliJsonInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        binary_dir = Path(os.environ["OPENOCEAN_KRAKENC_BINARY_DIR"])
        cls.binary = binary_dir / "OpenOcean-Krakenc.exe"
        if not cls.binary.exists():
            cls.binary = binary_dir / "OpenOcean-Krakenc"
        cls.source = Path(__file__).resolve().parents[1]
        cls.env = cls.source / "for_test" / "fixtures" / "two_profile_small.env"

    def run_cli(self, *arguments):
        completed = subprocess.run(
            [str(self.binary), *map(str, arguments)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(completed.stdout)

    def test_env_json_equivalence_and_field(self):
        with tempfile.TemporaryDirectory(prefix="ookc_cli_json_") as temp:
            temp = Path(temp)
            env_path = temp / "case.env"
            flp_path = temp / "case.flp"
            shutil.copy2(self.env, env_path)
            shutil.copy2(self.env.with_suffix(".flp"), flp_path)
            json_path = temp / "case_json.json"
            converted = self.run_cli("--convert-env", env_path, json_path)
            self.assertTrue(converted["passes"])
            self.assertTrue(json_path.exists())

            env_eigen = self.run_cli("--eigen", env_path)
            json_eigen = self.run_cli("--eigen", json_path)
            self.assertEqual(env_eigen["mode_count"], json_eigen["mode_count"])
            self.assertEqual(env_eigen["wavenumbers"], json_eigen["wavenumbers"])

            env_mod = temp / "env.mod"
            json_mod = temp / "json.mod"
            self.run_cli("--mod", env_path, env_mod)
            self.run_cli("--mod", json_path, json_mod)
            self.assertEqual(env_mod.read_bytes(), json_mod.read_bytes())

            env_field = self.run_cli("--field", env_path)
            json_field = self.run_cli("--field", json_path)
            self.assertTrue(env_field["passes"])
            self.assertTrue(json_field["passes"])
            env_shd = Path(env_field["path"])
            json_shd = Path(json_field["path"])
            self.assertTrue(env_shd.exists())
            self.assertTrue(json_shd.exists())
            self.assertEqual(env_shd.read_bytes(), json_shd.read_bytes())


if __name__ == "__main__":
    unittest.main()
