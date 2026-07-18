from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT))

from py_kraken.ook_data_model import ConfigModel, FieldData, ModeData
from py_kraken.ook_read import read_mod, read_shd


class ReaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        repository = Path(__file__).resolve().parents[3]
        executable = repository / "bin" / "OpenOceanKraken.exe"
        fixture = (repository / ".." / "test" / "MunkK.env").resolve()
        cls.temporary_directory = tempfile.TemporaryDirectory()
        cls.result_root = Path(cls.temporary_directory.name) / "munk"
        result = subprocess.run(
            [
                str(executable),
                str(fixture),
                "--output",
                str(cls.result_root),
                "--threads",
                "1",
                "--mod",
            ],
            text=True,
            capture_output=True,
            check=False,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_config_model_accepts_native_and_extension_options(self):
        config = ConfigModel.model_validate(
            {"input": "MunkK.env", "threads": 2, "project_option": True}
        )
        self.assertEqual(config.input, "MunkK.env")
        self.assertEqual(config.threads, 2)
        self.assertTrue(config.project_option)

    def test_read_shd_returns_source_range_depth_complex64(self):
        field = read_shd(self.result_root.with_suffix(".shd"))
        self.assertIsInstance(field, FieldData)
        self.assertEqual(field.values.dtype, np.complex64)
        self.assertEqual(field.values.ndim, 3)
        self.assertEqual(field.values.shape[0], field.source_depths.size)
        self.assertEqual(field.values.shape[1], field.receiver_ranges.size)
        self.assertEqual(field.values.shape[2], field.receiver_depths.size)
        self.assertGreater(np.abs(field.values).max(), 0.0)

    def test_read_mod_returns_profiles_modes_and_depths(self):
        modes = read_mod(self.result_root.with_suffix(".mod"))
        self.assertIsInstance(modes, ModeData)
        self.assertGreater(len(modes.profiles), 0)
        first = modes.profiles[0]
        self.assertGreater(first.wavenumbers.size, 0)
        self.assertEqual(first.mode_shapes.shape[0], first.wavenumbers.size)
        self.assertEqual(first.mode_shapes.shape[1], first.depth.size)

    def test_truncated_files_raise_value_error(self):
        for suffix, reader in ((".shd", read_shd), (".mod", read_mod)):
            source = self.result_root.with_suffix(suffix)
            truncated = Path(self.temporary_directory.name) / f"truncated{suffix}"
            truncated.write_bytes(source.read_bytes()[:32])
            with self.subTest(suffix=suffix):
                with self.assertRaises(ValueError):
                    reader(truncated)


if __name__ == "__main__":
    unittest.main()
