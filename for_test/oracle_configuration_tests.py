import hashlib
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from for_test.oracle_config import OracleConfigurationError, oracle_identity, resolve_oracle_binary


class OracleConfigurationTests(unittest.TestCase):
    def test_explicit_root_resolves_binary_and_reports_identity(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "krakenc.exe"
            payload = b"frozen-krakenc-oracle"
            executable.write_bytes(payload)

            resolved = resolve_oracle_binary("krakenc.exe", root)
            identity = oracle_identity(resolved)

            self.assertEqual(executable.resolve(), resolved)
            self.assertEqual(str(executable.resolve()), identity["path"])
            self.assertEqual(hashlib.sha256(payload).hexdigest(), identity["sha256"])

    def test_environment_root_is_used_without_legacy_fallback(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "field.exe"
            executable.write_bytes(b"field")
            with patch.dict(os.environ, {"OPENOCEANKRAKENC_FORTRAN_ROOT": str(root)}, clear=True):
                self.assertEqual(executable.resolve(), resolve_oracle_binary("field.exe"))

    def test_missing_configuration_fails_clearly(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(
                OracleConfigurationError,
                "OPENOCEANKRAKENC_FORTRAN_ROOT",
            ):
                resolve_oracle_binary("krakenc.exe")

    def test_missing_binary_names_configured_root(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            with self.assertRaisesRegex(OracleConfigurationError, str(root).replace("\\", "\\\\")):
                resolve_oracle_binary("krakenc.exe", root)


if __name__ == "__main__":
    unittest.main()
