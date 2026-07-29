import os
import subprocess
import tempfile
import unittest
from pathlib import Path


class RemediationCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        binary_dir = Path(
            os.environ["OPENOCEAN_KRAKENC_BINARY_DIR"])
        cls.binary = binary_dir / "OpenOcean-Krakenc.exe"
        if not cls.binary.exists():
            cls.binary = binary_dir / "OpenOcean-Krakenc"

    def run_failure(self, command, input_path, output_path=None):
        arguments = [str(self.binary), command, str(input_path)]
        if output_path is not None:
            arguments.append(str(output_path))
        completed = subprocess.run(
            arguments,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 1, completed.stderr)
        self.assertEqual(completed.stdout, "")
        return completed.stderr

    def run_both_failures(self, root, input_path):
        results = {}
        for command in ("--eigen", "--mod"):
            with self.subTest(command=command):
                output = (
                    root / "out.mod"
                    if command == "--mod"
                    else None
                )
                results[command] = self.run_failure(
                    command, input_path, output)
        return results

    def assert_error_contract(
            self, stderr, category, path, reason):
        self.assertIn(f"[{category}]", stderr)
        self.assertIn(str(path), stderr)
        self.assertIn(reason, stderr)

    def test_malformed_json(self):
        with tempfile.TemporaryDirectory(
                prefix="ookc_cli_error_") as temp:
            root = Path(temp)
            path = root / "malformed.json"
            path.write_text(
                '{"schemaVersion": 1,', encoding="utf-8")
            for stderr in self.run_both_failures(
                    root, path).values():
                self.assert_error_contract(
                    stderr, "parse_error", path, "line")
                self.assertIn("column", stderr)
                self.assertNotIn(
                    "parser byte/line-column position", stderr)

    def test_missing_required_field(self):
        with tempfile.TemporaryDirectory(
                prefix="ookc_cli_error_") as temp:
            root = Path(temp)
            path = root / "missing_field.json"
            path.write_text(
                '{"schemaVersion": 1, "Title": "bad"}',
                encoding="utf-8",
            )
            for stderr in self.run_both_failures(
                    root, path).values():
                self.assert_error_contract(
                    stderr, "schema_error", path, "freqinfo")

    def test_wrong_json_type(self):
        with tempfile.TemporaryDirectory(
                prefix="ookc_cli_error_") as temp:
            root = Path(temp)
            path = root / "wrong_type.json"
            path.write_text(
                '{"schemaVersion": "one"}',
                encoding="utf-8",
            )
            for stderr in self.run_both_failures(
                    root, path).values():
                self.assert_error_contract(
                    stderr, "schema_error", path, "number")

    def test_missing_file_for_both_commands(self):
        with tempfile.TemporaryDirectory(
                prefix="ookc_cli_error_") as temp:
            root = Path(temp)
            path = root / "does_not_exist.json"
            for stderr in self.run_both_failures(
                    root, path).values():
                self.assert_error_contract(
                    stderr,
                    "missing_file",
                    path,
                    "does not exist or is not a regular file",
                )
                self.assertNotIn("field=", stderr)


if __name__ == "__main__":
    unittest.main()
