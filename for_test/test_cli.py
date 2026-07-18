from pathlib import Path
import argparse
import subprocess
import tempfile
import unittest


class CliTests(unittest.TestCase):
    executable: Path
    fixture: Path

    def run_cli(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(self.executable), *args],
            text=True,
            capture_output=True,
            check=False,
            encoding="utf-8",
            errors="replace",
        )

    def test_help_is_successful(self):
        result = self.run_cli("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--output", result.stdout)
        self.assertIn("--mod-only", result.stdout)

    def test_missing_input_is_usage_error(self):
        result = self.run_cli()
        self.assertEqual(result.returncode, 2)
        self.assertIn("Usage:", result.stderr)

    def test_unknown_option_is_usage_error(self):
        result = self.run_cli(str(self.fixture), "--unknown")
        self.assertEqual(result.returncode, 2)
        self.assertIn("Unknown option", result.stderr)

    def test_mod_only_writes_mod(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "case"
            result = self.run_cli(
                str(self.fixture),
                "--output",
                str(root),
                "--threads",
                "1",
                "--mod-only",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(root.with_suffix(".mod").is_file())

    def test_velocity_writes_three_shd_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "case"
            result = self.run_cli(
                str(self.fixture),
                "--output",
                str(root),
                "--threads",
                "1",
                "--velocity",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            for suffix in ("_P.shd", "_V.shd", "_H.shd"):
                self.assertTrue(Path(f"{root}{suffix}").is_file())


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--fixture", required=True, type=Path)
    args, remaining = parser.parse_known_args()
    CliTests.executable = args.executable.resolve()
    CliTests.fixture = args.fixture.resolve()
    unittest.main(argv=[__file__, *remaining])


if __name__ == "__main__":
    main()
