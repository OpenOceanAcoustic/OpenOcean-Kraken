import json
import tempfile
import unittest
from pathlib import Path

from for_test.reference_pipeline import ReferencePipelineError, run_reference_pipeline


class ReferencePipelineTests(unittest.TestCase):
    def test_missing_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ReferencePipelineError, "missing required input"):
                run_reference_pipeline(
                    Path(tmp) / "missing",
                    Path("krakenc.exe"),
                    Path("field.exe"),
                    Path(tmp) / "artifacts",
                )

    def test_munk_kleaky_produces_mod_and_shd(self):
        project = Path(__file__).resolve().parents[1]
        workspace = project.parent
        result = run_reference_pipeline(
            workspace / "test" / "MunkKleaky",
            workspace / "krakenFortran" / "build_mingw" / "out" / "krakenc.exe",
            workspace / "krakenFortran" / "build_mingw" / "out" / "field.exe",
            project / "build" / "reference" / "MunkKleaky",
        )
        self.assertTrue(result["passes"], result)
        self.assertGreater(result["outputs"]["mod"]["bytes"], 0)
        self.assertGreater(result["outputs"]["shd"]["bytes"], 0)
        manifest = json.loads(
            (
                project
                / "build"
                / "reference"
                / "MunkKleaky"
                / "manifest.json"
            ).read_text(encoding="utf-8")
        )
        self.assertTrue(manifest["passes"])


if __name__ == "__main__":
    unittest.main()
