from contextlib import redirect_stderr, redirect_stdout
import io
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from tools import test_bio_kraken_equivalence as gate
from tools.test_bio_kraken_equivalence import (
    _load_readers,
    _read_reference_ranges,
    modal_errors,
    slope_errors,
    tl_errors,
)


class BiologicalGateTests(unittest.TestCase):
    def test_biological_fixture_ignore_rules_are_exact(self):
        repository_root = Path(__file__).resolve().parents[2]
        git = shutil.which("git")
        if git is None:
            self.skipTest("git executable is unavailable")
        worktree = subprocess.run(
            [git, "rev-parse", "--is-inside-work-tree"],
            cwd=repository_root,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
        if worktree.returncode != 0 or worktree.stdout.strip() != "true":
            self.skipTest("test source is not inside a Git worktree")

        def is_ignored(relative_path):
            completed = subprocess.run(
                [
                    git,
                    "check-ignore",
                    "--no-index",
                    "--quiet",
                    relative_path,
                ],
                cwd=repository_root,
                check=False,
                timeout=10,
            )
            self.assertIn(
                completed.returncode,
                (0, 1),
                f"git check-ignore failed for {relative_path}",
            )
            return completed.returncode == 0

        allowed = (
            "for_test/fixtures/biological/bio_uniform_off.env",
            "for_test/fixtures/biological/bio_uniform_off.flp",
            "for_test/fixtures/biological/bio_uniform_resonance.env",
            "for_test/fixtures/biological/bio_uniform_resonance.flp",
        )
        ignored = (
            "for_test/fixtures/biological/extra.env",
            "for_test/fixtures/biological/extra.flp",
            "for_test/fixtures/biological/bio_uniform_off.mod",
            "for_test/fixtures/biological/bio_uniform_off.shd",
        )
        for relative_path in allowed:
            with self.subTest(relative_path=relative_path):
                self.assertFalse(is_ignored(relative_path))
        for relative_path in ignored:
            with self.subTest(relative_path=relative_path):
                self.assertTrue(is_ignored(relative_path))

    def test_reference_gate_configuration_requires_python(self):
        repository_root = Path(__file__).resolve().parents[2]
        cmake = shutil.which("cmake")
        self.assertIsNotNone(cmake)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            build = root / "build"
            kraken = root / "kraken.exe"
            field = root / "field.exe"
            kraken.touch()
            field.touch()
            completed = subprocess.run(
                [
                    cmake,
                    "-S",
                    str(repository_root),
                    "-B",
                    str(build),
                    "-DBUILD_TESTING=ON",
                    "-DBUILD_AS_EXE=ON",
                    "-DBUILD_AS_PYTHON=OFF",
                    "-DOOK_ENABLE_KRAKEN_REFERENCE_TESTS=ON",
                    f"-DOOK_KRAKEN_EXE={kraken}",
                    f"-DOOK_FIELD_EXE={field}",
                    f"-DPython3_EXECUTABLE={root / 'missing-python.exe'}",
                ],
                check=False,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=120,
            )
        output = completed.stdout + completed.stderr
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(
            "Reference tests require a Python3 interpreter",
            output,
        )

    def test_fortran_range_record_is_decoded_as_float64(self):
        record_bytes = 32
        raw = bytearray(10 * record_bytes)
        struct.pack_into("<i", raw, 0, record_bytes // 4)
        struct.pack_into(
            "<7i",
            raw,
            2 * record_bytes,
            1,
            1,
            1,
            1,
            1,
            1,
            3,
        )
        struct.pack_into(
            "<3d", raw, 9 * record_bytes, 1000.0, 5000.0, 10000.0
        )

        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "reference.shd"
            path.write_bytes(raw)
            ranges = _read_reference_ranges(path)

        np.testing.assert_array_equal(
            ranges, np.array([1000.0, 5000.0, 10000.0])
        )

    def test_readers_load_without_plotting_package_initialization(self):
        read_mod, read_shd = _load_readers()

        self.assertTrue(callable(read_mod))
        self.assertTrue(callable(read_shd))

    def _run_mocked_gate(
        self,
        *,
        off_ook_wavenumbers=None,
        off_fortran_wavenumbers=None,
        off_ook_pressure_shape=(1, 3, 1),
        off_process_errors=(),
    ):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture_root = root / "fixtures"
            fixture_root.mkdir()
            for case_name in (
                "bio_uniform_off",
                "bio_uniform_resonance",
            ):
                for suffix in (".env", ".flp"):
                    (fixture_root / f"{case_name}{suffix}").touch()
            executables = []
            for name in ("ook.exe", "kraken.exe", "field.exe"):
                executable = root / name
                executable.touch()
                executables.append(executable)

            def prepare_case(**kwargs):
                case_name = kwargs["case_name"]
                paths = {}
                for label in (
                    "ook_mod",
                    "ook_shd",
                    "fortran_mod",
                    "fortran_shd",
                    "common_shd",
                ):
                    path = root / f"{case_name}_{label}.out"
                    path.touch()
                    paths[label] = path
                case_errors = (
                    list(off_process_errors)
                    if case_name == "bio_uniform_off"
                    else []
                )
                return paths, case_errors

            valid_off = np.array([1.0 - 1.0e-4j, 0.9 - 2.0e-4j])
            valid_on = np.array([1.0 - 2.0e-4j, 0.9 - 3.0e-4j])

            def read_mod(path):
                name = Path(path).name
                is_off = "bio_uniform_off" in name
                is_ook = "_ook_mod" in name
                if is_off and is_ook and off_ook_wavenumbers is not None:
                    wavenumbers = off_ook_wavenumbers
                elif (
                    is_off
                    and not is_ook
                    and off_fortran_wavenumbers is not None
                ):
                    wavenumbers = off_fortran_wavenumbers
                else:
                    wavenumbers = valid_off if is_off else valid_on
                return SimpleNamespace(
                    profiles=[
                        SimpleNamespace(wavenumbers=wavenumbers)
                    ]
                )

            def read_shd(path):
                invalid_off_ook = (
                    "bio_uniform_off_ook_shd" in Path(path).name
                )
                shape = (
                    off_ook_pressure_shape
                    if invalid_off_ook
                    else (1, 3, 1)
                )
                return SimpleNamespace(
                    values=np.ones(shape, dtype=np.complex128),
                    receiver_ranges=np.array([1000.0, 5000.0, 10000.0]),
                )

            args = SimpleNamespace(
                ook_exe=executables[0],
                kraken_exe=executables[1],
                field_exe=executables[2],
                fixture_root=fixture_root,
            )
            stderr = io.StringIO()
            stdout = io.StringIO()
            with (
                patch.object(
                    gate, "_prepare_and_run_case", side_effect=prepare_case
                ),
                patch.object(
                    gate, "_load_readers", return_value=(read_mod, read_shd)
                ),
                patch.object(
                    gate,
                    "_read_reference_ranges",
                    return_value=np.array([1000.0, 5000.0, 10000.0]),
                ),
                patch.object(gate, "parse_args", return_value=args),
                redirect_stderr(stderr),
                redirect_stdout(stdout),
            ):
                exit_code = gate.main()

        return exit_code, stdout.getvalue(), stderr.getvalue()

    def test_pressure_shape_mismatch_returns_readable_gate_error(self):
        exit_code, _, diagnostic = self._run_mocked_gate(
            off_ook_pressure_shape=(1, 2, 1)
        )

        self.assertEqual(exit_code, 1)
        self.assertIn("bio_uniform_off OOK pressure shape", diagnostic)
        self.assertNotIn("operands could not be broadcast", diagnostic)

    def test_mod_shape_mismatch_skips_modal_and_cross_case_math(self):
        for pressure_shape in ((1, 3, 1), (1, 2, 1)):
            with self.subTest(pressure_shape=pressure_shape):
                exit_code, stdout, stderr = self._run_mocked_gate(
                    off_ook_wavenumbers=np.ones(
                        (2, 2), dtype=np.complex128
                    ),
                    off_fortran_wavenumbers=np.ones(
                        (1, 4), dtype=np.complex128
                    ),
                    off_ook_pressure_shape=pressure_shape,
                    off_process_errors=(
                        "bio_uniform_off synthetic process error",
                    ),
                )

                self.assertEqual(exit_code, 1)
                self.assertIn(
                    "bio_uniform_off: OOK wavenumber shape",
                    stderr,
                )
                self.assertIn("wavenumber shape mismatch", stderr)
                if pressure_shape == (1, 2, 1):
                    self.assertIn(
                        "bio_uniform_off OOK pressure shape",
                        stderr,
                    )
                else:
                    self.assertNotIn(
                        "bio_uniform_off OOK pressure shape",
                        stderr,
                    )
                self.assertIn(
                    "bio_uniform_off synthetic process error",
                    stderr,
                )
                self.assertIn("modal comparison skipped", stdout)
                self.assertNotIn(
                    "operands could not be broadcast",
                    stderr,
                )

    def test_modal_threshold_breach_is_reported(self):
        errors = modal_errors(
            np.array([1.0 - 1.0e-4j]),
            np.array([1.0 - 1.0e-4j]),
        )
        self.assertEqual(errors, [])
        errors = modal_errors(
            np.array([1.0 + 2.0e-7 - 1.0e-4j]),
            np.array([1.0 - 1.0e-4j]),
        )
        self.assertTrue(any("real wavenumber" in item for item in errors))

    def test_off_mod_rejects_nonfinite_nonfirst_wavenumber_parts(self):
        valid = np.array([1.0 - 1.0e-4j, 0.9 - 2.0e-4j])
        cases = (
            ("OOK", "real", np.nan, -2.0e-4),
            ("OOK", "real", np.inf, -2.0e-4),
            ("OOK", "imaginary", 0.9, np.nan),
            ("OOK", "imaginary", 0.9, np.inf),
            ("Fortran", "real", np.nan, -2.0e-4),
            ("Fortran", "real", np.inf, -2.0e-4),
            ("Fortran", "imaginary", 0.9, np.nan),
            ("Fortran", "imaginary", 0.9, np.inf),
        )
        for side, component, real, imag in cases:
            with self.subTest(side=side, component=component, value=(real, imag)):
                corrupted = valid.copy()
                corrupted[1] = complex(real, imag)
                ook = corrupted if side == "OOK" else valid
                fortran = corrupted if side == "Fortran" else valid

                errors = modal_errors(ook, fortran)

                diagnostic = " ".join(errors)
                self.assertIn(side, diagnostic)
                self.assertIn(component, diagnostic)
                self.assertIn("non-finite", diagnostic)

    def test_numeric_gate_functions_reject_bad_shapes(self):
        modal_diagnostic = " ".join(
            modal_errors(
                np.array([[1.0 - 1.0e-4j]]),
                np.array([1.0 - 1.0e-4j]),
            )
        )
        self.assertIn("shape", modal_diagnostic)
        self.assertIn("OOK", modal_diagnostic)

        slope_diagnostic = " ".join(
            slope_errors(
                np.array([[1.0], [5.0], [10.0]]),
                np.array([1.0, 5.0, 10.0]),
            )
        )
        self.assertIn("shape", slope_diagnostic)
        self.assertIn("ranges", slope_diagnostic)

        tl_diagnostic = " ".join(
            tl_errors(
                np.array([[10.0], [10.1]]),
                np.array([10.0, 10.0]),
            )
        )
        self.assertIn("shape", tl_diagnostic)
        self.assertIn("OOK", tl_diagnostic)

    def test_slope_and_tl_reject_nonfinite_inputs_by_side(self):
        finite_ranges = np.array([1.0, 5.0, 10.0])
        finite_tl = np.array([1.0, 5.0, 10.0])
        for value in (np.nan, np.inf):
            with self.subTest(function="slope", side="ranges", value=value):
                ranges = finite_ranges.copy()
                ranges[1] = value
                diagnostic = " ".join(slope_errors(ranges, finite_tl))
                self.assertIn("ranges", diagnostic)
                self.assertIn("non-finite", diagnostic)
            with self.subTest(function="slope", side="delta TL", value=value):
                delta_tl = finite_tl.copy()
                delta_tl[1] = value
                diagnostic = " ".join(
                    slope_errors(finite_ranges, delta_tl)
                )
                self.assertIn("delta TL", diagnostic)
                self.assertIn("non-finite", diagnostic)
            for side in ("OOK", "Fortran"):
                with self.subTest(function="TL", side=side, value=value):
                    corrupted = finite_tl.copy()
                    corrupted[1] = value
                    ook = corrupted if side == "OOK" else finite_tl
                    fortran = corrupted if side == "Fortran" else finite_tl
                    diagnostic = " ".join(tl_errors(ook, fortran))
                    self.assertIn(side, diagnostic)
                    self.assertIn("non-finite", diagnostic)

    def test_tl_and_slope_threshold_breaches_are_reported(self):
        self.assertEqual(
            tl_errors(np.array([10.0, 10.1]), np.array([10.0, 10.0])),
            [],
        )
        self.assertTrue(
            tl_errors(np.array([10.0, 10.6]), np.array([10.0, 10.0]))
        )
        self.assertEqual(
            slope_errors(
                np.array([1.0, 5.0, 10.0]),
                np.array([1.0, 5.0, 10.0]),
            ),
            [],
        )
        self.assertTrue(
            slope_errors(
                np.array([1.0, 5.0, 10.0]),
                np.array([0.0, 0.0, 0.0]),
            )
        )


if __name__ == "__main__":
    unittest.main()
