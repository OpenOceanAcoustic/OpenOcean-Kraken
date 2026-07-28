"""Gate OOK biological attenuation against the KrakenFortran reference."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import types
from typing import Any, Callable

import numpy as np


MAX_REAL_K = 1.0e-7
MIN_IMAG_K = 1.0e-9
REL_IMAG_K = 1.0e-3
MAX_SLOPE_ERROR_DB_PER_KM = 0.02
MAX_P95_TL_DB = 0.20
MAX_SINGLE_MODE_TL_DB = 0.50
MAX_COMMON_FIELD_P95_DB = 0.05
EXPECTED_RANGES_KM = np.array([1.0, 5.0, 10.0])


def modal_errors(ook_k: np.ndarray, fortran_k: np.ndarray) -> list[str]:
    ook_k = np.asarray(ook_k)
    fortran_k = np.asarray(fortran_k)
    errors: list[str] = []
    for label, values in (("OOK", ook_k), ("Fortran", fortran_k)):
        if values.ndim != 1:
            errors.append(
                f"{label} wavenumber shape {values.shape} must be "
                "one-dimensional"
            )
        if values.size == 0:
            errors.append(f"{label} wavenumbers contain no modes")
        if not np.all(np.isfinite(values.real)):
            errors.append(
                f"{label} wavenumber real parts contain non-finite values"
            )
        if not np.all(np.isfinite(values.imag)):
            errors.append(
                f"{label} wavenumber imaginary parts contain "
                "non-finite values"
            )
    if ook_k.shape != fortran_k.shape:
        errors.append(
            f"wavenumber shape mismatch: OOK={ook_k.shape}, "
            f"Fortran={fortran_k.shape}"
        )
    if ook_k.ndim == 1 and fortran_k.ndim == 1:
        if ook_k.size != fortran_k.size:
            errors.append(
                f"mode count mismatch: OOK={ook_k.size}, "
                f"Fortran={fortran_k.size}"
            )
    if errors:
        return errors

    real_delta = np.abs(ook_k.real - fortran_k.real)
    if np.any(real_delta > MAX_REAL_K):
        errors.append(
            f"real wavenumber error {real_delta.max()} exceeds "
            f"{MAX_REAL_K}"
        )
    imag_delta = np.abs(ook_k.imag - fortran_k.imag)
    imag_limit = np.maximum(
        MIN_IMAG_K, REL_IMAG_K * np.abs(fortran_k.imag)
    )
    if np.any(imag_delta > imag_limit):
        errors.append("imaginary wavenumber error exceeds per-mode limit")
    return errors


def slope_errors(
    ranges_km: np.ndarray, delta_tl: np.ndarray
) -> list[str]:
    ranges_km = np.asarray(ranges_km, dtype=float)
    delta_tl = np.asarray(delta_tl, dtype=float)
    errors: list[str] = []
    for label, values in (
        ("ranges", ranges_km),
        ("delta TL", delta_tl),
    ):
        if values.ndim != 1:
            errors.append(
                f"{label} shape {values.shape} must be one-dimensional"
            )
        if values.size < 2:
            errors.append(f"{label} must contain at least two samples")
        if not np.all(np.isfinite(values)):
            errors.append(f"{label} contains non-finite values")
    if ranges_km.shape != delta_tl.shape:
        errors.append(
            f"slope input shape mismatch: ranges={ranges_km.shape}, "
            f"delta TL={delta_tl.shape}"
        )
    if errors:
        return errors

    slope = float(np.polyfit(ranges_km, delta_tl, 1)[0])
    if abs(slope - 1.0) > MAX_SLOPE_ERROR_DB_PER_KM:
        return [
            f"attenuation slope {slope} dB/km differs from "
            "1.0 dB/km by more than 0.02 dB/km"
        ]
    return []


def tl_errors(
    ook_tl: np.ndarray, fortran_tl: np.ndarray
) -> list[str]:
    ook_tl = np.asarray(ook_tl, dtype=float)
    fortran_tl = np.asarray(fortran_tl, dtype=float)
    errors: list[str] = []
    for label, values in (("OOK", ook_tl), ("Fortran", fortran_tl)):
        if values.ndim != 1:
            errors.append(
                f"{label} TL shape {values.shape} must be one-dimensional"
            )
        if values.size == 0:
            errors.append(f"{label} TL contains no samples")
        if not np.all(np.isfinite(values)):
            errors.append(f"{label} TL contains non-finite values")
    if ook_tl.shape != fortran_tl.shape:
        errors.append(
            f"TL shape mismatch: OOK={ook_tl.shape}, "
            f"Fortran={fortran_tl.shape}"
        )
    if errors:
        return errors

    delta = np.abs(ook_tl - fortran_tl)
    if float(np.percentile(delta, 95)) > MAX_P95_TL_DB:
        errors.append("P95 |delta TL| exceeds 0.20 dB")
    if float(delta.max()) > MAX_SINGLE_MODE_TL_DB:
        errors.append("single-mode max |delta TL| exceeds 0.50 dB")
    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ook-exe", required=True, type=Path)
    parser.add_argument("--kraken-exe", required=True, type=Path)
    parser.add_argument("--field-exe", required=True, type=Path)
    parser.add_argument("--fixture-root", required=True, type=Path)
    return parser.parse_args()


def _single_thread_env() -> dict[str, str]:
    environment = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = "1"
    return environment


def _run_process(
    label: str,
    command: list[Path | str],
    cwd: Path,
) -> tuple[str, list[str]]:
    absolute_command = [str(item) for item in command]
    try:
        completed = subprocess.run(
            absolute_command,
            cwd=str(cwd),
            env=_single_thread_env(),
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=300,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return "", [f"{label} could not run: {error}"]

    output = "\n".join((completed.stdout, completed.stderr))
    errors: list[str] = []
    if completed.returncode != 0:
        errors.append(f"{label} exited with code {completed.returncode}")
    lowered = output.lower()
    for marker in ("fatal error", "no modes"):
        if marker in lowered:
            errors.append(f"{label} output contains {marker}")
    return output, errors


def _diagnostic_errors(label: str, path: Path) -> tuple[str, list[str]]:
    if not path.is_file():
        return "", [f"{label} diagnostic file is missing: {path}"]
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        return "", [f"{label} diagnostic file cannot be read: {error}"]

    errors: list[str] = []
    lowered = text.lower()
    for marker in ("fatal error", "no modes"):
        if marker in lowered:
            errors.append(f"{label} diagnostic contains {marker}")
    return text, errors


def _copy_fixture_pair(
    fixture_root: Path, case_name: str, destination: Path
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for suffix in (".env", ".flp"):
        source = fixture_root / f"{case_name}{suffix}"
        shutil.copy2(source, destination / source.name)


def _read_reference_ranges(path: Path) -> np.ndarray:
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read reference SHD {path}: {error}") from error
    if len(raw) < 4:
        raise ValueError(f"reference SHD {path} is too short")

    record_words = struct.unpack_from("<i", raw, 0)[0]
    record_bytes = 4 * record_words
    if record_words <= 0 or record_bytes > len(raw):
        raise ValueError(
            f"reference SHD {path} has invalid record length {record_words}"
        )

    dimension_offset = 2 * record_bytes
    if dimension_offset + 28 > len(raw):
        raise ValueError(f"reference SHD {path} has no dimension record")
    dimensions = struct.unpack_from("<7i", raw, dimension_offset)
    range_count = dimensions[-1]
    range_offset = 9 * record_bytes
    range_bytes = range_count * np.dtype("<f8").itemsize
    if (
        range_count <= 0
        or range_bytes > record_bytes
        or range_offset + range_bytes > len(raw)
    ):
        raise ValueError(
            f"reference SHD {path} has invalid range count {range_count}"
        )
    return np.frombuffer(
        raw,
        dtype="<f8",
        count=range_count,
        offset=range_offset,
    ).copy()


def _load_readers() -> tuple[Callable[[Path], Any], Callable[[Path], Any]]:
    repository_root = Path(__file__).resolve().parents[1]
    package_root = (
        repository_root / "wrappers" / "py_kraken" / "py_kraken"
    )
    package_name = "_ook_biological_gate_readers"
    if package_name not in sys.modules:
        package = types.ModuleType(package_name)
        package.__path__ = [str(package_root)]
        package.__package__ = package_name
        sys.modules[package_name] = package

    for module_basename in ("ook_data_model", "ook_read"):
        module_name = f"{package_name}.{module_basename}"
        if module_name in sys.modules:
            continue
        module_path = package_root / f"{module_basename}.py"
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"cannot load reader module {module_path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        try:
            spec.loader.exec_module(module)
        except Exception:
            sys.modules.pop(module_name, None)
            raise

    reader_module = sys.modules[f"{package_name}.ook_read"]
    return reader_module.read_mod, reader_module.read_shd


def _pressure_errors(label: str, pressure: np.ndarray) -> list[str]:
    errors: list[str] = []
    expected_shape = (1, EXPECTED_RANGES_KM.size, 1)
    if pressure.shape != expected_shape:
        errors.append(
            f"{label} pressure shape {pressure.shape} does not match "
            f"{expected_shape}"
        )
    if pressure.size != EXPECTED_RANGES_KM.size:
        errors.append(
            f"{label} pressure sample count {pressure.size} does not match "
            f"{EXPECTED_RANGES_KM.size}"
        )
    if not np.all(np.isfinite(pressure)):
        errors.append(f"{label} pressure contains non-finite values")
    if np.any(np.abs(pressure) == 0.0):
        errors.append(f"{label} pressure contains zero values")
    return errors


def _tl(pressure: np.ndarray) -> np.ndarray:
    values = np.asarray(pressure, dtype=np.complex128).reshape(-1)
    return -20.0 * np.log10(np.abs(values))


def _prepare_and_run_case(
    *,
    case_name: str,
    fixture_root: Path,
    work_root: Path,
    ook_exe: Path,
    kraken_exe: Path,
    field_exe: Path,
) -> tuple[dict[str, Path], list[str]]:
    case_tag = (
        "off" if case_name == "bio_uniform_off" else "on"
    )
    case_root = work_root / case_tag
    ook_dir = case_root / "o"
    fortran_dir = case_root / "f"
    common_dir = case_root / "c"
    for directory in (ook_dir, fortran_dir, common_dir):
        _copy_fixture_pair(fixture_root, case_name, directory)

    errors: list[str] = []
    ook_root = (ook_dir / case_name).resolve()
    ook_env = ook_root.with_suffix(".env")
    _, process_errors = _run_process(
        f"OOK {case_name}",
        [
            ook_exe,
            ook_env,
            "--output",
            ook_root,
            "--threads",
            "1",
            "--mod",
        ],
        ook_dir.resolve(),
    )
    errors.extend(process_errors)

    fortran_root = (fortran_dir / case_name).resolve()
    _, process_errors = _run_process(
        f"KrakenFortran {case_name}",
        [kraken_exe, fortran_root],
        fortran_dir.resolve(),
    )
    errors.extend(process_errors)
    _, process_errors = _run_process(
        f"KrakenFortran field {case_name}",
        [field_exe, fortran_root],
        fortran_dir.resolve(),
    )
    errors.extend(process_errors)

    prt_text, diagnostic = _diagnostic_errors(
        f"KrakenFortran {case_name}", fortran_root.with_suffix(".prt")
    )
    errors.extend(diagnostic)
    _, diagnostic = _diagnostic_errors(
        f"KrakenFortran field {case_name}",
        fortran_root.with_suffix(".fprt"),
    )
    errors.extend(diagnostic)
    if case_name == "bio_uniform_resonance":
        if "Biological attenuation" not in prt_text:
            errors.append(
                "KrakenFortran resonance PRT does not contain "
                "Biological attenuation"
            )

    common_root = (common_dir / case_name).resolve()
    ook_mod = ook_root.with_suffix(".mod")
    if ook_mod.is_file():
        shutil.copy2(ook_mod, common_root.with_suffix(".mod"))
        _, process_errors = _run_process(
            f"KrakenFortran common field {case_name}",
            [field_exe, common_root],
            common_dir.resolve(),
        )
        errors.extend(process_errors)
        _, diagnostic = _diagnostic_errors(
            f"KrakenFortran common field {case_name}",
            common_root.with_suffix(".fprt"),
        )
        errors.extend(diagnostic)
    else:
        errors.append(f"OOK {case_name} MOD file is missing: {ook_mod}")

    paths = {
        "ook_mod": ook_mod,
        "ook_shd": ook_root.with_suffix(".shd"),
        "fortran_mod": fortran_root.with_suffix(".mod"),
        "fortran_shd": fortran_root.with_suffix(".shd"),
        "common_shd": common_root.with_suffix(".shd"),
    }
    for label, path in paths.items():
        if not path.is_file():
            errors.append(f"{case_name} {label} output is missing: {path}")
    return paths, errors


def run_gate(args: argparse.Namespace) -> list[str]:
    ook_exe = Path(args.ook_exe).resolve()
    kraken_exe = Path(args.kraken_exe).resolve()
    field_exe = Path(args.field_exe).resolve()
    fixture_root = Path(args.fixture_root).resolve()

    errors: list[str] = []
    for label, executable in (
        ("OOK executable", ook_exe),
        ("KrakenFortran executable", kraken_exe),
        ("KrakenFortran field executable", field_exe),
    ):
        if not executable.is_file():
            errors.append(f"{label} is missing: {executable}")
    if not fixture_root.is_dir():
        errors.append(f"fixture root is missing: {fixture_root}")

    case_names = ("bio_uniform_off", "bio_uniform_resonance")
    for case_name in case_names:
        for suffix in (".env", ".flp"):
            fixture = fixture_root / f"{case_name}{suffix}"
            if not fixture.is_file():
                errors.append(f"fixture is missing: {fixture}")
    if errors:
        return errors

    try:
        read_mod, read_shd = _load_readers()
    except Exception as error:
        return [f"cannot import OOK MOD/SHD readers: {error}"]

    case_data: dict[str, dict[str, np.ndarray]] = {}
    # KrakenFortran stores the case root in a fixed-width character buffer.
    # Keep the otherwise explicit absolute temporary roots comfortably short.
    with tempfile.TemporaryDirectory(prefix="ob-") as temp:
        work_root = Path(temp).resolve()
        for case_name in case_names:
            paths, case_errors = _prepare_and_run_case(
                case_name=case_name,
                fixture_root=fixture_root,
                work_root=work_root,
                ook_exe=ook_exe,
                kraken_exe=kraken_exe,
                field_exe=field_exe,
            )
            errors.extend(case_errors)
            if any(not path.is_file() for path in paths.values()):
                continue

            try:
                ook_modes = read_mod(paths["ook_mod"])
                fortran_modes = read_mod(paths["fortran_mod"])
                ook_field = read_shd(paths["ook_shd"])
                fortran_field = read_shd(paths["fortran_shd"])
                common_field = read_shd(paths["common_shd"])
                fortran_ranges_km = (
                    _read_reference_ranges(paths["fortran_shd"]) / 1000.0
                )
                common_ranges_km = (
                    _read_reference_ranges(paths["common_shd"]) / 1000.0
                )
            except (OSError, TypeError, ValueError) as error:
                errors.append(f"{case_name} output cannot be read: {error}")
                continue

            if len(ook_modes.profiles) != 1:
                errors.append(
                    f"{case_name} OOK MOD has {len(ook_modes.profiles)} "
                    "profiles; expected 1"
                )
                continue
            if len(fortran_modes.profiles) != 1:
                errors.append(
                    f"{case_name} Fortran MOD has "
                    f"{len(fortran_modes.profiles)} profiles; expected 1"
                )
                continue

            ook_k = ook_modes.profiles[0].wavenumbers
            fortran_k = fortran_modes.profiles[0].wavenumbers
            modal_case_errors = modal_errors(ook_k, fortran_k)
            errors.extend(
                f"{case_name}: {item}" for item in modal_case_errors
            )

            ook_pressure = ook_field.values
            fortran_pressure = fortran_field.values
            common_pressure = common_field.values
            pressure_errors: list[str] = []
            pressure_errors.extend(
                _pressure_errors(f"{case_name} OOK", ook_pressure)
            )
            pressure_errors.extend(
                _pressure_errors(f"{case_name} Fortran", fortran_pressure)
            )
            pressure_errors.extend(
                _pressure_errors(
                    f"{case_name} common field", common_pressure
                )
            )
            errors.extend(pressure_errors)
            range_sets = (
                (
                    "OOK",
                    np.asarray(ook_field.receiver_ranges, dtype=float)
                    / 1000.0,
                ),
                ("Fortran", fortran_ranges_km),
                ("common field", common_ranges_km),
            )
            for label, ranges_km in range_sets:
                if ranges_km.shape != EXPECTED_RANGES_KM.shape:
                    errors.append(
                        f"{case_name} {label} range shape "
                        f"{ranges_km.shape} does not match "
                        f"{EXPECTED_RANGES_KM.shape}"
                    )
                elif not np.array_equal(ranges_km, EXPECTED_RANGES_KM):
                    errors.append(
                        f"{case_name} {label} ranges "
                        f"{ranges_km.tolist()} do not match fixed FLP order "
                        f"{EXPECTED_RANGES_KM.tolist()}"
                    )

            common_p95 = float("nan")
            if common_pressure.shape != fortran_pressure.shape:
                errors.append(
                    f"{case_name} common field shape "
                    f"{common_pressure.shape} does not match Fortran "
                    f"{fortran_pressure.shape}"
                )
            elif (
                np.all(np.isfinite(common_pressure))
                and np.all(np.isfinite(fortran_pressure))
                and np.all(np.abs(common_pressure) > 0.0)
                and np.all(np.abs(fortran_pressure) > 0.0)
            ):
                common_delta = np.abs(
                    _tl(common_pressure) - _tl(fortran_pressure)
                )
                common_p95 = float(np.percentile(common_delta, 95))
                if common_p95 > MAX_COMMON_FIELD_P95_DB:
                    errors.append(
                        f"{case_name}: common-field P95 |delta TL| "
                        f"{common_p95} exceeds "
                        f"{MAX_COMMON_FIELD_P95_DB}"
                    )

            if modal_case_errors:
                modal_summary = (
                    "modal comparison skipped due to invalid MOD inputs"
                )
            else:
                max_real = float(
                    np.max(np.abs(ook_k.real - fortran_k.real))
                )
                max_imag = float(
                    np.max(np.abs(ook_k.imag - fortran_k.imag))
                )
                modal_summary = (
                    f"max |delta Re(k)|={max_real}; "
                    f"max |delta Im(k)|={max_imag}"
                )
            print(
                f"{case_name}: modes OOK={ook_k.size}, "
                f"Fortran={fortran_k.size}; {modal_summary}; "
                f"common-field P95 |delta TL|={common_p95}"
            )
            if not modal_case_errors and not pressure_errors:
                case_data[case_name] = {
                    "ook_k": ook_k,
                    "fortran_k": fortran_k,
                    "ook_pressure": ook_pressure,
                    "fortran_pressure": fortran_pressure,
                }

    if len(case_data) != len(case_names):
        return errors

    off = case_data["bio_uniform_off"]
    on = case_data["bio_uniform_resonance"]
    for implementation in ("ook", "fortran"):
        on_k = on[f"{implementation}_k"]
        off_k = off[f"{implementation}_k"]
        display = "OOK" if implementation == "ook" else "Fortran"
        if not np.all(on_k.imag < 0.0):
            errors.append(
                f"{display} resonance MOD does not have Im(k) < 0 "
                "for every mode"
            )
        if on_k.size != off_k.size:
            errors.append(
                f"{display} on/off mode count mismatch: "
                f"on={on_k.size}, off={off_k.size}"
            )
        elif not np.any(np.abs(on_k.imag - off_k.imag) > 1.0e-9):
            errors.append(
                f"{display} on/off imaginary wavenumbers do not differ "
                "by more than 1e-9"
            )

    ook_delta_tl = _tl(on["ook_pressure"]) - _tl(off["ook_pressure"])
    fortran_delta_tl = (
        _tl(on["fortran_pressure"]) - _tl(off["fortran_pressure"])
    )
    errors.extend(
        f"OOK: {item}"
        for item in slope_errors(EXPECTED_RANGES_KM, ook_delta_tl)
    )
    errors.extend(
        f"Fortran: {item}"
        for item in slope_errors(EXPECTED_RANGES_KM, fortran_delta_tl)
    )
    errors.extend(
        f"self-field biological delta: {item}"
        for item in tl_errors(ook_delta_tl, fortran_delta_tl)
    )

    ook_slope = float(
        np.polyfit(EXPECTED_RANGES_KM, ook_delta_tl, 1)[0]
    )
    fortran_slope = float(
        np.polyfit(EXPECTED_RANGES_KM, fortran_delta_tl, 1)[0]
    )
    own_delta = np.abs(ook_delta_tl - fortran_delta_tl)
    print(
        "OOK TL_on-TL_off="
        f"{ook_delta_tl.tolist()}; slope={ook_slope} dB/km"
    )
    print(
        "Fortran TL_on-TL_off="
        f"{fortran_delta_tl.tolist()}; slope={fortran_slope} dB/km"
    )
    print(
        "Self-field biological delta: "
        f"P95={float(np.percentile(own_delta, 95))} dB; "
        f"max={float(own_delta.max())} dB"
    )
    return errors


def main() -> int:
    args = parse_args()
    errors = run_gate(args)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("Biological attenuation equivalence gate passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
