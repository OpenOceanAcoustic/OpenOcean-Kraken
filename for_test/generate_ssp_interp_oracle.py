"""Generate the frozen ALG-022 oracle directly from the Fortran KrakenC."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

from mod_reader import read_first_wavenumber_set


PROFILE_ENV_NAMES = {
    "N": "ssp_interp_n.env",
    "C": "ssp_interp_c.env",
    "P": "ssp_interp_p.env",
    "S": "ssp_interp_s.env",
}
LOSSY_ENV_NAME = "lossy_gradient_interpolation.env"
SOURCE_NAMES = (
    "sspMod.f90",
    "pchipMod.f90",
    "splinec.f90",
    "AttenMod.f90",
)
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--fortran-source-root", required=True, type=Path)
    parser.add_argument("--fortran-build-root", required=True, type=Path)
    parser.add_argument("--gfortran", required=True, type=Path)
    parser.add_argument("--date", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def _required_file(path: Path, description: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{description} is not a file: {resolved}")
    return resolved


def _required_directory(path: Path, description: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_dir():
        raise FileNotFoundError(
            f"{description} is not a directory: {resolved}"
        )
    return resolved


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    value = digest.hexdigest()
    if not SHA256_PATTERN.fullmatch(value):
        raise RuntimeError(f"invalid SHA-256 generated for {path}")
    return value


def _file_record(path: Path) -> dict[str, Any]:
    resolved = _required_file(path, "provenance input")
    return {"path": str(resolved), "sha256": _sha256(resolved)}


def _runtime_environment(gfortran: Path) -> dict[str, str]:
    environment = os.environ.copy()
    existing_path = environment.get("PATH", "")
    environment["PATH"] = (
        str(gfortran.parent)
        if not existing_path
        else str(gfortran.parent) + os.pathsep + existing_path
    )
    return environment


def _command_record(
    command: list[str], completed: subprocess.CompletedProcess[str], cwd: Path
) -> dict[str, Any]:
    return {
        "command": command,
        "returncode": completed.returncode,
        "cwd": str(cwd.resolve()),
    }


def _run_profile(
    oracle_exe: Path,
    profile_type: str,
    environment: dict[str, str],
) -> tuple[list[dict[str, float]], dict[str, Any]]:
    with tempfile.TemporaryDirectory(
        prefix=f"ookc-alg022-profile-{profile_type.lower()}-"
    ) as run_name:
        run_dir = Path(run_name).resolve()
        command = [str(oracle_exe.resolve()), profile_type]
        completed = subprocess.run(
            command,
            cwd=run_dir,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
            env=environment,
        )
        record = _command_record(command, completed, run_dir)
        if completed.returncode != 0:
            raise RuntimeError(
                f"Fortran SSP oracle failed for {profile_type}: "
                f"return code {completed.returncode}; "
                f"stderr={completed.stderr.strip()}"
            )
        lines = [
            line.strip()
            for line in completed.stdout.splitlines()
            if line.strip()
        ]
        if len(lines) != 5:
            raise RuntimeError(
                f"Fortran SSP oracle returned {len(lines)} non-empty lines "
                f"for {profile_type}, expected 5"
            )
        samples: list[dict[str, float]] = []
        keys = ("depth", "cp_re", "cp_im", "cs_re", "cs_im", "rho")
        for index, line in enumerate(lines):
            fields = line.split()
            if len(fields) != 7 or fields[0] != profile_type:
                raise RuntimeError(
                    f"invalid Fortran SSP oracle row {index + 1} for "
                    f"{profile_type}: {line!r}"
                )
            try:
                values = [float(field) for field in fields[1:]]
            except ValueError as error:
                raise RuntimeError(
                    f"non-numeric Fortran SSP oracle row for {profile_type}: "
                    f"{line!r}"
                ) from error
            if not all(math.isfinite(value) for value in values):
                raise RuntimeError(
                    f"non-finite Fortran SSP oracle row for {profile_type}"
                )
            expected_depth = 50.0 + 12.5 * index
            if abs(values[0] - expected_depth) > 1.0e-12:
                raise RuntimeError(
                    f"unexpected depth for {profile_type} row {index + 1}: "
                    f"{values[0]}"
                )
            samples.append(dict(zip(keys, values, strict=True)))
        return samples, record


def _run_krakenc(
    krakenc: Path,
    env_path: Path,
    environment: dict[str, str],
) -> tuple[list[dict[str, float]], dict[str, Any]]:
    with tempfile.TemporaryDirectory(
        prefix=f"ookc-alg022-krakenc-{env_path.stem}-"
    ) as run_name:
        run_dir = Path(run_name).resolve()
        target_env = run_dir / env_path.name
        shutil.copy2(env_path, target_env)
        command = [str(krakenc), target_env.stem]
        completed = subprocess.run(
            command,
            cwd=run_dir,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
            env=environment,
        )
        record = _command_record(command, completed, run_dir)
        if completed.returncode != 0:
            raise RuntimeError(
                f"Fortran KrakenC failed for {env_path.name}: "
                f"return code {completed.returncode}; "
                f"stderr={completed.stderr.strip()}"
            )
        mod_path = run_dir / f"{target_env.stem}.mod"
        if not mod_path.is_file():
            raise RuntimeError(
                f"Fortran KrakenC did not create {mod_path.name}"
            )
        values = read_first_wavenumber_set(mod_path)
        if not values:
            raise RuntimeError(
                f"Fortran KrakenC returned no modes for {env_path.name}"
            )
        wavenumbers = []
        for index, value in enumerate(values):
            if not math.isfinite(value.real) or not math.isfinite(value.imag):
                raise RuntimeError(
                    f"Fortran KrakenC returned non-finite mode {index + 1} "
                    f"for {env_path.name}"
                )
            wavenumbers.append({"real": value.real, "imag": value.imag})
        return wavenumbers, record


def _validate_payload(payload: dict[str, Any]) -> None:
    if payload["schema"] != "OpenOcean-Krakenc.ssp-interpolation-oracle":
        raise RuntimeError("oracle schema mismatch before write")
    if payload["schema_version"] != 1:
        raise RuntimeError("oracle schema version mismatch before write")

    provenance = payload["provenance"]
    records = [provenance["krakenc"], provenance["libMisc"]]
    records.extend(provenance["sources"].values())
    records.extend(provenance["inputs"].values())
    records.append(provenance["driver"])
    for record in records:
        if not Path(record["path"]).is_absolute():
            raise RuntimeError("provenance path is not absolute")
        if not SHA256_PATTERN.fullmatch(record["sha256"]):
            raise RuntimeError("provenance SHA-256 is missing or malformed")

    for profile_type in PROFILE_ENV_NAMES:
        profile = payload["profiles"][profile_type]
        if len(profile["samples"]) != 5:
            raise RuntimeError(
                f"profile {profile_type} does not contain five samples"
            )
        if profile["mode_count"] != len(profile["wavenumbers"]):
            raise RuntimeError(
                f"profile {profile_type} mode count mismatch"
            )
        if profile["mode_count"] < 1:
            raise RuntimeError(f"profile {profile_type} has no modes")
    lossy = payload["lossy_gradient_interpolation"]
    if lossy["mode_count"] != len(lossy["wavenumbers"]):
        raise RuntimeError("lossy mode count mismatch")
    if lossy["mode_count"] < 1:
        raise RuntimeError("lossy profile has no modes")

    def validate_finite(value: Any, location: str) -> None:
        if isinstance(value, float) and not math.isfinite(value):
            raise RuntimeError(f"non-finite value at {location}")
        if isinstance(value, dict):
            for key, item in value.items():
                validate_finite(item, f"{location}.{key}")
        elif isinstance(value, list):
            for index, item in enumerate(value):
                validate_finite(item, f"{location}[{index}]")

    validate_finite(payload, "oracle")


def _atomic_write(payload: dict[str, Any], output: Path) -> None:
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=output.name + ".",
            suffix=".tmp",
            dir=output.parent,
            delete=False,
        ) as stream:
            temporary_path = Path(stream.name)
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        temporary_path.replace(output)
    except BaseException:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
        raise


def main() -> int:
    arguments = _arguments()
    project_root = _required_directory(
        arguments.project_root, "project root"
    )
    fortran_source_root = _required_directory(
        arguments.fortran_source_root, "Fortran source root"
    )
    fortran_build_root = _required_directory(
        arguments.fortran_build_root, "Fortran build root"
    )
    gfortran = _required_file(arguments.gfortran, "gfortran")
    generated_date = dt.date.fromisoformat(arguments.date).isoformat()

    module_directory = _required_directory(
        fortran_build_root / "fortran_prj" / "module",
        "Fortran module directory",
    )
    _required_file(module_directory / "sspmod.mod", "sspmod module")
    lib_misc = _required_file(
        fortran_build_root / "out" / "libMisc.a", "libMisc"
    )
    krakenc = _required_file(
        fortran_build_root / "out" / "krakenc.exe", "krakenc"
    )
    driver = _required_file(
        project_root / "for_test" / "fortran_ssp_oracle_driver.f90",
        "Fortran SSP oracle driver",
    )
    fixture_directory = _required_directory(
        project_root / "for_test" / "fixtures", "fixture directory"
    )
    env_paths = {
        profile_type: _required_file(
            fixture_directory / env_name, f"{profile_type} ENV"
        )
        for profile_type, env_name in PROFILE_ENV_NAMES.items()
    }
    lossy_env = _required_file(
        fixture_directory / LOSSY_ENV_NAME, "lossy ENV"
    )
    source_paths = {
        source_name: _required_file(
            fortran_source_root / "misc" / source_name,
            f"Fortran source {source_name}",
        )
        for source_name in SOURCE_NAMES
    }
    runtime_environment = _runtime_environment(gfortran)

    with tempfile.TemporaryDirectory(prefix="ookc-alg022-") as temp_name:
        compile_dir = Path(temp_name).resolve()
        oracle_exe = compile_dir / "fortran_ssp_oracle.exe"
        compile_command = [
            str(gfortran),
            "-I",
            str(module_directory),
            str(driver),
            str(lib_misc),
            "-o",
            str(oracle_exe),
        ]
        compile_completed = subprocess.run(
            compile_command,
            cwd=compile_dir,
            check=False,
            text=True,
            capture_output=True,
            timeout=180,
            env=runtime_environment,
        )
        if compile_completed.returncode != 0:
            raise RuntimeError(
                "Fortran SSP oracle compilation failed: "
                f"return code {compile_completed.returncode}; "
                f"stderr={compile_completed.stderr.strip()}"
            )
        if not oracle_exe.is_file():
            raise RuntimeError(
                "Fortran compiler returned success without producing "
                f"{oracle_exe}"
            )

        profiles: dict[str, Any] = {}
        profile_commands: dict[str, Any] = {}
        krakenc_commands: dict[str, Any] = {}
        for profile_type, env_path in env_paths.items():
            samples, profile_record = _run_profile(
                oracle_exe, profile_type, runtime_environment
            )
            wavenumbers, krakenc_record = _run_krakenc(
                krakenc, env_path, runtime_environment
            )
            profiles[profile_type] = {
                "samples": samples,
                "mode_count": len(wavenumbers),
                "wavenumbers": wavenumbers,
            }
            profile_commands[profile_type] = profile_record
            krakenc_commands[profile_type] = krakenc_record

        lossy_wavenumbers, lossy_command = _run_krakenc(
            krakenc, lossy_env, runtime_environment
        )
        krakenc_commands["lossy_gradient_interpolation"] = lossy_command

        input_paths = dict(env_paths)
        input_paths["lossy_gradient_interpolation"] = lossy_env
        payload = {
            "schema": "OpenOcean-Krakenc.ssp-interpolation-oracle",
            "schema_version": 1,
            "generated_date": generated_date,
            "provenance": {
                "krakenc": _file_record(krakenc),
                "libMisc": _file_record(lib_misc),
                "sources": {
                    name: _file_record(path)
                    for name, path in source_paths.items()
                },
                "inputs": {
                    name: _file_record(path)
                    for name, path in input_paths.items()
                },
                "driver": {
                    **_file_record(driver),
                    "compile_command": compile_command,
                    "compile_returncode": compile_completed.returncode,
                    "compile_cwd": str(compile_dir),
                },
            },
            "commands": {
                "profile_runs": profile_commands,
                "krakenc_runs": krakenc_commands,
            },
            "profiles": profiles,
            "lossy_gradient_interpolation": {
                "mode_count": len(lossy_wavenumbers),
                "wavenumbers": lossy_wavenumbers,
            },
        }
        _validate_payload(payload)
        _atomic_write(payload, arguments.output)

    print(f"wrote Fortran SSP oracle: {arguments.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
