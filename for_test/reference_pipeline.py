import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path


class ReferencePipelineError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path) -> dict:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": _sha256(path)}


def _reset_artifact_dir(path: Path) -> None:
    resolved = path.resolve()
    if resolved == Path(resolved.anchor) or len(resolved.parts) < 3:
        raise ReferencePipelineError(f"unsafe artifact directory: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)


def _run(executable: Path, case_name: str, cwd: Path, label: str) -> dict:
    command = [str(executable), case_name]
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=cwd,
            capture_output=True,
            timeout=180,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise ReferencePipelineError(f"{label} timed out after 180 seconds") from error
    elapsed = time.perf_counter() - started
    log_path = cwd / f"{label}.log"
    log_path.write_bytes(completed.stdout + b"\n--- STDERR ---\n" + completed.stderr)
    return {
        "command": command,
        "returncode": completed.returncode,
        "seconds": elapsed,
        "log": _file_record(log_path),
    }


def _write_manifest(path: Path, manifest: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    temporary.replace(path)


def run_reference_pipeline(
    case_root: Path, krakenc: Path, field: Path, artifact_dir: Path,
    run_field: bool = True,
) -> dict:
    case_root = case_root.resolve()
    krakenc = krakenc.resolve()
    field = field.resolve()
    artifact_dir = artifact_dir.resolve()
    required_inputs = [case_root.with_suffix(".env"), case_root.with_suffix(".flp")]
    missing = [str(path) for path in required_inputs if not path.is_file()]
    if missing:
        raise ReferencePipelineError("missing required input: " + ", ".join(missing))
    for executable in (krakenc, field):
        if not executable.is_file():
            raise ReferencePipelineError(f"missing executable: {executable}")

    _reset_artifact_dir(artifact_dir)
    copied_inputs = []
    allowed_suffixes = {".env", ".flp", ".sbp", ".brc", ".irc"}
    for source in case_root.parent.glob(case_root.name + ".*"):
        if source.is_file() and source.suffix.lower() in allowed_suffixes:
            target = artifact_dir / source.name
            shutil.copy2(source, target)
            copied_inputs.append(target)

    manifest = {
        "passes": False,
        "case_root": str(case_root),
        "executables": {
            "krakenc": _file_record(krakenc),
            "field": _file_record(field),
        },
        "inputs": {path.name: _file_record(path) for path in copied_inputs},
        "runs": {},
        "outputs": {},
    }

    case_name = case_root.name
    krakenc_run = _run(krakenc, case_name, artifact_dir, "krakenc")
    manifest["runs"]["krakenc"] = krakenc_run
    mod_path = artifact_dir / f"{case_name}.mod"
    if (
        krakenc_run["returncode"] != 0
        or not mod_path.is_file()
        or mod_path.stat().st_size == 0
    ):
        _write_manifest(artifact_dir / "manifest.json", manifest)
        raise ReferencePipelineError("krakenc failed to produce a nonempty MOD file")
    manifest["outputs"]["mod"] = _file_record(mod_path)

    if not run_field:
        manifest["passes"] = True
        _write_manifest(artifact_dir / "manifest.json", manifest)
        return manifest

    field_run = _run(field, case_name, artifact_dir, "field")
    manifest["runs"]["field"] = field_run
    shd_path = artifact_dir / f"{case_name}.shd"
    if (
        field_run["returncode"] != 0
        or not shd_path.is_file()
        or shd_path.stat().st_size == 0
    ):
        _write_manifest(artifact_dir / "manifest.json", manifest)
        raise ReferencePipelineError("field failed to produce a nonempty SHD file")
    manifest["outputs"]["shd"] = _file_record(shd_path)
    manifest["passes"] = True
    _write_manifest(artifact_dir / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-root", required=True, type=Path)
    parser.add_argument("--krakenc", required=True, type=Path)
    parser.add_argument("--field", required=True, type=Path)
    parser.add_argument("--artifact-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = run_reference_pipeline(
            args.case_root, args.krakenc, args.field, args.artifact_dir
        )
    except ReferencePipelineError as error:
        print(str(error))
        return 1
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
