import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path


class PipelineError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": sha256(path)}


def run(command, cwd: Path, label: str) -> dict:
    started = time.perf_counter()
    completed = subprocess.run(
        [str(value) for value in command], cwd=cwd, capture_output=True,
        timeout=600, check=False)
    elapsed = time.perf_counter() - started
    log = cwd / f"{label}.log"
    log.write_bytes(completed.stdout + b"\n--- STDERR ---\n" + completed.stderr)
    if completed.returncode != 0:
        raise PipelineError(
            f"{label} failed ({completed.returncode}); see {log}")
    return {"command": [str(value) for value in command],
            "returncode": completed.returncode, "seconds": elapsed,
            "log": file_record(log)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-root", required=True, type=Path)
    parser.add_argument("--cpp", required=True, type=Path)
    parser.add_argument("--artifact-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--repeat", action="store_true")
    args = parser.parse_args()

    case_root = args.case_root.resolve()
    executable = args.cpp.resolve()
    artifact_dir = args.artifact_dir.resolve()
    if args.threads < 1 or not executable.is_file():
        raise PipelineError("invalid executable or thread count")
    if artifact_dir == Path(artifact_dir.anchor) or len(artifact_dir.parts) < 3:
        raise PipelineError(f"unsafe artifact directory: {artifact_dir}")
    if artifact_dir.exists():
        shutil.rmtree(artifact_dir)
    artifact_dir.mkdir(parents=True)

    inputs = []
    allowed = {".env", ".flp", ".sbp", ".brc", ".irc"}
    for source in case_root.parent.glob(case_root.name + ".*"):
        if source.is_file() and source.suffix.lower() in allowed:
            target = artifact_dir / source.name
            shutil.copy2(source, target)
            inputs.append(target)
    env = artifact_dir / f"{case_root.name}.env"
    flp = artifact_dir / f"{case_root.name}.flp"
    if not env.is_file() or not flp.is_file():
        raise PipelineError("ENV and FLP inputs are required")

    manifest = {
        "passes": False,
        "case_root": str(case_root),
        "threads": args.threads,
        "executable": file_record(executable),
        "inputs": {path.name: file_record(path) for path in inputs},
        "runs": {}, "outputs": {}, "determinism": {},
    }
    names = [case_root.name, case_root.name + "_repeat"] if args.repeat else [case_root.name]
    try:
        for index, name in enumerate(names):
            root = artifact_dir / name
            if index > 0:
                for source in inputs:
                    shutil.copy2(source, root.with_suffix(source.suffix))
            manifest["runs"][f"mod_{index}"] = run(
                [executable, "--mod", env, root.with_suffix(".mod"),
                 "--threads", str(args.threads)], artifact_dir, f"mod_{index}")
            manifest["runs"][f"field_{index}"] = run(
                [executable, "--field", root], artifact_dir, f"field_{index}")
            manifest["outputs"][f"mod_{index}"] = file_record(root.with_suffix(".mod"))
            manifest["outputs"][f"shd_{index}"] = file_record(root.with_suffix(".shd"))
        if args.repeat:
            manifest["determinism"] = {
                "mod_byte_identical": manifest["outputs"]["mod_0"]["sha256"] ==
                                      manifest["outputs"]["mod_1"]["sha256"],
                "shd_byte_identical": manifest["outputs"]["shd_0"]["sha256"] ==
                                      manifest["outputs"]["shd_1"]["sha256"],
            }
            if not all(manifest["determinism"].values()):
                raise PipelineError("repeated outputs are not byte-identical")
        manifest["passes"] = True
    finally:
        (artifact_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"passes": True, "manifest": str(artifact_dir / "manifest.json")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
