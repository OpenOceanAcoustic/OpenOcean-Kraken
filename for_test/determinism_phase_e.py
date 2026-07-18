import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


DEFAULT_CASES = ("MunkKleaky", "solve3_mode_gain", "stepK_rd", "wedge")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--cases", nargs="*", default=list(DEFAULT_CASES))
    args = parser.parse_args()
    if args.repetitions < 2:
        parser.error("repetitions must be at least two")

    workspace = args.workspace.resolve()
    binary = args.binary.resolve()
    report = {
        "binary": str(binary),
        "binary_sha256": sha256(binary),
        "thread_counts": [1, 4],
        "repetitions": args.repetitions,
        "cases": {},
        "passes": False,
    }

    for case in args.cases:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            copied = []
            for source in (workspace / "test").glob(case + ".*"):
                if source.suffix.lower() in {".env", ".flp", ".sbp", ".brc", ".irc"}:
                    shutil.copy2(source, directory / source.name)
                    copied.append(source.name)
            env = directory / f"{case}.env"
            if not env.is_file():
                raise RuntimeError(f"missing test ENV for {case}")
            case_record = {"inputs": copied, "runs": {}, "passes": False}
            all_hashes = []
            for threads in (1, 4):
                records = []
                for repetition in range(args.repetitions):
                    output = directory / f"{case}_t{threads}_r{repetition}.mod"
                    completed = subprocess.run(
                        [str(binary), "--mod", str(env), str(output),
                         "--threads", str(threads)],
                        cwd=directory,
                        capture_output=True,
                        text=True,
                        timeout=600,
                        check=False,
                    )
                    if completed.returncode != 0 or not output.is_file():
                        raise RuntimeError(
                            f"{case} threads={threads} repetition={repetition} failed: "
                            f"{completed.stderr}")
                    digest = sha256(output)
                    all_hashes.append(digest)
                    records.append({
                        "repetition": repetition,
                        "sha256": digest,
                        "bytes": output.stat().st_size,
                        "cli": json.loads(completed.stdout),
                    })
                case_record["runs"][str(threads)] = records
            case_record["passes"] = len(set(all_hashes)) == 1
            report["cases"][case] = case_record

    report["passes"] = all(item["passes"] for item in report["cases"].values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "passes": report["passes"],
        "case_count": len(report["cases"]),
        "runs": len(report["cases"]) * 2 * args.repetitions,
    }))
    return 0 if report["passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
