import argparse
import hashlib
import json
import shutil
import statistics
import subprocess
import tempfile
import threading
import time
from pathlib import Path

import psutil


CASES = (
    "MunkKleaky",
    "elastic_fd_two_layer",
    "multilayer_elastic_stack",
    "multilayer_mud_sand",
    "neggradC_brc",
    "neggradC_irc",
    "stepK_rd",
    "wedge",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def measured_run(command, cwd, timeout=600):
    started = time.perf_counter()
    process = subprocess.Popen(
        [str(value) for value in command], cwd=cwd,
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    peak_rss = [0]
    ps_process = psutil.Process(process.pid)
    stop_sampling = threading.Event()

    def sample_memory():
        while not stop_sampling.is_set():
            try:
                peak_rss[0] = max(peak_rss[0], ps_process.memory_info().rss)
            except psutil.Error:
                return
            stop_sampling.wait(0.005)

    sampler = threading.Thread(target=sample_memory, daemon=True)
    sampler.start()
    try:
        try:
            _, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            _, stderr = process.communicate()
            raise
    finally:
        stop_sampling.set()
        sampler.join()
    elapsed = time.perf_counter() - started
    if process.returncode != 0:
        raise RuntimeError(
            f"command failed ({process.returncode}): {command}\n"
            f"{stderr.decode(errors='replace')}")
    return elapsed, peak_rss[0]


def benchmark(command_factory, cwd, warmups, repetitions):
    for iteration in range(warmups):
        measured_run(command_factory(f"warmup_{iteration}"), cwd)
    samples = []
    peaks = []
    for iteration in range(repetitions):
        seconds, peak = measured_run(command_factory(f"sample_{iteration}"), cwd)
        samples.append(seconds)
        peaks.append(peak)
    return {
        "seconds": samples,
        "median_seconds": statistics.median(samples),
        "peak_rss_bytes": max(peaks, default=0),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--fortran", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--heavy-min-fortran-seconds", type=float, default=0.5)
    parser.add_argument("--performance-min-fortran-seconds", type=float, default=0.1)
    args = parser.parse_args()
    if (args.threads < 1 or args.warmups < 0 or args.repetitions < 1 or
            args.heavy_min_fortran_seconds <= 0 or
            args.performance_min_fortran_seconds <= 0):
        parser.error("threads/repetitions/thresholds must be positive")

    workspace = args.workspace.resolve()
    cpp = args.cpp.resolve()
    fortran = args.fortran.resolve()
    report = {
        "rules": {
            "build": "Release",
            "cpp_threads": args.threads,
            "fortran_threads": 1,
            "warmups": args.warmups,
            "repetitions": args.repetitions,
            "clock": "time.perf_counter",
            "memory": "peak child-process RSS sampled every 5 ms",
        },
        "executables": {
            "cpp": {"path": str(cpp), "sha256": sha256(cpp)},
            "fortran": {"path": str(fortran), "sha256": sha256(fortran)},
        },
        "cases": {},
    }

    for case in CASES:
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            inputs = {}
            for source in (workspace / "test").glob(case + ".*"):
                if source.suffix.lower() in {".env", ".flp", ".sbp", ".brc", ".irc"}:
                    target = directory / source.name
                    shutil.copy2(source, target)
                    inputs[source.name] = sha256(source)

            env = directory / f"{case}.env"

            def cpp_command(label):
                output = directory / f"cpp_{label}.mod"
                output.unlink(missing_ok=True)
                return [cpp, "--mod", env, output,
                        "--threads", str(args.threads)]

            def fortran_command(label):
                del label
                (directory / f"{case}.mod").unlink(missing_ok=True)
                return [fortran, case]

            cpp_result = benchmark(
                cpp_command, directory, args.warmups, args.repetitions)
            fortran_result = benchmark(
                fortran_command, directory, args.warmups, args.repetitions)
            report["cases"][case] = {
                "inputs": inputs,
                "cpp": cpp_result,
                "fortran": fortran_result,
                "time_ratio_cpp_over_fortran":
                    cpp_result["median_seconds"] /
                    fortran_result["median_seconds"],
            }

    ratios = [entry["time_ratio_cpp_over_fortran"]
              for entry in report["cases"].values()]
    heavy_ratios = [entry["time_ratio_cpp_over_fortran"]
                    for entry in report["cases"].values()
                    if entry["fortran"]["median_seconds"] >=
                    args.heavy_min_fortran_seconds]
    performance_ratios = [entry["time_ratio_cpp_over_fortran"]
                          for entry in report["cases"].values()
                          if entry["fortran"]["median_seconds"] >=
                          args.performance_min_fortran_seconds]
    total_cpp_seconds = sum(entry["cpp"]["median_seconds"]
                            for entry in report["cases"].values())
    total_fortran_seconds = sum(entry["fortran"]["median_seconds"]
                                for entry in report["cases"].values())
    total_ratio = total_cpp_seconds / total_fortran_seconds
    report["summary"] = {
        "median_time_ratio_cpp_over_fortran": statistics.median(ratios),
        "maximum_time_ratio_cpp_over_fortran": max(ratios),
        "heavy_min_fortran_seconds": args.heavy_min_fortran_seconds,
        "maximum_heavy_time_ratio_cpp_over_fortran": max(heavy_ratios),
        "performance_min_fortran_seconds":
            args.performance_min_fortran_seconds,
        "performance_case_count": len(performance_ratios),
        "performance_median_time_ratio_cpp_over_fortran":
            statistics.median(performance_ratios),
        "total_cpp_median_seconds": total_cpp_seconds,
        "total_fortran_median_seconds": total_fortran_seconds,
        "total_time_ratio_cpp_over_fortran": total_ratio,
        "passes_single_thread_correctness_gate":
            args.threads != 1 or total_ratio <= 1.5,
        "passes_optimized_median_gate":
            statistics.median(performance_ratios) <= 1.0,
        "passes_optimized_heavy_gate": max(heavy_ratios) <= 1.25,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"]))


if __name__ == "__main__":
    main()
