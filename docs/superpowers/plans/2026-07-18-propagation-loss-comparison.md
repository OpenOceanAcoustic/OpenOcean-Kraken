# OOKc and KrakenC Propagation-Loss Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run OOKc and Fortran KrakenC/Field for every `test/*.env`, produce one four-panel propagation-loss comparison per source depth, and deliver reproducible metrics plus a Chinese physical/numerical analysis report.

**Architecture:** Keep numerical analysis pure and independently tested in `propagation_loss_analysis.py`; keep filesystem discovery, solver execution, manifests, and reporting in one CLI orchestrator. Reuse the existing SHD reader, stage every case in isolated output directories, reject grid mismatches instead of interpolating, and continue the batch after per-case failures.

**Tech Stack:** Python 3, `unittest`, NumPy, Matplotlib, existing C++17 OOKc Release executable, Fortran `krakenc.exe` and `field.exe`.

## Global Constraints

- Discover all 11 `test/*.env` inputs dynamically; do not edit any file under `test`.
- Render both source depths for `MunkKleaky`, yielding 12 source-depth figures when all calculations succeed.
- Use `TL = -20 log10(max(|p|, 1e-12))` and complex-pressure relative L2 as the primary agreement metric.
- Compute TL metrics only where either model is within `1e-8` of the joint field peak; report mask coverage.
- OOKc and KrakenC panels for a source depth must share axes and TL color limits.
- A failed case must not stop later cases; retain command, return code, elapsed time, stdout, stderr, and status.
- Preserve the user's pre-existing `.gitignore` modification.

---

### Task 1: Pure SHD reshaping and numerical metrics

**Files:**
- Create: `for_test/propagation_loss_analysis.py`
- Create: `for_test/test_propagation_loss_analysis.py`

**Interfaces:**
- Consumes: `for_test.shd_reader.ShadeFile`
- Produces: `PressureCube`, `ComparisonMetrics`, `shade_to_cube`, `compute_tl`, `valid_tl_mask`, `compute_metrics`, `assert_matching_grids`

- [ ] **Step 1: Write failing tests for pressure order, TL, masking, metrics, and grid rejection**

```python
import unittest
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from for_test.propagation_loss_analysis import (
    assert_matching_grids,
    compute_metrics,
    compute_tl,
    shade_to_cube,
    valid_tl_mask,
)
from for_test.shd_reader import ShadeFile


class PropagationLossAnalysisTests(unittest.TestCase):
    def test_shade_to_cube_preserves_source_depth_receiver_range_order(self):
        shade = ShadeFile(
            "synthetic", 50.0, [25.0, 250.0], [10.0, 20.0],
            [0.0, 1000.0, 2000.0],
            [complex(value, 0.0) for value in range(1, 13)],
        )
        cube = shade_to_cube(shade)
        self.assertEqual(cube.pressure.shape, (2, 2, 3))
        np.testing.assert_array_equal(cube.pressure[1, 0], [7, 8, 9])

    def test_compute_tl_applies_pressure_floor(self):
        actual = compute_tl(np.array([1.0, 0.1, 0.0]))
        np.testing.assert_allclose(actual, [0.0, 20.0, 240.0])

    def test_valid_tl_mask_uses_joint_peak(self):
        ookc = np.array([1.0, 1.0e-9, 0.0])
        krakenc = np.array([0.5, 0.0, 1.0e-7])
        np.testing.assert_array_equal(
            valid_tl_mask(ookc, krakenc), [True, False, True]
        )

    def test_compute_metrics_reports_complex_l2_and_db_statistics(self):
        reference = np.array([1.0 + 0.0j, 0.1 + 0.0j])
        actual = np.array([1.0 + 0.0j, 0.01 + 0.0j])
        metrics = compute_metrics(actual, reference)
        self.assertAlmostEqual(metrics.complex_relative_l2, 0.09 / np.sqrt(1.01))
        self.assertAlmostEqual(metrics.tl_mae_db, 10.0)
        self.assertAlmostEqual(metrics.tl_rmse_db, np.sqrt(200.0))
        self.assertAlmostEqual(metrics.tl_p95_db, 19.0)
        self.assertAlmostEqual(metrics.tl_max_db, 20.0)
        self.assertEqual(metrics.valid_fraction, 1.0)

    def test_assert_matching_grids_rejects_receiver_depth_difference(self):
        left = ShadeFile("a", 50.0, [25.0], [10.0], [1000.0], [1 + 0j])
        right = ShadeFile("b", 50.0, [25.0], [11.0], [1000.0], [1 + 0j])
        with self.assertRaisesRegex(ValueError, "receiver-depth grid mismatch"):
            assert_matching_grids(shade_to_cube(left), shade_to_cube(right))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run tests and verify RED**

Run: `python -m unittest for_test.test_propagation_loss_analysis -v`

Expected: FAIL with `ModuleNotFoundError: No module named 'for_test.propagation_loss_analysis'`.

- [ ] **Step 3: Implement the pure numerical API**

```python
from dataclasses import asdict, dataclass

import numpy as np

from for_test.shd_reader import ShadeFile

PRESSURE_FLOOR = 1.0e-12
VALID_RELATIVE_AMPLITUDE = 1.0e-8


@dataclass(frozen=True)
class PressureCube:
    title: str
    frequency: float
    source_depths: np.ndarray
    receiver_depths: np.ndarray
    ranges_metres: np.ndarray
    pressure: np.ndarray


@dataclass(frozen=True)
class ComparisonMetrics:
    complex_relative_l2: float
    tl_mae_db: float
    tl_rmse_db: float
    tl_p95_db: float
    tl_max_db: float
    valid_fraction: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


def shade_to_cube(shade: ShadeFile) -> PressureCube:
    shape = (len(shade.source_depths), len(shade.receiver_depths),
             len(shade.ranges_metres))
    pressure = np.asarray(shade.pressure, dtype=np.complex128)
    if pressure.size != int(np.prod(shape)):
        raise ValueError(f"pressure size {pressure.size} does not match grid {shape}")
    return PressureCube(
        shade.title, float(shade.frequency),
        np.asarray(shade.source_depths, dtype=float),
        np.asarray(shade.receiver_depths, dtype=float),
        np.asarray(shade.ranges_metres, dtype=float),
        pressure.reshape(shape),
    )


def compute_tl(pressure: np.ndarray) -> np.ndarray:
    return -20.0 * np.log10(np.maximum(np.abs(pressure), PRESSURE_FLOOR))


def valid_tl_mask(ookc: np.ndarray, krakenc: np.ndarray) -> np.ndarray:
    joint = np.maximum(np.abs(ookc), np.abs(krakenc))
    peak = float(np.max(joint, initial=0.0))
    return joint >= peak * VALID_RELATIVE_AMPLITUDE if peak > 0.0 else np.ones(joint.shape, bool)


def compute_metrics(ookc: np.ndarray, krakenc: np.ndarray) -> ComparisonMetrics:
    difference = ookc - krakenc
    denominator = float(np.linalg.norm(krakenc.ravel()))
    relative_l2 = float(np.linalg.norm(difference.ravel()) / denominator) if denominator else float("inf")
    mask = valid_tl_mask(ookc, krakenc)
    db_error = np.abs(compute_tl(ookc)[mask] - compute_tl(krakenc)[mask])
    return ComparisonMetrics(
        relative_l2,
        float(np.mean(db_error)),
        float(np.sqrt(np.mean(np.square(db_error)))),
        float(np.percentile(db_error, 95.0)),
        float(np.max(db_error)),
        float(np.mean(mask)),
    )


def assert_matching_grids(left: PressureCube, right: PressureCube) -> None:
    if not np.isclose(left.frequency, right.frequency):
        raise ValueError("frequency mismatch")
    for label, a, b in (
        ("source-depth", left.source_depths, right.source_depths),
        ("receiver-depth", left.receiver_depths, right.receiver_depths),
        ("range", left.ranges_metres, right.ranges_metres),
    ):
        if a.shape != b.shape or not np.allclose(a, b, rtol=0.0, atol=1.0e-8):
            raise ValueError(f"{label} grid mismatch")
    if left.pressure.shape != right.pressure.shape:
        raise ValueError("pressure grid mismatch")
```

- [ ] **Step 4: Run tests and verify GREEN**

Run: `python -m unittest for_test.test_propagation_loss_analysis -v`

Expected: 5 tests, all `ok`.

- [ ] **Step 5: Commit**

```powershell
git add for_test/propagation_loss_analysis.py for_test/test_propagation_loss_analysis.py
git commit -m "test: define propagation loss comparison metrics"
```

---

### Task 2: Four-panel plotting and overview

**Files:**
- Modify: `for_test/propagation_loss_analysis.py`
- Modify: `for_test/test_propagation_loss_analysis.py`

**Interfaces:**
- Consumes: two matching `PressureCube` objects and a source index
- Produces: `comparison_color_limits`, `nearest_receiver_index`, `plot_comparison`, `create_overview`

- [ ] **Step 1: Add failing tests for color limits, receiver selection, PNG output, and overview**

```python
    def test_comparison_color_limits_share_tl_range_and_keep_difference_nonzero(self):
        ookc = np.array([1.0, 0.1])
        krakenc = np.array([1.0, 0.1])
        tl_min, tl_max, difference_half_width = comparison_color_limits(ookc, krakenc)
        self.assertLess(tl_min, tl_max)
        self.assertEqual(difference_half_width, 0.1)

    def test_nearest_receiver_index_selects_source_depth(self):
        self.assertEqual(nearest_receiver_index(np.array([0.0, 20.0, 40.0]), 31.0), 2)

    def test_plot_comparison_and_overview_write_nonempty_png_files(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            cube = synthetic_cube()
            figure = directory / "comparison.png"
            plot_comparison("synthetic", cube, cube, 0, figure)
            self.assertGreater(figure.stat().st_size, 10_000)
            overview = directory / "overview.png"
            create_overview([figure], overview)
            self.assertGreater(overview.stat().st_size, 10_000)
```

Add imports for `tempfile`, the four plotting functions, and a `synthetic_cube()` fixture with one source, three receiver depths, four positive ranges, and nonzero complex pressure.

- [ ] **Step 2: Run focused tests and verify RED**

Run: `python -m unittest for_test.test_propagation_loss_analysis.PropagationLossAnalysisTests.test_plot_comparison_and_overview_write_nonempty_png_files -v`

Expected: FAIL because `plot_comparison` is not defined.

- [ ] **Step 3: Implement plotting with the noninteractive Agg backend**

Implement these exact signatures in `propagation_loss_analysis.py`:

```python
def comparison_color_limits(ookc: np.ndarray, krakenc: np.ndarray) -> tuple[float, float, float]:
    ookc_tl, krakenc_tl = compute_tl(ookc), compute_tl(krakenc)
    combined = np.concatenate([ookc_tl.ravel(), krakenc_tl.ravel()])
    finite = combined[np.isfinite(combined)]
    low, high = np.percentile(finite, [2.0, 98.0])
    if not high > low:
        high = low + 1.0
    absolute_difference = np.abs(ookc_tl - krakenc_tl)
    difference_half_width = max(0.1, float(np.percentile(absolute_difference, 98.0)))
    return float(low), float(high), difference_half_width


def nearest_receiver_index(receiver_depths: np.ndarray, source_depth: float) -> int:
    return int(np.argmin(np.abs(receiver_depths - source_depth)))


def plot_comparison(case_name: str, ookc: PressureCube, krakenc: PressureCube,
                    source_index: int, output_path: Path) -> ComparisonMetrics:
    assert_matching_grids(ookc, krakenc)
    ookc_pressure = ookc.pressure[source_index]
    krakenc_pressure = krakenc.pressure[source_index]
    ookc_tl = compute_tl(ookc_pressure)
    krakenc_tl = compute_tl(krakenc_pressure)
    difference = ookc_tl - krakenc_tl
    metrics = compute_metrics(ookc_pressure, krakenc_pressure)
    tl_min, tl_max, difference_half_width = comparison_color_limits(
        ookc_pressure, krakenc_pressure
    )
    ranges_km = ookc.ranges_metres / 1000.0
    depths = ookc.receiver_depths
    source_depth = float(ookc.source_depths[source_index])
    figure, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    for axis, values, title in (
        (axes[0, 0], ookc_tl, "OOKc transmission loss"),
        (axes[0, 1], krakenc_tl, "KrakenC transmission loss"),
    ):
        mesh = axis.pcolormesh(ranges_km, depths, values, shading="auto",
                               cmap="viridis", vmin=tl_min, vmax=tl_max)
        axis.invert_yaxis()
        axis.set(xlabel="Range (km)", ylabel="Receiver depth (m)", title=title)
        figure.colorbar(mesh, ax=axis, label="TL (dB)")
    difference_norm = TwoSlopeNorm(vmin=-difference_half_width, vcenter=0.0,
                                   vmax=difference_half_width)
    difference_mesh = axes[1, 0].pcolormesh(
        ranges_km, depths, difference, shading="auto", cmap="RdBu_r",
        norm=difference_norm
    )
    axes[1, 0].invert_yaxis()
    axes[1, 0].set(xlabel="Range (km)", ylabel="Receiver depth (m)",
                   title="OOKc - KrakenC TL difference")
    figure.colorbar(difference_mesh, ax=axes[1, 0], label="Difference (dB)")
    receiver_index = nearest_receiver_index(depths, source_depth)
    axes[1, 1].plot(ranges_km, ookc_tl[receiver_index], label="OOKc", linewidth=1.2)
    axes[1, 1].plot(ranges_km, krakenc_tl[receiver_index], label="KrakenC",
                    linewidth=1.0, linestyle="--")
    axes[1, 1].invert_yaxis()
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()
    axes[1, 1].set(
        xlabel="Range (km)", ylabel="TL (dB)",
        title=f"Range cut at receiver depth {depths[receiver_index]:.1f} m",
    )
    axes[1, 1].text(
        0.02, 0.02,
        f"relative L2 = {metrics.complex_relative_l2:.3e}\n"
        f"TL RMSE = {metrics.tl_rmse_db:.3f} dB\n"
        f"TL P95 = {metrics.tl_p95_db:.3f} dB",
        transform=axes[1, 1].transAxes, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"},
    )
    figure.suptitle(
        f"{case_name} | {ookc.frequency:g} Hz | source depth {source_depth:g} m"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return metrics


def create_overview(figure_paths: list[Path], output_path: Path) -> None:
    if not figure_paths:
        raise ValueError("at least one figure is required")
    columns = min(3, len(figure_paths))
    rows = (len(figure_paths) + columns - 1) // columns
    figure, axes = plt.subplots(rows, columns, figsize=(6 * columns, 4 * rows),
                                squeeze=False, constrained_layout=True)
    for axis, path in zip(axes.ravel(), figure_paths):
        axis.imshow(plt.imread(path))
        axis.set_title(path.stem, fontsize=9)
        axis.axis("off")
    for axis in axes.ravel()[len(figure_paths):]:
        axis.axis("off")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)


def failure_figures(case_name: str, source_count: int, message: str,
                    figures_directory: Path) -> list[Path]:
    figures_directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for source_index in range(source_count):
        suffix = f"source{source_index + 1}"
        path = figures_directory / f"{case_name}_{suffix}_failed.png"
        figure, axis = plt.subplots(figsize=(12, 8), constrained_layout=True)
        axis.axis("off")
        axis.text(0.5, 0.55, f"{case_name}\ncomparison unavailable",
                  ha="center", va="center", fontsize=22)
        axis.text(0.5, 0.42, message, ha="center", va="center", fontsize=13,
                  wrap=True)
        figure.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(figure)
        paths.append(path)
    return paths
```

The implementation must import `matplotlib`, call `matplotlib.use("Agg")` before importing `pyplot`, and type paths with `pathlib.Path`.

- [ ] **Step 4: Run all analysis tests and verify GREEN**

Run: `python -m unittest for_test.test_propagation_loss_analysis -v`

Expected: 8 tests, all `ok`, with no Matplotlib GUI window.

- [ ] **Step 5: Commit**

```powershell
git add for_test/propagation_loss_analysis.py for_test/test_propagation_loss_analysis.py
git commit -m "feat: render propagation loss comparison panels"
```

---

### Task 3: Case discovery, staging, command execution, and failure continuation

**Files:**
- Create: `for_test/plot_propagation_loss_comparison.py`
- Modify: `for_test/test_propagation_loss_analysis.py`

**Interfaces:**
- Produces: `CaseInput`, `discover_cases`, `stage_case_inputs`, `run_logged_command`, `safe_reset_output`, `run_case`
- `run_case(...) -> dict` always returns a status record; expected calculation failures are represented in the record, not raised out of the batch loop.

- [ ] **Step 1: Add failing tests for 11-case discovery, auxiliary mapping, safe output reset, and continued processing**

```python
class BatchOrchestrationTests(unittest.TestCase):
    def test_discover_cases_finds_all_workspace_env_files(self):
        workspace = Path(__file__).resolve().parents[1].parent
        cases = discover_cases(workspace / "test")
        self.assertEqual(len(cases), 11)
        self.assertEqual(sum(case.name == "MunkKleaky" for case in cases), 1)

    def test_stage_inputs_maps_missing_neggradk_reflection_files(self):
        workspace = Path(__file__).resolve().parents[1].parent
        case = next(case for case in discover_cases(workspace / "test")
                    if case.name == "neggradK_brc")
        with tempfile.TemporaryDirectory() as directory:
            record = stage_case_inputs(case, Path(directory))
            self.assertTrue((Path(directory) / "neggradK_brc.brc").is_file())
            self.assertEqual(record["neggradK_brc.brc"]["source"], "neggradC_brc.brc")

    def test_safe_reset_output_rejects_test_directory(self):
        workspace = Path(__file__).resolve().parents[1].parent
        with self.assertRaisesRegex(ValueError, "unsafe output directory"):
            safe_reset_output(workspace / "test", workspace)

    def test_batch_loop_records_failure_and_continues(self):
        records = process_cases(["first", "second"], lambda name: (
            {"case": name, "status": "failed"} if name == "first"
            else {"case": name, "status": "success"}
        ))
        self.assertEqual([record["status"] for record in records], ["failed", "success"])
```

- [ ] **Step 2: Run the new class and verify RED**

Run: `python -m unittest for_test.test_propagation_loss_analysis.BatchOrchestrationTests -v`

Expected: FAIL because the CLI module does not exist.

- [ ] **Step 3: Implement discovery, staging, safe reset, and logged execution**

Implement these public shapes:

```python
@dataclass(frozen=True)
class CaseInput:
    name: str
    root: Path
    env: Path
    flp: Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict:
    return {
        "path": str(path.resolve()), "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def discover_cases(test_dir: Path) -> list[CaseInput]:
    cases = []
    for env in sorted(test_dir.glob("*.env"), key=lambda path: path.name.lower()):
        root = env.with_suffix("")
        flp = root.with_suffix(".flp")
        if not flp.is_file():
            raise FileNotFoundError(f"missing FLP for {env.name}: {flp}")
        cases.append(CaseInput(env.stem, root, env, flp))
    if not cases:
        raise FileNotFoundError(f"no ENV inputs in {test_dir}")
    return cases


def stage_case_inputs(case: CaseInput, destination: Path) -> dict[str, dict]:
    destination.mkdir(parents=True, exist_ok=True)
    allowed = {".env", ".flp", ".sbp", ".brc", ".irc"}
    sources = {
        path.suffix.lower(): path
        for path in case.root.parent.glob(case.name + ".*")
        if path.is_file() and path.suffix.lower() in allowed
    }
    mapping = {
        ("neggradK_brc", ".brc"): case.root.parent / "neggradC_brc.brc",
        ("neggradK_irc", ".irc"): case.root.parent / "neggradC_irc.irc",
    }
    for key, source in mapping.items():
        if key[0] == case.name and key[1] not in sources:
            if not source.is_file():
                raise FileNotFoundError(f"missing mapped auxiliary input: {source}")
            sources[key[1]] = source
    records = {}
    for suffix, source in sorted(sources.items()):
        target = destination / f"{case.name}{suffix}"
        shutil.copy2(source, target)
        records[target.name] = {
            "source": source.name,
            "source_path": str(source.resolve()),
            "path": str(target.resolve()),
            "bytes": target.stat().st_size,
            "sha256": sha256(target),
        }
    return records


def safe_reset_output(output_dir: Path, workspace: Path) -> None:
    output = output_dir.resolve()
    workspace = workspace.resolve()
    forbidden = {Path(output.anchor), workspace, workspace / "test"}
    try:
        output.relative_to(workspace.parent)
    except ValueError as error:
        raise ValueError(f"unsafe output directory: {output}") from error
    if output in forbidden:
        raise ValueError(f"unsafe output directory: {output}")
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)


def run_logged_command(command: list[Path | str], cwd: Path, label: str,
                       timeout_seconds: int = 600) -> dict:
    values = [str(value) for value in command]
    started = time.perf_counter()
    timed_out = False
    try:
        completed = subprocess.run(values, cwd=cwd, capture_output=True,
                                   timeout=timeout_seconds, check=False)
        returncode = completed.returncode
        stdout, stderr = completed.stdout, completed.stderr
    except subprocess.TimeoutExpired as error:
        timed_out = True
        returncode = None
        stdout = error.stdout or b""
        stderr = error.stderr or b""
    elapsed = time.perf_counter() - started
    stdout_path = cwd / f"{label}.stdout.log"
    stderr_path = cwd / f"{label}.stderr.log"
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)
    return {
        "command": values, "cwd": str(cwd.resolve()),
        "returncode": returncode, "timed_out": timed_out, "seconds": elapsed,
        "stdout": file_record(stdout_path), "stderr": file_record(stderr_path),
    }


def process_cases(cases, runner) -> list[dict]:
    records = []
    for case in cases:
        try:
            records.append(runner(case))
        except Exception as error:
            records.append({"case": getattr(case, "name", str(case)),
                            "status": "analysis_failed", "error": str(error)})
    return records
```

- [ ] **Step 4: Run orchestration and analysis tests and verify GREEN**

Run: `python -m unittest for_test.test_propagation_loss_analysis -v`

Expected: 12 tests, all `ok`.

- [ ] **Step 5: Commit**

```powershell
git add for_test/plot_propagation_loss_comparison.py for_test/test_propagation_loss_analysis.py
git commit -m "feat: orchestrate isolated OOKc and KrakenC cases"
```

---

### Task 4: End-to-end case comparison, CSV, manifest, and Chinese analysis

**Files:**
- Modify: `for_test/plot_propagation_loss_comparison.py`
- Modify: `for_test/test_propagation_loss_analysis.py`

**Interfaces:**
- Consumes: executable paths and staged inputs
- Produces: `run_case`, `write_metrics_csv`, `write_manifest`, `write_analysis_report`, `main`

- [ ] **Step 1: Add failing report-structure and CLI-validation tests**

```python
    def test_report_writers_emit_required_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            rows = [{
                "case": "synthetic", "source_depth_m": 25.0,
                "complex_relative_l2": 0.0, "tl_mae_db": 0.0,
                "tl_rmse_db": 0.0, "tl_p95_db": 0.0,
                "tl_max_db": 0.0, "valid_fraction": 1.0,
                "max_difference_range_km": 1.0,
                "max_difference_depth_m": 20.0,
                "figure": "figures/synthetic_sd25m.png",
            }]
            cases = [{"case": "synthetic", "status": "success",
                      "source_count": 1}]
            write_metrics_csv(rows, output / "metrics.csv")
            write_manifest(cases, rows, output / "manifest.json")
            write_analysis_report(cases, rows, output / "analysis.md")
            self.assertIn("complex_relative_l2", (output / "metrics.csv").read_text())
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(manifest["case_count"], 1)
            self.assertIn("数值一致性", (output / "analysis.md").read_text(encoding="utf-8"))
            self.assertIn("声传播物理特征", (output / "analysis.md").read_text(encoding="utf-8"))
```

- [ ] **Step 2: Run the report test and verify RED**

Run: `python -m unittest for_test.test_propagation_loss_analysis.BatchOrchestrationTests.test_report_writers_emit_required_artifacts -v`

Expected: FAIL because the report writer functions are not defined.

- [ ] **Step 3: Implement end-to-end orchestration and reports**

`run_case` must:

1. stage identical inputs under `runs/<case>/ookc` and `runs/<case>/krakenc`;
2. run OOKc `--mod <env> <mod> --threads N`, then OOKc `--field <root>`;
3. run Fortran `krakenc.exe <case>`, then `field.exe <case>`;
4. classify the first failed stage without raising;
5. read both SHD files, assert grids match, and loop all source indices;
6. call `plot_comparison`, locate maximum valid absolute TL difference, and append one metrics row per source depth;
7. return a case record with source count, all command records, staged-input provenance, MOD/SHD/figure file records, and status.

Use this concrete implementation shape (with `failure_figures` rendering a white PNG containing the failure stage and retained-log path for each expected source depth):

```python
def expected_source_count(case_name: str) -> int:
    return 2 if case_name == "MunkKleaky" else 1


def run_case(case: CaseInput, ookc: Path, krakenc: Path, field: Path,
             output: Path, threads: int) -> dict:
    case_directory = output / "runs" / case.name
    ookc_directory = case_directory / "ookc"
    krakenc_directory = case_directory / "krakenc"
    ookc_inputs = stage_case_inputs(case, ookc_directory)
    krakenc_inputs = stage_case_inputs(case, krakenc_directory)
    record = {
        "case": case.name, "status": "analysis_failed",
        "source_count": expected_source_count(case.name),
        "inputs": {"ookc": ookc_inputs, "krakenc": krakenc_inputs},
        "runs": {}, "outputs": {}, "metrics": [], "figures": [],
    }
    ookc_root = ookc_directory / case.name
    krakenc_root = krakenc_directory / case.name
    record["runs"]["ookc_mod"] = run_logged_command(
        [ookc, "--mod", ookc_root.with_suffix(".env"),
         ookc_root.with_suffix(".mod"), "--threads", str(threads)],
        ookc_directory, "ookc_mod",
    )
    if record["runs"]["ookc_mod"]["returncode"] != 0:
        record["status"] = "ookc_failed"
        record["figures"].extend(file_record(path) for path in failure_figures(
            case.name, record["source_count"], "OOKc MOD failed", output / "figures"))
        return record
    record["runs"]["ookc_field"] = run_logged_command(
        [ookc, "--field", ookc_root], ookc_directory, "ookc_field"
    )
    if record["runs"]["ookc_field"]["returncode"] != 0:
        record["status"] = "ookc_failed"
        record["figures"].extend(file_record(path) for path in failure_figures(
            case.name, record["source_count"], "OOKc Field failed", output / "figures"))
        return record
    record["runs"]["krakenc"] = run_logged_command(
        [krakenc, case.name], krakenc_directory, "krakenc"
    )
    if record["runs"]["krakenc"]["returncode"] != 0:
        record["status"] = "krakenc_failed"
        record["figures"].extend(file_record(path) for path in failure_figures(
            case.name, record["source_count"], "KrakenC MOD failed", output / "figures"))
        return record
    record["runs"]["field"] = run_logged_command(
        [field, case.name], krakenc_directory, "field"
    )
    if record["runs"]["field"]["returncode"] != 0:
        record["status"] = "krakenc_failed"
        record["figures"].extend(file_record(path) for path in failure_figures(
            case.name, record["source_count"], "KrakenC Field failed", output / "figures"))
        return record
    ookc_shd = ookc_root.with_suffix(".shd")
    krakenc_shd = krakenc_root.with_suffix(".shd")
    try:
        ookc_cube = shade_to_cube(read_shade_file(ookc_shd))
        krakenc_cube = shade_to_cube(read_shade_file(krakenc_shd))
        assert_matching_grids(ookc_cube, krakenc_cube)
    except ValueError as error:
        record["status"] = "grid_mismatch"
        record["error"] = str(error)
        record["figures"].extend(file_record(path) for path in failure_figures(
            case.name, record["source_count"], str(error), output / "figures"))
        return record
    record["source_count"] = len(ookc_cube.source_depths)
    for source_index, source_depth in enumerate(ookc_cube.source_depths):
        depth_token = f"{source_depth:g}".replace(".", "p")
        figure = output / "figures" / f"{case.name}_sd{depth_token}m.png"
        metrics = plot_comparison(case.name, ookc_cube, krakenc_cube,
                                  source_index, figure)
        ookc_plane = ookc_cube.pressure[source_index]
        krakenc_plane = krakenc_cube.pressure[source_index]
        difference = np.abs(compute_tl(ookc_plane) - compute_tl(krakenc_plane))
        mask = valid_tl_mask(ookc_plane, krakenc_plane)
        ranked = np.where(mask, difference, -np.inf)
        depth_index, range_index = np.unravel_index(np.argmax(ranked), ranked.shape)
        row = {
            "case": case.name, "frequency_hz": ookc_cube.frequency,
            "source_depth_m": float(source_depth), **metrics.to_dict(),
            "max_difference_range_km": float(ookc_cube.ranges_metres[range_index] / 1000.0),
            "max_difference_depth_m": float(ookc_cube.receiver_depths[depth_index]),
            "figure": str(figure.relative_to(output)).replace("\\", "/"),
        }
        record["metrics"].append(row)
        record["figures"].append(file_record(figure))
    record["outputs"] = {
        "ookc_mod": file_record(ookc_root.with_suffix(".mod")),
        "ookc_shd": file_record(ookc_shd),
        "krakenc_mod": file_record(krakenc_root.with_suffix(".mod")),
        "krakenc_shd": file_record(krakenc_shd),
    }
    record["status"] = "success"
    return record
```

The report functions must use these signatures:

```python
def write_metrics_csv(rows: list[dict], path: Path) -> None:
    fields = [
        "case", "frequency_hz", "source_depth_m", "complex_relative_l2",
        "tl_mae_db", "tl_rmse_db", "tl_p95_db", "tl_max_db",
        "valid_fraction", "max_difference_range_km",
        "max_difference_depth_m", "figure",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(case_records: list[dict], rows: list[dict], path: Path) -> None:
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "case_count": len(case_records),
        "source_task_count": sum(record.get("source_count", 1)
                                 for record in case_records),
        "successful_source_task_count": len(rows),
        "cases": case_records,
        "metrics": rows,
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                    encoding="utf-8")


def write_analysis_report(case_records: list[dict], rows: list[dict], path: Path) -> None:
    status_lines = ["| 环境 | 状态 | 声源数 |", "|---|---:|---:|"]
    for record in case_records:
        status_lines.append(
            f"| {record['case']} | {record['status']} | {record.get('source_count', 1)} |"
        )
    metric_lines = [
        "| 环境 | 声源深度 (m) | 复压力相对 L2 | TL RMSE (dB) | TL P95 (dB) |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda value: value["complex_relative_l2"]):
        metric_lines.append(
            f"| {row['case']} | {row['source_depth_m']:.3g} | "
            f"{row['complex_relative_l2']:.3e} | {row['tl_rmse_db']:.3f} | "
            f"{row['tl_p95_db']:.3f} |"
        )
    physical_notes = {
        "MunkKleaky": "深海声道使声能在声速极小值附近反复折射；两个声源深度用于观察轴附近激发与浅层激发的差异。",
        "elastic_fd_two_layer": "水层下的有限厚弹性层允许纵波、横波和声学模态耦合，底部泄漏使高阶条纹更快衰减。",
        "multilayer_elastic_stack": "分级弹性层改变底部阻抗与模态截止，垂向条纹由多个固体界面共同塑造。",
        "multilayer_mud_sand": "软泥层先吸收并弱化高掠射角能量，下面较硬砂层产生更强反射，形成复合底损失。",
        "solve3_mode_gain": "两个距离剖面的声速改变会引起模式增减与投影，重点观察连接处相位和幅度连续性。",
        "stepK_rd": "阶跃地形在剖面变化处触发模式耦合，差异若沿变化距离扩散通常与耦合相位有关。",
        "wedge": "楔形波导深度持续改变，使模态逐步截止并交换能量，传播条纹随距离发生明显弯曲。",
    }
    detail_lines = []
    for row in rows:
        name = row["case"]
        note = physical_notes.get(
            name,
            "负声速梯度把声线向深水方向折射；外部或内部反射系数文件控制边界幅度和相位。",
        )
        detail_lines.extend([
            f"### {name}，声源深度 {row['source_depth_m']:.3g} m",
            "",
            f"复压力相对 L2 为 {row['complex_relative_l2']:.3e}；TL RMSE 为 "
            f"{row['tl_rmse_db']:.3f} dB。最大有效 TL 差异位于 "
            f"{row['max_difference_range_km']:.3f} km、"
            f"{row['max_difference_depth_m']:.3f} m。{note}",
            "",
        ])
    failed = [record for record in case_records if record["status"] != "success"]
    summary = (
        f"共处理 {len(case_records)} 个环境，获得 {len(rows)} 个成功声源场；"
        f"{len(failed)} 个环境存在失败或不可比较状态。"
    )
    content = "\n".join([
        "# OOKc 与 KrakenC 传播损失对比分析", "", summary, "",
        "## 运行状态", "", *status_lines, "", "## 数值一致性", "",
        "复压力相对 L2 使用全网格；TL 统计仅使用相对联合峰值不低于 `1e-8` 的网格。",
        "", *metric_lines, "", "## 声传播物理特征", "", *detail_lines,
        "## 方法与限制", "",
        "传播损失定义为 `-20 log10(max(|p|, 1e-12))` dB。深衰落零点附近的 dB 差值可能很大，因此以复压力相对 L2 为首要一致性指标。",
        "",
    ])
    path.write_text(content, encoding="utf-8")
```

`main` must parse `--test-dir`, `--ookc`, `--krakenc`, `--field`, `--output`, and `--threads`; validate every executable and directory before clearing output; write partial manifest/report even when some cases fail; create `overview.png`; print a one-line JSON summary; and return nonzero if any case lacks a successful comparison.

```python
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-dir", required=True, type=Path)
    parser.add_argument("--ookc", required=True, type=Path)
    parser.add_argument("--krakenc", required=True, type=Path)
    parser.add_argument("--field", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    test_dir = args.test_dir.resolve()
    executables = [args.ookc.resolve(), args.krakenc.resolve(), args.field.resolve()]
    if not test_dir.is_dir():
        parser.error(f"test directory does not exist: {test_dir}")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    for executable in executables:
        if not executable.is_file():
            parser.error(f"executable does not exist: {executable}")
    workspace = test_dir.parent
    output = args.output.resolve()
    safe_reset_output(output, workspace)
    (output / "figures").mkdir()
    cases = discover_cases(test_dir)
    records = process_cases(
        cases,
        lambda case: run_case(case, executables[0], executables[1],
                              executables[2], output, args.threads),
    )
    rows = [row for record in records for row in record.get("metrics", [])]
    write_metrics_csv(rows, output / "metrics.csv")
    write_manifest(records, rows, output / "manifest.json")
    write_analysis_report(records, rows, output / "analysis.md")
    figures = sorted((output / "figures").glob("*.png"))
    create_overview(figures, output / "overview.png")
    passes = all(record["status"] == "success" for record in records)
    print(json.dumps({
        "passes": passes, "case_count": len(records),
        "source_task_count": sum(record.get("source_count", 1) for record in records),
        "successful_source_task_count": len(rows), "output": str(output),
    }, ensure_ascii=False))
    return 0 if passes else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run all tests and verify GREEN**

Run: `python -m unittest for_test.test_propagation_loss_analysis -v`

Expected: 13 tests, all `ok`.

- [ ] **Step 5: Commit**

```powershell
git add for_test/plot_propagation_loss_comparison.py for_test/test_propagation_loss_analysis.py
git commit -m "feat: report full propagation loss comparison"
```

---

### Task 5: Smoke run, full 11-environment run, visual QA, and final evidence

**Files:**
- Generate: `E:/my_project/02_ookc_project/propagation_loss_comparison/**`
- Modify only if a reproduced defect requires it: the three new Python files from Tasks 1–4

**Interfaces:**
- Consumes: tested CLI plus frozen executable paths
- Produces: 12 figures/status tasks, overview, CSV, manifest, Chinese analysis, retained solver artifacts and logs

- [ ] **Step 1: Verify plotting dependencies and executable identities**

Run:

```powershell
python -c "import numpy, matplotlib; print(numpy.__version__, matplotlib.__version__)"
Get-FileHash build-release/OpenOcean-Krakenc.exe -Algorithm SHA256
Get-FileHash ../krakenFortran/build-acceptance-mingw/out/krakenc.exe -Algorithm SHA256
Get-FileHash ../krakenFortran/build-acceptance-mingw/out/field.exe -Algorithm SHA256
```

Expected: imports succeed and all three hashes are printed.

- [ ] **Step 2: Run a MunkKleaky smoke comparison in a temporary one-case input directory**

Copy only `MunkKleaky.env`, `.flp`, and `.sbp` into `build/propagation-smoke-input`, then run the CLI with output `build/propagation-smoke-output`.

Expected: exit 0, two figures, two CSV rows, one case with `status=success`.

- [ ] **Step 3: Run the complete batch**

Run:

```powershell
python for_test/plot_propagation_loss_comparison.py `
  --test-dir ../test `
  --ookc build-release/OpenOcean-Krakenc.exe `
  --krakenc ../krakenFortran/build-acceptance-mingw/out/krakenc.exe `
  --field ../krakenFortran/build-acceptance-mingw/out/field.exe `
  --output ../propagation_loss_comparison `
  --threads 4
```

Expected: 11 cases recorded and 12 source-depth figures/status tasks. If the CLI exits nonzero, inspect retained logs, write a failing regression test for any code defect, and follow RED-GREEN before rerunning. Do not hide a solver limitation; retain it as an explicit failed status and produce its placeholder/evidence.

- [ ] **Step 4: Programmatically verify the fresh output**

Run a read-only Python check that asserts:

```python
manifest["case_count"] == 11
manifest["source_task_count"] == 12
len(list((output / "figures").glob("*.png"))) == 12
(output / "overview.png").stat().st_size > 10_000
len(list(csv.DictReader((output / "metrics.csv").open()))) == manifest["successful_source_task_count"]
all(figure.stat().st_size > 10_000 for figure in figures)
```

- [ ] **Step 5: Render and inspect the overview plus representative full-resolution plots**

Use local image inspection for:

- `overview.png`;
- both `MunkKleaky` source depths;
- one elastic multilayer case;
- one reflection-coefficient case;
- `stepK_rd` and `wedge`.

Check titles, axes, inverted depth, common TL scales, centered difference scale, line panel, missing-data artifacts, and readable annotations. Correct visual defects with a failing structural test where practical, rerun all tests, then rerun the full batch.

- [ ] **Step 6: Fresh final verification**

Run:

```powershell
python -m unittest for_test.test_propagation_loss_analysis -v
git diff --check
git status --short
```

Read the complete test output, manifest, `metrics.csv`, and `analysis.md`. Verify each design requirement against the fresh artifacts before reporting completion.

- [ ] **Step 7: Commit final verified implementation**

```powershell
git add for_test/propagation_loss_analysis.py for_test/plot_propagation_loss_comparison.py for_test/test_propagation_loss_analysis.py
git commit -m "feat: compare OOKc and KrakenC propagation loss"
```
