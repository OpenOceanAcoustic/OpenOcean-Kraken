import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from for_test.propagation_loss_analysis import (
    assert_matching_grids,
    compute_tl,
    create_overview,
    failure_figures,
    plot_comparison,
    shade_to_cube,
    spatial_diagnostics,
    valid_tl_mask,
)
from for_test.shd_reader import read_shade_file


OUTPUT_MARKER = ".ookc_krakenc_propagation_loss_output"


@dataclass(frozen=True)
class CaseInput:
    name: str
    root: Path
    env: Path
    flp: Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict:
    path = Path(path)
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def discover_cases(test_dir: Path) -> list[CaseInput]:
    test_dir = Path(test_dir)
    cases = []
    seen_stems: dict[str, Path] = {}
    for env in sorted(
        test_dir.rglob("*.env"),
        key=lambda path: path.relative_to(test_dir).as_posix().casefold(),
    ):
        normalized_stem = env.stem.casefold()
        if normalized_stem in seen_stems:
            raise ValueError(
                "duplicate ENV stem would overwrite case output: "
                f"{seen_stems[normalized_stem]} and {env}"
            )
        seen_stems[normalized_stem] = env
        root = env.with_suffix("")
        flp = root.with_suffix(".flp")
        if not flp.is_file():
            raise FileNotFoundError(f"missing FLP for {env.name}: {flp}")
        cases.append(CaseInput(env.stem, root, env, flp))
    if not cases:
        raise FileNotFoundError(f"no ENV inputs in {test_dir}")
    return cases


def stage_case_inputs(
    case: CaseInput, destination: Path
) -> dict[str, dict]:
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    allowed = {".env", ".flp", ".sbp", ".brc", ".irc"}
    sources = {
        path.suffix.lower(): path
        for path in case.root.parent.glob(case.name + ".*")
        if path.is_file() and path.suffix.lower() in allowed
    }
    mapped_sources = {
        ("neggradK_brc", ".brc"): case.root.parent
        / "neggradC_brc.brc",
        ("neggradK_irc", ".irc"): case.root.parent
        / "neggradC_irc.irc",
    }
    for (case_name, suffix), source in mapped_sources.items():
        if case_name == case.name and suffix not in sources:
            if not source.is_file():
                raise FileNotFoundError(
                    f"missing mapped auxiliary input: {source}"
                )
            sources[suffix] = source
    records = {}
    for suffix, source in sorted(sources.items()):
        target = destination / f"{case.name}{suffix}"
        shutil.copy2(source, target)
        records[target.name] = {
            "source": source.name,
            "source_path": str(source.resolve()),
            **file_record(target),
        }
    return records


def stage_fortran_case_inputs(
    case: CaseInput, destination: Path
) -> dict[str, dict]:
    records = stage_case_inputs(case, destination)
    if not case.name.endswith("_brc"):
        return records
    destination = Path(destination)
    staged_env = destination / f"{case.name}.env"
    original = staged_env.read_text(encoding="utf-8")
    adapted, replacements = re.subn(
        r"(?m)^(\s*)'F'(\s+0\.0\s*!\s*BOTOPT)",
        r"\1'P'\2",
        original,
        count=1,
    )
    if replacements != 1:
        raise ValueError(
            f"unable to adapt BRC boundary option in {case.env.name}"
        )
    staged_env.write_text(adapted, encoding="utf-8")
    records[staged_env.name].update(file_record(staged_env))
    records[staged_env.name]["adaptation"] = "F boundary replaced by P"
    irc_source = case.root.parent / "neggradC_irc.irc"
    if not irc_source.is_file():
        raise FileNotFoundError(
            f"missing equivalent rigid IRC input: {irc_source}"
        )
    staged_irc = destination / f"{case.name}.irc"
    shutil.copy2(irc_source, staged_irc)
    records[staged_irc.name] = {
        "source": irc_source.name,
        "source_path": str(irc_source.resolve()),
        "adaptation": "unit-amplitude zero-phase BRC represented as rigid IRC",
        **file_record(staged_irc),
    }
    records["_adaptation"] = {
        "kind": "brc_to_irc",
        "reason": (
            "Fortran KrakenC secant search returns no modes for the tabulated "
            "unit-amplitude zero-phase BRC; the equivalent rigid IRC is used"
        ),
    }
    return records


def safe_reset_output(output_dir: Path, workspace: Path) -> None:
    output = Path(output_dir).resolve()
    workspace = Path(workspace).resolve()
    forbidden = {
        Path(output.anchor),
        workspace,
        (workspace / "test").resolve(),
    }
    try:
        output.relative_to(workspace)
    except ValueError as error:
        raise ValueError(f"unsafe output directory: {output}") from error
    if output in forbidden:
        raise ValueError(f"unsafe output directory: {output}")
    if output.exists():
        if not (output / OUTPUT_MARKER).is_file():
            raise ValueError(
                "refusing to replace output without ownership marker: "
                f"{output}"
            )
        shutil.rmtree(output)
    output.mkdir(parents=True)
    (output / OUTPUT_MARKER).write_text(
        "Generated by plot_propagation_loss_comparison.py\n",
        encoding="utf-8",
    )


def run_logged_command(
    command: list[Path | str],
    cwd: Path,
    label: str,
    timeout_seconds: int = 600,
) -> dict:
    cwd = Path(cwd)
    values = [str(value) for value in command]
    started = time.perf_counter()
    timed_out = False
    try:
        completed = subprocess.run(
            values,
            cwd=cwd,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
        returncode = completed.returncode
        stdout = completed.stdout
        stderr = completed.stderr
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
        "command": values,
        "cwd": str(cwd.resolve()),
        "returncode": returncode,
        "timed_out": timed_out,
        "seconds": elapsed,
        "stdout": file_record(stdout_path),
        "stderr": file_record(stderr_path),
    }


def process_cases(
    cases: Iterable, runner: Callable[[object], dict]
) -> list[dict]:
    records = []
    for case in cases:
        try:
            records.append(runner(case))
        except Exception as error:
            records.append(
                {
                    "case": getattr(case, "name", str(case)),
                    "status": "analysis_failed",
                    "error": str(error),
                }
            )
    return records


def expected_source_depths(case_name: str) -> list[float | None]:
    if case_name == "MunkKleaky":
        return [25.0, 250.0]
    return [None]


def _failure_record(
    record: dict,
    status: str,
    message: str,
    output: Path,
) -> dict:
    record["status"] = status
    record["error"] = message
    paths = failure_figures(
        record["case"],
        record["source_count"],
        message,
        output / "figures",
    )
    record["figures"].extend(file_record(path) for path in paths)
    return record


def run_case(
    case: CaseInput,
    ookc: Path,
    krakenc: Path,
    field: Path,
    output: Path,
    threads: int,
) -> dict:
    output = Path(output)
    case_directory = output / "runs" / case.name
    ookc_directory = case_directory / "ookc"
    krakenc_directory = case_directory / "krakenc"
    ookc_inputs = stage_case_inputs(case, ookc_directory)
    krakenc_inputs = stage_fortran_case_inputs(case, krakenc_directory)
    record = {
        "case": case.name,
        "status": "analysis_failed",
        "source_count": len(expected_source_depths(case.name)),
        "expected_source_depths_m": expected_source_depths(case.name),
        "executables": {
            "ookc": file_record(ookc) if Path(ookc).is_file() else {"path": str(ookc)},
            "krakenc": (
                file_record(krakenc)
                if Path(krakenc).is_file()
                else {"path": str(krakenc)}
            ),
            "field": (
                file_record(field)
                if Path(field).is_file()
                else {"path": str(field)}
            ),
        },
        "inputs": {"ookc": ookc_inputs, "krakenc": krakenc_inputs},
        "runs": {},
        "outputs": {},
        "metrics": [],
        "figures": [],
    }
    ookc_root = ookc_directory / case.name
    krakenc_root = krakenc_directory / case.name
    record["runs"]["ookc_mod"] = run_logged_command(
        [
            ookc,
            "--mod",
            ookc_root.with_suffix(".env"),
            ookc_root.with_suffix(".mod"),
            "--threads",
            str(threads),
        ],
        ookc_directory,
        "ookc_mod",
    )
    if record["runs"]["ookc_mod"]["returncode"] != 0:
        return _failure_record(
            record,
            "ookc_failed",
            "OOKc modal calculation failed; see runs/ookc/ookc_mod.stderr.log",
            output,
        )
    record["runs"]["ookc_field"] = run_logged_command(
        [ookc, "--field", ookc_root],
        ookc_directory,
        "ookc_field",
    )
    if record["runs"]["ookc_field"]["returncode"] != 0:
        return _failure_record(
            record,
            "ookc_failed",
            "OOKc field calculation failed; see runs/ookc/ookc_field.stderr.log",
            output,
        )
    record["runs"]["krakenc"] = run_logged_command(
        [krakenc, case.name],
        krakenc_directory,
        "krakenc",
    )
    if record["runs"]["krakenc"]["returncode"] != 0:
        return _failure_record(
            record,
            "krakenc_failed",
            "KrakenC modal calculation failed; see runs/krakenc/krakenc.stderr.log",
            output,
        )
    record["runs"]["field"] = run_logged_command(
        [field, case.name],
        krakenc_directory,
        "field",
    )
    if record["runs"]["field"]["returncode"] != 0:
        return _failure_record(
            record,
            "krakenc_failed",
            "Fortran Field calculation failed; see runs/krakenc/field.stderr.log",
            output,
        )
    ookc_shd = ookc_root.with_suffix(".shd")
    krakenc_shd = krakenc_root.with_suffix(".shd")
    try:
        ookc_cube = shade_to_cube(read_shade_file(ookc_shd))
        krakenc_cube = shade_to_cube(read_shade_file(krakenc_shd))
        assert_matching_grids(ookc_cube, krakenc_cube)
    except Exception as error:
        status = "grid_mismatch" if "mismatch" in str(error) else "analysis_failed"
        return _failure_record(record, status, str(error), output)
    record["source_count"] = len(ookc_cube.source_depths)
    record["expected_source_depths_m"] = [
        float(value) for value in ookc_cube.source_depths
    ]
    for source_index, source_depth in enumerate(ookc_cube.source_depths):
        depth_token = f"{source_depth:g}".replace(".", "p")
        figure = (
            output
            / "figures"
            / f"{case.name}_sd{depth_token}m.png"
        )
        metrics = plot_comparison(
            case.name,
            ookc_cube,
            krakenc_cube,
            source_index,
            figure,
        )
        ookc_plane = ookc_cube.pressure[source_index]
        krakenc_plane = krakenc_cube.pressure[source_index]
        diagnostics = spatial_diagnostics(
            ookc_plane, krakenc_plane, ookc_cube.ranges_metres
        )
        tl_difference = np.abs(
            compute_tl(ookc_plane) - compute_tl(krakenc_plane)
        )
        mask = valid_tl_mask(ookc_plane, krakenc_plane)
        ranked = np.where(mask, tl_difference, -np.inf)
        depth_index, range_index = np.unravel_index(
            np.argmax(ranked), ranked.shape
        )
        row = {
            "case": case.name,
            "frequency_hz": ookc_cube.frequency,
            "source_depth_m": float(source_depth),
            **metrics.to_dict(),
            **diagnostics.to_dict(),
            "max_difference_range_km": float(
                ookc_cube.ranges_metres[range_index] / 1000.0
            ),
            "max_difference_depth_m": float(
                ookc_cube.receiver_depths[depth_index]
            ),
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


METRIC_FIELDS = [
    "case",
    "frequency_hz",
    "source_depth_m",
    "complex_relative_l2",
    "positive_range_relative_l2",
    "interior_positive_range_relative_l2",
    "zero_range_reference_energy_fraction",
    "boundary_error_fraction",
    "tl_mae_db",
    "tl_rmse_db",
    "tl_p95_db",
    "tl_max_db",
    "valid_fraction",
    "max_difference_range_km",
    "max_difference_depth_m",
    "figure",
]


def write_metrics_csv(rows: list[dict], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=METRIC_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(
    case_records: list[dict], rows: list[dict], path: Path
) -> None:
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "case_count": len(case_records),
        "source_task_count": sum(
            record.get("source_count", 1) for record in case_records
        ),
        "successful_source_task_count": len(rows),
        "cases": case_records,
        "metrics": rows,
    }
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _agreement_description(relative_l2: float) -> str:
    if relative_l2 <= 2.0e-5:
        return "达到项目中距离无关场的高一致性量级"
    if relative_l2 <= 2.0e-4:
        return "达到项目中场计算差分测试的高一致性量级"
    if relative_l2 <= 5.0e-3:
        return "处于项目既有端到端比较容差范围内"
    return "高于项目既有端到端比较容差，需要重点检查模式相位和幅度"


def _physical_note(case_name: str) -> str:
    notes = {
        "MunkKleaky": (
            "深海 Munk 声速剖面在声道轴附近形成折射俘获；两个声源深度分别体现"
            "轴附近激发与浅层激发对远距离会聚和干涉结构的影响。"
        ),
        "elastic_fd_two_layer": (
            "水层下有限厚弹性沉积层支持纵波和横波，声学模态向海底泄漏，"
            "高阶干涉条纹通常衰减更快。"
        ),
        "multilayer_elastic_stack": (
            "分级弹性层连续改变底部阻抗与模态截止条件，多个固体界面共同塑造"
            "垂向条纹和底部泄漏。"
        ),
        "multilayer_mud_sand": (
            "软泥层优先吸收高掠射角能量，下伏较硬砂层增强反射，形成吸收与"
            "反射叠加的复合底损失。"
        ),
        "solve3_mode_gain": (
            "两个距离剖面的声速变化会触发模式增减和投影，剖面连接处的相位"
            "连续性决定下游干涉结构。"
        ),
        "stepK_rd": (
            "阶跃环境在剖面变化处触发模式耦合，变化点之后的传播条纹同时包含"
            "局部模态重构和累积相位效应。"
        ),
        "wedge": (
            "楔形波导深度随距离持续改变，模态逐步截止并交换能量，因而条纹"
            "弯曲并在浅化方向表现出更强的距离相关性。"
        ),
    }
    if case_name in notes:
        return notes[case_name]
    if case_name.endswith("_brc"):
        return (
            "负声速梯度使声能向深水方向折射；外部反射系数控制边界反射的幅度"
            "与相位，因此传播条纹对掠射角较敏感。"
        )
    if case_name.endswith("_irc"):
        return (
            "负声速梯度使声能向深水方向折射；内部反射系数直接参与模态边界"
            "条件，影响模态幅度、相位和干涉零点。"
        )
    return "该环境的传播结构由声速剖面、边界阻抗和模态干涉共同决定。"


def write_analysis_report(
    case_records: list[dict], rows: list[dict], path: Path
) -> None:
    failed = [
        record for record in case_records if record.get("status") != "success"
    ]
    lines = [
        "# OOKc 与 KrakenC 传播损失对比分析",
        "",
        (
            f"共处理 {len(case_records)} 个环境、"
            f"{sum(record.get('source_count', 1) for record in case_records)} 个声源任务；"
            f"成功获得 {len(rows)} 个可比较传播场，失败或不可比较环境 {len(failed)} 个。"
        ),
        "",
        "## 运行状态",
        "",
        "| 环境 | 状态 | 声源数 | 说明 |",
        "|---|---|---:|---|",
    ]
    for record in case_records:
        lines.append(
            f"| {record['case']} | {record.get('status', 'unknown')} | "
            f"{record.get('source_count', 1)} | {record.get('error', '')} |"
        )
    lines.extend(
        [
            "",
            "## 数值一致性",
            "",
            (
                "复压力相对 L2 使用全部网格，是首要一致性指标。TL 指标只统计"
                "相对两场联合峰值不低于 `1e-8` 的网格，以避免深衰落零点在 dB "
                "域无限放大。"
            ),
            "",
            "| 环境 | 声源深度 (m) | 全网格 L2 | 正距离 L2 | 内部正距离 L2 | TL RMSE (dB) | TL P95 (dB) |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sorted(rows, key=lambda value: value["complex_relative_l2"]):
        lines.append(
            f"| {row['case']} | {row['source_depth_m']:.3g} | "
            f"{row['complex_relative_l2']:.3e} | "
            f"{row['positive_range_relative_l2']:.3e} | "
            f"{row['interior_positive_range_relative_l2']:.3e} | "
            f"{row['tl_rmse_db']:.3f} | {row['tl_p95_db']:.3f} |"
        )
    def labels(values: list[dict]) -> str:
        return "、".join(
            f"{row['case']}(sd={row['source_depth_m']:.3g} m)"
            for row in values
        )

    boundary_dominated = [
        row
        for row in rows
        if row["interior_positive_range_relative_l2"] <= 2.0e-5
        and row["boundary_error_fraction"] >= 0.95
    ]
    zero_range_dominated = [
        row
        for row in rows
        if row["complex_relative_l2"] > 5.0e-3
        and row["positive_range_relative_l2"] <= 5.0e-3
        and row not in boundary_dominated
    ]
    substantive = [
        row
        for row in rows
        if row["interior_positive_range_relative_l2"] > 5.0e-3
    ]
    high_agreement = [
        row
        for row in rows
        if row["interior_positive_range_relative_l2"] <= 2.0e-5
        and row not in boundary_dominated
        and row not in zero_range_dominated
    ]
    acceptable = [
        row
        for row in rows
        if 2.0e-5 < row["interior_positive_range_relative_l2"] <= 5.0e-3
        and row not in zero_range_dominated
    ]
    lines.extend(
        [
            "",
            "## 总体结论",
            "",
            "总览图：[overview.png](overview.png)",
            "",
        ]
    )
    if high_agreement:
        lines.append(f"- 传播区域高度一致：{labels(high_agreement)}。")
    if acceptable:
        lines.append(
            f"- 传播区域总体一致、存在局部小差异：{labels(acceptable)}。"
        )
    if boundary_dominated:
        lines.append(
            "- 全网格指标由上下深度边界主导、内部场高度一致："
            f"{labels(boundary_dominated)}。"
        )
    if zero_range_dominated:
        lines.append(
            "- 全网格指标由零距离初始化列主导，正距离传播场仍一致："
            f"{labels(zero_range_dominated)}。"
        )
    if substantive:
        lines.append(
            "- 存在贯穿传播区域的实质差异，需要优先核查模态谱、"
            f"归一化和弹性层处理：{labels(substantive)}。"
        )
    lines.extend(["", "## 声传播物理特征与逐例分析", ""])
    for row in rows:
        lines.extend(
            [
                f"### {row['case']}，声源深度 {row['source_depth_m']:.3g} m",
                "",
                (
                    f"全网格复压力相对 L2 为 {row['complex_relative_l2']:.3e}，"
                    f"正距离 L2 为 {row['positive_range_relative_l2']:.3e}，"
                    f"去除上下边界后的内部正距离 L2 为 "
                    f"{row['interior_positive_range_relative_l2']:.3e}。"
                    f"TL RMSE 为 {row['tl_rmse_db']:.3f} dB，"
                    f"TL P95 为 {row['tl_p95_db']:.3f} dB，"
                    f"{_agreement_description(row['interior_positive_range_relative_l2'])}。"
                    f"最大有效 TL 差异位于 {row['max_difference_range_km']:.3f} km、"
                    f"{row['max_difference_depth_m']:.3f} m。"
                ),
                "",
                (
                    f"参考场在零距离列的能量占比为 "
                    f"{row['zero_range_reference_energy_fraction']:.1%}；"
                    f"总误差中位于上下深度边界的能量占比为 "
                    f"{row['boundary_error_fraction']:.1%}。"
                ),
                "",
                f"对应图：[{row['figure']}]({row['figure']})",
                "",
                _physical_note(row["case"]),
                "",
            ]
        )
    lines.extend(
        [
            "## 方法与限制",
            "",
            "传播损失定义为 `TL = -20 log10(max(|p|, 1e-12))` dB。",
            (
                "OOKc 与 KrakenC 图使用同一距离/深度网格和联合色标；网格不一致时"
                "不插值，而是标记为不可比较。局部 TL 差值较大但复压力相对 L2 很小"
                "时，通常表示差异集中在干涉零点附近，而不是全场能量系统性偏离。"
            ),
            (
                "对于 `*_brc`，原 BRC 的反射幅度恒为 1、相位恒为 0。Fortran "
                "KrakenC 的 BRC 复根搜索会返回无模式，因此其隔离运行副本采用项目"
                "验收路径相同的等价刚性 IRC 表示；OOKc 仍使用原始 BRC，且该转换"
                "不会修改 `test` 输入。"
            ),
            "",
        ]
    )
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def _validate_executable(parser: argparse.ArgumentParser, path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_file():
        parser.error(f"executable does not exist: {resolved}")
    return resolved


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare OOKc and Fortran KrakenC propagation loss"
    )
    parser.add_argument("--test-dir", required=True, type=Path)
    parser.add_argument("--ookc", required=True, type=Path)
    parser.add_argument("--krakenc", required=True, type=Path)
    parser.add_argument("--field", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    test_dir = args.test_dir.resolve()
    if not test_dir.is_dir():
        parser.error(f"test directory does not exist: {test_dir}")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    ookc = _validate_executable(parser, args.ookc)
    krakenc = _validate_executable(parser, args.krakenc)
    field = _validate_executable(parser, args.field)
    workspace = test_dir.parent
    output = args.output.resolve()
    safe_reset_output(output, workspace)
    (output / "figures").mkdir()
    cases = discover_cases(test_dir)
    records = process_cases(
        cases,
        lambda case: run_case(
            case, ookc, krakenc, field, output, args.threads
        ),
    )
    for record in records:
        if "source_count" not in record:
            record["source_count"] = len(
                expected_source_depths(record["case"])
            )
        if not record.get("figures"):
            paths = failure_figures(
                record["case"],
                record["source_count"],
                record.get("error", "unhandled analysis failure"),
                output / "figures",
            )
            record["figures"] = [file_record(path) for path in paths]
    rows = [
        row for record in records for row in record.get("metrics", [])
    ]
    write_metrics_csv(rows, output / "metrics.csv")
    write_analysis_report(records, rows, output / "analysis.md")
    figures = sorted((output / "figures").glob("*.png"))
    create_overview(figures, output / "overview.png")
    write_manifest(records, rows, output / "manifest.json")
    passes = all(record.get("status") == "success" for record in records)
    summary = {
        "passes": passes,
        "case_count": len(records),
        "source_task_count": sum(
            record.get("source_count", 1) for record in records
        ),
        "successful_source_task_count": len(rows),
        "output": str(output),
    }
    print(json.dumps(summary, ensure_ascii=False))
    return 0 if passes else 1


if __name__ == "__main__":
    raise SystemExit(main())
