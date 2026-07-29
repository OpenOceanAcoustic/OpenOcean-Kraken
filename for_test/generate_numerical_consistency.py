"""Generate OOKc/KrakenC numerical-consistency acceptance artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from for_test.numerical_consistency import (
    DEFAULT_THRESHOLDS,
    classify_slice,
)


GRADE_LABELS = {"pass": "通过", "watch": "关注", "fail": "不通过"}
ANOMALY_LABELS = {
    "zero_range": "r=0 零距离",
    "boundary": "上下边界",
    "interior_positive_range": "内部正距离声场",
    "local_interference_spike": "干涉零点/局部尖峰",
    "none": "未见显著异常区域",
}
PLOT_ANOMALY_LABELS = {
    "zero_range": "r=0 zero range",
    "boundary": "top/bottom boundary",
    "interior_positive_range": "interior positive range",
    "local_interference_spike": "local interference spike",
    "none": "no significant anomaly",
}


def build_records(metrics: list[dict]) -> list[dict]:
    records = []
    for metric in metrics:
        result = classify_slice(metric)
        records.append(
            {
                **metric,
                "grade": result.grade,
                "grade_zh": GRADE_LABELS[result.grade],
                "anomaly_regions": ";".join(result.anomaly_regions),
                "anomaly_regions_zh": "；".join(
                    ANOMALY_LABELS[value] for value in result.anomaly_regions
                ),
                "rationale": result.rationale,
            }
        )
    return records


def summary_payload(records: list[dict]) -> dict:
    counts = Counter(record["grade"] for record in records)
    grade_counts = {
        "pass": counts.get("pass", 0),
        "watch": counts.get("watch", 0),
        "fail": counts.get("fail", 0),
    }
    anomaly_counts = Counter(
        anomaly
        for record in records
        for anomaly in record["anomaly_regions"].split(";")
    )
    if grade_counts["fail"]:
        overall = "fail"
    elif grade_counts["watch"]:
        overall = "watch"
    else:
        overall = "pass"
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "thresholds": asdict(DEFAULT_THRESHOLDS),
        "case_count": len({record["case"] for record in records}),
        "slice_count": len(records),
        "grade_counts": grade_counts,
        "pass_rate": (
            grade_counts["pass"] / len(records) if records else 0.0
        ),
        "overall_grade": overall,
        "overall_grade_zh": GRADE_LABELS[overall],
        "anomaly_counts": dict(sorted(anomaly_counts.items())),
    }


def _write_csv(records: list[dict], path: Path) -> None:
    fields = list(records[0]) if records else ["case"]
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def _write_json(records: list[dict], summary: dict, path: Path) -> None:
    path.write_text(
        json.dumps({**summary, "records": records}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _record_labels(records: list[dict]) -> list[str]:
    return [
        f"{record['case']}\n{float(record['source_depth_m']):g}m"
        for record in records
    ]


def _write_plot(records: list[dict], summary: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    labels = _record_labels(records)
    x = np.arange(len(records))
    colours = {
        "pass": "#54a24b",
        "watch": "#f2cf5b",
        "fail": "#e45756",
    }
    point_colours = [colours[record["grade"]] for record in records]

    l2 = np.asarray(
        [record["interior_positive_range_relative_l2"] for record in records],
        dtype=float,
    )
    axes[0, 0].scatter(x, np.maximum(l2, 1.0e-16), c=point_colours, s=45)
    axes[0, 0].axhline(
        DEFAULT_THRESHOLDS.pass_interior_relative_l2,
        color="#54a24b",
        ls="--",
        label="pass limit 1e-3",
    )
    axes[0, 0].axhline(
        DEFAULT_THRESHOLDS.watch_interior_relative_l2,
        color="#e45756",
        ls=":",
        label="watch limit 5e-3",
    )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylabel("Interior positive-range relative L2")
    axes[0, 0].set_title("Complex-pressure consistency")
    axes[0, 0].legend()

    p95 = np.asarray([record["tl_p95_db"] for record in records], dtype=float)
    axes[0, 1].scatter(x, np.maximum(p95, 1.0e-12), c=point_colours, s=45)
    axes[0, 1].axhline(
        DEFAULT_THRESHOLDS.pass_tl_p95_db,
        color="#54a24b",
        ls="--",
        label="pass limit 0.05 dB",
    )
    axes[0, 1].axhline(
        DEFAULT_THRESHOLDS.watch_tl_p95_db,
        color="#e45756",
        ls=":",
        label="watch limit 0.5 dB",
    )
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_ylabel("TL P95 absolute difference (dB)")
    axes[0, 1].set_title("Transmission-loss consistency")
    axes[0, 1].legend()

    grade_names = ["pass", "watch", "fail"]
    axes[1, 0].bar(
        grade_names,
        [summary["grade_counts"][name] for name in grade_names],
        color=[colours[name] for name in grade_names],
    )
    axes[1, 0].set_ylabel("Source slices")
    axes[1, 0].set_title(
        f"Acceptance: {summary['overall_grade'].upper()} "
        f"(pass rate {summary['pass_rate']:.1%})"
    )

    anomaly_names = [
        "zero_range",
        "boundary",
        "interior_positive_range",
        "local_interference_spike",
        "none",
    ]
    axes[1, 1].barh(
        [PLOT_ANOMALY_LABELS[name] for name in anomaly_names],
        [summary["anomaly_counts"].get(name, 0) for name in anomaly_names],
        color="#4c78a8",
    )
    axes[1, 1].set_xlabel("Source slices (multi-label)")
    axes[1, 1].set_title("Attributed anomaly regions")

    for axis in axes[0, :]:
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
        axis.grid(alpha=0.25)
    for axis in axes[1, :]:
        axis.grid(axis="x" if axis is axes[1, 1] else "y", alpha=0.25)
    fig.suptitle("OOKc vs KrakenC numerical-consistency acceptance")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _write_markdown(records: list[dict], summary: dict, path: Path) -> None:
    lines = [
        "# OOKc 与 KrakenC 数值一致性验收",
        "",
        (
            f"覆盖 {summary['case_count']} 个环境、{summary['slice_count']} 个声源切片。"
            f"总体判定：**{summary['overall_grade_zh']}**；通过率 "
            f"**{summary['pass_rate']:.1%}**。"
        ),
        "",
        "## 验收判据",
        "",
        "- 通过：内部正距离相对 L2 ≤ `1e-3`，且 TL P95 ≤ `0.05 dB`。",
        "- 关注：未通过，但相对 L2 ≤ `5e-3`，且 TL P95 ≤ `0.5 dB`。",
        "- 不通过：超过关注阈值。异常区域标签用于解释误差来源，不改变等级。",
        "",
        "## 通过情况",
        "",
        f"- 通过：{summary['grade_counts']['pass']} 个切片。",
        f"- 关注：{summary['grade_counts']['watch']} 个切片。",
        f"- 不通过：{summary['grade_counts']['fail']} 个切片。",
        "",
        "## 逐切片误差与异常区域",
        "",
        "| 环境 | 声源深度(m) | 内部正距离相对L2 | TL RMSE(dB) | TL P95(dB) | 最大差位置(km,m) | 判定 | 异常区域 |",
        "|---|---:|---:|---:|---:|---|---|---|",
    ]
    for record in records:
        lines.append(
            f"| {record['case']} | {float(record['source_depth_m']):.3g} | "
            f"{float(record['interior_positive_range_relative_l2']):.3e} | "
            f"{float(record['tl_rmse_db']):.6f} | "
            f"{float(record['tl_p95_db']):.6f} | "
            f"{float(record['max_difference_range_km']):.3f}, "
            f"{float(record['max_difference_depth_m']):.3f} | "
            f"{record['grade_zh']} | {record['anomaly_regions_zh']} |"
        )

    non_pass = [record for record in records if record["grade"] != "pass"]
    lines.extend(["", "## 关注与不通过项", ""])
    if non_pass:
        for record in non_pass:
            lines.append(
                f"- {record['case']}(sd={float(record['source_depth_m']):g} m)："
                f"{record['grade_zh']}；{record['rationale']}；"
                f"异常区域：{record['anomaly_regions_zh']}。"
            )
    else:
        lines.append("- 无。全部切片满足“通过”阈值。")

    lines.extend(["", "## 异常区域说明", ""])
    for anomaly, label in ANOMALY_LABELS.items():
        affected = [
            record
            for record in records
            if anomaly in record["anomaly_regions"].split(";")
        ]
        if not affected:
            continue
        examples = "、".join(
            f"{record['case']}(sd={float(record['source_depth_m']):g} m)"
            for record in affected
        )
        lines.append(f"- {label}（{len(affected)} 个切片）：{examples}。")
    lines.extend(
        [
            "",
            "## 结论",
            "",
            (
                "等级由内部正距离复声压相对 L2 和 TL P95 共同决定；"
                "`r=0`、边界和局部尖峰均单独标注，因此不会用少量特殊网格点"
                "掩盖主体传播场的一致性。"
            ),
            "",
            "汇总图见 [numerical_consistency.png](numerical_consistency.png)。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_artifacts(manifest_path: Path, output_dir: Path) -> list[Path]:
    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    metrics = payload.get("metrics")
    if not isinstance(metrics, list) or not metrics:
        raise ValueError("manifest contains no propagation metrics")
    records = build_records(metrics)
    summary = summary_payload(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        output_dir / "numerical_consistency.csv",
        output_dir / "numerical_consistency.json",
        output_dir / "numerical_consistency.png",
        output_dir / "numerical_consistency.md",
    ]
    _write_csv(records, paths[0])
    _write_json(records, summary, paths[1])
    _write_plot(records, summary, paths[2])
    _write_markdown(records, summary, paths[3])
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Grade OOKc/KrakenC propagation numerical consistency"
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if not args.manifest.is_file():
        parser.error(f"manifest does not exist: {args.manifest}")
    generated = generate_artifacts(args.manifest, args.output_dir)
    payload = json.loads(generated[1].read_text(encoding="utf-8"))
    print(
        json.dumps(
            {
                "case_count": payload["case_count"],
                "slice_count": payload["slice_count"],
                "grade_counts": payload["grade_counts"],
                "overall_grade": payload["overall_grade"],
                "output_dir": str(args.output_dir.resolve()),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
