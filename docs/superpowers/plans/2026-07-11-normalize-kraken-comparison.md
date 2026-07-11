# Normalize 与 Kraken 计算结果对比测试 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 建立一个使用相同 ENV 运行 OOK 与 Fortran Kraken、比较归一化模态函数和复波数并保留环境与原始产物的可重复测试。

**Architecture:** 一个纯 Python 比较模块负责模态符号对齐、误差度量和近似归一化泛函检查；同一模块的运行层负责隔离执行两个真实求解器并保存日志与 JSON。`unittest` 先用合成 MOD 数据验证比较器会正确通过和失败，再运行专用声学环境的端到端 Kraken oracle 测试。

**Tech Stack:** Python 3 标准库、`unittest`、现有 `tools/read_mod_summary.py` MOD 解析器、OOK C++17 可执行文件、Fortran Kraken 可执行文件。

## Global Constraints

- 不修改 `src/algorithm/run.cpp` 或 Fortran Kraken 的生产算法。
- 数值基准固定为 `krakenFortran/build_mingw/out/kraken.exe`，不回退到 `krakenc.exe`。
- 两端使用同一个 `for_test/fixtures/normalize/normalize_reference.env`。
- 模态函数允许整模态乘以 `-1`，但不允许重新缩放。
- `Phi` 最大绝对误差或最大相对误差必须不超过 `5e-5`；近零阈值固定为 `1e-10`。
- 深度网格最大绝对误差不超过 `1e-6 m`；复波数最大绝对误差不超过 `1e-7 1/m`。
- 近似归一化泛函同时要求 OOK/Kraken 最大互差不超过 `5e-5`，且各自与 1 的最大偏差不超过 `5e-2`。后一个较宽容差只补偿 MOD 中单精度波数对临界模态导纳导数的放大，不替代严格的 OOK/Kraken 互差检查。
- 最新一次产物固定保留在 `for_test/artifacts/normalize_kraken/`，清理范围不得越出该目录。
- 不恢复或提交工作区中既有的无关删除项。

## File Structure

- Create: `for_test/normalize_kraken_compare.py` — 纯比较函数、归一化泛函检查、真实进程运行、产物写入和 CLI。
- Create: `for_test/test_normalize_kraken_compare.py` — 合成数据单元测试和真实端到端测试。
- Create: `for_test/fixtures/normalize/normalize_reference.env` — 单频、单声学层、1000 网格点和 1001 MOD 输出深度的固定环境。
- Generate: `for_test/artifacts/normalize_kraken/comparison.json` — 完整指标、阈值、命令、路径、哈希和结论。
- Generate: `for_test/artifacts/normalize_kraken/ook/normalize_reference.mod` — OOK 原始 MOD。
- Generate: `for_test/artifacts/normalize_kraken/kraken/normalize_reference.mod` — Kraken 原始 MOD。
- Generate: `for_test/artifacts/normalize_kraken/*.log` — 两端进程输出。
- Reference: `../tools/read_mod_summary.py` — 已验证的 MOD 二进制解析器。

---

### Task 1: 建立纯模态比较器

**Files:**
- Create: `for_test/test_normalize_kraken_compare.py`
- Create: `for_test/normalize_kraken_compare.py`

**Interfaces:**
- Consumes: `compare_modal_data(ook: dict, kraken: dict, thresholds: Thresholds) -> dict` 的两个 `parse_mod` 风格字典。
- Produces: `Thresholds` 数据类、`compare_modal_data`、`normalization_values` 和包含 `passes`、最差误差位置及逐模态摘要的字典。

- [ ] **Step 1: 写入符号对齐与差异检测的失败测试**

在 `for_test/test_normalize_kraken_compare.py` 中创建 `NormalizeComparatorTests`，使用以下工厂构造最小 MOD 数据：

```python
def modal_data(phi, k=(1.0, 0.0)):
    return {
        "nfreq": 1, "nmedia": 1, "ntot": len(phi), "nmesh": [len(phi) - 1],
        "freqs": [250.0], "z": [float(i) for i in range(len(phi))], "M": 1,
        "phi": [[{"real": float(v), "imag": 0.0} for v in phi]],
        "k": [{"real": float(k[0]), "imag": float(k[1])}],
        "halfspace": {"top": {"bc": "V"}, "bottom": {
            "bc": "A", "cp": {"real": 1590.0, "imag": 0.0}, "rho": 1.2,
        }},
    }

def test_sign_flipped_mode_passes(self):
    result = compare_modal_data(modal_data([0.0, 1.0, 0.0]), modal_data([0.0, -1.0, 0.0]))
    self.assertTrue(result["passes"])
    self.assertEqual(result["modes"][0]["sign"], -1)

def test_shape_difference_fails(self):
    result = compare_modal_data(modal_data([0.0, 1.0, 0.0]), modal_data([0.0, 0.5, 0.0]))
    self.assertFalse(result["passes"])

def test_nonfinite_value_fails(self):
    result = compare_modal_data(modal_data([0.0, float("nan"), 0.0]), modal_data([0.0, 1.0, 0.0]))
    self.assertFalse(result["passes"])
    self.assertIn("non-finite", " ".join(result["failures"]))
```

- [ ] **Step 2: 运行测试确认因模块缺失而失败**

Run:

```powershell
python -m unittest for_test/test_normalize_kraken_compare.py -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'normalize_kraken_compare'`。

- [ ] **Step 3: 实现最小比较器**

在 `for_test/normalize_kraken_compare.py` 中实现：

```python
@dataclass(frozen=True)
class Thresholds:
    phi_abs: float = 5.0e-5
    phi_rel: float = 5.0e-5
    near_zero: float = 1.0e-10
    z_abs: float = 1.0e-6
    k_abs: float = 1.0e-7
    normalization_pair: float = 5.0e-5
    normalization_one: float = 5.0e-2

def as_complex(value):
    return complex(value["real"], value["imag"])

def compare_modal_data(ook, kraken, thresholds=Thresholds()):
    failures = []
    for key in ("nfreq", "nmedia", "ntot", "M"):
        if ook[key] != kraken[key]:
            failures.append(f"{key} mismatch: OOK={ook[key]} Kraken={kraken[key]}")
    if failures:
        return {"passes": False, "failures": failures, "modes": []}
    if not all(math.isfinite(v) for data in (ook, kraken)
               for mode in data["phi"] for pair in mode for v in pair.values()):
        return {"passes": False, "failures": ["non-finite modal value"], "modes": []}
    z_max = max((abs(a - b) for a, b in zip(ook["z"], kraken["z"])), default=0.0)
    k_max = max((abs(as_complex(a) - as_complex(b)) for a, b in zip(ook["k"], kraken["k"])), default=0.0)
    modes = []
    phi_abs = phi_rel = 0.0
    for index, (left, right) in enumerate(zip(ook["phi"], kraken["phi"]), 1):
        a = [as_complex(v) for v in left]
        b = [as_complex(v) for v in right]
        sign = 1 if sum(abs(x-y)**2 for x, y in zip(a, b)) <= sum(abs(x+y)**2 for x, y in zip(a, b)) else -1
        abs_err = [abs(x-sign*y) for x, y in zip(a, b)]
        rel_err = [e/max(abs(x), abs(sign*y)) for e, x, y in zip(abs_err, a, b)
                   if max(abs(x), abs(sign*y)) >= thresholds.near_zero]
        row = {"mode": index, "sign": sign, "max_abs_phi": max(abs_err, default=0.0),
               "max_rel_phi": max(rel_err, default=0.0)}
        modes.append(row)
        phi_abs = max(phi_abs, row["max_abs_phi"])
        phi_rel = max(phi_rel, row["max_rel_phi"])
    if z_max > thresholds.z_abs:
        failures.append(f"z max abs {z_max} exceeds {thresholds.z_abs}")
    if k_max > thresholds.k_abs:
        failures.append(f"k max abs {k_max} exceeds {thresholds.k_abs}")
    if phi_abs > thresholds.phi_abs and phi_rel > thresholds.phi_rel:
        failures.append(f"Phi errors abs={phi_abs} rel={phi_rel}")
    return {"passes": not failures, "failures": failures, "max_abs_z": z_max,
            "max_abs_k": k_max, "max_abs_phi": phi_abs, "max_rel_phi": phi_rel,
            "modes": modes}
```

- [ ] **Step 4: 运行比较器单元测试确认通过**

Run: `python -m unittest for_test/test_normalize_kraken_compare.py -v`

Expected: 3 tests PASS。

- [ ] **Step 5: 提交纯比较器**

```powershell
git add for_test/normalize_kraken_compare.py for_test/test_normalize_kraken_compare.py
git commit -m "test: add Normalize modal comparator"
```

### Task 2: 添加固定环境和真实求解器运行层

**Files:**
- Modify: `for_test/test_normalize_kraken_compare.py`
- Modify: `for_test/normalize_kraken_compare.py`
- Create: `for_test/fixtures/normalize/normalize_reference.env`

**Interfaces:**
- Consumes: `run_comparison(workspace_root: Path, artifact_root: Path, timeout: int = 120) -> dict`。
- Produces: 隔离运行目录、两份 MOD、两份日志和完整 `comparison.json`。

- [ ] **Step 1: 写入端到端失败测试**

追加：

```python
class NormalizeKrakenIntegrationTests(unittest.TestCase):
    def test_reference_environment_matches_fortran_kraken(self):
        workspace = Path(__file__).resolve().parents[2]
        artifacts = Path(__file__).resolve().parent / "artifacts" / "normalize_kraken"
        result = run_comparison(workspace, artifacts, timeout=120)
        self.assertTrue(result["passes"], result["failures"])
        self.assertTrue((artifacts / "comparison.json").is_file())
        self.assertTrue((artifacts / "ook" / "normalize_reference.mod").is_file())
        self.assertTrue((artifacts / "kraken" / "normalize_reference.mod").is_file())
```

- [ ] **Step 2: 运行端到端测试确认因运行层或环境缺失而失败**

Run:

```powershell
python -m unittest for_test.test_normalize_kraken_compare.NormalizeKrakenIntegrationTests -v
```

Expected: FAIL with missing `run_comparison` or missing `normalize_reference.env`。

- [ ] **Step 3: 写入固定 ENV**

将 `calibK.env` 的物理参数固定为专用用例，并把输出深度改为完整 1001 点网格：

```text
'Normalize reference: vacuum top, acoustic half-space bottom'
250.0
1
'CVW'
1000  0.0  100.0
     0.0  1500.0 /
   100.0  1500.0 /
'A'  0.0
   100.0  1590.0  0.0  1.2 0.5 /
1400.0  20000
0.0
1
50.0 /
1001
0.0  100.0 /
```

- [ ] **Step 4: 实现安全产物目录、进程执行和 MOD 解析加载**

实现 `prepare_artifacts`，先解析目标路径并验证其位于 `OOK/for_test/artifacts` 下，之后才允许清理；实现 `run_command` 使用 `subprocess.run(..., timeout=timeout, capture_output=True, text=True, encoding="utf-8", errors="replace")`，并把命令、工作目录、返回码和输出写入日志。通过 `importlib.util.spec_from_file_location` 从 `workspace_root/tools/read_mod_summary.py` 加载 `parse_mod`。

OOK 命令固定为：

```python
[str(workspace_root / "OOK/bin/ook_case_runner.exe"), str(env_path),
 str(artifact_root / "ook/normalize_reference"), "1", "--mod"]
```

Kraken 命令固定为：

```python
[str(workspace_root / "krakenFortran/build_mingw/out/kraken.exe"), "normalize_reference"]
```

Kraken 的 `cwd` 固定为 `artifact_root / "kraken"`，且先把同一个 ENV 复制进去。

- [ ] **Step 5: 实现近似归一化泛函和 JSON 报告**

对每个模态按完整 1001 点 `z` 网格计算单层介质梯形积分；对声学下半空间，用 MOD 中的 `cp`、`rho` 和 `k.real**2`，按 `x1=0.9999999*x`、`x2=1.0000001*x` 对 `real(sqrt(x-omega**2/cp**2))/rho` 做相同有限差分，得到底导纳导数并加入末端项。把两端逐模态值、最大互差和各自最大 `|RN-1|` 写入比较结果，并纳入 Global Constraints 的通过条件。

`comparison.json` 还必须写入：ISO 时间、环境和两端可执行文件绝对路径及 SHA-256、实际命令、阈值、元数据、全局最差误差、前 20/中 20/后 20 模态摘要和 `passes`。

- [ ] **Step 6: 运行端到端测试并按实际红灯校准诊断而非放宽阈值**

Run:

```powershell
python -m unittest for_test.test_normalize_kraken_compare.NormalizeKrakenIntegrationTests -v
```

Expected: PASS。若失败，先检查命令、MOD 格式、模态匹配和归一化公式；不得把 `Phi=5e-5`、`z=1e-6`、`k=1e-7` 或归一化互差 `5e-5` 放宽。

- [ ] **Step 7: 提交环境和运行层**

```powershell
git add for_test/normalize_kraken_compare.py for_test/test_normalize_kraken_compare.py for_test/fixtures/normalize/normalize_reference.env
git commit -m "test: compare Normalize output with Kraken"
```

### Task 3: 完整验证和产物复核

**Files:**
- Verify: `for_test/artifacts/normalize_kraken/comparison.json`
- Verify: `for_test/artifacts/normalize_kraken/ook/normalize_reference.mod`
- Verify: `for_test/artifacts/normalize_kraken/kraken/normalize_reference.mod`

**Interfaces:**
- Consumes: Task 1 与 Task 2 的测试和运行器。
- Produces: 新鲜的测试证据和保留在固定目录中的最终产物。

- [ ] **Step 1: 运行全部 Normalize 测试**

Run: `python -m unittest for_test/test_normalize_kraken_compare.py -v`

Expected: 所有比较器与端到端测试 PASS，无 skipped。

- [ ] **Step 2: 运行现有 MOD 工具测试**

Run:

```powershell
$env:PYTHONPATH='E:\my_project\ook_project\tools'
python -m unittest E:\my_project\ook_project\tools\test_mod_tools.py -v
```

Expected: 4 tests PASS。

- [ ] **Step 3: 运行 OOK CTest 回归**

Run: `ctest --test-dir build_acceptance --output-on-failure`

Expected: `100% tests passed`。

- [ ] **Step 4: 复核 JSON 与原始文件**

Run:

```powershell
$r = Get-Content -Raw for_test/artifacts/normalize_kraken/comparison.json | ConvertFrom-Json
$r | Select-Object passes,failures,max_abs_z,max_abs_k,max_abs_phi,max_rel_phi
Get-Item for_test/artifacts/normalize_kraken/ook/normalize_reference.mod,
         for_test/artifacts/normalize_kraken/kraken/normalize_reference.mod,
         for_test/fixtures/normalize/normalize_reference.env |
  Select-Object FullName,Length,LastWriteTime
```

Expected: `passes=True`、`failures=[]`，三项文件均存在且非空。

- [ ] **Step 5: 检查最终差异和无关修改**

Run:

```powershell
git status --short
git diff --check HEAD~2..HEAD
```

Expected: 只包含本计划文件、比较器、测试和环境的预期提交；用户原有两个删除项仍保持为未提交删除。
