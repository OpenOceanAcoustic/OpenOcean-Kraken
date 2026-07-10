# OOK、Kraken 与 Krakenc 方案选型对比报告 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 基于当前工作区源码、构建配置和验收产物，生成一份中立、可追溯、离线可打开的 OOK、Kraken 与 Krakenc 中文方案选型 HTML 报告。

**Architecture:** 先建立可复核的证据清单并完成三程序运行探测，再以单文件 HTML 固化事实、限制和场景适配结论。独立 PowerShell 校验脚本检查报告结构、编码、离线性、证据标签和交互契约。

**Tech Stack:** C++17/CMake、Fortran/CMake、PowerShell 7/Windows PowerShell、HTML5、CSS3、原生 JavaScript。

## Global Constraints

- 主报告必须是单文件 HTML，CSS 与 JavaScript 内嵌，不依赖网络资源。
- 报告必须保持中立，不给出唯一推荐方案，也不计算综合总分。
- 报告同时覆盖数值可信度、功能覆盖、运行性能、集成与维护成本。
- 证据分为 A 级实测、B 级源码确认、C 级历史报告和未验证。
- 不支持的案例标记为“不适用”，缺少同条件证据的项目标记为“未验证”。
- 主报告保存到 `E:/my_project/ook_project/docs/07_ook_kraken_krakenc_selection_report.html`。

## File Structure

- Create: `E:/my_project/ook_project/tools/validate_ook_kraken_krakenc_report.ps1` — 验证报告编码、结构、离线性、证据标记和筛选脚本。
- Create: `E:/my_project/ook_project/docs/07_ook_kraken_krakenc_selection_report.html` — 最终单文件方案选型报告。
- Reference: `E:/my_project/ook_project/docs/06_ook_acceptance_with_velocity/04_artifacts/kraken_numeric_current_test/summary.json` — 当前 11 用例状态与历史耗时。
- Reference: `E:/my_project/ook_project/docs/06_ook_acceptance_with_velocity/02_numeric_regression_and_velocity_acceptance_report.html` — TL 与振速验收指标。
- Reference: `E:/my_project/ook_project/OOK/CMakeLists.txt` and `OOK/include/*.h` — OOK 构建与接口证据。
- Reference: `E:/my_project/ook_project/krakenFortran/kraken/CMakeLists.txt`, `kraken.f90`, `krakenc.f90` — 两个 Fortran 求解器的实现差异。

---

### Task 1: 建立当前证据基线

**Files:**
- Inspect: `E:/my_project/ook_project/OOK/include/OpenOceanKrakenInterface.h`
- Inspect: `E:/my_project/ook_project/OOK/include/OpenOceanKrakenParams.h`
- Inspect: `E:/my_project/ook_project/OOK/src/algorithm/run.cpp`
- Inspect: `E:/my_project/ook_project/krakenFortran/kraken/kraken.f90`
- Inspect: `E:/my_project/ook_project/krakenFortran/kraken/krakenc.f90`
- Inspect: `E:/my_project/ook_project/docs/06_ook_acceptance_with_velocity/04_artifacts/kraken_numeric_current_test/summary.json`

**Interfaces:**
- Consumes: 当前源码、已编译可执行文件和既有验收产物。
- Produces: 报告可直接使用的三类事实：源码能力、当前可运行状态、现有数值与耗时指标。

- [ ] **Step 1: 验证三个核心构建目标存在**

Run:

```powershell
$paths = @(
  'E:/my_project/ook_project/OOK/bin/ook_case_runner.exe',
  'E:/my_project/ook_project/krakenFortran/build_mingw/out/kraken.exe',
  'E:/my_project/ook_project/krakenFortran/build_mingw/out/krakenc.exe'
)
$paths | ForEach-Object { [pscustomobject]@{ Path = $_; Exists = Test-Path -LiteralPath $_ } }
```

Expected: 三行 `Exists=True`；若有缺失，在报告中标为“当前构建未验证”，不推断源代码能力缺失。

- [ ] **Step 2: 运行 OOK 接口测试**

Run:

```powershell
ctest --test-dir E:/my_project/ook_project/OOK/build_acceptance --output-on-failure
```

Expected: `100% tests passed`；若失败，记录失败测试名和返回码，并将当前构建状态标为 A 级失败事实。

- [ ] **Step 3: 复核运行时间数据口径**

Run:

```powershell
$p = 'E:/my_project/ook_project/docs/06_ook_acceptance_with_velocity/04_artifacts/kraken_numeric_current_test/summary.json'
$j = Get-Content -LiteralPath $p -Raw | ConvertFrom-Json
$j.cases | Select-Object case,ook_status,kraken_status,seconds_ook,seconds_kraken | Format-Table -AutoSize
```

Expected: 11 个用例均为 `ok`。将 `seconds_kraken` 标注为脚本“Kraken 优先、Krakenc 回退”的组合口径，不把它当作 Krakenc 独立性能数据。

- [ ] **Step 4: 对 Kraken 与 Krakenc 做隔离冒烟探测**

Run both executables in separate temporary directories copied from the `calibK` case; require exit code `0` and a generated `.mod`. Measure wall-clock time with `Measure-Command`. If either result is missing, retain its status as “未验证” instead of substituting the other solver’s result.

Expected: 获得两个独立的成功/失败状态与墙钟耗时；这些单用例数据只作为冒烟性能证据，不外推为完整性能排名。

- [ ] **Step 5: 汇总确定性事实**

Record these already established facts for the report:

```text
OOK: C++17 static core library + executable + interface test target; supports programmatic parameter/interface use and configurable thread pool.
Kraken: Fortran real-coefficient normal-mode solver; dynamic allocation for main modal arrays; command-line file-root workflow.
Krakenc: Fortran complex-coefficient variant; complex finite-difference matrices and PekRoot dependency; command-line file-root workflow.
Current OOK–Kraken evidence: 11/11 pressure/TL cases passed; maximum P95 |ΔTL| = 0.087538 dB.
Velocity evidence: current OOK pressure/horizontal-velocity equivalent TL median error is roughly 0.20–0.33 dB in the three focused cases; vertical velocity is present but lacks an equally independent oracle.
```

### Task 2: 编写报告契约校验器

**Files:**
- Create: `E:/my_project/ook_project/tools/validate_ook_kraken_krakenc_report.ps1`
- Test: `E:/my_project/ook_project/docs/07_ook_kraken_krakenc_selection_report.html`

**Interfaces:**
- Consumes: `-ReportPath <string>`，默认指向最终 HTML。
- Produces: 成功时输出 `REPORT_VALIDATION_OK` 并返回 0；任何契约失败时抛出异常并返回非零。

- [ ] **Step 1: 写入校验脚本**

```powershell
param(
  [string]$ReportPath = (Join-Path $PSScriptRoot '../docs/07_ook_kraken_krakenc_selection_report.html')
)

$resolved = [System.IO.Path]::GetFullPath($ReportPath)
if (-not (Test-Path -LiteralPath $resolved)) { throw "Report not found: $resolved" }
$html = [System.IO.File]::ReadAllText($resolved, [System.Text.Encoding]::UTF8)
$required = @(
  '<meta charset="utf-8">', 'OOK', 'Kraken', 'Krakenc',
  '数值可信度', '功能覆盖', '性能与工程化', '场景适配',
  'data-evidence=', 'data-filter=', '未验证', '@media print'
)
foreach ($token in $required) {
  if (-not $html.Contains($token)) { throw "Missing required token: $token" }
}
if ($html -match '<script[^>]+src=' -or $html -match '<link[^>]+href=["'']https?://') {
  throw 'External script or stylesheet dependency detected.'
}
if ($html -match 'æ•°|å¯¹æ¯”|ç»“è®º') { throw 'Likely UTF-8 mojibake detected.' }
$ids = [regex]::Matches($html, 'id="([^"]+)"') | ForEach-Object { $_.Groups[1].Value }
$hrefs = [regex]::Matches($html, 'href="#([^"]+)"') | ForEach-Object { $_.Groups[1].Value }
foreach ($href in $hrefs) {
  if ($href -notin $ids) { throw "Broken internal anchor: #$href" }
}
'REPORT_VALIDATION_OK'
```

- [ ] **Step 2: 在报告创建前运行校验器并确认失败**

Run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File E:/my_project/ook_project/tools/validate_ook_kraken_krakenc_report.ps1
```

Expected: FAIL with `Report not found` or a missing-token error.

### Task 3: 创建单文件 HTML 报告

**Files:**
- Create: `E:/my_project/ook_project/docs/07_ook_kraken_krakenc_selection_report.html`
- Test: `E:/my_project/ook_project/tools/validate_ook_kraken_krakenc_report.ps1`

**Interfaces:**
- Consumes: Task 1 的证据事实及 Task 2 的 HTML 契约。
- Produces: 不依赖网络、可筛选、可打印且证据可追溯的单文件中文报告。

- [ ] **Step 1: 写入语义化 HTML 骨架和内嵌样式**

Use this exact top-level structure and keep all CSS inside `<style>`:

```html
<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>OOK、Kraken 与 Krakenc 方案选型对比报告</title>
  <style>
    :root { color-scheme: light; --ink:#14213d; --paper:#f6f2ea; --card:#fff; --line:#d8d3c8; --accent:#0f766e; }
    * { box-sizing:border-box; }
    body { margin:0; color:var(--ink); background:var(--paper); font:16px/1.65 "Microsoft YaHei","PingFang SC",sans-serif; }
    header,nav,main { width:min(1180px,calc(100% - 32px)); margin-inline:auto; }
    section { margin:24px 0; padding:24px; background:var(--card); border:1px solid var(--line); border-radius:16px; }
    table { width:100%; border-collapse:collapse; }
    th,td { padding:10px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }
    button[aria-pressed="true"] { color:#fff; background:var(--accent); }
    @media (max-width:760px) { section { padding:16px; } .table-wrap { overflow-x:auto; } }
    @media print { .filters { display:none; } [data-evidence] { display:table-row !important; } }
  </style>
</head>
<body>
  <header id="top"><h1>OOK、Kraken 与 Krakenc 方案选型对比报告</h1></header>
  <nav aria-label="报告目录"><a href="#summary">执行摘要</a> · <a href="#matrix">综合矩阵</a> · <a href="#evidence">证据索引</a></nav>
  <main>
    <section id="summary"><h2>执行摘要</h2></section>
    <section id="profiles"><h2>项目画像</h2></section>
    <section id="features"><h2>算法与功能</h2></section>
    <section id="numerics"><h2>数值可信度</h2></section>
    <section id="engineering"><h2>性能与工程化</h2></section>
    <section id="matrix"><h2>综合对比矩阵</h2></section>
    <section id="scenarios"><h2>场景适配</h2></section>
    <section id="evidence"><h2>证据索引</h2></section>
  </main>
  <script>
    const buttons=[...document.querySelectorAll('[data-filter]')];
    buttons.forEach(button=>button.addEventListener('click',()=>{
      const filter=button.dataset.filter;
      buttons.forEach(item=>item.setAttribute('aria-pressed',String(item===button)));
      document.querySelectorAll('[data-evidence]').forEach(row=>{
        row.hidden=filter!=='all' && row.dataset.evidence!==filter;
      });
    }));
  </script>
</body>
</html>
```

The finished file must add complete Chinese prose, tables and evidence links below every section heading.

- [ ] **Step 2: 填入三者画像与算法差异**

State the verified architectural distinction precisely: OOK is a C++17 engineering implementation whose current algorithm path is benchmarked against Kraken; Kraken uses predominantly real finite-difference coefficients with complex modal wavenumbers; Krakenc uses complex finite-difference coefficients and the `PekRoot` path for lossy/complex problems. Avoid calling Krakenc “coupled mode” because the inspected source identifies it as a complex-coefficient normal-mode solver.

- [ ] **Step 3: 填入数值、性能和证据限制**

Include the 11/11 OOK–Kraken result, `0.087538 dB` maximum P95 TL difference, velocity evidence, historical per-case elapsed-time table, and the isolated `calibK` smoke probe. Explicitly label historical `seconds_kraken` as Kraken-first/fallback-composite timing and label the isolated timing as a single-case smoke result.

- [ ] **Step 4: 填入中立矩阵和场景适配表**

Every row must carry `data-evidence="verified"` or `data-evidence="unverified"`. Use “强/中/弱/未验证/不适用” only with a textual reason; do not add a total score, winner badge or single recommended product.

- [ ] **Step 5: 添加筛选交互和打印样式**

Implement buttons with `data-filter="all"`, `data-filter="verified"`, and `data-filter="unverified"`. The script must toggle rows by their `data-evidence` value, update `aria-pressed`, and leave all content visible when printing. Add `@media print` and responsive rules below 760 px.

- [ ] **Step 6: 运行自动校验**

Run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File E:/my_project/ook_project/tools/validate_ook_kraken_krakenc_report.ps1
```

Expected: `REPORT_VALIDATION_OK`.

- [ ] **Step 7: 进行浏览器视觉检查**

Open the HTML locally and inspect at desktop width and approximately 390 px width. Confirm there is no horizontal page overflow, comparison tables remain readable, filters work, anchors reach the correct section, and print preview hides controls while retaining all evidence rows.

- [ ] **Step 8: 复核最终文件与仓库状态**

Run:

```powershell
Get-Item E:/my_project/ook_project/docs/07_ook_kraken_krakenc_selection_report.html | Select-Object FullName,Length,LastWriteTime
git -C E:/my_project/ook_project/OOK status --short
```

Expected: HTML exists and is non-empty; OOK status contains no unintended source changes. The workspace root is not a Git repository, so the report and validator are delivered as workspace artifacts rather than committed OOK source changes.
