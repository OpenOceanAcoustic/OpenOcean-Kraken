# OOKc 与 KrakenC 全环境传播损失对比设计

## 1. 目标

对 `E:/my_project/02_ookc_project/test` 中全部 11 个 `.env` 环境分别运行：

- OOKc：`OpenOcean-Krakenc.exe --mod` 与 `--field`；
- KrakenC：Fortran `krakenc.exe` 与 `field.exe`。

从两条计算链产生的 SHD 文件读取复压力，绘制距离—深度传播损失场，并同时完成数值一致性分析和声传播物理特征分析。`MunkKleaky` 的 25 m、250 m 两个声源深度分别作图，其余环境各有一个声源深度，因此预期生成 12 张组合主图。

本任务只新增分析与测试代码，不修改 OOKc 或 Fortran 求解器算法，也不修改 `test` 中的原始输入文件。

## 2. 计算对象与输入处理

输入集合由 `test/*.env` 动态发现。每个 ENV 必须存在同名 FLP；同名的 `.sbp`、`.brc`、`.irc` 文件若存在，也复制到各自隔离的运行目录。

两个特殊映射只在生成的运行目录内执行，并记录在清单中：

- `neggradK_brc` 缺少同名 BRC，使用 `neggradC_brc.brc` 作为 `neggradK_brc.brc`；
- `neggradK_irc` 缺少同名 IRC，使用 `neggradC_irc.irc` 作为 `neggradK_irc.irc`。

这两个文件描述反射系数数据，与接收深度网格无关。映射不会回写 `test`。

每个环境在以下独立目录中运行，避免 MOD、SHD 和日志互相覆盖：

```text
propagation_loss_comparison/
  runs/<case>/ookc/
  runs/<case>/krakenc/
```

## 3. 代码结构

新增以下文件：

- `for_test/propagation_loss_analysis.py`：SHD 压力重排、TL 计算、有效网格掩码、误差指标、组合图绘制等纯分析逻辑；
- `for_test/plot_propagation_loss_comparison.py`：输入发现、两条求解链编排、隔离目录、日志与清单、CSV 和 Markdown 报告生成；
- `for_test/test_propagation_loss_analysis.py`：分析逻辑及批处理辅助行为的自动化测试。

复用现有组件：

- `for_test/shd_reader.py` 读取 Fortran direct-access SHD；
- `for_test/cpp_pipeline.py` 和 `for_test/reference_pipeline.py` 的运行、复制和清单模式，但批处理入口自行统一两条链的错误记录，以便单案例失败后继续运行其余案例。

命令行入口支持显式传入测试目录、三个可执行文件、输出目录和 OOKc 线程数。默认调用形式为：

```powershell
python for_test/plot_propagation_loss_comparison.py `
  --test-dir ../test `
  --ookc build-release/OpenOcean-Krakenc.exe `
  --krakenc ../krakenFortran/build-acceptance-mingw/out/krakenc.exe `
  --field ../krakenFortran/build-acceptance-mingw/out/field.exe `
  --output ../propagation_loss_comparison
```

## 4. 数据流

对每个环境依次执行：

1. 枚举并复制 ENV、FLP 和辅助文件到 OOKc、KrakenC 隔离目录；
2. 在 OOKc 目录执行 `--mod`，再执行 `--field`；
3. 在 KrakenC 目录依次执行 `krakenc.exe <case>` 和 `field.exe <case>`；
4. 读取两侧 SHD，并验证频率、声源深度、接收深度、距离和压力数组尺寸完全一致；
5. 将 `shd_reader` 的扁平压力按 `(source_depth, receiver_depth, range)` 重排；
6. 对每个声源深度计算 TL、差值和指标，生成一张组合图与一行 CSV；
7. 汇总案例状态、文件哈希、命令、耗时和输出文件，生成 JSON 清单与中文分析报告。

网格不一致时不插值，直接把该环境标记为不可比较，从而避免插值掩盖接口或文件格式问题。

## 5. 传播损失和误差口径

传播损失定义为：

```text
TL = -20 log10(max(|p|, 1e-12)) dB
```

其中 `p` 为 SHD 中的复压力。误差指标包括：

- 复压力相对 L2：`||p_ookc - p_krakenc||2 / ||p_krakenc||2`；
- TL 平均绝对误差 MAE；
- TL 均方根误差 RMSE；
- TL 绝对误差 P95；
- TL 最大绝对误差；
- TL 指标有效网格覆盖率。

为避免深衰落零点把 dB 指标无限放大，TL 汇总指标只使用满足下式的网格点：

```text
max(|p_ookc|, |p_krakenc|) >= max(global_peak_ookc, global_peak_krakenc) * 1e-8
```

复压力相对 L2 仍使用全部网格点，作为首要一致性指标。

## 6. 图形设计

每个环境和声源深度生成一张高分辨率 PNG，采用 2×2 布局：

1. OOKc 二维 TL 场；
2. KrakenC 二维 TL 场；
3. `TL_OOKc - TL_KrakenC` 二维差值场；
4. 最接近声源深度的接收深度处，两条距离—TL 曲线及关键误差指标。

绘图约束：

- 距离统一显示为 km，接收深度显示为 m，深度轴向下；
- OOKc 和 KrakenC 使用由两者联合有限值的第 2、98 百分位确定的相同 TL 色标；
- 差值图使用以 0 dB 为中心的对称发散色标，范围取绝对差值第 98 百分位，最小半宽为 0.1 dB；
- 标题包含环境名、频率和声源深度；
- 图片文件名包含环境名和声源深度，避免多声源覆盖；
- 图中文字使用英文以避免运行环境缺少中文字体，中文解释放在 Markdown 报告中。

若某环境无法完成两侧计算，仍生成同名失败占位图，注明失败阶段和日志位置，确保 12 个声源任务均有可见状态；如果在读取 SHD 前无法确定实际声源数，则至少按 FLP/已成功一侧 SHD 中可识别的声源数生成状态记录。

## 7. 对比分析报告

`analysis.md` 同时覆盖两类分析：

### 7.1 数值一致性

- 按复压力相对 L2 从小到大汇总所有声源场；
- 报告 TL MAE、RMSE、P95、最大差异和有效网格覆盖率；
- 定位差值最大的距离与深度；
- 结合项目既有端到端阈值解释一致性等级，但不把经验等级冒充新的验收标准；
- 区分广泛相位偏差与仅发生在深衰落零点附近的 dB 放大。

### 7.2 声传播物理特征

- `MunkKleaky`：深海声道、两个声源深度及远距离会聚/干涉特征；
- `elastic_fd_two_layer`、`multilayer_elastic_stack`、`multilayer_mud_sand`：弹性层、软泥和砂层导致的泄漏、底损失和垂向结构；
- `neggrad*_brc`、`neggrad*_irc`：负声速梯度和外部/内部反射系数边界影响；
- `solve3_mode_gain`：多剖面连接时的模式增减和相位连续性；
- `stepK_rd`、`wedge`：距离相关阶跃、楔形波导和耦合模传播特征。

每个环境至少包含一段基于实际图形和指标的结论，不只复述输入文件名称。

## 8. 错误处理与可追溯性

单个环境失败不终止批处理。每个阶段记录：

- 完整命令和工作目录；
- 返回码与耗时；
- 标准输出、标准错误日志；
- ENV、FLP、辅助文件、可执行文件、MOD、SHD 和图片的绝对路径、大小与 SHA-256；
- 辅助文件映射来源；
- 状态：`success`、`ookc_failed`、`krakenc_failed`、`grid_mismatch` 或 `analysis_failed`。

路径安全检查禁止把文件系统根目录或项目输入目录作为可清理的输出目录。批处理只重建显式输出目录下的 `runs`、`figures` 和汇总文件。

## 9. 测试策略与完成标准

实现遵循测试先行。自动化测试至少覆盖：

- 多声源 SHD 压力重排顺序；
- TL 下限与有限值；
- 有效网格掩码；
- 复压力 L2 与 TL 各项指标；
- 两幅 TL 图的联合色标和零差异色标下限；
- 最近声源深度的接收深度选择；
- BRC/IRC 辅助文件映射及来源记录；
- 网格不一致拒绝比较；
- 单案例失败后批处理继续；
- 图片、CSV、清单和报告的最小结构。

完成前执行：

1. 分析单元测试；
2. 一个小环境的端到端冒烟运行；
3. 全部 11 个环境的正式批处理；
4. 重新运行测试；
5. 核对清单中 11 个环境、12 个声源任务和对应图形/状态；
6. 生成图片缩略图总览，并逐图检查坐标、色标、标题、空白区域和异常伪影；
7. 逐项复核 `analysis.md` 是否覆盖数值一致性和物理特征。

最终交付目录包含：

```text
propagation_loss_comparison/
  figures/*.png
  runs/<case>/{ookc,krakenc}/...
  metrics.csv
  manifest.json
  analysis.md
  overview.png
```
