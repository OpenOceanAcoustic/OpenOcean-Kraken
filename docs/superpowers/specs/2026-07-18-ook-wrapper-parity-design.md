# OOK Wrapper 与接口完善设计

## 目标

参考 OpenOcean-Bellhop（OOB）的接口分层和使用体验，为 OpenOcean-Kraken（OOK）完整补齐 Python 与 MATLAB 工作流，并整理第三方依赖与 CMake 构建。实现必须适配 Kraken 的模态与声场能力，不能照搬 Bellhop 专属的射线、到达和波束接口，也不能改变现有 Kraken 核心计算结果。

## 范围

本次交付包括：

- pybind11 原生模块 `OpenOceanKraken`。
- Python 高层包 `wrappers/py_kraken`。
- MATLAB 包装层 `wrappers/m_kraken`。
- 可供 MATLAB 和终端调用的正式命令行接口。
- SHD、MOD、JSON 的读取、数据模型和绘图支持。
- 第三方依赖迁移及 CMake 构建选项。
- C++、CLI、Python、MATLAB 和回归测试。

本次不包括：

- Bellhop 的 ray、arrival、beam 或 angle 接口。
- Kraken 核心算法重构。
- 对 OOK 现有数值行为的主动修改。

## 总体架构

### C++ 核心层

保留现有 `OpenOceanKraken::Interface` 的计算行为。为绑定层增加结果访问所需的安全、带类型接口，使 Python 能取得压力、垂直速度、水平速度和模态结果，而不直接管理 OOK 的裸指针。

Python 默认取得拥有自身内存的 NumPy 副本。这样，在接口重新运行、调用 `clearResults()` 或 `free()` 后，已返回的数组仍然有效。原始裸地址不作为推荐公开接口。

### 构建与命令行层

第三方依赖统一迁移为：

- `third-party/eigen`
- `third-party/nlohmann-json`
- `third-party/pybind11`

迁移完成后删除根目录旧的 `eigen/` 和 `nlohmann/`，不保留重复副本。所有项目内 include 路径和 CMake target include 路径同步更新。

CMake 提供 `BUILD_AS_EXE`、`BUILD_AS_STATIC`、`BUILD_AS_SHARED`、`BUILD_AS_PYTHON` 和 CMake 标准的 `BUILD_TESTING` 独立构建选项。Python 模块使用当前 CMake 发现的解释器构建，不写死 Python 小版本或扩展文件后缀。

当前 `main.cpp` 中写死的本机输入输出路径改为正式 CLI。CLI 接受 ENV 或 JSON 输入，并支持：

```text
OpenOceanKraken <input.env|input.json>
  --output <result-root>
  --threads <N>
  --velocity
  --mod
  --mod-only
  --json
```

`--velocity` 输出压力、垂直速度和水平速度 SHD；`--mod` 在正常场计算后额外输出 MOD；`--mod-only` 只运行模态计算并输出 MOD；`--json` 导出规范化配置。无效参数、输入失败和输出失败均返回非零退出码，并将可操作的错误信息写到 stderr。

### Python 层

Python 包名为 `py_kraken`，原生模块名为 `OpenOceanKraken`，高层门面类为 `OpenOceanKraken_interface`，高层方法统一使用 `ook_` 前缀。

主要文件职责如下：

- `ook_data_model.py`：配置、声场和模态数据模型，以及原生枚举别名。
- `ook_interface.py`：面向用户的高层门面。
- `ook_read.py`：读取 SHD、MOD 和 JSON。
- `ook_plot.py`：压力/速度传播损失、模态波数和模态函数绘图。
- `demo.py`：ENV、JSON、内存结果和文件结果示例。
- 安装脚本、包元数据、`__init__.py` 和类型提示：完成可复现安装与 IDE 支持。

高层门面支持正常构造多个独立实例；同时提供 `get_instance()` 作为可选便捷入口。独立实例是默认推荐用法，避免多个算例共享状态。

核心高层 API 为：

```python
ook = OpenOceanKraken_interface(thread_num=8)

ook.ook_load_env("case.env")
ook.ook_load_json("case.json")
ook.ook_export_json("case.json")

ook.ook_set_frequency(50.0)
ook.ook_set_source_depth(np.array([1000.0]))
ook.ook_set_receiver_depth(start=0.0, end=5000.0, count=501)
ook.ook_set_receiver_range(start=0.0, end=100000.0, count=1001)
ook.ook_set_phase_speed(c_low, c_high)
ook.ook_set_ssp([range_independent_area])
ook.ook_set_run_mode(Run_Mode.MODE_B_Both)
ook.ook_set_velocity_enable(True)

ook.ook_run()
ook.ook_run_eigen()
ook.ook_run_field()

pressure = ook.ook_get_pressure()
vertical = ook.ook_get_vertical_velocity()
horizontal = ook.ook_get_horizontal_velocity()
modes = ook.ook_get_modes()

ook.ook_export_shd("result", data_type="pressure")
ook.ook_export_mod("result")
ook.ook_export_result("result")
```

pybind 原生层公开 OOK 实际使用的线程池、枚举、半空间、SSP 分层、剖面和反射系数值类型，以及 `Interface` 的输入、运行、清理、导出和安全结果访问方法。Eigen 向量和矩阵通过 pybind11 Eigen/NumPy 支持转换。

### Python 返回模型

`FieldData.values` 使用 `numpy.complex64`，维度顺序固定为 `[source, range, depth]`。模型同时携带标题、频率、声源深度、接收距离和接收深度坐标。单声源便捷访问可以返回二维视图，但底层标准模型始终保留声源维。

`ModeData` 携带剖面距离、深度网格、复波数、群速度和模态函数。由于不同剖面的模态数允许不同，数据按剖面存放为数组列表，不强制填充成一个矩形张量。

### MATLAB 层

MATLAB 不使用 MEX，采用 OOB 已验证的文件驱动方式：数据模型写入 JSON，调用正式 `OpenOceanKraken` CLI，再读取 MOD/SHD。

`krakenDataModel.m` 覆盖：

- 标题、线程数、频率和频率向量。
- 声源深度、接收深度和接收距离。
- 分层 SSP、剖面、上下半空间和衰减单位。
- 相速度范围、网格类型和源类型。
- 模态上限、剖面数和剖面距离。
- coherent/incoherent 与 adiabatic/coupled 运行选项。
- 反射系数、SBP 和速度场开关。

MATLAB 包含 SHD/MOD reader，以及压力场、速度场、模态波数和模态函数绘图函数。示例脚本覆盖 ENV 输入和数据模型生成 JSON 两种工作流。

## 数据流

Python 直接调用路径为：

```text
py_kraken high-level API
  -> OpenOceanKraken pybind module
  -> OpenOceanKraken::Interface
  -> owned NumPy result copies
  -> data model / plotting
```

MATLAB 调用路径为：

```text
krakenDataModel
  -> JSON
  -> OpenOceanKraken CLI
  -> MOD / SHD
  -> MATLAB readers
  -> arrays / plotting
```

## 错误处理

Python 使用明确的异常类型：

- 文件不存在抛出 `FileNotFoundError`。
- 数组维度、枚举、索引或数值范围错误抛出 `ValueError` 或 `IndexError`。
- 未运行就读取结果、速度未启用却读取速度、结果已释放等生命周期错误抛出 `RuntimeError`。
- C++ 标准异常由 pybind 转换为 Python 异常，不用空指针或静默失败表示错误。

MATLAB 检查 CLI 退出码、输出文件是否存在以及记录是否完整。失败异常必须包含执行命令和 CLI 错误输出。SHD/MOD reader 在记录长度、维度或文件截断不合法时立即报错。

## 兼容性

- 支持 CPython 3.9 及以上；本次在工作区现有 Python 3.12 环境完成实际验证。
- Windows 是本次完整验收平台；CMake 和源代码保持 Linux/macOS 可构建，不提交平台固定的 `.pyd` 文件名逻辑。
- MATLAB 通过 `system` 调用 CLI，不依赖 MEX ABI。
- 现有 C++ 公开接口保持兼容，新增接口不改变已有调用语义。

## 测试与验收

所有新增行为按测试驱动方式实现。

### 第三方依赖迁移

从全新构建目录运行配置与编译。代码和 CMake 中不得再引用根目录 `eigen/`、`nlohmann/`，根目录不得保留这两个依赖副本。

### C++ 与 CLI

保留全部现有 CTest。新增测试覆盖 CLI 参数解析、ENV/JSON 输入、线程数、压力/速度 SHD、MOD、输出目录、帮助信息和错误退出码。

### Python

测试覆盖：

- 原生模块导入、枚举和值类型。
- ENV/JSON 加载、运行、清理和释放生命周期。
- 压力、速度、模态数组的形状、dtype、所有权和数值一致性。
- SHD/MOD reader 的有效文件、截断文件和非法维度。
- 高层参数校验与异常类型。
- Matplotlib 无界面后端下的绘图函数。

### MATLAB

使用当前环境可用的 `matlab -batch` 实际执行：

- 数据模型写 JSON。
- 调用正式 CLI。
- 读取压力 SHD 与 MOD。
- 检查维度、有限值和关键数值。
- 使用无界面方式生成压力场图和模态图。

### 数值回归

运行 OOK 现有完整 CTest 套件。Python 和 MATLAB 读取的关键结果与同一 fixture 的 C++ 输出比较，确保 wrapper 和构建改造不改变核心计算结果。

## 完成标准

以下条件全部满足才视为完成：

- 三个第三方库位于 `third-party/`，旧依赖目录已移除。
- CMake 可以分别构建核心库、CLI、Python 模块和测试。
- Python 高层包可安装、导入、运行算例、读取内存结果并绘图。
- MATLAB 可写配置、调用 CLI、读取 MOD/SHD 并绘图。
- 文档和 demo 给出可直接执行的最小示例。
- 新增测试和现有回归测试全部通过。
