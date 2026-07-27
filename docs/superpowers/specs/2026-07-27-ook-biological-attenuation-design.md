# OOK 生物衰减兼容性设计

**日期：** 2026-07-27

**状态：** 已确认设计

**参考实现：** `krakenFortran`

**范围：** C++ 核心、ENV、JSON、公开 C++ API、Python、MATLAB、回归算例、Fortran 数值差分和 TL 验收的端到端支持。

## 1. 目标

为 OOK 增加 `B` 生物衰减选项，使其在有效输入范围内与 `krakenFortran` 保持相同的数值语义。

实现必须满足：

- 保持 Fortran 生物衰减公式及运算顺序；
- 保持深度闭区间判定和按输入顺序累加重叠层；
- 通过 `env_in_out.cpp` 和 `json_in_out.hpp` 完整传输生物层配置；
- 通过 C++、Python 和 MATLAB 对外暴露该配置；
- 在修改活动参数前拒绝格式错误或物理上未定义的输入；
- 保持现有弱损耗复声速换算和 SSP 插值行为；
- 从公式、模态波数和传播损失三个层面对照 `krakenFortran` 验证。

即使核心公式已经实现，只要 ENV 或 JSON 任一传输接口缺失，该功能仍视为未完成。

## 2. 已确认决策

以下决策已经明确确认：

1. 采用端到端实现，不做只修改核心公式的局部补丁。
2. 有效输入保持 Fortran 行为，非法输入在 OOK 中快速失败。
3. 生物衰减作为所有剖面共享的一套全局衰减配置。
4. 多剖面的衰减配置不一致时明确拒绝。
5. 扩展现有 `Atten_Mode`，不增加平行状态。
6. 接受 C++ ABI 变化，统一重新构建全部 OOK 原生二进制产物。
7. 不修改弱损耗复声速公式。
8. 不自动插入生物层上下边界节点。
9. `env_in_out.cpp` 和 `json_in_out.hpp` 是两个必须独立验收的交付项。

## 3. 当前兼容性缺口

`krakenFortran` 使用 Top Option 的第四个字符选择附加体积衰减：

- `T`：Thorp；
- `F`：Francois-Garrison；
- `B`：生物衰减。

OOK 当前只支持 `None`、`Thorpe` 和 `FrancGarr`。ENV 第四字符为 `B` 时会落入 `None`，并且后续生物层记录不会被消费。下一行随后会被错误解释为顶部半空间参数或介质头部。

OOK 当前也没有能够表达生物层的：

- JSON 模型；
- 公开参数类型；
- Python 绑定；
- MATLAB 配置入口。

## 4. 总体架构

### 4.1 单一所有权

继续由 `OOK_parameters::AttenUnit` 唯一持有衰减配置。

不得仿照 Fortran 引入模块级全局 `bio[]` 或 `NBioLayers`。这样可以保持 OOK 当前的对象参数传递方式，并维持：

- 线程安全；
- 可重入性；
- 不同求解器实例之间的状态隔离。

### 4.2 公开数据模型

在 `include/OpenOceanKrakenParams.h` 的现有衰减类型附近新增：

```cpp
inline constexpr std::size_t MaxBioLayers = 200;

struct BiologicalAttenuationLayer
{
    double Z1 = 0.0;  // m
    double Z2 = 0.0;  // m
    double f0 = 0.0;  // Hz
    double Q  = 0.0;  // 无量纲
    double a0 = 0.0;  // dB/km
};
```

新枚举值追加在末尾，不插入现有枚举值之间：

```cpp
enum class OceanAbsorptionModel
{
    None,
    Thorpe,
    FrancGarr,
    Biological
};
```

原子扩展现有配置：

```cpp
struct Atten_Mode
{
    AttenuationUnit attnUnit =
        AttenuationUnit::MODE_W_db_per_lambda;

    OceanAbsorptionModel absModel =
        OceanAbsorptionModel::None;

    std::vector<BiologicalAttenuationLayer> biologicalLayers;
};
```

生物层数量始终由 `biologicalLayers.size()` 得到，不再保存独立的 `NBioLayers` 字段。

### 4.3 多剖面全局语义

生物层数组属于全局 `Atten_Mode`，不属于单个 SSP 介质，也不下沉到 `Range_Independent_Area`。

多剖面输入的以下内容必须完全一致：

- 衰减单位；
- 吸收模型；
- 生物层数量；
- 每一层的参数；
- 生物层排列顺序。

发现差异时必须拒绝，不能继续沿用“只保留第一个剖面配置”的静默行为。

### 4.4 函数和数据流

保持现有主链路：

```text
OOK_parameters::AttenUnit
    → input_SSP
    → UpdateSSPLoss
    → CRCI
    → parseAttenuation
    → addOceanAbsorption
    → 复数 cp/cs
```

`Atten_Mode` 在损耗更新和衰减计算链路中统一通过 `const Atten_Mode&` 传递，避免复制内部 `vector`。

继续使用 `Interface::set_AttenUnit()` 作为 C++ 和 Python native 的主要配置入口，不增加可独立修改的第二套生物衰减状态。

## 5. ENV 传输接口

### 5.1 标准格式

Top Option 第四字符 `B` 启用生物衰减：

```text
'CVWB'
2
10.0 30.0 1000.0 5.0 0.04
40.0 60.0 1200.0 4.0 0.02
```

记录含义为：

```text
NBioLayers
Z1 Z2 f0 Q a0
... 每个生物层一行
```

单位为：

- `Z1`、`Z2`：米；
- `f0`：赫兹；
- `Q`：无量纲；
- `a0`：dB/km。

### 5.2 `env_in_out.cpp` 必须修改的内容

`src/module/env_in_out.cpp` 必须完成：

1. 在第四字符分支中显式增加 `case 'B'`。
2. 将 `absModel` 设为 `OceanAbsorptionModel::Biological`。
3. 在解析 Top Option 后立即读取层数和生物层记录。
4. 必须在读取 `A` 顶部半空间参数或内部介质头部之前完成上述读取。
5. 每个生物层必须严格解析为五个数值。
6. 先解析到临时 `Atten_Mode`。
7. 对完整临时配置执行统一校验。
8. 只有解析和校验全部成功后，才能写入 `params.AttenUnit`。
9. 错误信息必须包含生物层序号和错误字段。
10. 失败时通过现有 ENV 异常通道使 `from_env()` 返回 `false`。
11. 合并多剖面时比较完整衰减配置。

该文件当前只提供 ENV 读取功能，本功能不新增 ENV 写出接口。

### 5.3 ENV 兼容规则

- 缺少第四字符继续表示 `None`。
- `T` 和 `F` 保持现有 OOK 行为。
- `B` 必须读取紧随其后的生物层记录。
- 未知且非空的第四字符必须明确报错。
- `NBioLayers=0` 合法，对应无生物附加衰减。
- 层数超过 200 时拒绝。
- 缺行、字段不足、字段过多或数值非法时拒绝。

必须分别验证：

- Biological + 普通顶部边界；
- Biological + `A` 顶部半空间。

## 6. JSON 传输接口

### 6.1 标准结构

标准 JSON 表达为：

```json
{
  "AttenUnit": {
    "AttenuationUnit": "dB/lambda",
    "OceanAbsorptionModel": "Biological",
    "BiologicalLayers": [
      {
        "Z1": 10.0,
        "Z2": 30.0,
        "f0": 1000.0,
        "Q": 5.0,
        "a0": 0.04
      }
    ]
  }
}
```

保留现有 `"OceanAbsorptionModel"` 键名，以维持向后兼容。

### 6.2 `json_in_out.hpp` 必须修改的内容

`src/module/json_in_out.hpp` 必须完成：

1. 实现 `OceanAbsorptionModel::Biological` 与字符串 `"Biological"` 的双向转换。
2. 为 `BiologicalAttenuationLayer` 增加 `to_json/from_json`。
3. 扩展 `Atten_Mode` 的转换，加入 `"BiologicalLayers"`。
4. 模型为 Biological 时必须显式存在该数组。
5. 允许显式提供空数组。
6. 非 Biological 模型携带非空生物层数组时拒绝。
7. 层数由数组长度唯一派生。
8. 使用与 ENV 和 C++ setter 相同的结构校验。
9. 保持现有 `OOK_parameters` 顶层结构，不增加平行的顶层生物衰减键。

### 6.3 JSON 兼容规则

- 旧的非 Biological JSON 继续可以读取。
- 非 Biological 输出不写 `"BiologicalLayers"`，避免改变现有标准 JSON 和快照。
- 新 Biological JSON 被旧版 OOK 读取时应因未知 `"Biological"` 而安全失败。
- JSON 中不保存重复的层数。
- ENV `B` → 参数 → JSON → 参数，必须保持全部字段和层顺序不变。

## 7. 公开接口

### 7.1 C++

保留现有 setter：

```cpp
void Interface::set_AttenUnit(Atten_Mode mode);
```

该函数必须先校验传入对象，再修改活动参数。校验失败时抛出 `std::invalid_argument`，原配置保持不变。

内部高频调用链统一接收 `const Atten_Mode&`。

### 7.2 Python

pybind11 必须公开：

- `OceanAbsorptionModel.Biological`；
- `BiologicalAttenuationLayer`；
- 生物层的五个字段；
- `Atten_Mode.biologicalLayers`。

`.pyi` 文件必须同步描述上述接口。

Python 高层封装增加类型化便利入口：

```python
ook_set_biological_attenuation(
    layers,
    attenuation_unit=...
)
```

非法输入映射为 `ValueError`。

### 7.3 MATLAB

`krakenDataModel` 增加 Biological 配置方法，用于生成标准 JSON，并保持生物层输入顺序。

MATLAB 高层校验用于尽早提示；C++ JSON 导入仍是最终权威校验边界。

### 7.4 文档

主 README 和封装 README 必须说明：

- ENV `B` 参数块位置；
- 五个参数的单位；
- JSON 表达；
- Python native 和高层示例；
- MATLAB 示例；
- 闭区间和重叠层语义；
- `a0` 不是共振峰值，因为 `a(f0)=a0·Q²`。

## 8. 数值语义

### 8.1 基础材料衰减

首先使用现有 `N/M/m/F/W/Q/L` 逻辑，将用户输入的材料衰减转换为 Nepers/m。

### 8.2 生物衰减

按照输入顺序处理每个生物层：

\[
I_i(z)=
\begin{cases}
1,& Z_{1i}\le z\le Z_{2i}\\
0,& \text{其他}
\end{cases}
\]

\[
a_i(f,z)=
I_i(z)
\frac{a_{0i}}
{\left(1-\frac{f_{0i}^2}{f^2}\right)^2+\frac{1}{Q_i^2}}
\quad[\mathrm{dB/km}]
\]

每个命中的生物层都执行：

\[
\alpha_T\leftarrow\alpha_T+\frac{a_i}{8685.8896}
\quad[\mathrm{Np/m}]
\]

以下运算顺序是强制要求：

1. 按输入顺序遍历；
2. 使用不带 epsilon 的闭区间判断；
3. 计算当前层分母；
4. 计算当前层 `a`；
5. 将当前层结果除以 `8685.8896`；
6. 加入当前 `alphaT`；
7. 再进入下一层。

实现不得：

- 对生物层排序、合并或去重；
- 先汇总全部 dB 再统一换算；
- 在影响舍入顺序的情况下用预计算倒数替代除法；
- 将共振公式代数改写为 \(f^4\) 形式；
- 将 `a0` 误认为峰值衰减。

共振处：

\[
a(f_0)=a_0Q^2
\]

### 8.3 复声速

保持现有与 Fortran 对齐的弱损耗换算：

\[
\omega=2\pi f
\]

\[
c_I=\frac{\alpha_Tc^2}{\omega}
\]

\[
\widetilde c=c+i\,c_I
\]

本功能不引入所谓精确的逆波数改写。

### 8.4 SSP 与半空间

- 在已存储 SSP 控制点上计算生物衰减。
- P 波和 S 波的 `CRCI` 调用使用同一配置。
- 复声速更新后，再生成线性、PCHIP 或样条系数。
- 不自动插入 `Z1/Z2` 控制点。
- 重叠层全部累加。
- 相邻层共享端点时，该端点计入两层。
- `UpdateHSLoss()` 排除 Biological。
- 当前 `1e8` 半空间深度哨兵改为 `std::numeric_limits<double>::max()`。
- 最终声场继续通过复模态波数和 `exp(-ikr)` 获得衰减，不增加额外的最终压力乘子。

## 9. 校验与失败语义

### 9.1 结构规则

每个生物层必须满足：

```text
Z1、Z2、f0、Q、a0 全部有限
Z1 <= Z2
f0 > 0
Q > 0
a0 >= 0
```

模型级规则：

- `biologicalLayers.size() <= 200`；
- Biological + 空数组合法；
- 非 Biological + 非空数组非法；
- 重叠层合法；
- 共享端点合法；
- `Z1==Z2` 的零厚度层合法；
- `a0==0` 合法；
- 保持层顺序。

### 9.2 运行时规则

更新 SSP 损耗前检查：

- 频率有限且大于零；
- 所有分母有限且大于零；
- 每层贡献有限；
- 累计衰减有限；
- 最终复声速有限。

若换算后的复声速虚部大于实部，必须保留 Fortran 的致命错误语义，使本次求解失败。

流体中的 `cS==0` 仍然合法。

### 9.3 共享校验与原子修改

ENV、JSON、公开 setter 和运行时计算必须使用共享校验逻辑，不得分别复制规则。

所有修改状态的输入路径遵循：

```text
解析到临时对象
    → 校验完整对象
    → 一次性写入活动状态
```

任何失败都不能留下只填充了一部分的 Biological 配置。

### 9.4 错误映射

- ENV：输出诊断，`from_env()==false`；
- JSON：输出诊断，`from_json()==false`；
- C++ setter：`std::invalid_argument`；
- Python：`ValueError`；
- 非有限计算：`std::domain_error` 或 `std::overflow_error`；
- MATLAB：高层尽早提示，底层 C++ 继续作为最终失败边界。

## 10. 测试策略

### 10.1 公式单元测试

使用：

```text
Z1=20 m
Z2=40 m
f=f0=1000 Hz
Q=5
a0=0.04 dB/km
c=1500 m/s
基础衰减=0
```

预期共振衰减：

\[
a_B=1.0\ \mathrm{dB/km}
\]

\[
\alpha_B=\frac{1}{8685.8896}
=1.1512925515\times10^{-4}\ \mathrm{Np/m}
\]

测试必须覆盖：

- 下边界、层内和上边界；
- 紧邻层外的点；
- 重叠层累加；
- 基础衰减与 Biological 相加；
- 空生物层；
- 其他模型不受影响；
- 不同测试之间不存在全局状态残留；
- SSP 包含 Biological、半空间排除 Biological。

### 10.2 ENV 测试

增加正例：

- Biological + 普通顶部边界；
- Biological + `A` 顶部半空间。

增加负例：

- 缺少层数；
- 层数与行数不一致；
- 层数超过 200；
- 生物层行字段不完整；
- 非数值字段；
- 字段值非法；
- 未知非空第四字符；
- 多剖面衰减配置不一致。

### 10.3 JSON 测试

覆盖：

- Biological 枚举转换；
- 生物层对象转换；
- 完整 `Atten_Mode` 往返；
- ENV `B` → JSON → 参数往返；
- Biological 空数组；
- Biological 缺少数组；
- 字段类型错误；
- 非 Biological + 非空数组；
- 数组超过 200 层；
- None、Thorpe、FrancGarr 输出不变。

### 10.4 Python 与 MATLAB 测试

Python 自动测试必须验证：

- native 对象创建；
- 生物层数组赋值；
- 高层便利 setter；
- JSON 往返；
- `ValueError` 映射。

MATLAB 测试必须验证：

- 标准 JSON 生成；
- 字段和数组完整保留；
- OOK CLI 可以读取；
- 非法输入被拒绝。

### 10.5 Fortran 模态差分测试

增加固定算例：

```text
bio_uniform_off
bio_uniform_resonance
```

Biological-on 算例覆盖整个水层，以隔离生物衰减公式与 SSP 边界插值的影响。

使用 `kraken.exe` 作为主要真值实现。

验收条件：

- 模态数一致；
- `|ΔRe(k)| <= 1e-7 m^-1`；
- `|ΔIm(k)| <= max(1e-9, 1e-3·|Im(k_Fortran)|)`；
- Biological-on 时 `Im(k) < 0`；
- Biological-on 与 Biological-off 的波数虚部存在明确差异；
- Fortran 输出包含 `Biological attenuation`；
- 两边都不得报告 fatal error 或 no modes。

新增比较工具必须在超过阈值时返回非零退出码，不能只输出统计信息。

### 10.6 TL 端到端验收

在 1、5、10 km 比较 Biological-on 与 Biological-off。

验收条件：

- 每个实现的衰减斜率误差不超过 `0.02 dB/km`；
- OOK 对 Fortran 的 `P95 |ΔTL| <= 0.20 dB`；
- 单模态算例 `max |ΔTL| <= 0.50 dB`；
- 使用同一 Fortran `field.exe` 读取两份 MOD 时，`P95 <= 0.05 dB`；
- 所有压力样本有限且非零。

### 10.7 回归门禁

完成功能必须同时满足：

```text
公式单元测试
+ env_in_out 测试
+ json_in_out 测试
+ C++/Python/MATLAB 接口测试
+ Fortran 模态差分
+ TL 验收
+ 全部既有 CTest
```

## 11. 向后兼容与构建影响

### 11.1 源码和数据兼容

- 不含 `B` 的现有 ENV 继续支持。
- 现有非 Biological JSON 继续支持。
- MOD 和 SHD 格式不变。
- 正常使用保留字段的现有 C++ 源码保持源代码级兼容。
- 未知非空 ENV 第四字符从静默回退改为明确拒绝。

### 11.2 二进制兼容

增加 `std::vector` 会改变公开 `Atten_Mode` 及包含它的 `OOK_parameters` 的对象布局。

由于不存在必须保持兼容且无法重编译的外部 C++ 二进制客户端，本设计接受该 ABI 变化。

发布该功能时必须统一重新构建：

- OOK 核心库；
- OOK 可执行程序；
- 原生测试和 runner；
- Python `.pyd`；
- 仓库中的其他原生使用方。

纯 Python 文件、MATLAB `.m` 文件、ENV、JSON、MOD 和 SHD 不需要编译。

建议使用新的干净构建目录，避免新旧对象文件混用。

## 12. 非目标

本功能不负责：

- 修改 `krakenFortran`；
- 重构全部吸收模型；
- 修正当前 Francois-Garrison 环境参数读取；
- 改变弱损耗复声速近似；
- 向 SSP 网格自动插入生物层边界；
- 支持每个剖面独立的 Biological 配置；
- 增加 ENV 写出器；
- 修改 MOD 或 SHD 格式；
- 实现超出 Fortran 衰减表达式的生物散射或声速频散。

## 13. 组件范围

预计实施涉及：

- `include/OpenOceanKrakenParams.h`
- `src/algorithm/AttenMod.h`
- `src/algorithm/AttenMod.cpp`
- `src/algorithm/sspMod.h`
- `src/algorithm/sspMod.cpp`
- `src/module/env_in_out.cpp`
- `src/module/json_in_out.hpp`
- `src/module/OpenOceanKrakenInterface.cpp`
- `src/bind/OpenOceanKrakenBind.cpp`
- `wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi`
- `wrappers/py_kraken/py_kraken/ook_data_model.py`
- `wrappers/py_kraken/py_kraken/ook_interface.py`
- `wrappers/py_kraken/py_kraken/__init__.py`
- `wrappers/py_kraken/tests/`
- `wrappers/m_kraken/krakenDataModel.m`
- `wrappers/m_kraken/tests/`
- `for_test/test.cpp`
- `test_for_lcov/option_examples/`
- `test_for_lcov/bad_examples/`
- `test/`
- `tools/test_bio_kraken_equivalence.py`
- `README.md`
- `wrappers/py_kraken/README.md`
- `wrappers/m_kraken/README.md`

## 14. 完成定义

只有同时满足以下条件，功能才算完成：

1. OOK 公开参数模型可以表达 Biological。
2. ENV `B` 参数块在与 Fortran 一致的位置被完整消费。
3. JSON 可以完整往返模型和生物层数组。
4. C++、Python 和 MATLAB 都能配置该功能。
5. OOK 按批准的运算顺序计算 Fortran 公式。
6. 半空间排除 Biological。
7. 非法输入以一致方式原子失败。
8. 多剖面配置不一致时被拒绝。
9. 公式、模态和 TL 验收全部通过。
10. OOK 现有测试保持通过。
