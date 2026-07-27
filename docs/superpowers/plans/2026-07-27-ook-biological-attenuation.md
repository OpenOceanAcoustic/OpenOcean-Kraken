# OOK 生物衰减实施计划

> **面向代理执行者：** 必须使用 `superpowers:subagent-driven-development`（推荐）或 `superpowers:executing-plans`，严格按任务顺序实施；以下步骤使用复选框（`- [ ]`）跟踪。

**目标：** 为 OOK 增加与 `krakenFortran` 数值语义一致的 Biological 体积衰减，并让 C++、ENV、JSON、Python、MATLAB、模态波数和 TL 验收形成完整闭环。

**架构：** 继续由 `OOK_parameters::AttenUnit` 唯一持有全局衰减配置，在现有 `parseAttenuation → addOceanAbsorption → CRCI → UpdateSSPLoss` 链路内叠加 Biological，不在最终压力场增加第二个衰减通道。所有入口共用 `validateAttenuationMode()`，ENV/JSON/公开 setter 均采用“临时对象 → 完整校验 → 一次性提交”，半空间则构造只保留材料衰减的局部配置。

**技术栈：** C++17、CMake/Ninja、Eigen、nlohmann/json、pybind11、Python 3.12/NumPy/Pydantic 2/unittest、MATLAB function tests、KrakenFortran、CTest。

## 全局约束

- `MaxBioLayers` 固定为 `200`，层数只由 `biologicalLayers.size()` 派生，不新增重复计数字段。
- 数据模型固定为 `BiologicalAttenuationLayer{Z1, Z2, f0, Q, a0}`；单位依次为 m、m、Hz、无量纲、dB/km。
- `OceanAbsorptionModel::Biological` 只能追加在现有枚举末尾。
- Biological + 空层数组合法；非 Biological + 非空层数组非法。
- 每层必须满足：五字段全部有限、`Z1 <= Z2`、`f0 > 0`、`Q > 0`、`a0 >= 0`。
- 计算必须按输入顺序、闭区间逐层执行，不排序、不合并、不去重、不预先汇总 dB。
- 每层严格计算
  `a0 / ((1 - f0²/f²)² + 1/Q²)`，随后单独除以 `8685.8896` 并加入 `alphaT`。
- 共振处必须满足 `a(f0)=a0·Q²`；批准用例的结果是 `1.0 dB/km`。
- 保持 `cI = alphaT·c²/(2πf)`，不得改成精确逆波数公式。
- Biological 只作用于已存储 SSP 控制点；不得自动插入 `Z1/Z2`。
- P、S 波使用同一份 Biological 配置；流体 `cS==0` 继续合法。
- `UpdateHSLoss()` 必须排除 Biological，并将深度哨兵改为 `std::numeric_limits<double>::max()`。
- 多剖面 ENV 的衰减单位、模型、层数、每层字段及层顺序必须完全相同，否则拒绝。
- ENV `B` 块必须在顶部 `A` 半空间行或第一介质头之前被完整消费。
- JSON 的键固定为 `"OceanAbsorptionModel"` 和 `"BiologicalLayers"`；Biological 必须显式提供数组，空数组合法。
- 非 Biological JSON 输出不得写 `"BiologicalLayers"`，以保持现有快照。
- ENV、JSON、C++ setter 失败后，活动 `OOK_parameters` 不得发生部分修改。
- 不修改 `run.cpp`、`field.cpp`、MOD/SHD 格式、Francois-Garrison 环境参数行为或 ENV 写出能力。
- 公开 `Atten_Mode` 布局变化已获批准；发布前必须在新构建目录重编译核心库、CLI、共享库、runner、测试和 Python `.pyd`。
- 当前可验证工具链是 MinGW UCRT64 + Ninja；本功能不顺带处理现有 MSVC 编译选项问题。
- 现有用户修改 `wrappers/m_kraken/demo_env.m` 不属于本功能，所有提交都必须避开该文件。
- 外层 `../test`、`../test_for_lcov`、`../tools` 和 `../krakenFortran` 不属于 OOK Git 仓库；新增固定算例和 gate 工具必须存入 OOK 仓库内，Fortran 可执行文件只通过命令行参数引用。

---

## 文件结构

### 修改

- `include/OpenOceanKrakenParams.h`：公开层类型、枚举值和 `Atten_Mode` vector。
- `src/algorithm/AttenMod.h`：共享校验、比较、Biological 判定和新吸收函数签名。
- `src/algorithm/AttenMod.cpp`：共享校验、Fortran Biological 公式、有限性检查和复声速致命条件。
- `src/algorithm/sspMod.h`、`src/algorithm/sspMod.cpp`：`const Atten_Mode&` 传递、SSP 应用和半空间排除。
- `src/module/input_Freq.hpp`：默认化时整体清空 `Atten_Mode`。
- `src/module/input_SSP.hpp`：公开 setter 校验后提交。
- `src/module/env_in_out.cpp`：严格解析 `B` 块、校验、多剖面一致性。
- `src/module/json_in_out.hpp`：Biological 枚举、层对象、条件数组和校验。
- `src/module/OpenOceanKrakenInterface.cpp`：ENV/JSON candidate 原子导入。
- `src/bind/OpenOceanKrakenBind.cpp`：pybind11 Biological 类型和字段。
- `wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi`：原生类型存根。
- `wrappers/py_kraken/py_kraken/ook_data_model.py`：高层 Biological 层模型。
- `wrappers/py_kraken/py_kraken/ook_interface.py`：便利 setter。
- `wrappers/py_kraken/py_kraken/__init__.py`：导出高层层模型。
- `wrappers/py_kraken/tests/test_native.py`、`wrappers/py_kraken/tests/test_interface.py`：Python 回归。
- `wrappers/m_kraken/krakenDataModel.m`、`wrappers/m_kraken/tests/test_m_kraken.m`：MATLAB 配置和回归。
- `for_test/test.cpp`：C++ 公式、SSP、ENV、JSON、原子失败测试。
- `CMakeLists.txt`：注册 Python 工具单测和可选 Fortran/TL gate。
- `.gitignore`：只放行新增 Biological fixtures 和工具。
- `README.md`、`wrappers/py_kraken/README.md`、`wrappers/m_kraken/README.md`：端到端用法和单位。

### 新建

- `for_test/fixtures/biological/bio_uniform_off.env`
- `for_test/fixtures/biological/bio_uniform_off.flp`
- `for_test/fixtures/biological/bio_uniform_resonance.env`
- `for_test/fixtures/biological/bio_uniform_resonance.flp`
- `tools/__init__.py`
- `tools/test_bio_kraken_equivalence.py`
- `tools/tests/__init__.py`
- `tools/tests/test_bio_kraken_equivalence.py`

---

### 任务 1：公开数据模型、共享校验和原子 C++ setter

**文件：**

- 修改：`include/OpenOceanKrakenParams.h:86-110`
- 修改：`src/algorithm/AttenMod.h:11-29`
- 修改：`src/algorithm/AttenMod.cpp:3-5,117-125`
- 修改：`src/module/input_Freq.hpp:64-77`
- 修改：`src/module/input_SSP.hpp:142-145`
- 测试：`for_test/test.cpp:35-93,333-363,1905-1945`

**接口：**

- 产出：`MaxBioLayers`、`BiologicalAttenuationLayer`、`OceanAbsorptionModel::Biological`、`Atten_Mode::biologicalLayers`
- 产出：`void validateAttenuationMode(const Atten_Mode&)`
- 产出：`void validateAttenuationFrequency(double)`
- 产出：`bool attenuationModesEqual(const Atten_Mode&, const Atten_Mode&) noexcept`
- 产出：`bool isBiological(const Atten_Mode&)`
- 消费：后续 ENV、JSON、setter、SSP 和数值计算全部使用这些名称，不得另建平行校验器。

- [ ] **步骤 1：先写共享校验和 setter 原子性的失败测试**

在 `for_test/test.cpp` 增加 `using` 和测试函数，并在 `main()` 的衰减测试之后调用：

```cpp
using OpenOceanKraken::BiologicalAttenuationLayer;
using OpenOceanKraken::MaxBioLayers;
using OpenOceanKraken::attenuationModesEqual;
using OpenOceanKraken::validateAttenuationMode;

void test_biological_attenuation_validation(ook_test::TestRunner &test)
{
    Atten_Mode valid;
    valid.absModel = OceanAbsorptionModel::Biological;
    valid.biologicalLayers = {{20.0, 40.0, 1000.0, 5.0, 0.04}};
    validateAttenuationMode(valid);

    Atten_Mode empty = valid;
    empty.biologicalLayers.clear();
    validateAttenuationMode(empty);

    Atten_Mode copy = valid;
    test.require(attenuationModesEqual(valid, copy),
                 "identical Biological modes must compare equal");
    std::swap(copy.biologicalLayers[0].Z1, copy.biologicalLayers[0].Z2);
    test.require(!attenuationModesEqual(valid, copy),
                 "layer field or order changes must compare unequal");

    auto require_invalid = [&](Atten_Mode candidate, const std::string &name) {
        test.requireThrows(
            [&] { validateAttenuationMode(candidate); },
            name + " must be rejected");
    };

    Atten_Mode candidate = valid;
    candidate.biologicalLayers[0].Z1 =
        std::numeric_limits<double>::quiet_NaN();
    require_invalid(candidate, "non-finite Z1");
    candidate = valid;
    candidate.biologicalLayers[0].Z1 = 41.0;
    require_invalid(candidate, "Z1 greater than Z2");
    candidate = valid;
    candidate.biologicalLayers[0].f0 = 0.0;
    require_invalid(candidate, "non-positive f0");
    candidate = valid;
    candidate.biologicalLayers[0].Q = 0.0;
    require_invalid(candidate, "non-positive Q");
    candidate = valid;
    candidate.biologicalLayers[0].a0 = -0.01;
    require_invalid(candidate, "negative a0");
    candidate = valid;
    candidate.biologicalLayers.resize(MaxBioLayers + 1, valid.biologicalLayers[0]);
    require_invalid(candidate, "more than 200 layers");
    candidate = valid;
    candidate.absModel = OceanAbsorptionModel::Thorpe;
    require_invalid(candidate, "non-Biological model with layers");

    Interface iface;
    Atten_Mode original;
    original.absModel = OceanAbsorptionModel::Thorpe;
    iface.set_AttenUnit(original);
    test.requireThrows(
        [&] { iface.set_AttenUnit(candidate); },
        "invalid C++ setter input must throw");
    test.require(
        attenuationModesEqual(iface.getParams_const().AttenUnit, original),
        "failed C++ setter must preserve the original attenuation mode");
}
```

在 `main()` 中加入：

```cpp
test_biological_attenuation_validation(test);
```

- [ ] **步骤 2：运行测试，确认因模型和校验接口尚不存在而失败**

运行：

```powershell
cmake -S . -B build -G Ninja -DBUILD_TESTING=ON -DBUILD_AS_EXE=ON -DBUILD_AS_PYTHON=OFF
cmake --build build --target OpenOceanKraken_interface_tests
```

预期：编译失败，错误包含 `BiologicalAttenuationLayer`、`Biological` 或 `validateAttenuationMode` 未定义。

- [ ] **步骤 3：扩展公开参数模型**

在 `OpenOceanKrakenParams.h` 增加 `<cstddef>`，并将衰减模型定义改为：

```cpp
inline constexpr std::size_t MaxBioLayers = 200;

struct BiologicalAttenuationLayer
{
    double Z1 = 0.0;
    double Z2 = 0.0;
    double f0 = 0.0;
    double Q = 0.0;
    double a0 = 0.0;
};

enum class OceanAbsorptionModel
{
    None,
    Thorpe,
    FrancGarr,
    Biological
};

struct Atten_Mode
{
    AttenuationUnit attnUnit =
        AttenuationUnit::MODE_W_db_per_lambda;
    OceanAbsorptionModel absModel =
        OceanAbsorptionModel::None;
    std::vector<BiologicalAttenuationLayer> biologicalLayers;
};
```

- [ ] **步骤 4：声明并实现唯一的结构校验和精确比较**

在 `AttenMod.h` 声明：

```cpp
void validateAttenuationMode(const Atten_Mode &mode);
void validateAttenuationFrequency(double freq);
bool attenuationModesEqual(
    const Atten_Mode &lhs,
    const Atten_Mode &rhs) noexcept;
bool isBiological(const Atten_Mode &mode);
```

在 `AttenMod.cpp` 增加 `<stdexcept>`、`<string>`，然后实现。错误消息必须包含数组索引和字段名：

```cpp
namespace
{
    [[noreturn]] void invalidLayer(
        std::size_t index,
        const char *field,
        const char *requirement)
    {
        throw std::invalid_argument(
            "BiologicalLayers[" + std::to_string(index) + "]." +
            field + " " + requirement);
    }
}

void validateAttenuationMode(const Atten_Mode &mode)
{
    if (mode.biologicalLayers.size() > MaxBioLayers)
        throw std::invalid_argument(
            "BiologicalLayers must contain at most 200 layers.");

    if (!isBiological(mode) && !mode.biologicalLayers.empty())
        throw std::invalid_argument(
            "BiologicalLayers must be empty unless "
            "OceanAbsorptionModel is Biological.");

    for (std::size_t i = 0; i < mode.biologicalLayers.size(); ++i)
    {
        const auto &layer = mode.biologicalLayers[i];
        if (!std::isfinite(layer.Z1))
            invalidLayer(i, "Z1", "must be finite.");
        if (!std::isfinite(layer.Z2))
            invalidLayer(i, "Z2", "must be finite.");
        if (!std::isfinite(layer.f0))
            invalidLayer(i, "f0", "must be finite.");
        if (!std::isfinite(layer.Q))
            invalidLayer(i, "Q", "must be finite.");
        if (!std::isfinite(layer.a0))
            invalidLayer(i, "a0", "must be finite.");
        if (layer.Z1 > layer.Z2)
            invalidLayer(i, "Z1", "must be less than or equal to Z2.");
        if (layer.f0 <= 0.0)
            invalidLayer(i, "f0", "must be greater than zero.");
        if (layer.Q <= 0.0)
            invalidLayer(i, "Q", "must be greater than zero.");
        if (layer.a0 < 0.0)
            invalidLayer(i, "a0", "must be greater than or equal to zero.");
    }
}

void validateAttenuationFrequency(double freq)
{
    if (!std::isfinite(freq) || freq <= 0.0)
        throw std::domain_error(
            "Attenuation frequency must be finite and greater than zero.");
}

bool attenuationModesEqual(
    const Atten_Mode &lhs,
    const Atten_Mode &rhs) noexcept
{
    if (lhs.attnUnit != rhs.attnUnit ||
        lhs.absModel != rhs.absModel ||
        lhs.biologicalLayers.size() != rhs.biologicalLayers.size())
        return false;

    for (std::size_t i = 0; i < lhs.biologicalLayers.size(); ++i)
    {
        const auto &a = lhs.biologicalLayers[i];
        const auto &b = rhs.biologicalLayers[i];
        if (a.Z1 != b.Z1 || a.Z2 != b.Z2 || a.f0 != b.f0 ||
            a.Q != b.Q || a.a0 != b.a0)
            return false;
    }
    return true;
}

bool isBiological(const Atten_Mode &mode)
{
    return mode.absModel == OceanAbsorptionModel::Biological;
}
```

- [ ] **步骤 5：让默认化和 C++ setter 使用完整对象语义**

在 `input_Freq::setMunk()` 中用整体重置替换两个字段赋值：

```cpp
params.AttenUnit = Atten_Mode{};
```

在 `input_SSP::set_AttenUnit()` 中先校验再提交：

```cpp
void set_AttenUnit(OOK_parameters &params, Atten_Mode unit) const
{
    validateAttenuationMode(unit);
    params.AttenUnit = std::move(unit);
}
```

为 `input_SSP.hpp` 增加 `#include "AttenMod.h"` 和 `<utility>`。保持 `Interface::set_AttenUnit(Atten_Mode)` 的公开签名不变。

- [ ] **步骤 6：运行原生测试并确认通过**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests
ctest --test-dir build -R OpenOceanKraken_interface_tests --output-on-failure
```

预期：`OpenOceanKraken_interface_tests` 通过；默认 None/Thorpe/FrancGarr 行为不变。

- [ ] **步骤 7：提交任务 1**

```powershell
git add include/OpenOceanKrakenParams.h src/algorithm/AttenMod.h src/algorithm/AttenMod.cpp src/module/input_Freq.hpp src/module/input_SSP.hpp for_test/test.cpp
git commit -m "feat: add biological attenuation data model"
```

---

### 任务 2：Fortran Biological 公式、SSP 数据流和半空间排除

**文件：**

- 修改：`src/algorithm/AttenMod.h:11-16`
- 修改：`src/algorithm/AttenMod.cpp:5-80`
- 修改：`src/algorithm/sspMod.h:17-20`
- 修改：`src/algorithm/sspMod.cpp:477-574`
- 测试：`for_test/test.cpp:333-363,857-877,900-918`

**接口：**

- 修改产出：`double addOceanAbsorption(double alphaT, double z, double freq, const Atten_Mode&)`
- 修改产出：`void UpdateSSPLoss(double, double, int, SSP_Mode, const Atten_Mode&, SSPStructure&)`
- 修改产出：`void UpdateHSLoss(double&, double&, const Atten_Mode&, HSInfo&, HSInfo&)`
- 保持：`CRCI()` 的公开签名和弱损耗换算公式。

- [ ] **步骤 1：写公式、边界、重叠、复声速和半空间失败测试**

扩展现有 `test_attenuation_and_reflection_helpers()`：

```cpp
Atten_Mode biological;
biological.absModel = OceanAbsorptionModel::Biological;
biological.biologicalLayers = {{20.0, 40.0, 1000.0, 5.0, 0.04}};

constexpr double expected_np_per_m =
    1.0 / 8685.8896;
test.requireNear(
    addOceanAbsorption(0.0, 20.0, 1000.0, biological),
    expected_np_per_m, 1.0e-15, "lower boundary must be included");
test.requireNear(
    addOceanAbsorption(0.0, 30.0, 1000.0, biological),
    expected_np_per_m, 1.0e-15, "layer interior mismatch");
test.requireNear(
    addOceanAbsorption(0.0, 40.0, 1000.0, biological),
    expected_np_per_m, 1.0e-15, "upper boundary must be included");
test.requireNear(
    addOceanAbsorption(0.0, 19.999, 1000.0, biological),
    0.0, 1.0e-15, "point below layer must be excluded");
test.requireNear(
    addOceanAbsorption(0.0, 40.001, 1000.0, biological),
    0.0, 1.0e-15, "point above layer must be excluded");

Atten_Mode empty_biological = biological;
empty_biological.biologicalLayers.clear();
test.requireNear(
    addOceanAbsorption(0.25, 30.0, 1000.0, empty_biological),
    0.25, 0.0, "empty Biological layers must add no attenuation");

Atten_Mode zero_thickness = biological;
zero_thickness.biologicalLayers = {
    {30.0, 30.0, 1000.0, 5.0, 0.04},
    {20.0, 40.0, 1000.0, 5.0, 0.0}};
test.requireNear(
    addOceanAbsorption(0.0, 30.0, 1000.0, zero_thickness),
    expected_np_per_m, 1.0e-15,
    "zero-thickness layers are valid and a0 zero adds nothing");

Atten_Mode overlapping = biological;
overlapping.biologicalLayers.push_back(
    {30.0, 50.0, 1000.0, 5.0, 0.04});
test.requireNear(
    addOceanAbsorption(0.25, 30.0, 1000.0, overlapping),
    0.25 + 2.0 * expected_np_per_m,
    1.0e-15,
    "overlapping layers and shared endpoints must accumulate");

double bio_z = 30.0;
double bio_c = 1500.0;
double bio_alpha = 0.0;
double bio_freq = 1000.0;
double bio_freq0 = 1000.0;
double bio_beta = 1.0;
double bio_ft = 2000.0;
const auto bio_crci = CRCI(
    bio_z, bio_c, bio_alpha, bio_freq, bio_freq0,
    biological, bio_beta, bio_ft);
test.requireNear(
    std::imag(bio_crci),
    0.041227627617643738,
    1.0e-12,
    "weak-loss Biological complex sound speed mismatch");

test.requireThrows(
    [&] {
        double zero_frequency = 0.0;
        CRCI(bio_z, bio_c, bio_alpha, zero_frequency, bio_freq0,
             biological, bio_beta, bio_ft);
    },
    "zero attenuation frequency must fail");

Atten_Mode overflowing = biological;
overflowing.biologicalLayers[0].f0 =
    std::numeric_limits<double>::max();
test.requireThrows(
    [&] { addOceanAbsorption(0.0, 30.0, 1000.0, overflowing); },
    "non-finite Biological denominator must fail");
overflowing = biological;
overflowing.biologicalLayers[0].a0 =
    std::numeric_limits<double>::max();
test.requireThrows(
    [&] { addOceanAbsorption(0.0, 30.0, 1000.0, overflowing); },
    "non-finite Biological contribution must fail");

Atten_Mode fatal_loss = biological;
fatal_loss.biologicalLayers[0].a0 = 1.0e8;
test.requireThrows(
    [&] {
        CRCI(bio_z, bio_c, bio_alpha, bio_freq, bio_freq0,
             fatal_loss, bio_beta, bio_ft);
    },
    "imaginary sound speed greater than real sound speed must fail");

double fluid_shear_speed = 0.0;
const auto fluid_shear = CRCI(
    bio_z, fluid_shear_speed, bio_alpha, bio_freq, bio_freq0,
    biological, bio_beta, bio_ft);
test.requireNear(
    std::abs(fluid_shear), 0.0, 0.0,
    "fluid cS zero must remain valid");
```

更新原有 Thorpe/FrancGarr 调用，使其显式传 `z=0.0`。

在同一测试中加入独立实例检查，防止引入 Fortran 式全局状态：

```cpp
Interface biological_instance;
Interface default_instance;
biological_instance.set_AttenUnit(biological);
test.require(
    biological_instance.getParams_const().AttenUnit.absModel ==
        OceanAbsorptionModel::Biological,
    "first instance must retain Biological mode");
test.require(
    default_instance.getParams_const().AttenUnit.absModel ==
            OceanAbsorptionModel::None &&
        default_instance.getParams_const()
            .AttenUnit.biologicalLayers.empty(),
    "Biological state must not leak between Interface instances");
```

在现有 SSP/Half-space 测试旁增加：

```cpp
Atten_Mode bio_halfspace = biological;
HSInfo bio_top = top;
HSInfo bio_bottom = bot;
HSInfo baseline_top = top;
HSInfo baseline_bottom = bot;
Atten_Mode no_volume;

OpenOceanKraken::ssp::UpdateHSLoss(
    freq, freq0, bio_halfspace, bio_top, bio_bottom);
OpenOceanKraken::ssp::UpdateHSLoss(
    freq, freq0, no_volume, baseline_top, baseline_bottom);
test.requireNear(
    std::imag(bio_top.cp), std::imag(baseline_top.cp), 1.0e-15,
    "top half-space must exclude Biological attenuation");
test.requireNear(
    std::imag(bio_bottom.cp), std::imag(baseline_bottom.cp), 1.0e-15,
    "bottom half-space must exclude Biological attenuation");
```

把现有 `ssp_fixture()` 的 PCHIP/spline 测试扩展为：

```cpp
auto bio_ssp = ssp_fixture();
auto baseline_ssp = ssp_fixture();
const Eigen::Index stored_points = bio_ssp.z.size();
Atten_Mode whole_water_bio;
whole_water_bio.absModel = OceanAbsorptionModel::Biological;
whole_water_bio.biologicalLayers = {
    {0.0, 100.0, 50.0, 5.0, 0.04}};
Atten_Mode no_volume_ssp;

OpenOceanKraken::ssp::UpdateSSPLoss(
    50.0, 50.0, baseline_ssp.NMedia,
    SSP_Mode::MODE_P_cPCHIP, no_volume_ssp, baseline_ssp);
OpenOceanKraken::ssp::UpdateSSPLoss(
    50.0, 50.0, bio_ssp.NMedia,
    SSP_Mode::MODE_P_cPCHIP, whole_water_bio, bio_ssp);
test.require(
    bio_ssp.z.size() == stored_points,
    "Biological attenuation must not insert SSP boundary nodes");
test.require(
    std::imag(bio_ssp.cp(1)) > std::imag(baseline_ssp.cp(1)) &&
        std::imag(bio_ssp.cs(1)) > std::imag(baseline_ssp.cs(1)),
    "P and non-zero S speed must use the same Biological mode");
test.require(
    bio_ssp.cpCoef.allFinite() && bio_ssp.csCoef.allFinite(),
    "PCHIP coefficients must be rebuilt after Biological loss");

OpenOceanKraken::ssp::UpdateSSPLoss(
    50.0, 50.0, bio_ssp.NMedia,
    SSP_Mode::MODE_S_cCubic, whole_water_bio, bio_ssp);
test.require(
    bio_ssp.cpSpline.allFinite() && bio_ssp.csSpline.allFinite(),
    "spline coefficients must be rebuilt after Biological loss");
```

- [ ] **步骤 2：运行测试，确认新签名和数值行为尚未实现**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests
```

预期：编译因 `addOceanAbsorption` 参数数量不匹配而失败；修改声明但未实现时，共振或半空间断言失败。

- [ ] **步骤 3：按 Fortran 源表达式实现 Biological**

把声明和定义改为：

```cpp
double addOceanAbsorption(
    double alphaT,
    double z,
    double freq,
    const Atten_Mode &mode);
```

在 Thorpe/FrancGarr 之后增加独立 Biological 分支，保持以下语句顺序：

```cpp
else if (isBiological(mode))
{
    for (const auto &layer : mode.biologicalLayers)
    {
        if (z >= layer.Z1 && z <= layer.Z2)
        {
            const double denominator =
                std::pow(
                    1.0 -
                        std::pow(layer.f0, 2) /
                            std::pow(freq, 2),
                    2) +
                1.0 / std::pow(layer.Q, 2);
            if (!std::isfinite(denominator) || denominator <= 0.0)
                throw std::domain_error(
                    "Biological attenuation denominator must be "
                    "finite and greater than zero.");

            double contribution = layer.a0 / denominator;
            if (!std::isfinite(contribution))
                throw std::overflow_error(
                    "Biological attenuation contribution is not finite.");

            contribution = contribution / 8685.8896;
            alphaT = alphaT + contribution;
            if (!std::isfinite(alphaT))
                throw std::overflow_error(
                    "Accumulated attenuation is not finite.");
        }
    }
}
```

不得将 `std::pow(layer.f0, 2) / std::pow(freq, 2)` 改写成 `f0/f` 的预计算平方，也不得预计算 `1/8685.8896`。

- [ ] **步骤 4：补齐 CRCI 运行时校验和 Fortran 致命条件**

`CRCI()` 按以下顺序更新：

```cpp
validateAttenuationFrequency(freq);
double alphaT =
    parseAttenuation(freq, freq0, alpha, ft, beta, c, AttenUnit);
if (!std::isfinite(alphaT))
    throw std::overflow_error(
        "Base material attenuation is not finite.");

alphaT = addOceanAbsorption(alphaT, z, freq, AttenUnit);
const double omega = 2.0 * pi * freq;
const double cImag = alphaT * c * c / omega;
if (!std::isfinite(c) || !std::isfinite(cImag))
    throw std::overflow_error(
        "Complex sound speed is not finite.");
if (cImag > c)
    throw std::domain_error(
        "The complex sound speed has an imaginary part "
        "greater than its real part.");
return std::complex<double>(c, cImag);
```

`c==0` 时 `cImag==0`，不得将其当作错误。

- [ ] **步骤 5：按 const 引用更新 SSP，并在半空间移除 Biological**

在 `sspMod.h/.cpp` 将两个函数的 `Atten_Mode` 参数改为 `const Atten_Mode&`，并在 `sspMod.cpp` 增加 `<limits>`。`UpdateSSPLoss()` 在任何 SSP 数组写入之前调用：

```cpp
validateAttenuationMode(AttenUnit);
validateAttenuationFrequency(freq);
```

保留当前“先逐控制点更新 P/S 复声速，再生成 PCHIP/spline 系数”的顺序。

`UpdateHSLoss()` 使用：

```cpp
validateAttenuationMode(AttenUnit);
Atten_Mode halfspaceMode = AttenUnit;
if (isBiological(halfspaceMode))
{
    halfspaceMode.absModel = OceanAbsorptionModel::None;
    halfspaceMode.biologicalLayers.clear();
}
const double huge = std::numeric_limits<double>::max();
```

所有半空间 `CRCI()` 调用传 `halfspaceMode`。Thorpe 和 FrancGarr 半空间行为保持不变。

- [ ] **步骤 6：运行原生回归**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests
ctest --test-dir build -R "OpenOceanKraken_(interface|lcov_fixture)_tests" --output-on-failure
```

预期：公式共振值、闭区间、重叠层、复声速、SSP 和半空间测试通过；既有衰减模型通过。

- [ ] **步骤 7：提交任务 2**

```powershell
git add src/algorithm/AttenMod.h src/algorithm/AttenMod.cpp src/algorithm/sspMod.h src/algorithm/sspMod.cpp for_test/test.cpp
git commit -m "feat: apply biological attenuation to SSP"
```

---

### 任务 3：ENV `B` 传输、多剖面一致性和原子导入

**文件：**

- 修改：`src/module/env_in_out.cpp:16-73,565-751,1048-1115`
- 修改：`src/module/OpenOceanKrakenInterface.cpp:594-604`
- 测试：`for_test/test.cpp:1240-1400,1580-1614,1713-1773`

**接口：**

- 消费：`validateAttenuationMode()`、`attenuationModesEqual()`
- 产出：Top Option 第四字符 `B` 的严格读取。
- 产出：ENV 多剖面共享配置的一致性拒绝。
- 保持：`bool Interface::from_env(const std::string&)`；不增加 ENV writer。

- [ ] **步骤 1：写临时 ENV/FLP 工厂及完整正负测试**

在 `for_test/test.cpp` 增加：

```cpp
std::string biological_env_text(
    const std::string &top_option,
    const std::string &biological_block)
{
    std::ostringstream out;
    out << "'Biological ENV test'\n"
        << "1000.0\n"
        << "1\n"
        << "'" << top_option << "'\n"
        << biological_block;
    if (top_option.size() > 1 && top_option[1] == 'A')
        out << "0.0 1450.0 0.0 1.0 0.0 0.0 /\n";
    out << "1000 0.0 100.0\n"
        << "0.0 1500.0 /\n"
        << "100.0 1500.0 /\n"
        << "'A' 0.0\n"
        << "100.0 1600.0 0.0 1.5 0.0 /\n"
        << "1400.0 20000.0\n"
        << "10.0\n"
        << "1\n"
        << "50.0 /\n"
        << "1\n"
        << "50.0 /\n";
    return out.str();
}

const char *single_profile_flp =
    "/,\n"
    "'RA'\n"
    "1,\n"
    "1\n"
    "0.0 /\n"
    "3\n"
    "1.0 5.0 10.0 /\n"
    "1\n"
    "50.0 /\n"
    "1\n"
    "50.0 /\n"
    "1\n"
    "0.0 /\n";

void write_text(const fs::path &path, const std::string &text)
{
    std::ofstream stream(path);
    stream << text;
    if (!stream)
        throw std::runtime_error("cannot write test fixture " + path.string());
}
```

新增 `test_biological_env_transport()`，用唯一临时目录写同名 `.env/.flp`，覆盖：

```cpp
struct EnvVariant
{
    std::string name;
    std::string option;
    std::string block;
    bool expected;
};

const std::vector<EnvVariant> variants = {
    {"normal_top", "CVWB", "1\n20 40 1000 5 0.04\n", true},
    {"top_halfspace", "CAWB", "1\n20 40 1000 5 0.04\n", true},
    {"empty_layers", "CVWB", "0\n", true},
    {"missing_count", "CVWB", "", false},
    {"negative_count", "CVWB", "-1\n", false},
    {"too_many", "CVWB", "201\n", false},
    {"count_mismatch", "CVWB", "2\n20 40 1000 5 0.04\n", false},
    {"four_fields", "CVWB", "1\n20 40 1000 5\n", false},
    {"six_fields", "CVWB", "1\n20 40 1000 5 0.04 9\n", false},
    {"non_numeric", "CVWB", "1\n20 40 bad 5 0.04\n", false},
    {"invalid_depths", "CVWB", "1\n40 20 1000 5 0.04\n", false},
    {"unknown_option", "CVWX", "", false},
};
```

对成功项检查模型、层字段和 A 半空间的 `alphaR==1450`；对失败项先给接口设置 Thorpe，再确认 `from_env()==false` 且 Thorpe 原配置未变。

另外直接调用 `read_env_file()` 写两个连续 profile：

- 两个完全相同的 `CVWB` 配置必须成功；
- 第二个 profile 分别改变衰减单位、模型、层数、任一五字段或交换层顺序时必须失败。

- [ ] **步骤 2：运行测试并确认当前 `B` 未消费或未知字符静默回退**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests
ctest --test-dir build -R OpenOceanKraken_interface_tests --output-on-failure
```

预期：第一个 `CVWB` 正例失败，或被错误解析为介质头；`CVWX` 当前不会按新要求拒绝。

- [ ] **步骤 3：增加严格整数/浮点读取辅助函数**

在 `env_in_out.cpp` 的 token 辅助函数附近增加：

```cpp
static long long parse_strict_integer(
    const std::string &token,
    const std::string &context)
{
    std::size_t consumed = 0;
    long long value = 0;
    try
    {
        value = std::stoll(token, &consumed);
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(context + " must be an integer.");
    }
    if (consumed != token.size())
        throw std::runtime_error(context + " must be an integer.");
    return value;
}

static double parse_strict_double(
    const std::string &token,
    const std::string &context)
{
    std::size_t consumed = 0;
    double value = 0.0;
    try
    {
        value = std::stod(token, &consumed);
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(context + " must be numeric.");
    }
    if (consumed != token.size())
        throw std::runtime_error(context + " must be numeric.");
    return value;
}
```

新增完整层块读取器：

```cpp
static std::vector<BiologicalAttenuationLayer>
parse_biological_layers(
    const std::vector<std::string> &lines,
    std::size_t &line_index)
{
    if (line_index >= lines.size())
        throw std::runtime_error(
            "Biological attenuation layer count is missing.");

    const auto count_tokens = fortran_tokens(lines[line_index++]);
    if (count_tokens.size() != 1)
        throw std::runtime_error(
            "Biological attenuation layer count must contain one integer.");

    const long long count = parse_strict_integer(
        count_tokens[0], "NBioLayers");
    if (count < 0 || count > static_cast<long long>(MaxBioLayers))
        throw std::runtime_error(
            "NBioLayers must be between 0 and 200.");

    std::vector<BiologicalAttenuationLayer> layers;
    layers.reserve(static_cast<std::size_t>(count));
    for (long long i = 0; i < count; ++i)
    {
        if (line_index >= lines.size())
            throw std::runtime_error(
                "BiologicalLayers[" + std::to_string(i) +
                "] record is missing.");
        const auto tokens = fortran_tokens(lines[line_index++]);
        if (tokens.size() != 5)
            throw std::runtime_error(
                "BiologicalLayers[" + std::to_string(i) +
                "] must contain exactly five fields.");
        const std::string prefix =
            "BiologicalLayers[" + std::to_string(i) + "].";
        layers.push_back({
            parse_strict_double(tokens[0], prefix + "Z1"),
            parse_strict_double(tokens[1], prefix + "Z2"),
            parse_strict_double(tokens[2], prefix + "f0"),
            parse_strict_double(tokens[3], prefix + "Q"),
            parse_strict_double(tokens[4], prefix + "a0")});
    }
    return layers;
}
```

- [ ] **步骤 4：在顶部半空间前解析临时 `Atten_Mode`**

`parse_single_env()` 进入 Top Option 前创建：

```cpp
Atten_Mode candidateAttenuation;
```

第三字符只写 `candidateAttenuation.attnUnit`；第四字符按以下规则写模型：

```cpp
case 'T':
    candidateAttenuation.absModel = OceanAbsorptionModel::Thorpe;
    break;
case 'F':
    candidateAttenuation.absModel = OceanAbsorptionModel::FrancGarr;
    break;
case 'B':
    candidateAttenuation.absModel = OceanAbsorptionModel::Biological;
    candidateAttenuation.biologicalLayers =
        parse_biological_layers(lines, line_idx);
    break;
default:
    throw std::runtime_error(
        std::string("Unknown top option letter in fourth position: ") +
        option4);
```

缺少第四字符时仍设 None。Top Option 解析完成后、进入 `if (HSTop.BC == A)` 之前：

```cpp
validateAttenuationMode(candidateAttenuation);
params.AttenUnit = std::move(candidateAttenuation);
```

- [ ] **步骤 5：拒绝多剖面配置不一致**

在 `read_env_file()` 追加第二及后续剖面之前增加：

```cpp
if (!params.sspInput.empty() &&
    !attenuationModesEqual(params.AttenUnit, single_params.AttenUnit))
{
    std::cerr
        << "ENV profiles must use identical attenuation unit, "
           "absorption model, Biological layer values and layer order."
        << std::endl;
    return false;
}
```

第一个剖面仍完整赋给 `params`，随后清空并统一追加 `sspInput`。

- [ ] **步骤 6：让 `Interface::from_env()` 原子提交**

替换直接写活动参数的实现，并保证构建 SSP 时的异常也映射为 `false`：

```cpp
bool Interface::from_env(const std::string &envPath)
{
    try
    {
        OOK_parameters candidate;
        if (!read_env_file(envPath, candidate))
            return false;
        if (!read_flp_file(envPath, candidate))
            return false;

        this->impl->INPUT_SSP.set_SSP(
            candidate, candidate.sspInput);
        this->getParams() = std::move(candidate);
        markDirty(DirtyKind::All);
        return true;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: Failed to parse ENV file "
                  << envPath << ": " << e.what() << std::endl;
        return false;
    }
}
```

- [ ] **步骤 7：运行 ENV 独立验收**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests OpenOceanKraken
ctest --test-dir build -R "OpenOceanKraken_(interface|lcov_fixture|cli)_tests" --output-on-failure
```

预期：普通顶部、A 顶部、空层、所有负例、等价多剖面和不一致多剖面均符合断言；现有无第四字符/T/F ENV 保持通过。

- [ ] **步骤 8：提交任务 3**

```powershell
git add src/module/env_in_out.cpp src/module/OpenOceanKrakenInterface.cpp for_test/test.cpp
git commit -m "feat: parse biological attenuation from ENV"
```

---

### 任务 4：JSON 双向传输和原子导入

**文件：**

- 修改：`src/module/json_in_out.hpp:179-236,912-1003`
- 修改：`src/module/OpenOceanKrakenInterface.cpp:559-583`
- 测试：`for_test/test.cpp:1473-1577`

**接口：**

- 产出：`"Biological"` 枚举字符串。
- 产出：`BiologicalAttenuationLayer` 的 `to_json/from_json`。
- 产出：`Atten_Mode` 的条件 `"BiologicalLayers"`。
- 保持：顶层 `"AttenUnit"` 结构和 `Interface::from_json()` 返回 `bool`。

- [ ] **步骤 1：先扩展现有 JSON round-trip 失败测试**

在 `test_json_roundtrip_from_params()` 中构造两层 Biological，调用 setter 后检查：

```cpp
Atten_Mode biological;
biological.absModel = OceanAbsorptionModel::Biological;
biological.biologicalLayers = {
    {10.0, 30.0, 1000.0, 5.0, 0.04},
    {40.0, 60.0, 1200.0, 4.0, 0.02}};
params_iface.set_AttenUnit(biological);

const auto biological_json =
    OpenOcean_json::parse(params_iface.to_json_string());
test.require(
    biological_json["AttenUnit"]["OceanAbsorptionModel"] ==
        "Biological",
    "JSON must emit Biological model");
test.require(
    biological_json["AttenUnit"]["BiologicalLayers"].size() == 2,
    "JSON must preserve Biological layer count");
test.requireNear(
    biological_json["AttenUnit"]["BiologicalLayers"][1]["f0"]
        .get<double>(),
    1200.0, 0.0, "JSON must preserve layer order and fields");
```

往返后使用 `attenuationModesEqual()` 比较完整配置。

复用任务 3 的有效 `CVWB` 临时文件，增加 ENV → JSON → 参数往返：

```cpp
const fs::path bio_roundtrip_root =
    fs::temp_directory_path() / "ook_biological_env_json_roundtrip";
fs::create_directories(bio_roundtrip_root);
const fs::path valid_bio_env = bio_roundtrip_root / "case.env";
write_text(
    valid_bio_env,
    biological_env_text(
        "CVWB", "1\n20 40 1000 5 0.04\n"));
write_text(bio_roundtrip_root / "case.flp", single_profile_flp);

Interface env_source;
test.require(
    env_source.from_env(valid_bio_env.string()),
    "Biological ENV must load before JSON round-trip");
const fs::path env_json =
    fs::temp_directory_path() / "ook_biological_env_roundtrip.json";
test.require(
    env_source.to_json(env_json.string()),
    "Biological ENV must export to JSON");
Interface env_json_copy;
test.require(
    env_json_copy.from_json(env_json.string()),
    "Biological ENV-derived JSON must reload");
test.require(
    attenuationModesEqual(
        env_source.getParams_const().AttenUnit,
        env_json_copy.getParams_const().AttenUnit),
    "ENV to JSON must preserve every layer field and order");
fs::remove(env_json);
fs::remove_all(bio_roundtrip_root);
```

补充空数组、缺数组、数组类型错误、字段缺失、字段类型错误、201 层、非 Biological + 非空数组；复用现有 `expect_bad_json()`。每个失败项先给目标接口设置 Thorpe，并确认失败后仍是 Thorpe。

最后把模式改回 Thorpe，确认输出中：

```cpp
test.require(
    !non_biological_json["AttenUnit"].contains("BiologicalLayers"),
    "non-Biological JSON output must omit BiologicalLayers");
```

- [ ] **步骤 2：运行测试并确认 Biological JSON 尚未实现**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests
ctest --test-dir build -R OpenOceanKraken_interface_tests --output-on-failure
```

预期：JSON 输出仍把 Biological 当作未知枚举，或缺少 `BiologicalLayers`。

- [ ] **步骤 3：实现枚举和层对象转换**

在 `OceanAbsorptionModel` 双向转换中增加 `"Biological"`。

紧接其后实现：

```cpp
void to_json(
    OpenOcean_json &out,
    const BiologicalAttenuationLayer &layer)
{
    out = OpenOcean_json{
        {"Z1", layer.Z1},
        {"Z2", layer.Z2},
        {"f0", layer.f0},
        {"Q", layer.Q},
        {"a0", layer.a0}};
}

void from_json(
    const OpenOcean_json &in,
    BiologicalAttenuationLayer &layer)
{
    if (!in.is_object())
        throw std::runtime_error(
            "Each BiologicalLayers entry must be a JSON object.");
    BiologicalAttenuationLayer candidate;
    candidate.Z1 = in.at("Z1").get<double>();
    candidate.Z2 = in.at("Z2").get<double>();
    candidate.f0 = in.at("f0").get<double>();
    candidate.Q = in.at("Q").get<double>();
    candidate.a0 = in.at("a0").get<double>();
    layer = candidate;
}
```

- [ ] **步骤 4：实现条件数组和共享校验**

替换 `Atten_Mode` 转换：

```cpp
void to_json(OpenOcean_json &out, const Atten_Mode &mode)
{
    validateAttenuationMode(mode);
    out = OpenOcean_json{
        {"AttenuationUnit", mode.attnUnit},
        {"OceanAbsorptionModel", mode.absModel}};
    if (isBiological(mode))
        out["BiologicalLayers"] = mode.biologicalLayers;
}

void from_json(const OpenOcean_json &in, Atten_Mode &mode)
{
    if (!in.is_object())
        throw std::runtime_error(
            "Atten_Mode must be a JSON object.");

    Atten_Mode candidate;
    candidate.attnUnit =
        in.at("AttenuationUnit").get<AttenuationUnit>();
    candidate.absModel =
        in.at("OceanAbsorptionModel").get<OceanAbsorptionModel>();

    if (isBiological(candidate))
    {
        if (!in.contains("BiologicalLayers"))
            throw std::runtime_error(
                "Biological AttenUnit requires BiologicalLayers.");
        candidate.biologicalLayers =
            in.at("BiologicalLayers")
                .get<std::vector<BiologicalAttenuationLayer>>();
    }
    else if (in.contains("BiologicalLayers"))
    {
        candidate.biologicalLayers =
            in.at("BiologicalLayers")
                .get<std::vector<BiologicalAttenuationLayer>>();
    }

    validateAttenuationMode(candidate);
    mode = std::move(candidate);
}
```

非 Biological + 显式空数组可读，但规范化输出会省略该键。

- [ ] **步骤 5：让 `Interface::from_json()` 原子提交**

改为：

```cpp
bool Interface::from_json(const std::string &jsonPath)
{
    std::ifstream ifs(jsonPath);
    if (!ifs.is_open())
    {
        std::cerr << "Error: Failed to open JSON file "
                  << jsonPath << std::endl;
        return false;
    }
    try
    {
        OpenOcean_json json;
        ifs >> json;
        OOK_parameters candidate = json.get<OOK_parameters>();
        this->impl->INPUT_SSP.set_SSP(
            candidate, candidate.sspInput);
        this->getParams() = std::move(candidate);
        markDirty(DirtyKind::All);
        return true;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: Failed to parse JSON file "
                  << jsonPath << ": " << e.what() << std::endl;
        return false;
    }
}
```

- [ ] **步骤 6：运行 JSON 独立验收**

运行：

```powershell
cmake --build build --target OpenOceanKraken_interface_tests
ctest --test-dir build -R "OpenOceanKraken_(interface|snapshot)_tests" --output-on-failure
```

预期：Biological 两层和空层往返通过；所有错误 JSON 返回 `false` 且旧配置未变；None/Thorpe/FrancGarr JSON 快照不增加新键。

- [ ] **步骤 7：提交任务 4**

```powershell
git add src/module/json_in_out.hpp src/module/OpenOceanKrakenInterface.cpp for_test/test.cpp
git commit -m "feat: serialize biological attenuation in JSON"
```

---

### 任务 5：Python native 绑定、高层类型和便利入口

**文件：**

- 修改：`src/bind/OpenOceanKrakenBind.cpp:84-131`
- 修改：`wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi:18-65`
- 修改：`wrappers/py_kraken/py_kraken/ook_data_model.py:1-20,140`
- 修改：`wrappers/py_kraken/py_kraken/ook_interface.py:1-15,171-179`
- 修改：`wrappers/py_kraken/py_kraken/__init__.py:1-30`
- 测试：`wrappers/py_kraken/tests/test_native.py`
- 测试：`wrappers/py_kraken/tests/test_interface.py`

**接口：**

- 产出：`native.OceanAbsorptionModel.Biological`
- 产出：`native.BiologicalAttenuationLayer`
- 产出：`native.Atten_Mode.biologicalLayers`
- 产出：`BiologicalAttenuationLayerModel`
- 产出：`OpenOceanKraken_interface.ook_set_biological_attenuation()`

- [ ] **步骤 1：先写 native 和高层失败测试**

在 `test_native.py` 增加：

```python
def test_biological_attenuation_binding_and_value_error(self):
    layer = native.BiologicalAttenuationLayer()
    layer.Z1 = 10.0
    layer.Z2 = 30.0
    layer.f0 = 1000.0
    layer.Q = 5.0
    layer.a0 = 0.04

    mode = native.Atten_Mode()
    mode.absModel = native.OceanAbsorptionModel.Biological
    mode.biologicalLayers = [layer]

    interface = native.Interface()
    interface.set_AttenUnit(mode)
    payload = json.loads(interface.to_json_string())
    self.assertEqual(
        payload["AttenUnit"]["OceanAbsorptionModel"], "Biological")
    self.assertEqual(
        payload["AttenUnit"]["BiologicalLayers"][0]["f0"], 1000.0)

    layer.Q = 0.0
    mode.biologicalLayers = [layer]
    with self.assertRaises(ValueError):
        interface.set_AttenUnit(mode)
```

在 `test_interface.py` 增加：

```python
def test_high_level_biological_attenuation_preserves_order(self):
    interface = OpenOceanKraken_interface(thread_num=1)
    interface.ook_set_biological_attenuation([
        {"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04},
        {"Z1": 40.0, "Z2": 60.0, "f0": 1200.0, "Q": 4.0, "a0": 0.02},
    ])
    payload = json.loads(interface._interface.to_json_string())
    self.assertEqual(
        [item["f0"] for item in
         payload["AttenUnit"]["BiologicalLayers"]],
        [1000.0, 1200.0])

    interface.ook_set_biological_attenuation([])
    payload = json.loads(interface._interface.to_json_string())
    self.assertEqual(payload["AttenUnit"]["BiologicalLayers"], [])

    with self.assertRaises(ValueError):
        interface.ook_set_biological_attenuation([
            {"Z1": 30.0, "Z2": 10.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04}
        ])
```

为两个文件增加 `import json`。

- [ ] **步骤 2：在新的 ABI 构建目录确认测试先失败**

运行：

```powershell
$bioPython = (Get-Command python).Source
cmake -S . -B build_bio -G Ninja `
  -DCMAKE_CXX_COMPILER=D:/program/msys64/ucrt64/bin/g++.exe `
  -DBUILD_DEBUG=OFF `
  -DBUILD_TESTING=ON `
  -DBUILD_AS_EXE=ON `
  -DBUILD_AS_STATIC=ON `
  -DBUILD_AS_PYTHON=ON `
  "-DPython3_EXECUTABLE=$bioPython"
cmake --build build_bio --target OpenOceanKraken_python
ctest --test-dir build_bio -R OpenOceanKraken_python_tests --output-on-failure
```

预期：Python 因 Biological 符号或高层方法不存在而失败。

- [ ] **步骤 3：扩展 pybind11 和 `.pyi`**

在枚举末尾追加：

```cpp
.value("Biological", OceanAbsorptionModel::Biological)
```

在 `Atten_Mode` 绑定前增加：

```cpp
py::class_<BiologicalAttenuationLayer>(
    module, "BiologicalAttenuationLayer")
    .def(py::init<>())
    .def_readwrite("Z1", &BiologicalAttenuationLayer::Z1)
    .def_readwrite("Z2", &BiologicalAttenuationLayer::Z2)
    .def_readwrite("f0", &BiologicalAttenuationLayer::f0)
    .def_readwrite("Q", &BiologicalAttenuationLayer::Q)
    .def_readwrite("a0", &BiologicalAttenuationLayer::a0);
```

`Atten_Mode` 增加：

```cpp
.def_readwrite(
    "biologicalLayers",
    &Atten_Mode::biologicalLayers);
```

`.pyi` 增加对应枚举、类和：

```python
class BiologicalAttenuationLayer:
    Z1: float
    Z2: float
    f0: float
    Q: float
    a0: float
    def __init__(self) -> None: ...

class Atten_Mode:
    attnUnit: AttenuationUnit
    absModel: OceanAbsorptionModel
    biologicalLayers: list[BiologicalAttenuationLayer]
    def __init__(self) -> None: ...
```

- [ ] **步骤 4：增加 Pydantic 层模型**

`ook_data_model.py` 增加 `import math`，并定义：

```python
class BiologicalAttenuationLayerModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Z1: float
    Z2: float
    f0: float
    Q: float
    a0: float

    @model_validator(mode="after")
    def validate_layer(self) -> "BiologicalAttenuationLayerModel":
        for name in ("Z1", "Z2", "f0", "Q", "a0"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        if self.Z1 > self.Z2:
            raise ValueError("Z1 must be less than or equal to Z2")
        if self.f0 <= 0.0:
            raise ValueError("f0 must be positive")
        if self.Q <= 0.0:
            raise ValueError("Q must be positive")
        if self.a0 < 0.0:
            raise ValueError("a0 must be non-negative")
        return self
```

把类名加入 `ook_data_model.__all__` 和包级 `__init__.__all__`。

- [ ] **步骤 5：实现高层便利 setter**

在 `ook_interface.py` 导入该模型并加入：

```python
def ook_set_biological_attenuation(
    self,
    layers: list[
        BiologicalAttenuationLayerModel | dict[str, float]
    ],
    attenuation_unit: Any = None,
) -> None:
    if not isinstance(layers, list):
        raise ValueError("layers must be a list")
    if len(layers) > 200:
        raise ValueError("layers must contain at most 200 entries")

    parsed = [
        BiologicalAttenuationLayerModel.model_validate(layer)
        for layer in layers
    ]
    native_layers = []
    for item in parsed:
        layer = native.BiologicalAttenuationLayer()
        layer.Z1 = item.Z1
        layer.Z2 = item.Z2
        layer.f0 = item.f0
        layer.Q = item.Q
        layer.a0 = item.a0
        native_layers.append(layer)

    mode = native.Atten_Mode()
    mode.attnUnit = (
        native.AttenuationUnit.MODE_W_db_per_lambda
        if attenuation_unit is None
        else attenuation_unit
    )
    mode.absModel = native.OceanAbsorptionModel.Biological
    mode.biologicalLayers = native_layers
    self._interface.set_AttenUnit(mode)
```

空列表合法；最终 C++ setter 仍是权威校验边界。

- [ ] **步骤 6：重建 `.pyd` 并运行 Python 回归**

运行：

```powershell
cmake --build build_bio --target OpenOceanKraken_python
ctest --test-dir build_bio -R OpenOceanKraken_python_tests --output-on-failure
```

预期：native 创建、vector 赋值、高层顺序、JSON 和 `ValueError` 全部通过。

- [ ] **步骤 7：提交任务 5**

```powershell
git add src/bind/OpenOceanKrakenBind.cpp wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi wrappers/py_kraken/py_kraken/ook_data_model.py wrappers/py_kraken/py_kraken/ook_interface.py wrappers/py_kraken/py_kraken/__init__.py wrappers/py_kraken/tests/test_native.py wrappers/py_kraken/tests/test_interface.py
git commit -m "feat: expose biological attenuation to Python"
```

---

### 任务 6：MATLAB 标准 JSON 配置入口

**文件：**

- 修改：`wrappers/m_kraken/krakenDataModel.m:38-198`
- 测试：`wrappers/m_kraken/tests/test_m_kraken.m:28-54`

**接口：**

- 产出：`krakenDataModel.setBiologicalAttenuation(layers, attenuationUnit)`
- 输入：`layers` 为 `N×5` 数值矩阵，列顺序固定为 `[Z1 Z2 f0 Q a0]`。
- 输出：标准 `AttenUnit` JSON，保持行顺序。

- [ ] **步骤 1：先写 MATLAB JSON、CLI 和非法输入测试**

在 `test_m_kraken.m` 增加：

```matlab
function testBiologicalAttenuationJsonAndCli(testCase)
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
exe = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'MunkK.env');
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>

baseRoot = fullfile(work, 'base');
command = sprintf('"%s" "%s" --output "%s" --threads 1 --mod-only --json', ...
    exe, fixture, baseRoot);
[status, output] = system(command);
verifyEqual(testCase, status, 0, output);

model = krakenDataModel(exe);
model.loadJson([baseRoot '.json']);
model.setBiologicalAttenuation([
    10.0 30.0 1000.0 5.0 0.04
    40.0 60.0 1200.0 4.0 0.02
]);
jsonPath = model.Write(fullfile(work, 'biological.json'));
payload = jsondecode(fileread(jsonPath));
verifyEqual(testCase, payload.AttenUnit.OceanAbsorptionModel, 'Biological');
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(1).f0, 1000.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(2).f0, 1200.0);

resultRoot = fullfile(work, 'validated');
command = sprintf('"%s" "%s" --output "%s" --threads 1 --mod-only', ...
    exe, jsonPath, resultRoot);
[status, output] = system(command);
verifyEqual(testCase, status, 0, output);
verifyTrue(testCase, isfile([resultRoot '.mod']));

model.setBiologicalAttenuation(zeros(0, 5));
emptyJson = model.Write(fullfile(work, 'biological_empty.json'));
emptyPayload = jsondecode(fileread(emptyJson));
verifyTrue(testCase, isempty(emptyPayload.AttenUnit.BiologicalLayers));

verifyError(testCase, ...
    @() model.setBiologicalAttenuation([30 10 1000 5 0.04]), ...
    'OpenOceanKraken:InvalidBiologicalLayer');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation([10 30 1000 0 0.04]), ...
    'OpenOceanKraken:InvalidBiologicalLayer');
end
```

- [ ] **步骤 2：运行 MATLAB 测试并确认方法不存在**

运行：

```powershell
matlab -batch "restoredefaultpath; addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
```

预期：失败信息指出 `setBiologicalAttenuation` 未定义。

- [ ] **步骤 3：实现 MATLAB 配置方法**

在公开 `methods` 块加入：

```matlab
function setBiologicalAttenuation(obj, layers, attenuationUnit)
    if nargin < 3 || isempty(attenuationUnit)
        attenuationUnit = 'dB/lambda';
    end
    attenuationUnit = validatestring(char(attenuationUnit), {
        'dB/m/kHz', 'Params_Lose', 'dB/m', ...
        'Nepers/m', 'Quality_Factor', 'dB/lambda'});

    if isempty(layers)
        layers = zeros(0, 5);
    end
    if ~isnumeric(layers) || ~isreal(layers) || ...
            ndims(layers) ~= 2 || size(layers, 2) ~= 5 || ...
            size(layers, 1) > 200 || any(~isfinite(layers), 'all')
        error('OpenOceanKraken:InvalidBiologicalLayer', ...
            'layers must be a finite real N-by-5 matrix with N <= 200.');
    end
    if any(layers(:, 1) > layers(:, 2)) || ...
            any(layers(:, 3) <= 0) || ...
            any(layers(:, 4) <= 0) || ...
            any(layers(:, 5) < 0)
        error('OpenOceanKraken:InvalidBiologicalLayer', ...
            'Require Z1 <= Z2, f0 > 0, Q > 0, and a0 >= 0.');
    end

    encodedLayers = cell(1, size(layers, 1));
    for index = 1:size(layers, 1)
        encodedLayers{index} = struct( ...
            'Z1', layers(index, 1), ...
            'Z2', layers(index, 2), ...
            'f0', layers(index, 3), ...
            'Q', layers(index, 4), ...
            'a0', layers(index, 5));
    end
    obj.Attenuation = struct( ...
        'AttenuationUnit', attenuationUnit, ...
        'OceanAbsorptionModel', 'Biological', ...
        'BiologicalLayers', {encodedLayers});
end
```

`applyConfig()` 和 `makeConfig()` 已整体传输 `Attenuation`，不得另建顶层键。

- [ ] **步骤 4：运行 MATLAB 和 CLI 验收**

运行：

```powershell
cmake --build build_bio --target OpenOceanKraken
matlab -batch "restoredefaultpath; addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
```

预期：MATLAB 两层顺序、标准 JSON、OOK CLI 读取和非法输入拒绝全部通过。

- [ ] **步骤 5：提交任务 6，明确避开用户的 demo 文件**

```powershell
git add wrappers/m_kraken/krakenDataModel.m wrappers/m_kraken/tests/test_m_kraken.m
git commit -m "feat: configure biological attenuation from MATLAB"
```

---

### 任务 7：固定算例、Fortran 模态差分和 TL 非零 gate

**文件：**

- 修改：`.gitignore:27,45`
- 修改：`CMakeLists.txt:107-177`
- 新建：`for_test/fixtures/biological/bio_uniform_off.env`
- 新建：`for_test/fixtures/biological/bio_uniform_off.flp`
- 新建：`for_test/fixtures/biological/bio_uniform_resonance.env`
- 新建：`for_test/fixtures/biological/bio_uniform_resonance.flp`
- 新建：`tools/__init__.py`
- 新建：`tools/test_bio_kraken_equivalence.py`
- 新建：`tools/tests/__init__.py`
- 新建：`tools/tests/test_bio_kraken_equivalence.py`

**接口：**

- 产出：`python tools/test_bio_kraken_equivalence.py --ook-exe E:/my_project/01_ook_project/OOK/bin/OpenOceanKraken.exe --kraken-exe E:/my_project/01_ook_project/krakenFortran/build_mingw/out/kraken.exe --field-exe E:/my_project/01_ook_project/krakenFortran/build_mingw/out/field.exe --fixture-root E:/my_project/01_ook_project/OOK/for_test/fixtures/biological`
- 产出：超任一阈值返回退出码 `1`，全部通过返回 `0`。
- 产出：CTest 选项 `OOK_ENABLE_KRAKEN_REFERENCE_TESTS`、缓存路径 `OOK_KRAKEN_EXE`、`OOK_FIELD_EXE`。

- [ ] **步骤 1：先写 gate 纯函数的失败单测**

`tools/tests/test_bio_kraken_equivalence.py` 测试以下公开函数：

```python
import unittest
import numpy as np

from tools.test_bio_kraken_equivalence import (
    modal_errors,
    slope_errors,
    tl_errors,
)


class BiologicalGateTests(unittest.TestCase):
    def test_modal_threshold_breach_is_reported(self):
        errors = modal_errors(
            np.array([1.0 - 1.0e-4j]),
            np.array([1.0 - 1.0e-4j]),
        )
        self.assertEqual(errors, [])
        errors = modal_errors(
            np.array([1.0 + 2.0e-7 - 1.0e-4j]),
            np.array([1.0 - 1.0e-4j]),
        )
        self.assertTrue(any("real wavenumber" in item for item in errors))

    def test_tl_and_slope_threshold_breaches_are_reported(self):
        self.assertEqual(
            tl_errors(np.array([10.0, 10.1]), np.array([10.0, 10.0])),
            [],
        )
        self.assertTrue(tl_errors(
            np.array([10.0, 10.6]), np.array([10.0, 10.0])))
        self.assertEqual(
            slope_errors(np.array([1.0, 5.0, 10.0]),
                         np.array([1.0, 5.0, 10.0])),
            [],
        )
        self.assertTrue(slope_errors(
            np.array([1.0, 5.0, 10.0]),
            np.array([0.0, 0.0, 0.0])))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **步骤 2：运行工具单测并确认模块不存在**

运行：

```powershell
$env:PYTHONPATH = (Resolve-Path .)
python -m unittest tools.tests.test_bio_kraken_equivalence -v
```

预期：导入 `tools.test_bio_kraken_equivalence` 失败。

- [ ] **步骤 3：把固定算例和工具路径纳入 OOK Git**

在 `.gitignore` 末尾追加精确例外：

```gitignore
!/for_test/fixtures/
/for_test/fixtures/*
!/for_test/fixtures/biological/
!/for_test/fixtures/biological/**

!/tools/
/tools/*
!/tools/__init__.py
!/tools/test_bio_kraken_equivalence.py
!/tools/tests/
/tools/tests/*
!/tools/tests/__init__.py
!/tools/tests/test_bio_kraken_equivalence.py
```

把 `tools/__init__.py` 写为 `"""OOK test tools."""`，把
`tools/tests/__init__.py` 写为 `"""Tests for OOK tools."""`，确保
`python -m unittest tools.tests...` 在干净检出中可导入。

`bio_uniform_off.env` 使用：

```text
'Biological uniform off'
1000.0
1
'CVW'
1000 0.0 100.0
0.0 1500.0 /
100.0 1500.0 /
'A' 0.0
100.0 1600.0 0.0 1.5 0.0 /
1400.0 20000.0
10.0
1
50.0 /
1
50.0 /
```

`bio_uniform_resonance.env` 仅在 Top Option 后增加 Biological 块：

```text
'Biological uniform resonance'
1000.0
1
'CVWB'
1
0.0 100.0 1000.0 5.0 0.04
1000 0.0 100.0
0.0 1500.0 /
100.0 1500.0 /
'A' 0.0
100.0 1600.0 0.0 1.5 0.0 /
1400.0 20000.0
10.0
1
50.0 /
1
50.0 /
```

两个 `.flp` 内容完全相同，但各自保留同名文件：

```text
/,
'RA'
1,
1
0.0 /
3
1.0 5.0 10.0 /
1
50.0 /
1
50.0 /
1
0.0 /
```

- [ ] **步骤 4：实现严格阈值函数**

`tools/test_bio_kraken_equivalence.py` 定义：

```python
MAX_REAL_K = 1.0e-7
MIN_IMAG_K = 1.0e-9
REL_IMAG_K = 1.0e-3
MAX_SLOPE_ERROR_DB_PER_KM = 0.02
MAX_P95_TL_DB = 0.20
MAX_SINGLE_MODE_TL_DB = 0.50
MAX_COMMON_FIELD_P95_DB = 0.05
EXPECTED_RANGES_KM = np.array([1.0, 5.0, 10.0])


def modal_errors(ook_k: np.ndarray, fortran_k: np.ndarray) -> list[str]:
    errors: list[str] = []
    if ook_k.size != fortran_k.size:
        return [f"mode count mismatch: OOK={ook_k.size}, "
                f"Fortran={fortran_k.size}"]
    real_delta = np.abs(ook_k.real - fortran_k.real)
    if np.any(real_delta > MAX_REAL_K):
        errors.append(
            f"real wavenumber error {real_delta.max()} exceeds "
            f"{MAX_REAL_K}")
    imag_delta = np.abs(ook_k.imag - fortran_k.imag)
    imag_limit = np.maximum(
        MIN_IMAG_K, REL_IMAG_K * np.abs(fortran_k.imag))
    if np.any(imag_delta > imag_limit):
        errors.append("imaginary wavenumber error exceeds per-mode limit")
    return errors


def slope_errors(ranges_km: np.ndarray, delta_tl: np.ndarray) -> list[str]:
    slope = float(np.polyfit(ranges_km, delta_tl, 1)[0])
    if abs(slope - 1.0) > MAX_SLOPE_ERROR_DB_PER_KM:
        return [f"attenuation slope {slope} dB/km differs from "
                "1.0 dB/km by more than 0.02 dB/km"]
    return []


def tl_errors(ook_tl: np.ndarray, fortran_tl: np.ndarray) -> list[str]:
    delta = np.abs(ook_tl - fortran_tl)
    errors: list[str] = []
    if float(np.percentile(delta, 95)) > MAX_P95_TL_DB:
        errors.append("P95 |delta TL| exceeds 0.20 dB")
    if float(delta.max()) > MAX_SINGLE_MODE_TL_DB:
        errors.append("single-mode max |delta TL| exceeds 0.50 dB")
    return errors
```

- [ ] **步骤 5：实现可执行编排和所有验收条件**

主程序必须：

1. 用临时目录分别复制 off/on `.env/.flp`；
2. 对每个 case 运行 OOK：
   `OpenOceanKraken.exe <env> --output <root> --threads 1 --mod`；
3. 对每个 case 在独立目录运行：
   `kraken.exe <case-base>`，随后 `field.exe <case-base>`；
4. 再复制 OOK `.mod` 给同一个 Fortran `field.exe` 生成 common-field SHD；
5. 通过 `wrappers/py_kraken/py_kraken/ook_read.py` 的 `read_mod/read_shd` 读取结果；
6. 检查所有进程返回码，日志不得含 `fatal error` 或 `no modes`；
7. on case 的 Fortran `.prt` 必须含 `Biological attenuation`；
8. 检查模态数、Re/Im 阈值、on 时全部 `Im(k)<0`、on/off 虚部存在大于 `1e-9` 的差异；
9. 检查所有压力有限且非零；
10. 对 OOK 与 Fortran 分别计算 `TL_on-TL_off` 在 1/5/10 km 的斜率；
11. 检查 OOK/Fortran 自有 field 的 P95/max；
12. 检查同一 Fortran field 读取两份 MOD 的 P95；
13. 收集全部错误并返回 `1`，无错误返回 `0`。

主入口必须采用：

```python
def main() -> int:
    args = parse_args()
    errors = run_gate(args)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("Biological attenuation equivalence gate passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

所有 subprocess 使用 `check=False`、捕获 stdout/stderr、固定单线程环境，并把 executable、fixture 和输出路径作为显式参数，不依赖 PATH 猜测。

- [ ] **步骤 6：在 CTest 注册纯工具测试和可选参考实现 gate**

在 `CMakeLists.txt` 的 Python3 测试区加入：

```cmake
add_test(
    NAME OpenOceanKraken_biological_tool_tests
    COMMAND ${Python3_EXECUTABLE} -m unittest
        tools.tests.test_bio_kraken_equivalence -v
)
set_tests_properties(
    OpenOceanKraken_biological_tool_tests
    PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}"
)

option(
    OOK_ENABLE_KRAKEN_REFERENCE_TESTS
    "Enable KrakenFortran Biological modal/TL acceptance gate"
    OFF
)
set(OOK_KRAKEN_EXE "" CACHE FILEPATH "Path to KrakenFortran kraken.exe")
set(OOK_FIELD_EXE "" CACHE FILEPATH "Path to KrakenFortran field.exe")

if(OOK_ENABLE_KRAKEN_REFERENCE_TESTS)
    if(NOT TARGET OpenOceanKraken)
        message(FATAL_ERROR
            "Reference tests require BUILD_AS_EXE=ON")
    endif()
    if(NOT EXISTS "${OOK_KRAKEN_EXE}" OR NOT EXISTS "${OOK_FIELD_EXE}")
        message(FATAL_ERROR
            "Reference tests require valid OOK_KRAKEN_EXE and OOK_FIELD_EXE")
    endif()
    add_test(
        NAME OpenOceanKraken_biological_equivalence_tests
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_bio_kraken_equivalence.py
            --ook-exe $<TARGET_FILE:OpenOceanKraken>
            --kraken-exe ${OOK_KRAKEN_EXE}
            --field-exe ${OOK_FIELD_EXE}
            --fixture-root
                ${CMAKE_CURRENT_SOURCE_DIR}/for_test/fixtures/biological
    )
    set_tests_properties(
        OpenOceanKraken_biological_equivalence_tests
        PROPERTIES LABELS "biological;integration"
    )
endif()
```

- [ ] **步骤 7：运行工具单测和真实 Fortran/TL gate**

运行：

```powershell
$env:PYTHONPATH = (Resolve-Path .)
python -m unittest tools.tests.test_bio_kraken_equivalence -v

$bioPython = (Get-Command python).Source
cmake -S . -B build_bio -G Ninja `
  -DBUILD_DEBUG=OFF `
  -DBUILD_TESTING=ON `
  -DBUILD_AS_EXE=ON `
  -DBUILD_AS_PYTHON=ON `
  -DOOK_ENABLE_KRAKEN_REFERENCE_TESTS=ON `
  -DOOK_KRAKEN_EXE=E:/my_project/01_ook_project/krakenFortran/build_mingw/out/kraken.exe `
  -DOOK_FIELD_EXE=E:/my_project/01_ook_project/krakenFortran/build_mingw/out/field.exe `
  "-DPython3_EXECUTABLE=$bioPython"
cmake --build build_bio --parallel
ctest --test-dir build_bio -R "biological_(tool|equivalence)_tests" --output-on-failure
```

预期：工具单测通过；真实 gate 满足规范全部模态和 TL 阈值。手工把测试阈值降到不可能值时，脚本必须返回非零；随后恢复常量并重跑通过。

- [ ] **步骤 8：提交任务 7**

```powershell
git add .gitignore CMakeLists.txt
git add for_test/fixtures/biological/bio_uniform_off.env for_test/fixtures/biological/bio_uniform_off.flp for_test/fixtures/biological/bio_uniform_resonance.env for_test/fixtures/biological/bio_uniform_resonance.flp
git add tools/__init__.py tools/test_bio_kraken_equivalence.py tools/tests/__init__.py tools/tests/test_bio_kraken_equivalence.py
git commit -m "test: gate biological attenuation against KrakenFortran"
```

---

### 任务 8：文档、干净全量重编译和完成门禁

**文件：**

- 修改：`README.md:5-51`
- 修改：`wrappers/py_kraken/README.md:20-63`
- 修改：`wrappers/m_kraken/README.md:23-40`

**接口：**

- 文档化：ENV `B` 块、JSON、C++/Python/MATLAB、闭区间、重叠层、`a0·Q²`。
- 验证：全部 CTest、Python、MATLAB、Fortran 模态、TL 和干净 ABI 重编译。

- [ ] **步骤 1：补主 README 的 ENV/JSON 和 ABI 章节**

主 README 增加以下标准示例：

```text
'CVWB'
2
10.0 30.0 1000.0 5.0 0.04
40.0 60.0 1200.0 4.0 0.02
```

以及：

```json
{
  "AttenUnit": {
    "AttenuationUnit": "dB/lambda",
    "OceanAbsorptionModel": "Biological",
    "BiologicalLayers": [
      {"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04}
    ]
  }
}
```

文字明确：

- `B` 块紧跟 Top Option，并位于顶部 `A` 半空间/第一介质之前；
- `[Z1,Z2]` 是闭区间，重叠层和共享端点逐层累加；
- `a0` 单位是 dB/km，且 `a(f0)=a0·Q²`；
- Biological 不作用于半空间；
- 此版本改变 C++ ABI，必须从同一提交干净重编译全部原生产物。

- [ ] **步骤 2：补 Python README**

同时给 native 和 facade 示例。高层示例固定为：

```python
from py_kraken import OpenOceanKraken_interface

ook = OpenOceanKraken_interface(thread_num=1)
ook.ook_set_biological_attenuation([
    {"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04},
    {"Z1": 40.0, "Z2": 60.0, "f0": 1200.0, "Q": 4.0, "a0": 0.02},
])
```

说明非法输入由 Python 提前映射为 `ValueError`，C++ setter 仍是最终权威边界。

- [ ] **步骤 3：补 MATLAB README**

加入：

```matlab
model.setBiologicalAttenuation([
    10.0 30.0 1000.0 5.0 0.04
    40.0 60.0 1200.0 4.0 0.02
], 'dB/lambda');
jsonPath = model.Write('tmp/biological_case.json');
```

说明列顺序、单位、空 `0×5` 矩阵语义，以及最终由 OOK CLI JSON 导入校验。

- [ ] **步骤 4：在全新目录重编译所有 ABI 相关目标**

运行：

```powershell
$bioPython = (Get-Command python).Source
cmake -S . -B build_bio_release -G Ninja `
  -DCMAKE_CXX_COMPILER=D:/program/msys64/ucrt64/bin/g++.exe `
  -DBUILD_DEBUG=OFF `
  -DBUILD_TESTING=ON `
  -DBUILD_AS_EXE=ON `
  -DBUILD_AS_STATIC=ON `
  -DBUILD_AS_SHARED=ON `
  -DBUILD_AS_PYTHON=ON `
  -DOOK_ENABLE_KRAKEN_REFERENCE_TESTS=ON `
  -DOOK_KRAKEN_EXE=E:/my_project/01_ook_project/krakenFortran/build_mingw/out/kraken.exe `
  -DOOK_FIELD_EXE=E:/my_project/01_ook_project/krakenFortran/build_mingw/out/field.exe `
  "-DPython3_EXECUTABLE=$bioPython"
cmake --build build_bio_release --parallel
```

必须确认该次构建实际生成/重链接：

```text
OpenOceanKraken_core
OpenOceanKraken
OpenOceanKraken_shared
OpenOceanKraken_interface_tests
ook_case_runner
OpenOceanKraken_snapshot_tests
OpenOceanKraken_python
```

- [ ] **步骤 5：运行完整自动门禁**

运行：

```powershell
ctest --test-dir build_bio_release --output-on-failure

$env:PYTHONPATH = (Resolve-Path wrappers/py_kraken)
python -m unittest discover -s wrappers/py_kraken/tests -v

matlab -batch "restoredefaultpath; addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
```

预期同时通过：

```text
公式单元测试
+ env_in_out 测试
+ json_in_out 测试
+ C++/Python/MATLAB 接口测试
+ Fortran 模态差分
+ TL 验收
+ 全部既有 CTest
```

- [ ] **步骤 6：执行最终静态和工作树检查**

运行：

```powershell
git diff --check
git status --short
git diff --name-only
```

预期：

- `git diff --check` 无输出；
- 任务文件均已提交；
- 唯一剩余的既有用户修改是 `wrappers/m_kraken/demo_env.m`；
- 不存在 `.pyd`、`.dll`、`.mod`、`.shd` 或构建目录误入 Git。

- [ ] **步骤 7：提交文档**

```powershell
git add README.md wrappers/py_kraken/README.md wrappers/m_kraken/README.md
git commit -m "docs: document biological attenuation workflow"
```

---

## 完成判定

执行者只有在以下事实同时成立时才能宣布完成：

1. `Atten_Mode` 是唯一 Biological 配置所有者，没有模块级全局层数组。
2. ENV `B` 和 JSON `"BiologicalLayers"` 分别有独立成功/失败测试。
3. ENV `B` 块在顶部半空间或介质行之前被完整消费。
4. JSON Biological 空数组可往返，非 Biological 输出快照不变。
5. C++ setter、ENV、JSON 失败均保持旧活动参数。
6. 多剖面不一致配置被拒绝。
7. 公式的闭区间、顺序、逐层除法和重叠累加与 Fortran 一致。
8. SSP P/S 应用 Biological，半空间不应用。
9. Python `.pyd` 来自新 ABI 的干净重编译，Python/MATLAB 高层入口均通过。
10. 模态数、复波数和 TL 全部满足已批准阈值，gate 超阈值时返回非零。
11. 既有 CTest 全部通过。
12. 用户的 `wrappers/m_kraken/demo_env.m` 修改未被覆盖或提交。
