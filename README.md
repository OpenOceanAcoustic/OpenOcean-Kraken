# OpenOcean-Krakenc

`OpenOcean-Krakenc` 是原版 `krakenFortran/krakenc` 与其 `field.f90` 声场链路的纯 C++17 实现。工程组织和公共接口参考已验证的 `OpenOcean-Kraken`，数值验收以原版 Fortran `krakenc.exe + field.exe` 为基准。

## 已实现能力

- KRAKENC ENV 解析：单剖面与连续多剖面、声学层、弹性层、真空/刚性/半空间边界、BRC/IRC。
- 复色散函数、确定性复根搜索、根排除、网格加密与 Richardson 外推。
- 复三对角逆迭代、复模态归一化、确定性相位和群速度。
- 弹性 compound-matrix 传播以及复杂边界阻抗。
- Fortran direct-access 兼容 MOD 读写，支持一个文件内的多环境剖面。
- `field.f90` 对应声场：距离无关、绝热模态（AD）和耦合模态（CM）。
- 压力、水平振速、垂直振速；相干/非相干叠加；`R`、`X`、`S` 源类型；接收距离偏移和 SBP 波束图。
- Fortran 兼容 SHD 输出、命令行入口和可嵌入 C++ `KernelInterface`。
- 单线程与多线程 MOD/Field；多线程 MOD 输出经过 SHA-256 字节一致性验证。

当前按一次调用一个频率工作。多频批处理可由调用方逐频执行；一个 MOD 内的原生多频记录尚未实现。

## 构建与测试

```powershell
cmake -S . -B build-release -G Ninja `
  -DCMAKE_CXX_COMPILER=D:/program/mingw64/bin/g++.exe `
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON `
  -DOPENOCEANKRAKENC_EIGEN_INCLUDE=D:/deps/eigen `
  -DOPENOCEANKRAKENC_JSON_INCLUDE=D:/deps/nlohmann-json/include
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

仓库内置的契约、I/O 和小型数值夹具默认直接运行。完整参考数据套件通过
`OPENOCEANKRAKENC_TEST_ROOT` 显式指定；目录缺少任一必需案例时，依赖外部
数据的测试会标记为 disabled，而不会误报为数值回归。Fortran 对比另需设置
`OPENOCEANKRAKENC_FORTRAN_ROOT`。

Release 默认在编译器支持时启用 IPO/LTO。可用 `-DOPENOCEANKRAKENC_IPO=OFF` 关闭；`OPENOCEANKRAKENC_NATIVE_OPTIMIZATION` 默认关闭。

严格告警构建：

```powershell
cmake -S . -B build-warnings -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_CXX_COMPILER=D:/program/mingw64/bin/g++.exe `
  -DOPENOCEANKRAKENC_IPO=OFF `
  "-DCMAKE_CXX_FLAGS=-Wall -Wextra -Wpedantic -Werror"
cmake --build build-warnings --parallel
```

## 命令行

```powershell
# 复波数 JSON
build-release\OpenOcean-Krakenc.exe --eigen ..\test\MunkKleaky.env

# ENV -> MOD；多剖面时线程可并行求解
build-release\OpenOcean-Krakenc.exe --mod `
  ..\test\stepK_rd.env ..\test\stepK_rd.mod --threads 4

# 读取 <root>.mod/.flp，并在存在时读取同名 ENV/SBP/BRC/IRC，输出 SHD
build-release\OpenOcean-Krakenc.exe --field ..\test\stepK_rd
```

`--field <root>` 保留原有 MOD + FLP 回读方式；`--field <case.env|case.json>` 直接在内存中执行本征值与 Field 链路，不会隐式生成临时 MOD。公共接口的 `runField()` 在本征结果缺失时先在内存中求解，`run()` 则按 `RunMode` 调度完整链路。

## C++ 接口

```cpp
#include "OpenOceanKrakencKernelInterface.h"

OpenOceanKrakenc::KernelInterface solver;
solver.setNumThreads(4);
auto &params = solver.getParams();
params.envPath = R"(..\test\stepK_rd.env)";
solver.run();

const auto &eigen = solver.getOutput().eigen;
const auto &pressure = solver.getOutput().pressure;
const auto &horizontal = solver.getOutput().horizontalVelocity;
const auto &vertical = solver.getOutput().verticalVelocity;
```

### ENV、JSON 与 setter 统一入口

```cpp
OpenOceanKrakenc::KernelInterface solver;

solver.from_env("case.env");       // ENV + 同名 FLP；按需读取 SBP/BRC/IRC/TRC
solver.from_json("case.json");     // JSON schema v1

solver.set_Freq(50.0);              // 也可完全使用 setter 配置
solver.set_Sz(25.0, 25.0, 1);
solver.run();
solver.export_result("result");
```

三种入口最终都转换为同一个 `OOKC_parameters`，再进入同一套纯 C++ Krakenc + Field 数值核心。解析失败采用事务语义，不覆盖上一次有效参数。`get_u/get_v/get_h` 分别对应压力、垂直振速和水平振速，平铺顺序为 source → receiver-depth → range。

`KernelInterface::from_env()` 与 OpenOcean-Kraken 模板一致，面向完整 ENV + FLP 运行配置；同名 FLP 缺失时返回 `false` 且不修改原参数。仅求本征值时使用 CLI `--eigen case.env` 或底层 `read_env_file()`。纯 MOD + FLP 的 `runField()` 不要求 SSP，但仍会执行多频、吸收模型和网格能力检查。公共 `RProf`/剖面 `Range` 使用 km，接收距离 `Rr` 使用 m。

CLI 的本征、MOD 和 Field 均支持 ENV/JSON 双输入：

```powershell
build-release\OpenOcean-Krakenc.exe --convert-env case.env case.json
build-release\OpenOcean-Krakenc.exe --eigen case.env
build-release\OpenOcean-Krakenc.exe --eigen case.json
build-release\OpenOcean-Krakenc.exe --mod case.env case_env.mod --threads 4
build-release\OpenOcean-Krakenc.exe --mod case.json case.mod --threads 4
build-release\OpenOcean-Krakenc.exe --field case.env
build-release\OpenOcean-Krakenc.exe --field case.json
```

`export_result("root")` 生成 `root_P.mod` 和 `root_P.shd`；启用 `set_Velocity_enable(true)` 后同时生成 `root_V.shd` 与 `root_H.shd`。`export_mod()` 直接使用已经计算的 `ModeZ/PhiMode/k`，不会重新求解。

### 公共能力矩阵

| 能力 | JSON 往返 | 求解 |
|---|---:|---:|
| SSP N / C | 支持 | 支持 |
| SSP P / S / A | 支持 | 明确拒绝 |
| Boundary V / R / A / F / P | 支持 | 支持 |
| Boundary G | 支持 | 明确拒绝 |
| 顶部 TRC、底部 BRC / IRC | 支持 | 支持 |
| Attenuation F / L / M / m / N / Q / W | 支持 | 支持 |
| Thorpe / FrancGarr | 支持 | 明确拒绝 |
| Source R / X / S | 支持 | 支持 |
| Range-independent / AD / CM | 支持 | 支持 |
| Coherent / Incoherent | 支持 | 支持 |
| Rectangular grid | 支持 | 支持 |
| Irregular grid | 支持 | Field 前明确拒绝 |
| 单频 | 支持 | 支持 |
| 多频 `freqvec` | 支持 | 单次求解前明确拒绝 |

未实现能力会在进入数值核心或写文件前抛出 `std::invalid_argument`，不会静默降级。

`flpPath`、`modPath` 和 `shdPath` 留空时由 `envPath` 的文件根推导。

## 验证摘要

- Release 与 Debug 均为 20/20 CTest 通过，包含 8 个代表案例的 ENV/JSON 等价矩阵。
- `MunkKleaky`：329 个复波数全量对照；最大绝对误差约 `7.43e-9 1/m`，最大相对误差约 `5.16e-8`。
- 复模态形状、群速度、三个弹性案例、BRC/IRC 均有 Fortran 回归。
- 距离无关 Field 的复压力相对 L2 阈值为 `2e-5`。
- stepK/wedge 的 AD、CM 以及纯 C++ `ENV -> MOD -> Field -> SHD` 链路已通过；端到端压力相对 L2 阈值为 `5e-3`。
- 四线程重复生成的代表性 MOD/SHD 与自身重复结果字节一致。
- 严格 `-Wall -Wextra -Wpedantic -Werror` 构建通过。
- 本轮复测的单线程八案例累计耗时比 C++/Fortran 为 `1.196x`；四线程累计比为 `0.486x`，四线程性能案例中位比为 `0.685x`。
- 最终只读代码复审未发现 Critical 或 Important。

详细记录见 [实现状态](../docs/01_Initial_writing_of_Krakenc/OpenOcean-Krakenc_Final_Implementation_Status.md)、[阶段五验证与性能报告](../docs/01_Initial_writing_of_Krakenc/OpenOcean-Krakenc_Phase5_Verification_and_Performance_Report.md)、[API 兼容实施与复审报告](../docs/02/OpenOcean-Krakenc_API_Compatibility_Verification_Report.md) 和 [HTML 复审看板](../docs/02/OpenOcean-Krakenc_API_Compatibility_Review.html)。

## 已知参考差异与限制

- `multilayer_elastic_stack` 的 Fortran 粗/细网格尾部配对不完整；C++ 保留 26 个有限根，自动测试严格比较前 19 个稳定配对模态，不用截断或零值伪装成功。
- 原版 `neggradC_brc` 的常数 `R=1` 在阻抗换算中出现奇异除法；C++ 安全等价为刚性边界，并用等价 IRC/Fortran 结果验证。
- 原版 Fortran Field 在 wedge 的 AD 分支会异常退出，因此 wedge 使用 CM 完成 Fortran 对照；stepK 同时覆盖 AD 与 CM。
- 尚未实现一个 MOD 文件内的多频记录、Field 3D 和非模态声场分支。
- 参考 Fortran 源码声明 GPLv3。正式分发前仍需完成许可证兼容性、署名和源码提供义务检查。

本轮未执行 Git 提交、合并、推送或 PR。
