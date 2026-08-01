# OpenOceanKraken（OOK）

OpenOceanKraken 是并行水声简正模计算程序。本仓库提供统一的 C++ 核心、生产命令行、Python `py_kraken` wrapper 和 MATLAB `m_kraken` wrapper。Eigen、nlohmann/json 与 pybind11 均位于 `third-party/`，所有前端复用同一个 `OpenOceanKraken_core`。

## 构建

Windows + MinGW/Ninja 示例：

```powershell
cmake -S . -B build_final -G Ninja `
  -DCMAKE_CXX_COMPILER=D:/program/mingw64/bin/g++.exe `
  -DBUILD_TESTING=ON -DBUILD_AS_EXE=ON -DBUILD_AS_PYTHON=ON `
  -DPython3_EXECUTABLE=(Get-Command python).Source
cmake --build build_final --parallel
ctest --test-dir build_final --output-on-failure
```

可用构建选项：

- `BUILD_AS_EXE`：构建 `OpenOceanKraken` CLI，默认开启。
- `BUILD_AS_STATIC`：输出核心静态库，默认开启。
- `BUILD_AS_SHARED`：额外输出共享库，默认关闭。
- `BUILD_AS_PYTHON`：构建原生 Python 模块 `OpenOceanKraken`，默认关闭。

也可使用 Visual Studio 或其他 C++17 编译器；省略上例中的生成器和编译器参数即可使用 CMake 默认工具链。

## CLI

运行以下命令可以查看完整的中文帮助：

```powershell
bin\OpenOceanKraken.exe -h
bin\OpenOceanKraken.exe --help
```

基本调用格式：

```text
OpenOceanKraken.exe <文件名> [-t 线程数] [-v] [-j] [-out 输出目录] [-time]
```

参数说明：

- `<文件名>`：输入配置文件，支持 `.env` 或 `.json`；省略后缀时默认使用 `.env`。
- `-t <线程数>`：指定正整数线程数，等价于 `--threads <N>`；未指定时使用硬件线程数。
- `-v`：启用压力、垂直振速和水平振速计算，等价于 `--velocity`。它是开关参数，后面不需要填写 `0` 或 `1`。
- `-j`：运行前导出规范化 JSON 配置，等价于 `--json`。
- `-out <输出目录>`：将结果写入指定目录，结果根名称沿用输入文件名。
- `-time`、`--time`：打印本次计算耗时。
- `--output <路径>`：直接指定不含扩展名的结果根路径，不能和 `-out` 同时使用。
- `--mod`：在声场结果之外额外导出 MOD 文件。
- `--mod-only`：只计算本征模态并导出 MOD 文件。
- `-h`、`--help`：显示命令行帮助。

常用示例：

```powershell
# 默认调用
bin\OpenOceanKraken.exe ..\test\MunkK.env

# 输入文件省略 .env 后缀
bin\OpenOceanKraken.exe ..\test\MunkK

# 使用 8 个线程
bin\OpenOceanKraken.exe ..\test\MunkK.env -t 8

# 启用压力和振速计算
bin\OpenOceanKraken.exe ..\test\MunkK.env -v

# 运行前导出 JSON
bin\OpenOceanKraken.exe ..\test\MunkK.env -j

# 将结果写入 output 目录
bin\OpenOceanKraken.exe ..\test\MunkK.env -out output

# 打印计算耗时
bin\OpenOceanKraken.exe ..\test\MunkK.env -time

# 组合使用
bin\OpenOceanKraken.exe ..\test\MunkK.env -t 8 -v -j -out output -time

# 使用完整长选项，并额外导出 MOD
bin\OpenOceanKraken.exe ..\test\MunkK.env --output tmp\munk --threads 4 --mod

# 从 JSON 加载参数
bin\OpenOceanKraken.exe tmp\munk.json --output tmp\munk_json --threads 4
```

普通压力结果为 `<root>.shd`；启用振速后为 `<root>_P.shd`、`<root>_V.shd`、`<root>_H.shd`；模态结果为 `<root>.mod`。

OOB 的 `-M` 内存上限和 `-m` 内存报告依赖其专用内存管理模块。OOK 当前没有对应接口，使用这两个参数时会明确提示不支持。

## Biological 体积衰减

ENV 的 Top Option 第四个字符为 `B` 时，其后紧跟 Biological 层数和逐层参数。例如：

```text
'CVWB'
2
10.0 30.0 1000.0 5.0 0.04
40.0 60.0 1200.0 4.0 0.02
```

`B` 块必须紧跟 Top Option，并位于可选的顶部 `A` 半空间记录或第一条介质记录之前。每一层固定为 `Z1 Z2 f0 Q a0`：`Z1`、`Z2` 的单位为 m，`f0` 的单位为 Hz，`Q` 为无量纲品质因数，`a0` 的单位为 dB/km。

也可把以下 `AttenUnit` 配置片段合并到完整 OOK JSON 的标准顶层：

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

Biological 层遵循以下语义：

- 深度范围 `[Z1,Z2]` 是闭区间；共享端点和彼此重叠的层按输入顺序逐层累加，不合并或去重。
- 在共振频率处 `a(f0)=a0·Q²`，因此 `a0` 不是共振峰值本身。
- Biological 附加项作用于 SSP 介质中的 P 波和 S 波衰减计算；顶部、底部半空间仍保留各自的材料衰减，但不应用 Biological 附加项。
- Biological 的空 `BiologicalLayers` 数组可以导入、导出并保持为空；非 Biological 配置的 JSON 输出不会增加 `BiologicalLayers`。

C++ 可通过 `Atten_Mode::biologicalLayers` 配置层，并由 `KernelInterface::set_AttenUnit` 完成最终校验。此版本为公开 C++ 对象增加了布局成员，因而改变了二进制 ABI：核心库、CLI、共享库、Python `.pyd` 以及所有 C++ 二进制客户端必须从同一提交整体干净重编译，不能混用旧对象文件或旧动态库。`.mod`、`.shd` 的存储格式没有改变；变化只发生在内存对象布局和中间计算能力。格式不变不表示数值不变，启用 Biological 后输出的模态和声场数值会反映新增衰减。

## Wrapper

- Python：见 [`wrappers/py_kraken/README.md`](wrappers/py_kraken/README.md)。
- MATLAB：见 [`wrappers/m_kraken/README.md`](wrappers/m_kraken/README.md)。

两套读取器统一把场数据表示为 `[source, range, depth]`。Python 数组为拥有自身内存的 `numpy.complex64`，清理或释放原生接口后，已经取得的结果仍然有效。OOK 是简正模模型，不提供 Bellhop 专属的 ray、arrival、beam 或 angle API。
