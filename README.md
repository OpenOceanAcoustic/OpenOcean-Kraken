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

```text
OpenOceanKraken <input.env|input.json> [options]
  --output <root>  输出根路径（不含扩展名）
  --threads <N>    正整数线程数
  --velocity       输出压力、垂直速度、水平速度
  --mod             额外输出 MOD
  --mod-only        仅计算并输出 MOD
  --json            输出规范化 JSON
```

示例：

```powershell
bin\OpenOceanKraken.exe ..\test\MunkK.env --output tmp\munk --threads 4 --mod
bin\OpenOceanKraken.exe tmp\munk.json --output tmp\munk_json --threads 4
```

普通压力结果为 `<root>.shd`；启用振速后为 `<root>_P.shd`、`<root>_V.shd`、`<root>_H.shd`；模态结果为 `<root>.mod`。

## Wrapper

- Python：见 [`wrappers/py_kraken/README.md`](wrappers/py_kraken/README.md)。
- MATLAB：见 [`wrappers/m_kraken/README.md`](wrappers/m_kraken/README.md)。

两套读取器统一把场数据表示为 `[source, range, depth]`。Python 数组为拥有自身内存的 `numpy.complex64`，清理或释放原生接口后，已经取得的结果仍然有效。OOK 是简正模模型，不提供 Bellhop 专属的 ray、arrival、beam 或 angle API。
