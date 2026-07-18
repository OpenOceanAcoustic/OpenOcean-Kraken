# OOK Wrapper Parity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 OOK 增加可复现构建的 Python 原生绑定与高层包、文件驱动的 MATLAB wrapper、正式 CLI、SHD/MOD 读取绘图能力，并把第三方依赖统一迁移到 `third-party/`。

**Architecture:** 保持 `OpenOceanKraken_core` 的数值算法不变，在 C++ 接口上增加拥有内存的结果快照，再通过 pybind11 转换为 NumPy。Python 直接调用绑定；MATLAB 写 JSON、调用正式 CLI、读取 MOD/SHD。所有新行为先写失败测试，再做最小实现。

**Tech Stack:** C++17、CMake 3.10+、Eigen、nlohmann/json、pybind11、CPython 3.9+、NumPy、Pydantic 2、Matplotlib、MATLAB `-batch`。

## Global Constraints

- 不修改 Kraken 核心计算公式或现有数值容差。
- 不迁移 Bellhop 的 ray、arrival、beam、angle 接口。
- Python 包名为 `py_kraken`，原生模块名为 `OpenOceanKraken`，门面类名为 `OpenOceanKraken_interface`，高层方法使用 `ook_` 前缀。
- `FieldData.values` 的固定维度顺序为 `[source, range, depth]`，dtype 为 `numpy.complex64`。
- Python 默认取得拥有自身内存的结果副本，不把裸指针作为推荐 API。
- MATLAB 只通过 JSON/ENV、CLI、MOD/SHD 工作，不新增 MEX。
- 根目录 `eigen/`、`nlohmann/` 在迁移后删除，只保留 `third-party/` 下的版本。
- Windows 是完整验收平台；CMake 不写死 Python 小版本或 `.pyd` 文件名。

---

## File Structure

### C++ and build files

- `.gitignore`：只忽略根目录历史依赖路径，允许跟踪 `third-party/`。
- `CMakeLists.txt`：构建选项、第三方 target、CLI、pybind 和测试注册。
- `main.cpp`：正式 CLI 入口。
- `include/OpenOceanKrakenInterface.h`：拥有内存的 `FieldSnapshot`、`ModeProfileSnapshot` 与 getter 声明。
- `src/module/OpenOceanKrakenInterface.cpp`：结果深拷贝实现。
- `src/module/cli_options.hpp`、`src/module/cli_options.cpp`：CLI 参数模型和纯解析逻辑。
- `src/bind/OpenOceanKrakenBind.cpp`：pybind11 模块。
- `third-party/eigen/`、`third-party/nlohmann-json/`、`third-party/pybind11/`：本地依赖。

### Python files

- `wrappers/py_kraken/pyproject.toml`：Python 包元数据和依赖。
- `wrappers/py_kraken/py_kraken/__init__.py`：稳定导出面。
- `wrappers/py_kraken/py_kraken/ook_data_model.py`：Pydantic/NumPy 数据模型。
- `wrappers/py_kraken/py_kraken/ook_read.py`：SHD、MOD、JSON reader。
- `wrappers/py_kraken/py_kraken/ook_interface.py`：高层门面。
- `wrappers/py_kraken/py_kraken/ook_plot.py`：声场和模态绘图。
- `wrappers/py_kraken/py_kraken/demo.py`：可运行示例。
- `wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi`：原生 API 类型提示。
- `wrappers/py_kraken/install.py`：在活动虚拟环境中构建、安装、验证。
- `wrappers/py_kraken/tests/`：原生、reader、门面和绘图测试。

### MATLAB files

- `wrappers/m_kraken/krakenDataModel.m`：配置、JSON、CLI 与结果门面。
- `wrappers/m_kraken/read_shd.m`：SHD reader。
- `wrappers/m_kraken/read_mod.m`：MOD reader。
- `wrappers/m_kraken/plotshd.m`：压力/速度场图。
- `wrappers/m_kraken/plotmode.m`：模态图。
- `wrappers/m_kraken/demo_env.m`、`wrappers/m_kraken/demo_json.m`：示例。
- `wrappers/m_kraken/tests/test_m_kraken.m`：MATLAB batch 验收。

---

### Task 1: Migrate dependencies and establish build targets

**Files:**
- Create: `tools/tests/test_dependency_layout.py`
- Create: `third-party/eigen/`
- Create: `third-party/nlohmann-json/nlohmann/`
- Create: `third-party/pybind11/`
- Modify: `.gitignore`
- Modify: `CMakeLists.txt`
- Remove: `eigen/`
- Remove: `nlohmann/`

**Interfaces:**
- Consumes: current header includes such as `<Eigen/Dense>` and `<nlohmann/json.hpp>`.
- Produces: CMake targets `OpenOceanKraken_core`, `OpenOceanKraken`, `OpenOceanKraken_shared`, and build option `BUILD_AS_PYTHON`.

- [ ] **Step 1: Write the failing dependency-layout test**

```python
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


class DependencyLayoutTests(unittest.TestCase):
    def test_dependencies_live_only_under_third_party(self):
        self.assertFalse((ROOT / "eigen").exists())
        self.assertFalse((ROOT / "nlohmann").exists())
        self.assertTrue((ROOT / "third-party/eigen/Eigen/Core").is_file())
        self.assertTrue((ROOT / "third-party/nlohmann-json/nlohmann/json.hpp").is_file())
        self.assertTrue((ROOT / "third-party/pybind11/CMakeLists.txt").is_file())

    def test_cmake_exposes_required_build_options(self):
        text = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        for option in ("BUILD_AS_EXE", "BUILD_AS_STATIC", "BUILD_AS_SHARED", "BUILD_AS_PYTHON"):
            self.assertIn(f"option({option}", text)
        self.assertIn("third-party/eigen", text)
        self.assertIn("third-party/nlohmann-json", text)
        self.assertIn("third-party/pybind11", text)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the test and verify the layout test fails**

Run: `python -m unittest tools.tests.test_dependency_layout -v`

Expected: FAIL because root `eigen/` and `nlohmann/` exist and `third-party/` does not.

- [ ] **Step 3: Move the dependencies and make ignore rules root-specific**

Run these PowerShell operations after resolving each source and destination under the repository root:

```powershell
New-Item -ItemType Directory -Path third-party -Force
Move-Item -LiteralPath eigen -Destination third-party/eigen
New-Item -ItemType Directory -Path third-party/nlohmann-json -Force
Move-Item -LiteralPath nlohmann -Destination third-party/nlohmann-json/nlohmann
Copy-Item -Recurse -LiteralPath ../OpenOcean-Bellhop-main/third-party/pybind11-master -Destination third-party/pybind11
```

Change the dependency ignore entries to root-only entries:

```gitignore
/eigen/
/nlohmann/
/pybind11-master/
```

- [ ] **Step 4: Refactor CMake around explicit options and third-party paths**

Add these options and dependency variables while retaining the existing core source list and tests:

```cmake
option(BUILD_AS_EXE "Build the OpenOceanKraken CLI" ON)
option(BUILD_AS_STATIC "Build the OpenOceanKraken static library" ON)
option(BUILD_AS_SHARED "Build the OpenOceanKraken shared library" OFF)
option(BUILD_AS_PYTHON "Build the OpenOceanKraken Python module" OFF)

project(OpenOceanKraken LANGUAGES CXX)

set(OOK_EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third-party/eigen")
set(OOK_NLOHMANN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third-party/nlohmann-json")
set(OOK_PYBIND11_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third-party/pybind11")

set(OOK_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/src/module
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util
    ${OOK_EIGEN_DIR}
    ${OOK_NLOHMANN_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)
```

Use these target guards, keeping `OpenOceanKraken_core` as the single implementation linked by all front ends:

```cmake
add_library(OpenOceanKraken_core STATIC ${CORE_SOURCES})
target_include_directories(OpenOceanKraken_core PUBLIC ${OOK_INCLUDE_DIRS})

if(BUILD_AS_EXE)
    add_executable(OpenOceanKraken main.cpp)
    target_link_libraries(OpenOceanKraken PRIVATE OpenOceanKraken_core)
endif()

if(BUILD_AS_SHARED)
    add_library(OpenOceanKraken_shared SHARED ${CORE_SOURCES})
    target_include_directories(OpenOceanKraken_shared PUBLIC ${OOK_INCLUDE_DIRS})
endif()

if(BUILD_AS_PYTHON)
    add_subdirectory(${OOK_PYBIND11_DIR} EXCLUDE_FROM_ALL)
endif()
```

- [ ] **Step 5: Run layout and clean build tests**

Run:

```powershell
python -m unittest tools.tests.test_dependency_layout -v
cmake -S . -B build_wrappers -DBUILD_TESTING=ON -DBUILD_AS_EXE=ON -DBUILD_AS_PYTHON=OFF
cmake --build build_wrappers --parallel
```

Expected: dependency tests PASS; CMake configures without old dependency paths; CLI and existing tests compile.

- [ ] **Step 6: Commit the dependency migration**

```powershell
git add .gitignore CMakeLists.txt tools/tests/test_dependency_layout.py third-party
git commit -m "build: consolidate OOK third-party dependencies"
```

---

### Task 2: Replace the hard-coded executable with a tested CLI

**Files:**
- Create: `src/module/cli_options.hpp`
- Create: `src/module/cli_options.cpp`
- Create: `for_test/test_cli.py`
- Modify: `main.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: `OpenOceanKraken::Interface::{from_env,from_json,run,runEigen,export_shd,export_mod,export_result,to_json}`.
- Produces: `OpenOceanKraken::cli::Options parse_options(const std::vector<std::string>&)` and the documented CLI flags.

- [ ] **Step 1: Write black-box CLI tests**

Create `for_test/test_cli.py` with cases that run the passed executable:

```python
from pathlib import Path
import argparse
import subprocess
import tempfile
import unittest


class CliTests(unittest.TestCase):
    executable: Path
    fixture: Path

    def run_cli(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(self.executable), *args],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_help_is_successful(self):
        result = self.run_cli("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--output", result.stdout)
        self.assertIn("--mod-only", result.stdout)

    def test_missing_input_is_usage_error(self):
        result = self.run_cli()
        self.assertEqual(result.returncode, 2)
        self.assertIn("Usage:", result.stderr)

    def test_mod_only_writes_mod(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "case"
            result = self.run_cli(str(self.fixture), "--output", str(root), "--threads", "1", "--mod-only")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(root.with_suffix(".mod").is_file())

    def test_velocity_writes_three_shd_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "case"
            result = self.run_cli(str(self.fixture), "--output", str(root), "--threads", "1", "--velocity")
            self.assertEqual(result.returncode, 0, result.stderr)
            for suffix in ("_P.shd", "_V.shd", "_H.shd"):
                self.assertTrue(Path(f"{root}{suffix}").is_file())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--fixture", required=True, type=Path)
    args, remaining = parser.parse_known_args()
    CliTests.executable = args.executable.resolve()
    CliTests.fixture = args.fixture.resolve()
    unittest.main(argv=[__file__, *remaining])


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the CLI test and verify failure**

Run: `python for_test/test_cli.py --executable bin/OpenOcean-Kraken.exe --fixture ../test/solve3_mode_loss.env -v`

Expected: FAIL because the current `main.cpp` ignores command-line arguments and uses a hard-coded path.

- [ ] **Step 3: Implement the pure option parser**

Define the exact model:

```cpp
namespace OpenOceanKraken::cli
{
    struct Options
    {
        std::string input_path;
        std::string output_root;
        int threads = 0;
        bool velocity = false;
        bool export_mod = false;
        bool mod_only = false;
        bool export_json = false;
        bool help = false;
    };

    Options parse_options(const std::vector<std::string> &args);
    std::string usage();
}
```

`parse_options` accepts exactly `--output`, `--threads`, `--velocity`, `--mod`, `--mod-only`, `--json`, and `--help`; it rejects unknown flags, missing values, non-positive thread counts, and more than one input path with `std::invalid_argument`.

- [ ] **Step 4: Implement the CLI execution flow**

`main.cpp` must:

```cpp
int main(int argc, char **argv)
{
    try
    {
        std::vector<std::string> args(argv + 1, argv + argc);
        const auto options = OpenOceanKraken::cli::parse_options(args);
        if (options.help)
        {
            std::cout << OpenOceanKraken::cli::usage();
            return 0;
        }
        if (options.input_path.empty())
        {
            std::cerr << OpenOceanKraken::cli::usage();
            return 2;
        }
        return run_kraken(options);
    }
    catch (const std::invalid_argument &error)
    {
        std::cerr << "OpenOceanKraken: " << error.what() << '\n';
        return 2;
    }
    catch (const std::exception &error)
    {
        std::cerr << "OpenOceanKraken failed: " << error.what() << '\n';
        return 1;
    }
}
```

`run_kraken` selects `from_env` or `from_json` from the lower-case extension, creates the output directory, uses hardware concurrency when `threads == 0`, and follows the flag semantics in the design specification.

- [ ] **Step 5: Register and run CLI tests**

Add a CTest named `OpenOceanKraken_cli_tests` using `$<TARGET_FILE:OpenOceanKraken>` and `../test/solve3_mode_loss.env`.

Run:

```powershell
cmake --build build_wrappers --parallel
ctest --test-dir build_wrappers -R "OpenOceanKraken_cli_tests" --output-on-failure
```

Expected: all CLI tests PASS.

- [ ] **Step 6: Commit the CLI**

```powershell
git add CMakeLists.txt main.cpp src/module/cli_options.hpp src/module/cli_options.cpp
git add -f for_test/test_cli.py
git commit -m "feat: add reusable OpenOceanKraken CLI"
```

---

### Task 3: Add owned C++ snapshots and missing modal controls

**Files:**
- Modify: `include/OpenOceanKrakenInterface.h`
- Modify: `src/module/OpenOceanKrakenInterface.cpp`
- Create: `for_test/test_result_snapshots.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: `OOK_parameters`, `OOK_output`, `EigenParams`, and existing raw getters.
- Produces: `FieldSnapshot`, `ModeProfileSnapshot`, result-copy getters, and setters for profile ranges, mode limit, coherence type, and mode type.

- [ ] **Step 1: Write a compile-time and runtime snapshot test**

The test loads `solve3_mode_loss.env`, runs with one thread, and asserts:

```cpp
const auto pressure = interface.getPressureCopy();
require(pressure.source_count == static_cast<std::size_t>(interface.getParams_const().Pos.NSz));
require(pressure.range_count == static_cast<std::size_t>(interface.getParams_const().Pos.NRr));
require(pressure.depth_count == static_cast<std::size_t>(interface.getParams_const().Pos.NRz_per_range));
require(pressure.values.size() == pressure.source_count * pressure.range_count * pressure.depth_count);
require(pressure.source_depths.size() == pressure.source_count);
require(pressure.receiver_ranges.size() == pressure.range_count);

const auto modes = interface.getModesCopy();
require(modes.size() == static_cast<std::size_t>(interface.getParams_const().NProf));
require(!modes.front().wavenumbers.empty());
require(modes.front().mode_shapes.rows() == static_cast<Eigen::Index>(modes.front().wavenumbers.size()));

const auto saved = pressure.values.front();
interface.clearResults();
require(pressure.values.front() == saved);

interface.set_RProf(0.0, 10000.0, 3);
interface.set_MLimit(64);
interface.set_CoherenceType(OpenOceanKraken::CoherenceType::Coherent);
interface.set_ModeType(OpenOceanKraken::ModeType::Adiabatic);
require(interface.getParams_const().NProf == 3);
require(interface.getParams_const().MLimit == 64);
```

- [ ] **Step 2: Build and verify the test fails to compile**

Run: `cmake --build build_wrappers --target OpenOceanKraken_snapshot_tests --parallel`

Expected: compilation fails because snapshot types and getters do not exist.

- [ ] **Step 3: Add exact snapshot value types**

Add public value types in namespace `OpenOceanKraken`:

```cpp
struct FieldSnapshot
{
    std::string title;
    double frequency = 0.0;
    std::size_t source_count = 0;
    std::size_t range_count = 0;
    std::size_t depth_count = 0;
    std::vector<double> source_depths;
    std::vector<double> receiver_ranges;
    std::vector<double> receiver_depths;
    std::vector<std::complex<float>> values;
};

struct ModeProfileSnapshot
{
    double profile_range = 0.0;
    std::vector<double> depth;
    std::vector<std::complex<double>> wavenumbers;
    std::vector<double> group_velocity;
    Eigen::MatrixXcd mode_shapes;
};
```

Add these `Interface` methods:

```cpp
FieldSnapshot getPressureCopy() const;
FieldSnapshot getVerticalVelocityCopy() const;
FieldSnapshot getHorizontalVelocityCopy() const;
std::vector<ModeProfileSnapshot> getModesCopy() const;
void set_RProf(const Eigen::VectorXd &ranges);
void set_RProf(double start, double end, int count);
void set_MLimit(int limit);
void set_CoherenceType(CoherenceType type);
void set_ModeType(ModeType type);
```

- [ ] **Step 4: Implement deep-copy getters**

Use one private helper that validates the output pointer, copies coordinates from `params->Pos`, and copies exactly `NSz * NRr * NRz_per_range` complex values in existing source/range/depth storage order. `getModesCopy()` copies only the first `M` wavenumbers, group velocities, and mode-shape rows for each of the `NProf` profiles.

The modal setters validate non-empty, finite, non-decreasing profile ranges; `count > 0`; `limit > 0`; update `NProf`, `RProf`, `MLimit`, `coherenceType`, or `modeType`; and mark eigen and field results dirty through the existing invalidation mechanism.

- [ ] **Step 5: Run snapshot and existing interface tests**

Run:

```powershell
cmake --build build_wrappers --parallel
ctest --test-dir build_wrappers -R "snapshot|interface" --output-on-failure
```

Expected: snapshot test PASS and existing interface tests remain PASS.

- [ ] **Step 6: Commit snapshot APIs**

```powershell
git add include/OpenOceanKrakenInterface.h src/module/OpenOceanKrakenInterface.cpp CMakeLists.txt
git add -f for_test/test_result_snapshots.cpp
git commit -m "feat: expose owned Kraken result snapshots"
```

---

### Task 4: Build the pybind11 native module

**Files:**
- Modify: `src/bind/OpenOceanKrakenBind.cpp`
- Modify: `CMakeLists.txt`
- Create: `wrappers/py_kraken/tests/test_native.py`

**Interfaces:**
- Consumes: snapshot APIs from Task 3 and existing `Interface` setters/exporters.
- Produces: importable native module `OpenOceanKraken` with `ThreadPool`, enums, SSP value types, `Interface`, and NumPy result methods.

- [ ] **Step 1: Write native import and shape tests**

```python
from pathlib import Path
import sys
import unittest
import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[1] / "py_kraken"
sys.path.insert(0, str(PACKAGE_DIR))
import OpenOceanKraken as native


class NativeBindingTests(unittest.TestCase):
    def test_required_symbols_are_exported(self):
        for name in ("Interface", "ThreadPool", "Run_Mode", "SSP_Mode", "Media_Mode", "BC_Mode"):
            self.assertTrue(hasattr(native, name), name)

    def test_pressure_is_owned_complex64_source_range_depth_array(self):
        fixture = Path(__file__).resolve().parents[3] / "test/solve3_mode_loss.env"
        interface = native.Interface()
        pool = native.ThreadPool(1)
        interface.setThreadPool(pool)
        interface.setNumThreads(1)
        self.assertTrue(interface.from_env(str(fixture)))
        interface.run()
        pressure = interface.get_pressure()
        self.assertEqual(pressure.dtype, np.complex64)
        self.assertEqual(pressure.ndim, 3)
        saved = pressure.copy()
        interface.clearResults()
        np.testing.assert_array_equal(pressure, saved)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Configure Python build and verify the test fails**

Run:

```powershell
cmake -S . -B build_python -DBUILD_TESTING=ON -DBUILD_AS_PYTHON=ON -DPython3_EXECUTABLE=(Get-Command python).Source
cmake --build build_python --parallel
python -m unittest discover -s wrappers/py_kraken/tests -p test_native.py -v
```

Expected: import fails because `OpenOceanKrakenBind.cpp` is empty and the target is not implemented.

- [ ] **Step 3: Bind enums and value types**

Use `PYBIND11_MODULE(OpenOceanKraken, module)` and bind `ThreadPool`; `SSP_Mode`, `Media_Mode`, `AttenuationUnit`, `OceanAbsorptionModel`, `BC_Mode`, `Source_Mode`, `Grid_Mode`, `Run_Mode`, `CoherenceType`, `ModeType`; and `Atten_Mode`, `HSInfo`, `SSPLayer`, `Range_Independent_Area`, `ReflectionCoef`, `FieldSnapshot`, `ModeProfileSnapshot` with their writable or read-only value fields.

- [ ] **Step 4: Bind Interface with safe lifetimes and NumPy results**

Bind both constructors, using `py::keep_alive<1, 2>()` for the thread-pool constructor and `setThreadPool`. Bind the existing loading, serialization, setter, run, clear, free and export methods. Expose snapshot results with these names:

```cpp
.def("get_pressure", [](const Interface &self) { return field_to_numpy(self.getPressureCopy()); })
.def("get_vertical_velocity", [](const Interface &self) { return field_to_numpy(self.getVerticalVelocityCopy()); })
.def("get_horizontal_velocity", [](const Interface &self) { return field_to_numpy(self.getHorizontalVelocityCopy()); })
.def("get_pressure_snapshot", &Interface::getPressureCopy)
.def("get_vertical_velocity_snapshot", &Interface::getVerticalVelocityCopy)
.def("get_horizontal_velocity_snapshot", &Interface::getHorizontalVelocityCopy)
.def("get_modes", &Interface::getModesCopy)
```

`field_to_numpy` allocates a new `py::array_t<std::complex<float>>` with shape `{source_count, range_count, depth_count}` and copies all values. Do not expose `get_u_AllSources()` or other pointer addresses as normal Python methods.

- [ ] **Step 5: Configure output path and run native tests**

Set the module `LIBRARY_OUTPUT_DIRECTORY` and `RUNTIME_OUTPUT_DIRECTORY` to `wrappers/py_kraken/py_kraken`, including per-configuration variants for multi-config generators.

Run:

```powershell
cmake --build build_python --parallel
python -m unittest discover -s wrappers/py_kraken/tests -p test_native.py -v
```

Expected: native tests PASS.

- [ ] **Step 6: Commit the native module**

```powershell
git add CMakeLists.txt src/bind/OpenOceanKrakenBind.cpp wrappers/py_kraken/tests/test_native.py
git commit -m "feat: add OpenOceanKraken Python bindings"
```

---

### Task 5: Add Python data models and SHD/MOD readers

**Files:**
- Create: `wrappers/py_kraken/pyproject.toml`
- Create: `wrappers/py_kraken/py_kraken/__init__.py`
- Create: `wrappers/py_kraken/py_kraken/ook_data_model.py`
- Create: `wrappers/py_kraken/py_kraken/ook_read.py`
- Create: `wrappers/py_kraken/tests/test_readers.py`

**Interfaces:**
- Consumes: OOK JSON schema and OOK-exported MOD/SHD binary layout.
- Produces: `ConfigModel`, `FieldData`, `ModeProfileData`, `ModeData`, `read_json`, `read_shd`, `read_mod`.

- [ ] **Step 1: Write reader and model tests against generated fixtures**

Tests must assert:

```python
field = read_shd(generated_pressure_shd)
self.assertEqual(field.values.dtype, np.complex64)
self.assertEqual(field.values.shape, (len(field.source_depths), len(field.receiver_ranges), len(field.receiver_depths)))
self.assertTrue(np.isfinite(field.values).all())

modes = read_mod(generated_mod)
self.assertGreater(len(modes.profiles), 0)
self.assertGreater(modes.profiles[0].wavenumbers.size, 0)
self.assertEqual(modes.profiles[0].mode_shapes.shape[0], modes.profiles[0].wavenumbers.size)

with self.assertRaises(ValueError):
    read_shd(truncated_shd)
with self.assertRaises(ValueError):
    read_mod(truncated_mod)
```

Generate valid fixtures in `setUpClass` by calling the Task 2 CLI on `test/solve3_mode_loss.env`; create truncated fixtures by writing the first 32 bytes of each valid file to a temporary path.

- [ ] **Step 2: Run reader tests and verify import failure**

Run: `python -m unittest discover -s wrappers/py_kraken/tests -p test_readers.py -v`

Expected: FAIL because the package models and readers do not exist.

- [ ] **Step 3: Implement exact Pydantic models**

```python
class ConfigModel(BaseModel):
    model_config = ConfigDict(extra="allow")


class FieldData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    title: str
    frequency: float
    source_depths: np.ndarray
    receiver_ranges: np.ndarray
    receiver_depths: np.ndarray
    values: np.ndarray


class ModeProfileData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    profile_range: float
    depth: np.ndarray
    wavenumbers: np.ndarray
    group_velocity: np.ndarray
    mode_shapes: np.ndarray


class ModeData(BaseModel):
    profiles: list[ModeProfileData]
```

Validators require one-dimensional coordinate arrays, `values.shape == (NSz, NRr, NRz)`, complex field values, one-dimensional wavenumbers, and `mode_shapes.shape == (M, Nz)`.

- [ ] **Step 4: Implement strict readers**

Implement these exact entry points:

```python
def read_json(path: str | Path) -> ConfigModel:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    return ConfigModel.model_validate_json(source.read_text(encoding="utf-8"))


def read_shd(path: str | Path) -> FieldData:
    """Read Kraken direct-access SHD records and return source/range/depth order."""


def read_mod(path: str | Path) -> ModeData:
    """Read every Kraken MOD profile, including z, k, optional VG, and phi."""
```

For `read_shd`, port the direct-access record offsets from OOB `read_shd_bin.m`, reshape the file payload from record order into `(NSz, NRr, NRz)`, and then build `FieldData`. For `read_mod`, port the complete record traversal from `tools/read_mod_summary.py::parse_mod`; convert real/imaginary float pairs to `np.complex64`, return an empty `float64` group-velocity vector when the MOD record set does not contain VG, and advance the profile base record by the exact number of header, mode-shape, and wavenumber records consumed. Before each unpack, verify `offset + requested_bytes <= len(raw)` and raise `ValueError(f"{path} is truncated at record {record}")` when false.

- [ ] **Step 5: Add package metadata and run tests**

Use this dependency declaration:

```toml
[project]
name = "py-kraken"
version = "0.1.0"
requires-python = ">=3.9"
dependencies = ["numpy>=1.23", "pydantic>=2.0", "matplotlib>=3.6"]

[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[tool.setuptools]
packages = ["py_kraken"]
```

Run: `python -m unittest discover -s wrappers/py_kraken/tests -p test_readers.py -v`

Expected: all reader tests PASS.

- [ ] **Step 6: Commit Python models and readers**

```powershell
git add wrappers/py_kraken/pyproject.toml wrappers/py_kraken/py_kraken wrappers/py_kraken/tests/test_readers.py
git commit -m "feat: add Kraken Python result models and readers"
```

---

### Task 6: Add the Python high-level facade, plotting, typing, install flow, and demo

**Files:**
- Create: `wrappers/py_kraken/py_kraken/ook_interface.py`
- Create: `wrappers/py_kraken/py_kraken/ook_plot.py`
- Create: `wrappers/py_kraken/py_kraken/demo.py`
- Create: `wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi`
- Create: `wrappers/py_kraken/install.py`
- Modify: `wrappers/py_kraken/py_kraken/__init__.py`
- Create: `wrappers/py_kraken/tests/test_interface.py`
- Create: `wrappers/py_kraken/tests/test_plot.py`

**Interfaces:**
- Consumes: native module and data models/readers from Tasks 4–5.
- Produces: `OpenOceanKraken_interface` and all approved `ook_` methods.

- [ ] **Step 1: Write facade lifecycle and validation tests**

Test independent instances, optional singleton behavior, missing files, vector/range overload validation, pressure/velocity/mode results, and owned-array lifetime:

```python
first = OpenOceanKraken_interface(thread_num=1)
second = OpenOceanKraken_interface(thread_num=1)
self.assertIsNot(first, second)
self.assertIs(OpenOceanKraken_interface.get_instance(1), OpenOceanKraken_interface.get_instance(1))

with self.assertRaises(FileNotFoundError):
    first.ook_load_env("missing.env")
with self.assertRaises(ValueError):
    first.ook_set_receiver_range(start=0.0, end=1000.0, count=None)

first.ook_load_env(str(self.fixture))
first.ook_set_velocity_enable(True)
first.ook_run()
self.assertEqual(first.ook_get_pressure().values.ndim, 3)
self.assertEqual(first.ook_get_vertical_velocity().values.ndim, 3)
self.assertGreater(len(first.ook_get_modes().profiles), 0)
```

- [ ] **Step 2: Write headless plotting tests**

Use `matplotlib.use("Agg")`. Assert `plot_field(field)` returns an `Axes`, uses range in km on x, depth on y with inverted direction, and `plot_modes(modes, profile_index=0)` returns a `Figure` with at least two axes.

- [ ] **Step 3: Run tests and verify missing-module failures**

Run: `python -m unittest discover -s wrappers/py_kraken/tests -p "test_interface.py" -v`

Run: `python -m unittest discover -s wrappers/py_kraken/tests -p "test_plot.py" -v`

Expected: both fail because facade and plotting modules do not exist.

- [ ] **Step 4: Implement the high-level facade**

Implement exact methods:

```python
ook_load_env(path: str | Path) -> None
ook_load_json(path: str | Path) -> None
ook_export_json(path: str | Path) -> None
ook_get_config() -> ConfigModel
ook_set_frequency(frequency: float) -> None
ook_set_frequency_vector(frequencies: np.ndarray) -> None
ook_set_title(title: str) -> None
ook_set_source_depth(values=None, *, start=None, end=None, count=None) -> None
ook_set_receiver_depth(values=None, *, start=None, end=None, count=None) -> None
ook_set_receiver_range(values=None, *, start=None, end=None, count=None) -> None
ook_set_profile_range(values=None, *, start=None, end=None, count=None) -> None
ook_set_phase_speed(c_low: float, c_high: float) -> None
ook_set_mode_limit(limit: int) -> None
ook_set_ssp(areas: list[Range_Independent_Area]) -> None
ook_set_attenuation(mode: Atten_Mode) -> None
ook_set_grid_type(mode: Grid_Mode) -> None
ook_set_source_type(mode: Source_Mode) -> None
ook_set_run_mode(mode: Run_Mode) -> None
ook_set_coherence_type(mode: CoherenceType) -> None
ook_set_mode_type(mode: ModeType) -> None
ook_set_rmax(rmax: float) -> None
ook_set_reflection_top(coefficients: np.ndarray) -> None
ook_set_reflection_bottom(coefficients: np.ndarray) -> None
ook_set_sbp(pattern: np.ndarray, angles: np.ndarray) -> None
ook_set_velocity_enable(enabled: bool) -> None
ook_run() -> None
ook_run_eigen() -> None
ook_run_field() -> None
ook_clear() -> None
ook_free() -> None
ook_get_pressure() -> FieldData
ook_get_vertical_velocity() -> FieldData
ook_get_horizontal_velocity() -> FieldData
ook_get_modes() -> ModeData
ook_export_shd(path: str | Path, data_type: str = "pressure") -> None
ook_export_mod(path: str | Path) -> None
ook_export_result(path: str | Path) -> None
```

Keep the `ThreadPool` object as an instance field for the entire native interface lifetime.

- [ ] **Step 5: Implement plotting, typing, install, and demo**

`plot_field` computes transmission loss as `-20 * log10(max(abs(values[source_index]), tiny))`, plots range/depth using `pcolormesh`, and labels velocity plots with the selected component. `plot_modes` plots real/imaginary wavenumbers and selected normalized mode shapes. The demo loads `test/solve3_mode_loss.env`, runs one thread, prints result shapes, and saves pressure/mode PNGs without requiring interactive display.

The install script configures `BUILD_AS_PYTHON=ON` with the active interpreter, builds, performs `pip install -e`, then verifies imports of `py_kraken`, `OpenOceanKraken_interface`, and the embedded native module.

- [ ] **Step 6: Run the complete Python suite**

Run:

```powershell
python -m pip install -e wrappers/py_kraken --no-deps
python -m unittest discover -s wrappers/py_kraken/tests -v
python -m py_kraken.demo --fixture ../test/solve3_mode_loss.env --output-dir tmp/py_kraken_demo
```

Expected: all tests PASS; demo writes pressure and mode PNG files.

- [ ] **Step 7: Commit the Python high-level package**

```powershell
git add wrappers/py_kraken
git commit -m "feat: add py_kraken high-level workflow"
```

---

### Task 7: Add MATLAB data model, readers, plots, and batch tests

**Files:**
- Create: `wrappers/m_kraken/krakenDataModel.m`
- Create: `wrappers/m_kraken/read_shd.m`
- Create: `wrappers/m_kraken/read_mod.m`
- Create: `wrappers/m_kraken/plotshd.m`
- Create: `wrappers/m_kraken/plotmode.m`
- Create: `wrappers/m_kraken/demo_env.m`
- Create: `wrappers/m_kraken/demo_json.m`
- Create: `wrappers/m_kraken/tests/test_m_kraken.m`

**Interfaces:**
- Consumes: Task 2 CLI and its MOD/SHD output.
- Produces: `krakenDataModel.Write`, `run`, `getPressure`, `getVerticalVelocity`, `getHorizontalVelocity`, `getModes`, `plotPressure`, `plotModes`.

- [ ] **Step 1: Write a MATLAB batch test**

```matlab
function tests = test_m_kraken
tests = functiontests(localfunctions);
end

function testEnvCliAndReaders(testCase)
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
exe = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'solve3_mode_loss.env');
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's'));
model = krakenDataModel(exe);
model.NumThreads = 1;
model.loadEnv(fixture);
model.run(work, 'Velocity', true, 'Mod', true);
pressure = model.getPressure(work);
modes = model.getModes(work);
verifyEqual(testCase, ndims(pressure.values), 3);
verifyTrue(testCase, all(isfinite(real(pressure.values)), 'all'));
verifyGreaterThan(testCase, numel(modes.profiles), 0);
verifyGreaterThan(testCase, numel(modes.profiles(1).wavenumbers), 0);
end
```

Add this second test for JSON round-trip and invisible plots:

```matlab
function testJsonRoundTripAndPlots(testCase)
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
exe = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'solve3_mode_loss.env');
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's'));

model = krakenDataModel(exe);
model.NumThreads = 1;
model.loadEnv(fixture);
resultRoot = model.run(work, 'Json', true, 'Mod', true);

copy = krakenDataModel(exe);
copy.NumThreads = 1;
copy.loadJson([resultRoot '.json']);
jsonPath = copy.Write(work);
verifyTrue(testCase, isfile(jsonPath));

pressure = model.getPressure(work);
modes = model.getModes(work);
pressureFigure = figure('Visible', 'off');
modeFigure = figure('Visible', 'off');
figureCleanup = onCleanup(@() close([pressureFigure, modeFigure]));
model.plotPressure(work, 'Figure', pressureFigure);
model.plotModes(work, 'Figure', modeFigure);
verifyGreaterThan(testCase, numel(findall(pressureFigure, 'Type', 'axes')), 0);
verifyGreaterThan(testCase, numel(findall(modeFigure, 'Type', 'axes')), 1);
end
```

- [ ] **Step 2: Run MATLAB test and verify missing-class failure**

Run:

```powershell
matlab -batch "addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
```

Expected: FAIL because `krakenDataModel` and readers do not exist.

- [ ] **Step 3: Implement strict MATLAB readers**

`read_shd` returns a struct with fields `title`, `frequency`, `source_depths`, `receiver_ranges`, `receiver_depths`, and complex `values` ordered `[source, range, depth]`. `read_mod` returns `profiles`, where each profile has `profile_range`, `depth`, `wavenumbers`, `group_velocity`, and `mode_shapes`. Both readers check every `fread` count and throw `OpenOceanKraken:InvalidSHD` or `OpenOceanKraken:InvalidMOD` on truncation or invalid dimensions.

- [ ] **Step 4: Implement `krakenDataModel`**

Use a handle class with these public properties:

```matlab
Title
NumThreads
Frequency
FrequencyVector
SourceDepth
ReceiverDepth
ReceiverRange
ProfileRange
PhaseSpeed
SSP
Attenuation
GridType
SourceType
RunMode
CoherenceType
ModeType
MLimit
VelocityEnabled
ReflectionCoef
SBP
```

`Write(cachePath)` emits the existing OOK JSON schema with top-level keys `Title`, `freqinfo`, `AttenUnit`, `Pos`, `MLimit`, `NProf`, `RProf`, `hasModePos`, `sspInput`, `ReflectionCoef`, `SBP`, `is_Velocity`, `cLow`, `cHigh`, `Rmax`, `SourceType`, `RunMode`, `CoherenceType`, and `ModeType`. `loadEnv(path)` records an ENV input without converting it; `loadJson(path)` decodes the same schema into the public properties. `resultRoot = run(cachePath, Name, Value)` invokes the configured executable with quoted paths, checks `status == 0`, checks all requested files, and returns the output root. Result and plot methods call the strict readers.

- [ ] **Step 5: Implement plots and demos**

`plotshd` plots `-20*log10(max(abs(field), realmin))`, range in km, inverted depth axis, and a colorbar. `plotmode` creates wavenumber and selected mode-shape panels. `demo_env.m` runs the checked-in ENV fixture; `demo_json.m` constructs a valid minimal range-independent acoustic profile, writes JSON, runs, and plots.

- [ ] **Step 6: Run MATLAB tests**

Run:

```powershell
matlab -batch "addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
```

Expected: MATLAB reports all tests passed with no visible figures.

- [ ] **Step 7: Commit MATLAB wrappers**

```powershell
git add wrappers/m_kraken
git commit -m "feat: add m_kraken file-driven workflow"
```

---

### Task 8: Document, verify, and audit the complete integration

**Files:**
- Create: `README.md`
- Create: `wrappers/py_kraken/README.md`
- Create: `wrappers/m_kraken/README.md`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: all prior tasks.
- Produces: documented build/install/run commands and a single complete CTest/Python/MATLAB verification path.

- [ ] **Step 1: Write documentation smoke assertions**

Extend `tools/tests/test_dependency_layout.py`:

```python
def test_wrapper_documentation_names_real_commands(self):
    root_readme = (ROOT / "README.md").read_text(encoding="utf-8")
    python_readme = (ROOT / "wrappers/py_kraken/README.md").read_text(encoding="utf-8")
    matlab_readme = (ROOT / "wrappers/m_kraken/README.md").read_text(encoding="utf-8")
    self.assertIn("-DBUILD_AS_PYTHON=ON", root_readme)
    self.assertIn("OpenOceanKraken_interface", python_readme)
    self.assertIn("krakenDataModel", matlab_readme)
```

- [ ] **Step 2: Run the smoke assertion and verify missing docs fail**

Run: `python -m unittest tools.tests.test_dependency_layout -v`

Expected: FAIL because the three README files do not all exist.

- [ ] **Step 3: Write exact user documentation**

Document:

- Clean CMake commands for CLI/tests and `BUILD_AS_PYTHON=ON`.
- CLI flags with ENV and JSON examples.
- Editable Python install, minimal direct-run example, file reader example, and plotting example.
- MATLAB path setup, `krakenDataModel` example, direct ENV example, and headless test command.
- Result shape `[source, range, depth]`, units, output suffixes, ownership, and lifecycle rules.

- [ ] **Step 4: Register all executable test groups in CTest**

Ensure CTest contains these names when their prerequisites are enabled:

```text
OpenOceanKraken_interface_tests
OpenOceanKraken_snapshot_tests
OpenOceanKraken_cli_tests
OpenOceanKraken_python_tests
```

MATLAB remains an explicit batch command because MATLAB availability is environment-specific.

- [ ] **Step 5: Run clean full verification**

Run:

```powershell
cmake -S . -B build_final -DBUILD_TESTING=ON -DBUILD_AS_EXE=ON -DBUILD_AS_PYTHON=ON -DPython3_EXECUTABLE=(Get-Command python).Source
cmake --build build_final --parallel
ctest --test-dir build_final --output-on-failure
python -m unittest discover -s wrappers/py_kraken/tests -v
matlab -batch "addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
git diff --check
git status --short
```

Expected: configuration and build succeed; all CTest, Python, and MATLAB tests pass; `git diff --check` emits no errors; status contains only intended documentation changes before commit.

- [ ] **Step 6: Commit documentation and final CTest registration**

```powershell
git add README.md wrappers/py_kraken/README.md wrappers/m_kraken/README.md CMakeLists.txt tools/tests/test_dependency_layout.py
git commit -m "docs: document OOK wrapper workflows"
```

- [ ] **Step 7: Record final evidence**

Run:

```powershell
git status --short
git log -8 --oneline
```

Expected: clean working tree and a visible commit for each independently tested task.
