# py_kraken

`py_kraken` 包含三层接口：原生模块 `OpenOceanKraken`、Python facade `OpenOceanKraken_interface`，以及独立 SHD/MOD 读取和绘图函数。需要 Python 3.9+、NumPy、Pydantic 2 和 Matplotlib。

## 构建与安装

从 OOK 根目录执行：

```powershell
python wrappers\py_kraken\install.py `
  --compiler D:\program\mingw64\bin\g++.exe
```

或者已有 `.pyd` 后仅安装 Python 层：

```powershell
python -m pip install -e wrappers\py_kraken --no-deps
```

## 高层接口

```python
from py_kraken import OpenOceanKraken_interface

ook = OpenOceanKraken_interface(thread_num=2)
ook.ook_load_env("../test/MunkK.env")
ook.ook_set_velocity_enable(True)
ook.ook_run()

pressure = ook.ook_get_pressure()
modes = ook.ook_get_modes()
print(pressure.values.shape)  # [source, range, depth]
print(len(modes.profiles))
```

普通构造函数始终返回独立对象；只有显式调用 `OpenOceanKraken_interface.get_instance(thread_num)` 才使用可选单例。接口保留线程池直至原生对象生命周期结束。所有 `ook_get_*` 数组都是独立副本。

facade 提供 ENV/JSON 加载与导出、频率/深度/距离/SSP/衰减/模态控制、运行/清理/释放、压力/速度/模态读取以及 SHD/MOD 导出。向量设置支持 `values=[...]` 或完整的 `start=..., end=..., count=...`，不可混用。

### Biological 衰减

高层 facade 接受由字典或 `BiologicalAttenuationLayerModel` 组成的列表：

```python
from py_kraken import OpenOceanKraken_interface

ook = OpenOceanKraken_interface(thread_num=1)
ook.ook_set_biological_attenuation([
    {"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04},
    {"Z1": 40.0, "Z2": 60.0, "f0": 1200.0, "Q": 4.0, "a0": 0.02},
])
ook.ook_export_json("tmp/biological_case.json")
```

原生接口的等价配置如下：

```python
from py_kraken import native

layer = native.BiologicalAttenuationLayer()
layer.Z1 = 10.0
layer.Z2 = 30.0
layer.f0 = 1000.0
layer.Q = 5.0
layer.a0 = 0.04

mode = native.Atten_Mode()
mode.attnUnit = native.AttenuationUnit.MODE_W_db_per_lambda
mode.absModel = native.OceanAbsorptionModel.Biological
mode.biologicalLayers = [layer]

interface = native.Interface()
interface.set_AttenUnit(mode)
assert interface.to_json("biological_native.json")
```

字典键固定为 `Z1`、`Z2`、`f0`、`Q`、`a0`，单位依次为 m、m、Hz、无量纲和 dB/km。facade 未指定 `attenuation_unit` 时默认使用 `dB/lambda`；显式指定时应传入 native `AttenuationUnit` 枚举，而不是 JSON/MATLAB 使用的字符串。向 facade 传入空列表 `[]` 会清除 Biological 层，同时保留合法的 Biological 空数组配置。Python 会把层列表或层字段的类型错误、非有限数、`Z1 > Z2`、`f0 <= 0`、`Q <= 0`、`a0 < 0` 或层数超限等明显非法输入提前映射为 `ValueError`；错误的 `attenuation_unit` 类型不属于这一映射契约。C++ `set_AttenUnit` 仍是所有入口共享的最终权威校验边界。

native 与 facade 最终都写入标准顶层 `AttenUnit` JSON，其中模型名为 `Biological`、层数组为 `BiologicalLayers`。该输出及空数组、非法输入和失败回滚由真实 OOK `json_in_out` 导入路径验证。

## 文件读取与绘图

```python
from py_kraken import plot_field, plot_modes, read_mod, read_shd

field = read_shd("tmp/munk.shd")
modes = read_mod("tmp/munk.mod")
axes = plot_field(field, source_index=0)
figure = plot_modes(modes, profile_index=0)
```

SHD 数据为 `complex64` 的 `[source, range, depth]`。MOD reader 会读取文件中的全部环境剖面；MOD 格式本身不保存剖面距离，因此文件读取结果的 `profile_range` 为 `None`，原生内存结果则保留实际距离。

运行无界面 demo：

```powershell
python -m py_kraken.demo --fixture ..\test\MunkK.env --output-dir tmp\py_kraken_demo
```

## 测试

```powershell
python -m unittest discover -s wrappers\py_kraken\tests -v
```
