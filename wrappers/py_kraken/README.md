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
