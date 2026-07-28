# m_kraken

`m_kraken` 使用 OOK CLI 驱动计算，不依赖 MEX。`krakenDataModel` 管理 ENV/JSON 输入、CLI 参数和结果根路径；`read_shd`、`read_mod` 严格检查固定记录并返回 MATLAB struct；`plotshd`、`plotmode` 提供压力和模态绘图。

## 使用 ENV

```matlab
addpath('wrappers/m_kraken');
model = krakenDataModel(fullfile(pwd, 'bin', 'OpenOceanKraken.exe'));
model.NumThreads = 2;
model.loadEnv(fullfile(fileparts(pwd), 'test', 'MunkK.env'));
root = model.run('tmp/matlab_case', 'Velocity', true, 'Mod', true);

pressure = model.getPressure(root);
vertical = model.getVerticalVelocity(root);
modes = model.getModes(root);
model.plotPressure(root);
model.plotModes(root);
```

`pressure.values`、速度场的维度均为 `[source, range, depth]`。普通压力文件是 `.shd`；速度结果使用 `_P.shd`、`_V.shd`、`_H.shd`；模态文件是 `.mod`。

## JSON

```matlab
model.loadJson('case.json');
jsonPath = model.Write('tmp/json_cache');
model.loadJson(jsonPath);
root = model.run('tmp/json_result', 'Mod', true);
```

`Write` 保留 OOK schema，并显式保持单元素向量/SSP/layer 为 JSON 数组，避免 MATLAB JSON round-trip 改变输入形状。可运行 `demo_env` 和 `demo_json` 查看完整的无界面 PNG 输出流程。

## Biological 衰减

```matlab
model.loadJson('case.json');  % 先加载完整的基准配置
model.setBiologicalAttenuation([
    10.0 30.0 1000.0 5.0 0.04
    40.0 60.0 1200.0 4.0 0.02
], 'dB/lambda');
jsonPath = model.Write('tmp/biological_case.json');
```

矩阵列顺序固定为 `[Z1 Z2 f0 Q a0]`，单位依次为 m、m、Hz、无量纲和 dB/km。衰减单位只接受以下六个完整、区分大小写的规范字符串：`dB/m/kHz`、`Params_Lose`、`dB/m`、`Nepers/m`、`Quality_Factor`、`dB/lambda`。

合法空层必须是数值、实数的 `0×5` 矩阵，推荐用 `zeros(0, 5)` 构造；`[]` 的形状是 `0×0`，不会被当作空层配置。`setBiologicalAttenuation` 准备包含 `OceanAbsorptionModel` 和 `BiologicalLayers` 的衰减配置，随后由 `Write` 写入标准顶层 `AttenUnit`。输出的 JSON 最终由真实 OOK CLI 的 `json_in_out` 导入路径校验。

## 测试

```powershell
matlab -batch "restoredefaultpath; addpath('wrappers/m_kraken'); addpath('wrappers/m_kraken/tests'); results=runtests('wrappers/m_kraken/tests'); assertSuccess(results)"
```

这里使用 `restoredefaultpath` 是为了防止用户路径中的同名 `runtests.m` 遮蔽 MATLAB 内置测试框架；如环境没有该问题可以省略。
