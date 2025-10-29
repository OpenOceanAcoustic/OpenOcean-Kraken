# OpenOcean-Kraken

OpenOcean-Kraken是一个海洋声学模型SDK，用于海洋声学计算。

## 目录结构

```
OpenOcean-Kraken/
├── include/           # 头文件目录
│   └── kkc_params.h   # 主要参数定义
├── src/               # 源代码目录
│   ├── set_pekeris.cpp # Pekeris波导设置
│   └── set_pekeris.h   # Pekeris波导头文件
├── main.cpp           # 主程序入口
├── CMakeLists.txt     # CMake构建配置
└── README.md          # 项目说明文档
```

## 构建说明

### 依赖项

- CMake 3.10或更高版本
- C++17兼容编译器
- Eigen3库
- nlohmann_json库

### 构建步骤

#### 1. 克隆仓库（如果适用）

```bash
git clone <仓库URL>
cd OpenOcean-Kraken
```

#### 2. 创建构建目录

```bash
mkdir build
cd build
```

#### 3. 配置项目（默认Release模式）

```bash
cmake ..
```

#### 4. 编译项目

```bash
cmake --build .
```

#### 5. 运行可执行文件

```bash
./OpenOcean-Kraken
```

## 启用Debug模式

要启用调试模式，可以在配置阶段设置BUILD_DEBUG选项：

```bash
cmake .. -DBUILD_DEBUG=ON
```

这将：
- 设置构建类型为Debug
- 添加调试信息标志
- 禁用优化（-O0）
- 定义DEBUG宏
- 启用详细的警告信息

## 平台特定说明

### Windows (MSVC)
- Debug模式：使用`/Od`（禁用优化）和`/Zi`（生成调试信息）
- Release模式：使用`/O2`（最大优化）

### Linux/macOS (GCC/Clang)
- Debug模式：使用`-O0 -g`（禁用优化，生成调试信息）
- Release模式：使用`-O3`（最高优化级别）

## 许可证

请参考项目根目录下的LICENSE文件。
