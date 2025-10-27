cmake_minimum_required(VERSION 3.10)

# 项目名称和版本
project(OpenOceanKraken VERSION 1.0)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置包含目录
include_directories(include)
include_directories(src/algorithm)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 定义源文件
set(SOURCES
    main.cpp
    src/algorithm/sspMod.cpp
    src/algorithm/AttenMod.cpp
)

# 创建可执行文件
add_executable(${PROJECT_NAME} ${SOURCES})