cmake_minimum_required(VERSION 3.10)
project(OpenOcean-Kraken)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Release 模式开启 O3 优化
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
# 添加debug选项
option(BUILD_DEBUG "Build with debug information" ON)

if(BUILD_DEBUG)
    message(STATUS "Building in debug mode")
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_definitions(DEBUG=1)
else()
    message(STATUS "Building in release mode")
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------
# 输出目录统一到 /bin
# -------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# 包含目录
include_directories(include)

# 查找Eigen库
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/module
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util
    ${CMAKE_CURRENT_SOURCE_DIR}/for_test
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 添加源文件
file(GLOB SOURCES "src/*.cpp" "src/algorithm/*.cpp" "src/module/*.cpp" "src/util/*.cpp" "src/bind/*.cpp" "for_test/*.cpp" "main.cpp")

# 创建可执行文件
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC -static-libgcc -static-libstdc++)

