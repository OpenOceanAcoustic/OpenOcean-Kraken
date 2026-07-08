cmake_minimum_required(VERSION 3.10)
project(OpenOcean-Kraken)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

option(BUILD_DEBUG "Build with debug information" ON)

if(BUILD_DEBUG)
    message(STATUS "Building in debug mode")
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_definitions(DEBUG=1)
else()
    message(STATUS "Building in release mode")
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

set(OOK_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/src/module
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}
)

file(GLOB CORE_SOURCES
    "src/*.cpp"
    "src/algorithm/*.cpp"
    "src/module/*.cpp"
    "src/util/*.cpp"
)

add_library(OpenOceanKraken_core STATIC ${CORE_SOURCES})
target_include_directories(OpenOceanKraken_core PUBLIC ${OOK_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenOceanKraken_core -static-libgcc -static-libstdc++)

include(CTest)
if(BUILD_TESTING)
    add_executable(OpenOceanKraken_interface_tests for_test/test.cpp)
    target_include_directories(OpenOceanKraken_interface_tests PRIVATE ${OOK_INCLUDE_DIRS})
    target_link_libraries(OpenOceanKraken_interface_tests PRIVATE OpenOceanKraken_core -static-libgcc -static-libstdc++)
    add_test(
        NAME OpenOceanKraken_interface_tests
        COMMAND OpenOceanKraken_interface_tests ${CMAKE_CURRENT_SOURCE_DIR}/../test
    )
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../test_for_lcov")
        add_test(
            NAME OpenOceanKraken_lcov_fixture_tests
            COMMAND OpenOceanKraken_interface_tests ${CMAKE_CURRENT_SOURCE_DIR}/../test_for_lcov
        )
    endif()

    add_executable(ook_case_runner for_test/ook_case_runner.cpp)
    target_include_directories(ook_case_runner PRIVATE ${OOK_INCLUDE_DIRS})
    target_link_libraries(ook_case_runner PRIVATE OpenOceanKraken_core -static-libgcc -static-libstdc++)
endif()
