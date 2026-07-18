cmake_minimum_required(VERSION 3.15)
project(OpenOceanKraken LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

option(BUILD_DEBUG "Build with debug information" ON)
option(BUILD_AS_EXE "Build the OpenOceanKraken CLI" ON)
option(BUILD_AS_STATIC "Publish the OpenOceanKraken static library" ON)
option(BUILD_AS_SHARED "Build the OpenOceanKraken shared library" OFF)
option(BUILD_AS_PYTHON "Build the OpenOceanKraken Python module" OFF)

if(BUILD_DEBUG)
    message(STATUS "Building in debug mode")
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_definitions(DEBUG=1)
else()
    message(STATUS "Building in release mode")
    set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif()

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

file(GLOB CORE_SOURCES
    "src/*.cpp"
    "src/algorithm/*.cpp"
    "src/module/*.cpp"
    "src/util/*.cpp"
)

add_library(OpenOceanKraken_core STATIC ${CORE_SOURCES})
target_include_directories(OpenOceanKraken_core PUBLIC ${OOK_INCLUDE_DIRS})
set_target_properties(OpenOceanKraken_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(BUILD_AS_STATIC)
    set_target_properties(OpenOceanKraken_core PROPERTIES OUTPUT_NAME OpenOceanKraken_static)
endif()

if(BUILD_AS_EXE)
    add_executable(OpenOceanKraken main.cpp)
    target_link_libraries(OpenOceanKraken PRIVATE OpenOceanKraken_core)
    if(NOT MSVC)
        target_link_options(OpenOceanKraken PRIVATE -static-libgcc -static-libstdc++)
    endif()
endif()

if(BUILD_AS_SHARED)
    add_library(OpenOceanKraken_shared SHARED ${CORE_SOURCES})
    target_include_directories(OpenOceanKraken_shared PUBLIC ${OOK_INCLUDE_DIRS})
    set_target_properties(OpenOceanKraken_shared PROPERTIES OUTPUT_NAME OpenOceanKraken_shared)
endif()

if(BUILD_AS_PYTHON)
    add_subdirectory(${OOK_PYBIND11_DIR} EXCLUDE_FROM_ALL)
    pybind11_add_module(OpenOceanKraken_python src/bind/OpenOceanKrakenBind.cpp)
    target_include_directories(OpenOceanKraken_python PRIVATE ${OOK_INCLUDE_DIRS})
    target_link_libraries(OpenOceanKraken_python PRIVATE OpenOceanKraken_core)
    set(OOK_PYTHON_PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/py_kraken/py_kraken")
    set_target_properties(OpenOceanKraken_python PROPERTIES
        OUTPUT_NAME OpenOceanKraken
        LIBRARY_OUTPUT_DIRECTORY "${OOK_PYTHON_PACKAGE_DIR}"
        RUNTIME_OUTPUT_DIRECTORY "${OOK_PYTHON_PACKAGE_DIR}"
    )
    foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        set_target_properties(OpenOceanKraken_python PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY_${config} "${OOK_PYTHON_PACKAGE_DIR}"
            RUNTIME_OUTPUT_DIRECTORY_${config} "${OOK_PYTHON_PACKAGE_DIR}"
        )
    endforeach()
    if(NOT MSVC)
        target_link_options(OpenOceanKraken_python PRIVATE -static-libgcc -static-libstdc++)
        get_filename_component(OOK_COMPILER_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
        set(OOK_WINPTHREAD_DLL "${OOK_COMPILER_BIN_DIR}/libwinpthread-1.dll")
        if(WIN32 AND EXISTS "${OOK_WINPTHREAD_DLL}")
            add_custom_command(TARGET OpenOceanKraken_python POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${OOK_WINPTHREAD_DLL}"
                    "${OOK_PYTHON_PACKAGE_DIR}/libwinpthread-1.dll"
            )
        endif()
    endif()
endif()

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

    add_executable(OpenOceanKraken_snapshot_tests for_test/test_result_snapshots.cpp)
    target_include_directories(OpenOceanKraken_snapshot_tests PRIVATE ${OOK_INCLUDE_DIRS})
    target_link_libraries(OpenOceanKraken_snapshot_tests PRIVATE OpenOceanKraken_core -static-libgcc -static-libstdc++)
    add_test(
        NAME OpenOceanKraken_snapshot_tests
        COMMAND OpenOceanKraken_snapshot_tests ${CMAKE_CURRENT_SOURCE_DIR}/../test/MunkK.env
    )

    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../test/solve3_mode_loss.env")
        add_test(
            NAME OpenOceanKraken_case_runner_output_tests
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/for_test/test_case_runner_outputs.py
                --runner $<TARGET_FILE:ook_case_runner>
                --fixture ${CMAKE_CURRENT_SOURCE_DIR}/../test/solve3_mode_loss.env
        )
        add_test(
            NAME OpenOceanKraken_runner_validation_tests
            COMMAND ${Python3_EXECUTABLE} -m unittest tools.tests.test_runner_validation -v
        )
        set_tests_properties(
            OpenOceanKraken_runner_validation_tests
            PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
        )
        if(TARGET OpenOceanKraken AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../test/MunkK.env")
            add_test(
                NAME OpenOceanKraken_cli_tests
                COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/for_test/test_cli.py
                    --executable $<TARGET_FILE:OpenOceanKraken>
                    --fixture ${CMAKE_CURRENT_SOURCE_DIR}/../test/MunkK.env
                    -v
            )
        endif()
    endif()

    if(BUILD_AS_PYTHON AND Python3_Interpreter_FOUND)
        add_test(
            NAME OpenOceanKraken_python_tests
            COMMAND ${Python3_EXECUTABLE} -m unittest discover
                -s ${CMAKE_CURRENT_SOURCE_DIR}/wrappers/py_kraken/tests
                -v
        )
        set_tests_properties(
            OpenOceanKraken_python_tests
            PROPERTIES
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/wrappers/py_kraken"
        )
    endif()
endif()
