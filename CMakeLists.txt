cmake_minimum_required(VERSION 3.20)
project(OpenOcean-Krakenc VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckIPOSupported)
option(OPENOCEANKRAKENC_IPO "Enable interprocedural optimization for Release builds" ON)
if(OPENOCEANKRAKENC_IPO)
    check_ipo_supported(RESULT OPENOCEANKRAKENC_IPO_SUPPORTED OUTPUT OPENOCEANKRAKENC_IPO_ERROR)
    if(OPENOCEANKRAKENC_IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    endif()
endif()

option(BUILD_TESTING "Build OpenOcean-Krakenc tests" ON)
option(OPENOCEANKRAKENC_NATIVE_OPTIMIZATION
       "Tune Release numerical kernels for the current CPU" OFF)

add_library(OpenOceanKrakenc_core STATIC
    src/OpenOceanKrakencEnvironment.cpp
    src/OpenOceanKrakencParams.cpp
    src/algorithm/AcousticCase.cpp
    src/algorithm/ComplexDispersion.cpp
    src/algorithm/ComplexNumerics.cpp
    src/algorithm/ComplexMatrixBuilder.cpp
    src/algorithm/ComplexRootFinder.cpp
    src/algorithm/KrakencSolver.cpp
    src/algorithm/InverseIterationComplex.cpp
    src/algorithm/ComplexModeSolver.cpp
    src/algorithm/ComplexModeNormalization.cpp
    src/algorithm/ElasticCompound.cpp
    src/algorithm/ReflectionBoundary.cpp
    src/algorithm/ModeFileWriter.cpp
    src/algorithm/FieldSolver.cpp
    src/algorithm/CoupledFieldSolver.cpp
    src/module/ParameterAdapters.cpp
    src/module/env_in_out.cpp
    src/module/json_in_out.cpp
    src/module/ResultExport.cpp
    src/module/OpenOceanKrakencInterface.cpp
)
target_include_directories(OpenOceanKrakenc_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)
if(OPENOCEANKRAKENC_NATIVE_OPTIMIZATION AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(OpenOceanKrakenc_core PUBLIC
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-mtune=native>
    )
endif()

add_executable(OpenOcean-Krakenc main.cpp)
target_link_libraries(OpenOcean-Krakenc PRIVATE OpenOceanKrakenc_core)
if(MINGW)
    target_link_options(OpenOcean-Krakenc PRIVATE -static-libgcc -static-libstdc++)
endif()

add_executable(OpenOceanKrakenc_acoustic_case_runner
    for_test/krakenc_acoustic_case_runner.cpp
)
target_link_libraries(OpenOceanKrakenc_acoustic_case_runner PRIVATE OpenOceanKrakenc_core)
if(MINGW)
    target_link_options(OpenOceanKrakenc_acoustic_case_runner PRIVATE -static-libgcc -static-libstdc++)
endif()

add_executable(OpenOceanKrakenc_phase3_munk_runner
    for_test/phase3_munk_runner.cpp
)
target_link_libraries(OpenOceanKrakenc_phase3_munk_runner PRIVATE OpenOceanKrakenc_core)
if(MINGW)
    target_link_options(OpenOceanKrakenc_phase3_munk_runner PRIVATE -static-libgcc -static-libstdc++)
endif()

include(CTest)
find_package(Python3 COMPONENTS Interpreter)
if(BUILD_TESTING)
    add_executable(OpenOceanKrakenc_api_compatibility_tests
        for_test/api_compatibility_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_api_compatibility_tests PRIVATE OpenOceanKrakenc_core)
    if(MINGW)
        target_link_options(OpenOceanKrakenc_api_compatibility_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_api_compatibility_tests
             COMMAND OpenOceanKrakenc_api_compatibility_tests)

    add_executable(OpenOceanKrakenc_parameter_adapter_tests
        for_test/parameter_adapter_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_parameter_adapter_tests PRIVATE OpenOceanKrakenc_core)
    if(MINGW)
        target_link_options(OpenOceanKrakenc_parameter_adapter_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_parameter_adapter_tests
             COMMAND OpenOceanKrakenc_parameter_adapter_tests)

    add_executable(OpenOceanKrakenc_env_io_tests for_test/env_io_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_env_io_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_env_io_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_env_io_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_env_io_tests COMMAND OpenOceanKrakenc_env_io_tests)

    add_executable(OpenOceanKrakenc_json_io_tests for_test/json_io_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_json_io_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_json_io_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_json_io_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_json_io_tests COMMAND OpenOceanKrakenc_json_io_tests)

    add_executable(OpenOceanKrakenc_interface_in_memory_tests
        for_test/interface_in_memory_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_interface_in_memory_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_interface_in_memory_tests PRIVATE
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_interface_in_memory_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_interface_in_memory_tests
             COMMAND OpenOceanKrakenc_interface_in_memory_tests)

    add_executable(OpenOceanKrakenc_export_compatibility_tests
        for_test/export_compatibility_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_export_compatibility_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_export_compatibility_tests PRIVATE
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_export_compatibility_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_export_compatibility_tests
             COMMAND OpenOceanKrakenc_export_compatibility_tests)

    add_executable(OpenOceanKrakenc_contract_tests for_test/contract_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_contract_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_contract_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_contract_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_contract_tests COMMAND OpenOceanKrakenc_contract_tests)

    add_executable(OpenOceanKrakenc_phase2_unit_tests for_test/phase2_unit_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_phase2_unit_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_phase2_unit_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_phase2_unit_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_phase2_unit_tests COMMAND OpenOceanKrakenc_phase2_unit_tests)

    add_executable(OpenOceanKrakenc_phase3_unit_tests for_test/phase3_unit_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_phase3_unit_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_phase3_unit_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_phase3_unit_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_phase3_unit_tests COMMAND OpenOceanKrakenc_phase3_unit_tests)

    add_executable(OpenOceanKrakenc_phase4_unit_tests for_test/phase4_unit_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_phase4_unit_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_phase4_unit_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_phase4_unit_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_phase4_unit_tests COMMAND OpenOceanKrakenc_phase4_unit_tests)

    add_test(NAME OpenOceanKrakenc_cli_version COMMAND OpenOcean-Krakenc --version)
    set_tests_properties(OpenOceanKrakenc_cli_version PROPERTIES
        PASS_REGULAR_EXPRESSION "OpenOcean-Krakenc 0.1.0"
    )
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
    add_test(
        NAME OpenOceanKrakenc_cli_json_input_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_cli_json_input.py -v
    )
    set_tests_properties(OpenOceanKrakenc_cli_json_input_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 240
    )
    add_test(
        NAME OpenOceanKrakenc_env_json_equivalence_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_env_json_equivalence.py -v
    )
    set_tests_properties(OpenOceanKrakenc_env_json_equivalence_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 1800
    )
    add_test(
        NAME OpenOceanKrakenc_reference_pipeline_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_reference_pipeline.py -v
    )
    set_tests_properties(OpenOceanKrakenc_reference_pipeline_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        TIMEOUT 240
    )
    add_test(
        NAME OpenOceanKrakenc_phase2_munk_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase2_munk.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase2_munk_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 360
    )
    add_test(
        NAME OpenOceanKrakenc_phase3_munk_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase3_munk.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase3_munk_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 420
    )
    add_test(
        NAME OpenOceanKrakenc_phase3_elastic_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase3_elastic.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase3_elastic_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 900
    )
    add_test(
        NAME OpenOceanKrakenc_phase3_reflection_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase3_reflection.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase3_reflection_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 900
    )
    add_test(
        NAME OpenOceanKrakenc_phase4_field_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase4_field.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase4_field_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 900
    )
    add_test(
        NAME OpenOceanKrakenc_phase4_range_dependent_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase4_range_dependent.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase4_range_dependent_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        TIMEOUT 900
    )
endif()
