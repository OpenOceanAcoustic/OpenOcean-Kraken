cmake_minimum_required(VERSION 3.20)
project(OpenOcean-Krakenc VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckCXXSourceCompiles)

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
set(OPENOCEANKRAKENC_FORTRAN_ROOT "" CACHE PATH
    "Directory containing the frozen krakenc and field Fortran executables")
option(OPENOCEANKRAKENC_NATIVE_OPTIMIZATION
       "Tune Release numerical kernels for the current CPU" OFF)
option(OPENOCEANKRAKENC_ENABLE_ASAN
       "Build OpenOcean-Krakenc with AddressSanitizer" OFF)
option(OPENOCEANKRAKENC_ENABLE_UBSAN
       "Build OpenOcean-Krakenc with UndefinedBehaviorSanitizer" OFF)

find_package(Threads REQUIRED)

add_library(OpenOceanKrakenc_core STATIC
    src/OpenOceanKrakencEnvironment.cpp
    src/OpenOceanKrakencParams.cpp
    src/algorithm/AcousticCase.cpp
    src/algorithm/ComplexDispersion.cpp
    src/algorithm/ComplexNumerics.cpp
    src/algorithm/SspInterpolation.cpp
    src/algorithm/ComplexMatrixBuilder.cpp
    src/algorithm/ComplexRootFinder.cpp
    src/algorithm/KrakencSolver.cpp
    src/algorithm/InverseIterationComplex.cpp
    src/algorithm/ComplexModeSolver.cpp
    src/algorithm/ComplexModeNormalization.cpp
    src/algorithm/ElasticCompound.cpp
    src/algorithm/ReflectionBoundary.cpp
    src/algorithm/ScatteringLoss.cpp
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
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/nlohmann-json>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(OpenOceanKrakenc_core PUBLIC Threads::Threads)
set_target_properties(OpenOceanKrakenc_core PROPERTIES EXPORT_NAME core)
if(OPENOCEANKRAKENC_ENABLE_ASAN OR OPENOCEANKRAKENC_ENABLE_UBSAN)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(OPENOCEANKRAKENC_SANITIZERS "")
        if(OPENOCEANKRAKENC_ENABLE_ASAN)
            list(APPEND OPENOCEANKRAKENC_SANITIZERS address)
        endif()
        if(OPENOCEANKRAKENC_ENABLE_UBSAN)
            list(APPEND OPENOCEANKRAKENC_SANITIZERS undefined)
        endif()
        list(JOIN OPENOCEANKRAKENC_SANITIZERS "," OPENOCEANKRAKENC_SANITIZER_LIST)
        set(OPENOCEANKRAKENC_REQUIRED_FLAGS_SAVED "${CMAKE_REQUIRED_FLAGS}")
        set(OPENOCEANKRAKENC_REQUIRED_LINK_OPTIONS_SAVED
            "${CMAKE_REQUIRED_LINK_OPTIONS}")
        set(CMAKE_REQUIRED_FLAGS
            "${CMAKE_REQUIRED_FLAGS} -fsanitize=${OPENOCEANKRAKENC_SANITIZER_LIST} -fno-omit-frame-pointer")
        list(APPEND CMAKE_REQUIRED_LINK_OPTIONS
            "-fsanitize=${OPENOCEANKRAKENC_SANITIZER_LIST}")
        string(REPLACE "," "_" OPENOCEANKRAKENC_SANITIZER_CACHE_SUFFIX
            "${OPENOCEANKRAKENC_SANITIZER_LIST}")
        set(OPENOCEANKRAKENC_SANITIZER_PROBE
            "OPENOCEANKRAKENC_SANITIZER_${OPENOCEANKRAKENC_SANITIZER_CACHE_SUFFIX}_AVAILABLE")
        check_cxx_source_compiles("int main() { return 0; }"
            ${OPENOCEANKRAKENC_SANITIZER_PROBE})
        set(CMAKE_REQUIRED_FLAGS "${OPENOCEANKRAKENC_REQUIRED_FLAGS_SAVED}")
        set(CMAKE_REQUIRED_LINK_OPTIONS
            "${OPENOCEANKRAKENC_REQUIRED_LINK_OPTIONS_SAVED}")
        if(NOT ${OPENOCEANKRAKENC_SANITIZER_PROBE})
            message(FATAL_ERROR
                "Requested sanitizer runtime is unavailable for ${CMAKE_CXX_COMPILER}; install the ASan/UBSan runtime or select a supported Clang/GCC toolchain")
        endif()
        target_compile_options(OpenOceanKrakenc_core PUBLIC
            "-fsanitize=${OPENOCEANKRAKENC_SANITIZER_LIST}"
            -fno-omit-frame-pointer
        )
        target_link_options(OpenOceanKrakenc_core PUBLIC
            "-fsanitize=${OPENOCEANKRAKENC_SANITIZER_LIST}"
        )
    elseif(MSVC)
        if(OPENOCEANKRAKENC_ENABLE_UBSAN)
            message(FATAL_ERROR "UBSan is not supported by the configured MSVC toolchain")
        endif()
        target_compile_options(OpenOceanKrakenc_core PUBLIC /fsanitize=address /Oy-)
        target_link_options(OpenOceanKrakenc_core PUBLIC /fsanitize=address)
    else()
        message(FATAL_ERROR
            "Sanitizers are not configured for compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()
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

add_executable(OpenOceanKrakenc_velocity_case_runner
    for_test/krakenc_velocity_case_runner.cpp
)
target_link_libraries(OpenOceanKrakenc_velocity_case_runner PRIVATE OpenOceanKrakenc_core)
if(MINGW)
    target_link_options(OpenOceanKrakenc_velocity_case_runner PRIVATE -static-libgcc -static-libstdc++)
endif()

install(TARGETS OpenOceanKrakenc_core
    EXPORT OpenOceanKrakencTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY third_party/eigen/Eigen
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT OpenOceanKrakencTargets
    FILE OpenOceanKrakencTargets.cmake
    NAMESPACE OpenOceanKrakenc::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OpenOceanKrakenc
)
configure_package_config_file(
    cmake/OpenOceanKrakencConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/OpenOceanKrakencConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OpenOceanKrakenc
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/OpenOceanKrakencConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/OpenOceanKrakencConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/OpenOceanKrakencConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OpenOceanKrakenc
)

add_executable(OpenOceanKrakenc_root_audit_runner
    for_test/root_audit_runner.cpp
)
target_link_libraries(OpenOceanKrakenc_root_audit_runner PRIVATE OpenOceanKrakenc_core)
if(MINGW)
    target_link_options(OpenOceanKrakenc_root_audit_runner PRIVATE -static-libgcc -static-libstdc++)
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
    add_executable(OpenOceanKrakenc_remediation_ssp_interpolation_tests
        for_test/remediation_ssp_interpolation_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_remediation_ssp_interpolation_tests
        PRIVATE OpenOceanKrakenc_core
    )
    target_compile_definitions(OpenOceanKrakenc_remediation_ssp_interpolation_tests PRIVATE
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_remediation_ssp_interpolation_tests PRIVATE
            -static-libgcc -static-libstdc++
        )
    endif()
    add_test(
        NAME OpenOceanKrakenc_remediation_ssp_interpolation_alg_022
        COMMAND OpenOceanKrakenc_remediation_ssp_interpolation_tests ALG-022
    )
    set_tests_properties(
        OpenOceanKrakenc_remediation_ssp_interpolation_alg_022
        PROPERTIES
            LABELS "remediation;ssp-interpolation;ALG-022"
            TIMEOUT 60
    )

    add_executable(OpenOceanKrakenc_remediation_runtime_tests
        for_test/remediation_runtime_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_remediation_runtime_tests
        PRIVATE OpenOceanKrakenc_core
    )
    target_compile_definitions(OpenOceanKrakenc_remediation_runtime_tests PRIVATE
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
        NAME OpenOceanKrakenc_remediation_runtime_alg_004
        COMMAND OpenOceanKrakenc_remediation_runtime_tests ALG-004
    )
    set_tests_properties(OpenOceanKrakenc_remediation_runtime_alg_004 PROPERTIES
        LABELS "remediation;runtime;ALG-004"
        TIMEOUT 5
    )
    add_test(
        NAME OpenOceanKrakenc_remediation_runtime_alg_007
        COMMAND OpenOceanKrakenc_remediation_runtime_tests ALG-007
    )
    set_tests_properties(OpenOceanKrakenc_remediation_runtime_alg_007 PROPERTIES
        LABELS "remediation;runtime;ALG-007"
        TIMEOUT 5
    )
    foreach(OPENOCEANKRAKENC_ALG019_SCENARIO IN ITEMS
            completion id timed)
        add_test(
            NAME OpenOceanKrakenc_remediation_runtime_alg_019_${OPENOCEANKRAKENC_ALG019_SCENARIO}
            COMMAND OpenOceanKrakenc_remediation_runtime_tests
                    ALG-019-${OPENOCEANKRAKENC_ALG019_SCENARIO}
        )
        set_tests_properties(
            OpenOceanKrakenc_remediation_runtime_alg_019_${OPENOCEANKRAKENC_ALG019_SCENARIO}
            PROPERTIES
                LABELS "remediation;runtime;ALG-019"
                TIMEOUT 1
        )
    endforeach()
    if(MINGW)
        target_link_options(OpenOceanKrakenc_remediation_runtime_tests PRIVATE
            -static-libgcc -static-libstdc++
        )
    endif()

    add_executable(OpenOceanKrakenc_remediation_root_mode_tests
        for_test/remediation_root_mode_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_remediation_root_mode_tests
        PRIVATE OpenOceanKrakenc_core
    )
    add_test(
        NAME OpenOceanKrakenc_remediation_root_mode_alg_001
        COMMAND OpenOceanKrakenc_remediation_root_mode_tests ALG-001
    )
    set_tests_properties(
        OpenOceanKrakenc_remediation_root_mode_alg_001
        PROPERTIES LABELS "remediation;root-mode;ALG-001"
    )
    add_test(
        NAME OpenOceanKrakenc_remediation_root_mode_alg_002
        COMMAND OpenOceanKrakenc_remediation_root_mode_tests ALG-002
    )
    set_tests_properties(
        OpenOceanKrakenc_remediation_root_mode_alg_002
        PROPERTIES LABELS "remediation;root-mode;ALG-002"
    )
    add_test(
        NAME OpenOceanKrakenc_remediation_root_mode_alg_006
        COMMAND OpenOceanKrakenc_remediation_root_mode_tests ALG-006
    )
    set_tests_properties(
        OpenOceanKrakenc_remediation_root_mode_alg_006
        PROPERTIES LABELS "remediation;root-mode;ALG-006"
    )
    add_test(
        NAME OpenOceanKrakenc_remediation_root_mode_alg_009
        COMMAND OpenOceanKrakenc_remediation_root_mode_tests ALG-009
    )
    set_tests_properties(
        OpenOceanKrakenc_remediation_root_mode_alg_009
        PROPERTIES LABELS "remediation;root-mode;ALG-009"
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_remediation_root_mode_tests PRIVATE
            -static-libgcc -static-libstdc++
        )
    endif()

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

    add_executable(OpenOceanKrakenc_phase_d_interface_tests
        for_test/phase_d_interface_tests.cpp
    )
    target_link_libraries(OpenOceanKrakenc_phase_d_interface_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_phase_d_interface_tests PRIVATE
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_phase_d_interface_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_phase_d_interface_tests
        COMMAND OpenOceanKrakenc_phase_d_interface_tests)

    add_executable(OpenOceanKrakenc_phase2_unit_tests for_test/phase2_unit_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_phase2_unit_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_phase2_unit_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(MINGW)
        target_link_options(OpenOceanKrakenc_phase2_unit_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME OpenOceanKrakenc_phase2_unit_tests COMMAND OpenOceanKrakenc_phase2_unit_tests)

    add_executable(OpenOceanKrakenc_phase3_unit_tests for_test/phase3_unit_tests.cpp)
    target_link_libraries(OpenOceanKrakenc_phase3_unit_tests PRIVATE OpenOceanKrakenc_core)
    target_compile_definitions(OpenOceanKrakenc_phase3_unit_tests PRIVATE
        OPENOCEANKRAKENC_WORKSPACE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
        OPENOCEANKRAKENC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
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
        NAME OpenOceanKrakenc_ssp_interpolation_regression_tests
        COMMAND ${Python3_EXECUTABLE}
                -m unittest for_test.test_ssp_interpolation_regression -v
    )
    set_tests_properties(
        OpenOceanKrakenc_ssp_interpolation_regression_tests
        PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            LABELS "remediation;ssp-interpolation;ALG-022"
            TIMEOUT 30
    )

    add_test(
        NAME OpenOceanKrakenc_remediation_cli_tests
        COMMAND ${Python3_EXECUTABLE}
                -m unittest for_test/test_remediation_cli.py -v
    )
    set_tests_properties(OpenOceanKrakenc_remediation_cli_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT
            "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>"
        LABELS "remediation;runtime;ALG-021"
        TIMEOUT 30
    )

    add_test(
        NAME OpenOceanKrakenc_numerical_consistency_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_numerical_consistency.py -v
    )
    set_tests_properties(OpenOceanKrakenc_numerical_consistency_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        TIMEOUT 120
    )

    add_test(
        NAME OpenOceanKrakenc_velocity_analysis_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_velocity_analysis.py -v
    )
    set_tests_properties(OpenOceanKrakenc_velocity_analysis_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_VELOCITY_RUNNER=$<TARGET_FILE:OpenOceanKrakenc_velocity_case_runner>"
        TIMEOUT 300
    )

    add_test(
        NAME OpenOceanKrakenc_oracle_configuration_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/oracle_configuration_tests.py -v
    )
    set_tests_properties(OpenOceanKrakenc_oracle_configuration_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        TIMEOUT 60
    )

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
        ENVIRONMENT "OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 240
    )
    add_test(
        NAME OpenOceanKrakenc_phase2_munk_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase2_munk.py -v
    )
    add_test(
        NAME OpenOceanKrakenc_phase_b_physics_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase_b_physics.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase_b_physics_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 1800
    )
    set_tests_properties(OpenOceanKrakenc_phase2_munk_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 360
    )
    add_test(
        NAME OpenOceanKrakenc_phase3_munk_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase3_munk.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase3_munk_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 420
    )
    add_test(
        NAME OpenOceanKrakenc_phase3_elastic_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase3_elastic.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase3_elastic_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 900
    )
    add_test(
        NAME OpenOceanKrakenc_phase_c_physics_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase_c_physics.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase_c_physics_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 1200
    )
    add_test(
        NAME OpenOceanKrakenc_phase3_reflection_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase3_reflection.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase3_reflection_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 900
    )
    add_test(
        NAME OpenOceanKrakenc_phase4_field_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase4_field.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase4_field_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 900
    )
    add_test(
        NAME OpenOceanKrakenc_phase4_range_dependent_tests
        COMMAND ${Python3_EXECUTABLE} -m unittest for_test/test_phase4_range_dependent.py -v
    )
    set_tests_properties(OpenOceanKrakenc_phase4_range_dependent_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "OPENOCEAN_KRAKENC_BINARY_DIR=$<TARGET_FILE_DIR:OpenOcean-Krakenc>;OPENOCEANKRAKENC_FORTRAN_ROOT=${OPENOCEANKRAKENC_FORTRAN_ROOT}"
        TIMEOUT 900
    )
endif()
