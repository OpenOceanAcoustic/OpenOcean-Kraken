cmake_minimum_required(VERSION 3.20)
project(OpenOceanKrakencConsumer LANGUAGES CXX)

find_package(OpenOceanKrakenc CONFIG REQUIRED)

add_executable(openocean_krakenc_consumer main.cpp)
target_link_libraries(openocean_krakenc_consumer PRIVATE OpenOceanKrakenc::core)
target_compile_features(openocean_krakenc_consumer PRIVATE cxx_std_17)
if(MINGW)
    target_link_options(openocean_krakenc_consumer PRIVATE
        -static-libgcc -static-libstdc++)
endif()
