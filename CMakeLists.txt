cmake_minimum_required(VERSION 3.14)

project(quad-wavelet-tree)

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3 -march=native -DNDEBUG -flto")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
cmake_print_variables(CMAKE_CXX_FLAGS)

set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS YES)

# interface definitions
add_library(quad-wavelet-tree INTERFACE)
target_include_directories(quad-wavelet-tree INTERFACE include)

add_executable(test_qvector tests/test_qvector.cpp)
target_include_directories(test_qvector PRIVATE include tests)

add_executable(test_qwm tests/test_qwm.cpp)
target_include_directories(test_qwm PRIVATE include tests)

add_executable(create_dataset tests/create_dataset.cpp)
target_include_directories(create_dataset PRIVATE include tests)