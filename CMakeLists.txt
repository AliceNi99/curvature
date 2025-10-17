cmake_minimum_required(VERSION 3.15)
project(curvature LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization and CPU-specific flags
add_compile_options(-O3 -march=native -ffast-math)

add_executable(curvature
    src/main.cpp
    src/mesh.cpp
 "src/stl_reader.cpp")

find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
    message(STATUS "Using system Eigen")
    target_include_directories(curvature PRIVATE ${EIGEN3_INCLUDE_DIR})
elseif(EXISTS "${CMAKE_SOURCE_DIR}/external/Eigen/CMakeLists.txt" OR IS_DIRECTORY "${CMAKE_SOURCE_DIR}/external/Eigen")
    message(STATUS "Using local Eigen in external/Eigen")
    add_definitions(-DEIGEN_NO_DEBUG)
    target_include_directories(curvature PRIVATE ${CMAKE_SOURCE_DIR}/external)
else()
    message(FATAL_ERROR "Eigen3 not found! Install it system-wide or add it to external/Eigen.")
endif()

target_include_directories(curvature PRIVATE include)



