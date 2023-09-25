cmake_minimum_required(VERSION 3.9)
project(3D-CORRIDOR-GPU)

set(CMAKE_BUILD_TYPE "Release")

find_package(CUDA QUIET REQUIRED)

# Pass options to NVCC
set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3
    )
cuda_select_nvcc_arch_flags(ARCH_FLAGS)
cuda_add_executable(pip point_in_polyhedron.cu)
target_include_directories(pip PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


