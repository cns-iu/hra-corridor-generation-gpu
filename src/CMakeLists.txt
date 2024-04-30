cmake_minimum_required(VERSION 3.9)
project(3D-CORRIDOR-GPU CUDA CXX)
set(CGAL_DIR "$ENV{HOME}/CGAL-5.5.3")

set(CMAKE_C_COMPILER /usr/bin/gcc-10)
set(CMAKE_CXX_COMPILER /usr/bin/g++-10)

# find packages
find_package(CUDA QUIET REQUIRED)
find_package(CGAL REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost 1.45.0 COMPONENTS filesystem)
find_package(cpprestsdk REQUIRED)
find_package(assimp REQUIRED)


set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 14)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED OFF)  
set(Boost_USE_STATIC_RUNTIME OFF) 

# Pass options to NVCC
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -O3)

cuda_select_nvcc_arch_flags(ARCH_FLAGS)
# cuda_add_executable(pip point_in_polyhedron.cu)
# target_include_directories(pip PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# static library geometry 
file(GLOB GEOMETRY_SRC_FILES geometry/*.cpp)
add_library(geometry ${GEOMETRY_SRC_FILES})
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry)
target_link_libraries(geometry PUBLIC CGAL::CGAL)
target_link_libraries(geometry PUBLIC Eigen3::Eigen)
target_link_libraries(geometry PUBLIC assimp)

# static library corridor_gpu
set(CUDA_SEPARABLE_COMPILATION ON)
file(GLOB CORRIDOR_GPU_SRC_FILES corridor_gpu/*.cu)
cuda_add_library(corridor_gpu ${CORRIDOR_GPU_SRC_FILES})
set_property(TARGET corridor_gpu PROPERTY CUDA_ARCHITECTURES native)
target_include_directories(corridor_gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/corridor_gpu)
target_link_libraries(corridor_gpu ${CUDA_LIBRARIES})

#static library corridor_cpu
file(GLOB CORRIDOR_CPU_SRC_FILES corridor_cpu/*.cpp)
add_library(corridor_cpu ${CORRIDOR_CPU_SRC_FILES})
target_include_directories(corridor_cpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/corridor_cpu)
target_link_libraries(corridor_cpu geometry)


# cuda_add_executable(main corridor.cu main.cu)
# target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_link_libraries(main geometry)

# basic corridor generation
cuda_add_executable(main main.cpp)
target_link_libraries(main geometry)
target_link_libraries(main corridor_gpu)
target_link_libraries(main corridor_cpu)

# GPU corridor generation
cuda_add_executable(gpu_corridor_api server2.cpp)
target_link_libraries(gpu_corridor_api geometry)
target_link_libraries(gpu_corridor_api corridor_gpu)
target_link_libraries(gpu_corridor_api corridor_cpu)
target_link_libraries(gpu_corridor_api cpprestsdk::cpprest)


# verification
cuda_add_executable(comparison comparison.cpp)
target_include_directories(comparison PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/index)
target_link_libraries(comparison geometry)
target_link_libraries(comparison corridor_gpu)
target_link_libraries(comparison corridor_cpu)
# target_link_libraries(comparison cpprestsdk::cpprest)

# test aabb gpu
cuda_add_executable(test_aabbtree_gpu test_aabbtree.cpp)
# target_include_directories(test_aabbtree_gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/index)
target_link_libraries(test_aabbtree_gpu geometry)
target_link_libraries(test_aabbtree_gpu corridor_gpu)
target_link_libraries(test_aabbtree_gpu corridor_cpu)

if(Boost_FOUND)
    target_link_libraries(main ${Boost_LIBRARIES})
    target_link_libraries(gpu_corridor_api ${Boost_LIBRARIES})
    target_link_libraries(comparison ${Boost_LIBRARIES})
    target_link_libraries(test_aabbtree_gpu ${Boost_LIBRARIES})
endif()