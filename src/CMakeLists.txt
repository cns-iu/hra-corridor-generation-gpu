cmake_minimum_required(VERSION 3.9)
project(3D-CORRIDOR-GPU)

# find packages
find_package(CUDA QUIET REQUIRED)
find_package(CGAL REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost 1.45.0 COMPONENTS filesystem)

set(CMAKE_BUILD_TYPE "Release")
# set(CGAL_DIR "$ENV{HOME}/CGAL-5.5.2")
set(CMAKE_CXX_STANDARD 14)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED OFF)  
set(Boost_USE_STATIC_RUNTIME OFF) 

# Pass options to NVCC
set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3
    )

cuda_select_nvcc_arch_flags(ARCH_FLAGS)
# cuda_add_executable(pip point_in_polyhedron.cu)
# target_include_directories(pip PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


# file(GLOB GEOMETRY_SRC_FILES geometry/*.cpp)
# add_library(geometry ${GEOMETRY_SRC_FILES})
# target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry)
# target_link_libraries(geometry PUBLIC CGAL::CGAL)
# target_link_libraries(geometry PUBLIC Eigen3::Eigen)

cuda_add_executable(main corridor.cu main.cu)
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_link_libraries(main geometry)

if(Boost_FOUND)
    target_link_libraries(main ${Boost_LIBRARIES})
endif()