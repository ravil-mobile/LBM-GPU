cmake_minimum_required(VERSION 2.8)

project(lbm)

set(CUDA_HOST_COMPILER /usr/bin/gcc-6)
find_package(CUDA REQUIRED)

<<<<<<< HEAD
# HASAN: for some reason CUDA_HOST_COMPILER is not working
# TEMPORARY WORKAROUND - use --compiler-binder to specify gcc executable
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -O2;--compiler-bindir /usr/bin/gcc-6; -std=c++11; -arch=sm_50; -Xptxas; -v)
=======
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -O2; -std=c++11; -arch=sm_20; -Xptxas -v; --maxrregcount=32)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -O3; -std=c++11; -arch=sm_20; -Xptxas -v)
>>>>>>> 5dab4f8... Last commit in the 'ravil' branch

add_definitions(-DGRAPHICS -DDEBUG)
#add_definitions(-DDEBUG)
file(GLOB SRCS "./src/*.cpp" "./src/*.cu")

#include_directories("${CUDA_INCLUDE_DIRS}")
include_directories(${CMAKE_SOURCE_DIR}/headers)

cuda_add_executable(lbm ${SRCS})
target_link_libraries(lbm -lm -lglfw -lGL -lGLEW 
                      ${CUDA_LIBRARIES}
                      ${CUDA_cublas_LIBRARY})
# Add CUBLAS to target
