cmake_minimum_required(VERSION 2.8)

project(lbm)
set(CMAKE_C_COMPILER /usr/local/cuda-8.0/bin//nvcc)
set(CMAKE_C_FLAGS "-g")
set (CMAKE_CXX_STANDARD 11)

set(SRCS ./src/main.cpp
         ./src/init.cpp
         ./src/parameters.cpp
         ./src/stub.cpp
         ./src/collision.cpp
         ./src/stream.cpp
         ./src/scanning.cpp
         ./src/helper.cpp
         ./src/boundary_conditions.cpp)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/headers)

add_executable(lbm ${SRCS})
target_link_libraries(lbm -lm)
