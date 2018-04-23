cmake_minimum_required(VERSION 2.8)

project(lbm)
set(CMAKE_C_COMPILER /usr/local/cuda-8.0/bin//nvcc)

add_definitions(-DGRAPHICS -DDEBUG)
set(CMAKE_CXX_FLAGS "-g") 
set (CMAKE_CXX_STANDARD 11)

file(GLOB SRCS "./src/*.cpp")

include_directories(${CMAKE_SOURCE_DIR}/headers)

add_executable(lbm ${SRCS})
target_link_libraries(lbm -lm)
