cmake_minimum_required(VERSION 2.8)

project(lbm)
add_executable(lbm main.c stub.c parameters.c)
target_link_libraries(lbm -lm)
