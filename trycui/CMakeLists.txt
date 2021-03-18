cmake_minimum_required(VERSION 3.16)
project(trycui)

set(CMAKE_CXX_STANDARD 14)

add_executable(trycui main.cpp)

target_link_libraries(trycui levmar lapack blas m)