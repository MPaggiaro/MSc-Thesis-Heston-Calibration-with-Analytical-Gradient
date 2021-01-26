cmake_minimum_required(VERSION 3.16)
project(lin)

set(CMAKE_CXX_STANDARD 14)

add_executable(lin DIS.cpp Faddeeva.cc Faddeeva.hh)

target_link_libraries(lin levmar lapack blas m)