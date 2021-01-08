cmake_minimum_required(VERSION 3.16)
project(ex1)

set(CMAKE_CXX_STANDARD 14)

add_executable(ex1 main.cpp EuropeanOption.cpp EuropeanOption.h)
