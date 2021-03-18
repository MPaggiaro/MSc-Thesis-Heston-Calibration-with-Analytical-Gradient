cmake_minimum_required(VERSION 3.16)
project(validation)

set(CMAKE_CXX_STANDARD 14)

add_executable(validation validation.cpp EuropeanOption.cpp EuropeanOption.h
        Calibration.cpp Calibration.h SPXOption.cpp
        SPXOption.h VIXOption.cpp VIXOption.h Faddeeva.cc Faddeeva.hh doublecomparison.cpp doublecomparison.h)

target_link_libraries(validation levmar lapack blas m)

