cmake_minimum_required(VERSION 3.16)
project(marketCalibration)

set(CMAKE_CXX_STANDARD 14)

add_executable(main main.cpp EuropeanOption.cpp EuropeanOption.h
        Calibration.cpp Calibration.h SPXOption.cpp
        SPXOption.h VIXOption.cpp VIXOption.h Faddeeva.cc Faddeeva.hh utility.cpp utility.h)

target_link_libraries(main levmar lapack blas m)

