cmake_minimum_required(VERSION 3.5)
project(GrossPitaevskiiFDM)
include_directories(.)
include_directories(./src)


set(CMAKE_CXX_FLAGS "-std=c++17")

add_executable(GrossPitaevskiiFDM
        main.cpp
        )

