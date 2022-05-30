cmake_minimum_required(VERSION 3.5)
project(CQOrbit)
include_directories(.)
include_directories(OrbitFinder_Files)

set(CMAKE_CXX_FLAGS "-std=c++17")

add_executable(GrossPitaevskiiFDM
        src/main.cpp
        )

