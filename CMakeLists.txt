cmake_minimum_required(VERSION 3.5)

SET(proj_name "GrossPitaevskiiFDM")
project(${proj_name} LANGUAGES CXX CUDA)
enable_language(CUDA)

set(CMAKE_CXX_FLAGS "-std=c++17")

add_subdirectory(src)
add_subdirectory(test)
