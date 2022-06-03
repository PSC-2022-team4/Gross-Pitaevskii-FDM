cmake_minimum_required(VERSION 3.5)
project(GrossPitaevskiiFDM LANGUAGES CXX CUDA)
enable_language(CUDA)

include_directories(.)
include_directories(./src)


set(CMAKE_CXX_FLAGS "-std=c++17")

add_executable(GrossPitaevskiiFDM
        src/main.cpp
        )

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(GrossPitaevskiiFDM ${MPI_C_LIBRARIES})
