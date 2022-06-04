cmake_minimum_required(VERSION 3.5)
project(GrossPitaevskiiFDM LANGUAGES CXX CUDA)
enable_language(CUDA)

include_directories(.)
include_directories(./src)


set(CMAKE_CXX_FLAGS "-std=c++17")

add_executable(GrossPitaevskiiFDM
        src/main.cpp
        src/utils.h
        src/utils.cpp
        src/domain/base_domain.h
        src/domain/base_domain.cpp
        src/domain/polar_domain.h
        src/domain/polar_domain.cpp
        src/domain/rectangular_domain.h
        src/domain/rectangular_domain.cpp
        src/initial_condition/initial_condition.h
        src/initial_condition/initial_condition.cpp
        src/serial_solver/base_serial_solver.h
        src/serial_solver/base_serial_solver.cpp
        src/serial_solver/forward_euler/forward_euler_rectangular_solver.h
        src/serial_solver/forward_euler/forward_euler_rectangular_solver.cpp
        src/serial_solver/crank_nicolson/crank_nicolson_rectangular_solver.h
        src/serial_solver/crank_nicolson/crank_nicolson_rectangular_solver.cpp
        )

#MPI linking
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(GrossPitaevskiiFDM ${MPI_C_LIBRARIES})


#Armadillo linking
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(GrossPitaevskiiFDM ${ARMADILLO_LIBRARIES})

target_link_libraries(GrossPitaevskiiFDM stdc++fs)