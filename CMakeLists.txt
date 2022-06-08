cmake_minimum_required(VERSION 3.5)

SET(proj_name "GrossPitaevskiiFDM")
project(${proj_name} LANGUAGES CXX CUDA)
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
        src/domain/rect_domain.h
        src/domain/rect_domain.cpp
        src/potential/base_potential.h
        src/potential/base_potential.cpp
        src/potential/harmonic_potential.h
        src/potential/harmonic_potential.cpp
        src/initial_condition/initial_condition.h
        src/initial_condition/initial_condition.cpp
        src/solver/base_solver.h
        src/solver/base_solver.cpp
        src/solver/serial_solver/forward_euler/fe_rect_solver.h
        src/solver/serial_solver/forward_euler/fe_rect_solver.cpp
        src/solver/parallel_solver/forward_euler/fe_rect_psolver.cuh
        src/solver/parallel_solver/forward_euler/fe_rect_psolver.cu
        src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh
        src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cu
        src/test/test_psolver/test_cn_rect_psolver.cu
        src/test/test_psolver/test_cn_rect_psolver.cuh
        )

#MPI linking
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(GrossPitaevskiiFDM ${MPI_LIBRARIES})


# GTest linking
find_package(GTest REQUIRED)
include_directories(${GTest_INCLUDE_PATH})
target_link_libraries(GrossPitaevskiiFDM ${GTest_LIBRARIES})


#Armadillo linking
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(GrossPitaevskiiFDM ${ARMADILLO_LIBRARIES})

target_link_libraries(GrossPitaevskiiFDM stdc++fs)

#OpenACC linking
find_package(OpenACC REQUIRED)
# set(OpenACC_CXX_FLAGS "-acc=gpu -acc=noautopar -Minfo=accel -gpu=cc${CUDA_ARCH}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenACC_CXX_FLAGS} ${OpenACC_CXX_FLAGS_CUSTOM}")
include_directories(${OPENACC_INCLUDE_DIRS})
target_link_libraries(GrossPitaevskiiFDM ${OPENACC_LIBRARIES})

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
target_link_libraries(GrossPitaevskiiFDM ${CUDA_LIBRARIES})

