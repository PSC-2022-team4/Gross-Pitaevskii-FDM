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
        #src/domain/base_domain.h
        #src/domain/base_domain.cpp
        src/domain/base_domain_copy.h
        src/domain/base_domain_copy.cpp
        #src/domain/polar_domain.h
        #src/domain/polar_domain.cpp
        #src/domain/rect_domain.h
        #src/domain/rect_domain.cpp
        src/domain/rect_domain_copy.h
        src/domain/rect_domain_copy.cpp
        src/initial_condition/initial_condition.h
        src/initial_condition/initial_condition.cpp
        src/serial_solver/base_serial_solver.h
        src/serial_solver/base_serial_solver.cpp
        src/serial_solver/forward_euler/fe_rect_solver.h
        src/serial_solver/forward_euler/fe_rect_solver.cpp
        #src/serial_solver/crank_nicolson/cn_rect_solver.h
        #src/serial_solver/crank_nicolson/cn_rect_solver.cpp
        src/serial_solver/crank_nicolson/cn_rect_solver_copy.h
        src/serial_solver/crank_nicolson/cn_rect_solver_copy.cpp
        src/parallel_solver/forward_euler/fe_rect_psolver.cuh
        src/parallel_solver/forward_euler/fe_rect_psolver.cu
        src/parallel_solver/crank_nicolson/cn_rect_psolver.cuh
        src/parallel_solver/crank_nicolson/cn_rect_psolver.cu
        )

#MPI linking
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(GrossPitaevskiiFDM ${MPI_LIBRARIES})

#if(MPI_COMPILE_FLAGS)
#  set_target_properties(GrossPitaevskiiFDM PROPERTIES
#    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#endif()

#if(MPI_LINK_FLAGS)
#  set_target_properties(GrossPitaevskiiFDM PROPERTIES
#    LINK_FLAGS "${MPI_LINK_FLAGS}")
#endif()

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

