cmake_minimum_required(VERSION 3.5)


SET(BINARY ${CMAKE_PROJECT_NAME}_tst)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp *.cu *cuh)
SET(SOURCES ${TEST_SOURCES})

find_package(MPI REQUIRED)
find_package(CUDA REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories(${CUDA_INCLUDE_PATH})
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
add_executable(${BINARY} ${TEST_SOURCES})
add_test(NAME ${BINARY} COMMAND ${BINARY})
target_link_libraries(${BINARY} ${CMAKE_PROJECT_NAME}_lib)
target_link_libraries(${BINARY} gtest)
target_link_libraries(${BINARY} ${CUDA_LIBRARIES})
target_link_libraries(${BINARY} ${MPI_LIBRARIES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(${BINARY} -lnvToolsExt)