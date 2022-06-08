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
