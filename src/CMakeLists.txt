cmake_minimum_required(VERSION 3.5)
enable_language(CUDA)

SET( BINARY ${CMAKE_PROJECT_NAME} )
SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )
SET( LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )

file(GLOB_RECURSE SRC_FILES LIST_DIRECTORIES false *.h *.cpp *.cu *cuh)
ADD_EXECUTABLE( ${BINARY}_run  ${SRC_FILES} )
ADD_LIBRARY( ${BINARY}_lib STATIC ${SRC_FILES} )

#MPI linking
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${BINARY}_run ${MPI_LIBRARIES})
target_link_libraries(${BINARY}_lib ${MPI_LIBRARIES})


#Armadillo linking
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${BINARY}_run ${ARMADILLO_LIBRARIES})
target_link_libraries(${BINARY}_lib ${ARMADILLO_LIBRARIES})

# stl filesystem linking
target_link_libraries(${BINARY}_run stdc++fs)
target_link_libraries(${BINARY}_lib stdc++fs)

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
target_link_libraries(${BINARY}_run ${CUDA_LIBRARIES})
target_link_libraries(${BINARY}_lib ${CUDA_LIBRARIES})


target_link_libraries(${BINARY}_run -lnvToolsExt)
target_link_libraries(${BINARY}_lib -lnvToolsExt)