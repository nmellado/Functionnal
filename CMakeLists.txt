project(Functionnal)
cmake_minimum_required (VERSION 3.2)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set (ROOT ${CMAKE_CURRENT_SOURCE_DIR})

set (FUNCTIONNAL_DIRECTORY ${ROOT}/Functionnal)

file (GLOB_RECURSE HEADER ${FUNCTIONNAL_DIRECTORY}/*.h)
file (GLOB_RECURSE IMPL ${FUNCTIONNAL_DIRECTORY}/*.hpp)

set (GNUPLOT ${ROOT}/trace.plot )

set (CMAKE_CXX_STANDARD 11)

list( APPEND CMAKE_CXX_FLAGS "-fPIC -Wall -Wpedantic -Wextra")

add_custom_target( functionnal_copy_gnuplot_files
    COMMAND ${CMAKE_COMMAND} -E copy ${GNUPLOT} ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying gnuplot script to build dir" VERBATIM
    DEPENDS ${GNUPLOT}
)

set (DEMOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/demos)

add_library( functionnal INTERFACE )
target_sources( functionnal INTERFACE ${HEADER} ${IMPL} )

add_executable (bezier ${DEMOS_ROOT}/bezier.cpp )
add_dependencies( bezier functionnal functionnal_copy_gnuplot_files )
target_include_directories ( bezier PUBLIC ${FUNCTIONNAL_DIRECTORY} )

