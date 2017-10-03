project(Functionnal)
cmake_minimum_required (VERSION 2.8)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set (ROOT ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${ROOT})
file (GLOB_RECURSE HEADER ${ROOT}/Functionnal/*.h)
file (GLOB_RECURSE IMPL ${ROOT}/Functionnal/*.hpp)

set (GNUPLOT ${ROOT}/trace.plot )

list( APPEND CMAKE_CXX_FLAGS "-std=c++11 -fPIC -Wall -Wpedantic -Wextra")

add_custom_target( install_res
    COMMAND ${CMAKE_COMMAND} -E copy ${GNUPLOT} ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying gnuplot script to build dir" VERBATIM
    DEPENDS ${GNUPLOT}
)

set (DEMOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/demos)

add_executable (bezier ${DEMOS_ROOT}/bezier.cpp ${HEADER} ${IMPL})
add_dependencies( bezier install_res)

