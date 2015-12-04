project(Interpo-Apps)
cmake_minimum_required (VERSION 2.8)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set (ROOT ${CMAKE_CURRENT_SOURCE_DIR})

set (HEADER  ${ROOT}/common.h
             ${ROOT}/functional.h
             ${ROOT}/constrainedBezier.h
             ${ROOT}/functionalEvalFunc.h )
set (IMPL    ${ROOT}/functional.hpp )
set (SRC     ${ROOT}/main.cpp )


set (GNUPLOT ${ROOT}/trace.plot )

list( APPEND CMAKE_CXX_FLAGS "-std=c++11 -fPIC -Wall -Wpedantic -Wextra")

add_custom_target( install_res
    COMMAND ${CMAKE_COMMAND} -E copy ${GNUPLOT} ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying gnuplot script to build dir" VERBATIM
    #FOREACH (RES IN LISTS ${RES_FILES})
    #COMMAND message(STATUS Installing ${RES_FILES})
        #COMMAND ${CMAKE_COMMAND} -E copy ${RES} ${CMAKE_CURRENT_BINARY_DIR}
    #ENDFOREACH (RES)
)

add_executable (bezier ${SRC} ${HEADER} ${IMPL} ${GNUPLOT})
add_dependencies( bezier install_res)

