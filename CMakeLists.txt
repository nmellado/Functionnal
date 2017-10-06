project(Functionnal)
cmake_minimum_required (VERSION 3.2)

################################################################################
###### OPTIONS
option (FUNCTIONNAL_COMPILE_DEMOS "Compile demo applications" TRUE)

###### END OPTIONS
######
###### DO NOT EDIT BELOW THIS LINE
################################################################################


find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


set (FUNCTIONNAL_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB_RECURSE HEADER ${FUNCTIONNAL_DIRECTORY}/Functionnal/*.h)
file (GLOB_RECURSE IMPL ${FUNCTIONNAL_DIRECTORY}/Functionnal/*.hpp)

set (CMAKE_CXX_STANDARD 11)

add_library( functionnal INTERFACE )
target_sources( functionnal INTERFACE ${HEADER} ${IMPL} )
target_include_directories ( functionnal INTERFACE ${FUNCTIONNAL_DIRECTORY} )

install (DIRECTORY ${FUNCTIONNAL_DIRECTORY}/Functionnal/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

if ( FUNCTIONNAL_COMPILE_DEMOS )

    set (GNUPLOT ${FUNCTIONNAL_DIRECTORY}/trace.plot )

    add_custom_target( functionnal_copy_gnuplot_files
        COMMAND ${CMAKE_COMMAND} -E copy ${GNUPLOT} ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Copying gnuplot script to build dir" VERBATIM
        DEPENDS ${GNUPLOT}
    )

    set (DEMOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/demos)

    add_executable ( bezier ${DEMOS_ROOT}/bezier.cpp )
    add_dependencies( bezier functionnal functionnal_copy_gnuplot_files )
    target_include_directories ( bezier PUBLIC ${FUNCTIONNAL_DIRECTORY} )

    add_executable (cwise ${DEMOS_ROOT}/cwise.cpp )
    add_dependencies( cwise functionnal)
    target_include_directories ( cwise PUBLIC ${FUNCTIONNAL_DIRECTORY} )

    add_executable (polynomial ${DEMOS_ROOT}/polynomial.cpp )
    add_dependencies( polynomial functionnal)
    target_include_directories ( polynomial PUBLIC ${FUNCTIONNAL_DIRECTORY} )

endif ( FUNCTIONNAL_COMPILE_DEMOS )
