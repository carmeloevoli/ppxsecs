project(ppgamma C CXX)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_C_COMPILER "gcc" )
set( CMAKE_C_FLAGS "-O3" )
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-O3 -std=c++11" )
set( PPGAMMA_INCLUDES )
set( PPGAMMA_LIBRARIES )

# include
include_directories( "${CMAKE_SOURCE_DIR}/include" )

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
    set(HAVE_LIBGSL 1)
    list(APPEND PPGAMMA_INCLUDES ${GSL_INCLUDE_DIR})
    list(APPEND PPGAMMA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
    message(FATAL_ERROR "Require GSL: set GSL_HOME")
endif(GSL_FOUND)

# cparamlib (provided)
add_subdirectory(libs/cparamlib)
list(APPEND PPGAMMA_LIBRARIES cparamlib)
list(APPEND PPGAMMA_INCLUDES libs/cparamlib/include)

# add source files
set( PPGAMMA_SRCS
	src/AAfrag.cpp
	src/kamae.cpp
	src/kelner-aharonian.cpp
	)

# make library
include_directories (include ${PPGAMMA_INCLUDES})
add_library( PPGAMMA_LIB ${PPGAMMA_SRCS} ${PPGAMMA_INCLUDES} )

# make executables
add_executable(ppgamma src/main.cpp)
target_link_libraries (ppgamma PPGAMMA_LIB ${PPGAMMA_LIBRARIES})

# show summary
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_C_COMPILER ...... = ${CMAKE_C_COMPILER}")
message (STATUS "CMAKE_C_FLAGS ......... = ${CMAKE_C_FLAGS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "PPGAMMA_INCLUDES ...... = ${PPGAMMA_INCLUDES}")
message (STATUS "PPGAMMA_LIBRARIES ..... = ${PPGAMMA_LIBRARIES}")
