project (HODLR3Dlib)
set(CMAKE_C_COMPILER "/usr/local/bin/gcc-10")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-10")
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)
find_package(OpenMP)

# Provide the input file that you want to solve for:
# This can also be passed as a command line argument:
if(NOT INPUT_FILE)
  #set (INPUT_FILE "examples/testHODLR3Dsolve.cpp")
  set (INPUT_FILE "examples/testHODLR3D.cpp")
endif()
# Here set the name for the output executable:
# This can also be passed as a command line argument:
if(NOT OUTPUT_EXECUTABLE)
    set (OUTPUT_EXECUTABLE "exec1")
endif()
# Change this accordingly to point to the Eigen root folder
# By default, it will pick up the environment variable EIGEN_PATH
set(EIGEN_PATH "$ENV{EIGEN_PATH}")

################ FINDING EIGEN ####################
if(NOT EIGEN_PATH)
    message ( FATAL_ERROR "Please point the environment variable EIGEN_PATH to the root directory of your Eigen3 installation.")
endif()

################ COMPILERS: APPROPRIATELY USE INTEL OR GNU ####################

# Toggle this to 1 when Intel Compilers are available:
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++17 -O3 -ipo -qopenmp -xHost -funroll-loops -g -liomp5 -lpthread -ldl")
else()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++17 -g -fopenmp -lpthread -ldl")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++17 -O3 -funroll-loops -g -fopenmp -lpthread -ldl")
endif ()

# Flags for testing coverage information
# Need by coveralls:
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS} --coverage")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS} --coverage")

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("include"
                    "${EIGEN_PATH}"
                   )

################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################

set(HODLR3Dlib_SRCS
    ${PROJECT_SOURCE_DIR}/src/kernel.cpp
    ${PROJECT_SOURCE_DIR}/src/ACA.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR3DTree.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR3D.cpp
    ${PROJECT_SOURCE_DIR}/src/gmres.cpp
   )

if(INPUT_FILE)
  add_executable(${OUTPUT_EXECUTABLE}
                 ${INPUT_FILE}
                 ${HODLR3Dlib_SRCS}
                )
endif()

################ SUMMARY OF BUILD INFORMATION ################

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("Input File       : " ${INPUT_FILE})
message("CXX Compiler     : " ${CMAKE_CXX_COMPILER})
message("CXX_FLAGS        : " ${CMAKE_CXX_FLAGS})
message("Eigen dir        : " ${EIGEN_PATH})
# message("MKL dir          : " ${MKLROOT})
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("")
