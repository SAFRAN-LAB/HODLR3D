set(CMAKE_C_COMPILER "/path/to/bin/gcc--v10")
set(CMAKE_CXX_COMPILER "/path/to/bin/g++-v10")
project(HODLR3D)

cmake_minimum_required (VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(EIGEN_PATH "/path/to/eigen3")
set(HOME_PATH "/path/to/MPI/HODLR3D")
set(MPI_INC "/path/to/open-mpi/4.1.1_2/include")
set(MPI_LIB "/path/to/open-mpi/4.1.1_2/lib")

include_directories("${EIGEN_PATH}" "${HOME_PATH}" "${MPI_INC}")

set(DTYPE_FLAG "-DUSE_DOUBLE") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} -std=c++17 -O3 -g -fopenmp -L${MPI_LIB} -lmpi -ldl -lpthread")

# add the executable
add_executable(hodlr3d testFMM3DRAMeff2.cpp)