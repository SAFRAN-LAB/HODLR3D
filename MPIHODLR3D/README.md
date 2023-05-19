# Documentation - HODLR3D (MPI version)

This section explains how to replicate the results of scaling parallel HODLR3D based on MPI.

To do so, follow these steps:

1. The Eigen library can be downloaded from its [website](https://eigen.tuxfamily.org/index.php?title=Main_Page).
2. The code is tested with OpenMPI 4.1.1, which can be downloaded from its [website](https://www.open-mpi.org/software/ompi/v4.1/).
3. Set the following variables in the "CMakeLists.txt" file:
    - CMAKE_C_COMPILER: GCC version greater than GCC9
    - CMAKE_CXX_COMPILER: GCC version greater than GCC9
    - EIGEN_PATH: Set the path for the Eigen library.
    - HOME_PATH: Provide the path where the .cpp file is located.
    - MPI_INC: Provide the path to the Open MPI ‘include’ directory - path/to/open-mpi/4.1.x_x/include
    - MPI_LIB: Provide the path to the Open MPI ‘library’ directory - path/to/open-mpi/4.1.x_x/lib

Note: The code has been tested with other MPI wrapper compilers as well. The HODLR3D code has been tested with mpicxx, Intel-based mpi wrapper compilers. The scaling does not affect due to changes in the compiler.

**Sample CMakeLists.txt:**

```makefile
set(CMAKE_C_COMPILER "/path/to/bin/gcc-v10")
set(CMAKE_CXX_COMPILER "/path/to/bin/g++-v10")
project(HODLR3D)

cmake_minimum_required (VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(EIGEN_PATH "/path/to/eigen3")
set(HOME_PATH "/path/to/HODLR3D")
set(MPI_INC "/path/to/open-mpi/4.1.1_2/include")
set(MPI_LIB "/path/to/open-mpi/4.1.1_2/lib")
```

**Installation and Building**

Follow these steps:

1. Create a build directory.
2. Use `cmake /path/to/CMakeLists.txt`.
3. Use the `make` command to build the project. This will create an executable called "hodlr3d".
4. For a cluster with multiple nodes, decide the number of MPI processes to run in parallel. As described in the article, choose the number of MPI processes to be a power of 8. Ensure that you have the necessary nodes available.
5. Use the command `mpiexec.hydra -np 64 -genv I_MPI_PIN=1 -genv I_MPI_FABRICS=shm:ofi -hostfile $PBS_NODEFILE ./hodlr3d 50 10 1 6 1 > H3_1_50.txt` to run the program. Replace "64" with the number of MPI processes you want to run. The output will be stored in the file "H3_1_50.txt".

The inputs to the executable `./hodlr3d x1 x2 x3 x4 x5` are mandatory. Each input is explained below:

- x1 → cubeRootN → determines the system size, N, which is calculated as N = pow(cubeRootN, 3).
- x2 → nParticlesInLeafAlong1D → determines the maximum number of particles in a leaf node, calculated as pow(nParticlesInLeafAlong1D, 3).
- x3 → L is the half-side length of the cube and represents the computational domain.
- x4 → TOL_POW is the tolerance set for the ACA routine.
- x5 → Qchoice is used to select the kernel you want to use. For various choices, refer to the “kernel.hpp” file.

**Sample Installation and Building:**

```bash
user@computer HODLR3D$ mkdir build && cd build
user@computer build$ cmake ..
user@computer build$ make
user@computer build$ mpiexec.hydra -np 2 -genv I_MPI_PIN=1 -genv I_MPI_FABRICS=shm:ofi -hostfile $PBS_NODEFILE ./hodlr3d 50 10 1 6 1 > H3D_2_50.txt
```

Sample output file "H3D_2_50.txt".
```text
MPI Code with 2 processors..
Tree formed.. with 3 levels
System setting - 0,6
MPI Process Information set to tree
Target Level1
Scheduled ...
MPI Code with 2 processors..
Tree formed.. with 3 levels
System setting - 0,6
Scheduled ...
++++ Time to find Low-rank basis ++++
(Avg,Max) = 30.625,30.6883
Initialised
Initialised
++++ Time to Communicate among process ++++
(Avg,Max) = 0.0303617,0.0601609
++++ Time to generate entries for HODLR3D ++++
(Avg,Max) = 4.63271,4.65746
++++ Time to matrix-vector product ++++
(Avg,Max) = 0.278291,0.283
Error in sol..5.99558e-07
```