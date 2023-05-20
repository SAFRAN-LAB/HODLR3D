***************
Reproducibility
***************

This part of the documentation helps in reproducing the results shown in the article titled **HODLR3D: Hierarchical matrices for**  :math:`N`-**body problems in three dimensions**, authored by **V A Kandappan, Vaishnavi Gujjula, Sivaram Ambikasaran**.
The code is available as an open source library and can be found `here <https://github.com/SAFRAN-LAB/HODLR3D>`_.

Numerical rank for different kernels
------------------------------------

To obtain Figure 5 (numerical rank vs N) and Figure 7 (Plot of singular values :math:`\sigma_{i}` normalized with the first singular value versus index :math:`i`) of the article, the following steps are to be followed

#. Run the file ``Matlab_files/get_singular_values.m``

   * The input to the file is 'choice' that decides what kernel is to be used. Enter 1 for :math:`1/r`, 2 for :math:`1/r^4`, 3 for :math:`cos(r)/r`.

   * The outputs of the file are svd_f (singular values of the face sharing interaction), svd_e (edge sharing interaction), svd_v (vertex sharing interaction), svd_w (well-separated interaction), and N (the matrix size).

   * These outputs are written to file "output_file_%d_%d.mat". The first argument is 'N' (the size of the matrix) and the second argument is 'choice'.

   * For example, if 'choice' 1 is inputed, it means that the singular values of various off-diagonal interactions will be computed for the kernel :math:`1/r`.

#. Run the file ``Matlab_files/getRanks.m``

   * For the singular values loaded from a given file, the code outputs the numerical rank for a given tolerance of 'tolr'.

   * For example, in line 12 of the code, if 'data = load('output_file_125_1.mat');', it means it reads the singular values of various off-diagonal interactions with the size of the matrix set to 125 and the kernel set to choice 1.


Numerical benchmarks of HODLR3D matrix-vector product in comparison with those of HODLR and :math:`\mathcal{H}` matrix-vector products for the kernel :math:`\frac{1}{r}`
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
To reproduce the results illustrated in Figure 9 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, :math:`\mathcal{H}`, and HODLR.

1. cubeRootN = vary between 20 and 150 with a step size of 10
2. nParticlesInLeafAlong1D = 6
3. L = 1.0
4. TOL_POW = 7
5. Qchoice = 7

HODLR3D
^^^^^^^

Key in the file ```examples/testHODLR3D.cpp`` as input under ``INPUT_FILE`` in ``HODLR3Dlib/CMakeLists.txt``. Here you also set the name of the output executable, say ``testHODLR3D``, under ``OUTPUT_EXECUTABLE_NAME``.
Compile and build the executable as described in :ref:`Testing`.

For example, run the following command::

   ./testHODLR3D 20 6 1.0 7 7

:math:`\mathcal{H}`
^^^^^^^^^^^^^^^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3D.mk clean

Then make the file::

	 make -f Makefile3D.mk

Run the generated executable as, for instance,::

   ./testH 20 6 1.0 7 7

HODLR
^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3D.mk clean

Then make the file::

	 make -f Makefile3D.mk

Run the generated executable as, for instance,::

   ./testHODLR 20 6 1.0 7 7


Numerical benchmarks of HODLR3D matrix-vector product in comparison with those of HODLR and :math:`\mathcal{H}` matrix-vector products for the kernel :math:`\frac{1}{r^4}`
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the results illustrated in Figure 10 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, :math:`\mathcal{H}`, and HODLR.

1. cubeRootN = vary between 20 and 150 with a step size of 10
2. nParticlesInLeafAlong1D = 6
3. L = 1.0
4. TOL_POW = 7
5. Qchoice = 8

HODLR3D
^^^^^^^

Key in the file ```examples/testHODLR3D.cpp`` as input under ``INPUT_FILE`` in ``HODLR3Dlib/CMakeLists.txt``. Here you also set the name of the output executable, say ``testHODLR3D``, under ``OUTPUT_EXECUTABLE_NAME``.
Compile and build the executable as described in :ref:`Testing`.

For example, run the following command::

   ./testHODLR3D 20 6 1.0 7 8

:math:`\mathcal{H}`
^^^^^^^^^^^^^^^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3D.mk clean

Then make the file::

	 make -f Makefile3D.mk

Run the generated executable as, for instance,::

 ./testH 20 6 1.0 7 8

HODLR
^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3D.mk clean

Then make the file::

	 make -f Makefile3D.mk

Run the generated executable as, for instance,::

 ./testHODLR 20 6 1.0 7 8


Numerical benchmarks of HODLR3D matrix-vector product in comparison with those of HODLR and :math:`\mathcal{H}` matrix-vector products for the kernel :math:`\frac{cos(r)}{r}`
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the results illustrated in Figure 11 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, :math:`\mathcal{H}`, and HODLR.

1. cubeRootN = vary between 20 and 150 with a step size of 10
2. nParticlesInLeafAlong1D = 6
3. L = 1.0
4. TOL_POW = 7
5. Qchoice = 13

HODLR3D
^^^^^^^

Key in the file ```examples/testHODLR3D.cpp`` as input under ``INPUT_FILE`` in ``HODLR3Dlib/CMakeLists.txt``. Here you also set the name of the output executable, say ``testHODLR3D``, under ``OUTPUT_EXECUTABLE_NAME``.
Compile and build the executable as described in :ref:`Testing`.

For example, run the following command::

   ./testHODLR3D 20 6 1.0 7 13

:math:`\mathcal{H}`
^^^^^^^^^^^^^^^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3D.mk clean

Then make the file::

	 make -f Makefile3D.mk

Run the generated executable as, for instance,::

 ./testH 20 6 1.0 7 13

HODLR
^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3D.mk clean

Then make the file::

	 make -f Makefile3D.mk

Run the generated executable as, for instance,::

 ./testHODLR 20 6 1.0 7 13


Numerical benchmarks of the HODLR3D accelerated iterative solver for the integral equation in comparison with those of HODLR and :math:`\mathcal{H}`
----------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the results illustrated in Figure 12 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, :math:`\mathcal{H}`, and HODLR.

1. cubeRootN = vary between 20 and 150 with a step size of 10
2. nParticlesInLeafAlong1D = 6
3. L = 1.0
4. TOL_POW = 7
5. Qchoice = 16

HODLR3D
^^^^^^^

Key in the file ```examples/testHODLR3Dsolve.cpp`` as input under ``INPUT_FILE`` in ``HODLR3Dlib/CMakeLists.txt``. Here you also set the name of the output executable, say ``testHODLR3Dsolve``, under ``OUTPUT_EXECUTABLE_NAME``.
Compile and build the executable as described in :ref:`Testing`.

For example, run the following command::

   ./testHODLR3Dsolve 20 6 1.0 7 16


:math:`\mathcal{H}`
^^^^^^^^^^^^^^^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3Dsolve.mk clean

Then make the file::

	 make -f Makefile3Dsolve.mk

Run the generated executable as, for instance,::

./testH 20 6 1.0 7 16

HODLR
^^^^^

Clean using make clean before running the code, i.e.,::

	 make -f Makefile3Dsolve.mk clean

Then make the file::

	 make -f Makefile3Dsolve.mk

Run the generated executable as, for instance,::

./testHODLR 20 6 1.0 7 16


Numerical benchmarks of parallel HODLR3D matrix-vector product using MPI
------------------------------------------------------------------------
To reproduce the results illustrated in Table 4 of the article, follow the instructions given here.

1. The Eigen library can be downloaded from its `website <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_.
2. The code is tested with OpenMPI 4.1.1, which can be downloaded from its `website <https://www.open-mpi.org/software/ompi/v4.1/>`_.
3. Set the following variables in the "CMakeLists.txt" file
    - CMAKE_C_COMPILER: GCC version greater than GCC9
    - CMAKE_CXX_COMPILER: GCC version greater than GCC9
    - EIGEN_PATH: Set the path for the Eigen library.
    - HOME_PATH: Provide the path where the .cpp file is located.
    - MPI_INC: Provide the path to the Open MPI ‘include’ directory - path/to/open-mpi/4.1.x_x/include
    - MPI_LIB: Provide the path to the Open MPI ‘library’ directory - path/to/open-mpi/4.1.x_x/lib

Note: The code has been tested with other MPI wrapper compilers as well. The HODLR3D code has been tested with mpicxx, Intel-based mpi wrapper compilers. The scaling does not affect due to changes in the compiler.

Sample CMakeLists.txt::


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

Sample Installation and Building::


  user@computer HODLR3D$ mkdir build && cd build
  user@computer build$ cmake ..
  user@computer build$ make
  user@computer build$ mpiexec.hydra -np 2 -genv I_MPI_PIN=1 -genv I_MPI_FABRICS=shm:ofi -hostfile $PBS_NODEFILE ./hodlr3d 50 10 1 6 1 > H3D_2_50.txt


Sample output file "H3D_2_50.txt"::

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
