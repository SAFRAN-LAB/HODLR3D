The folder HODLR3D/H contains the files related to the H matrix that is considered in the article. The purpose of the code in this folder is to make a comparative study of HODLR3D with HODLR matrix. The HODLR matrix is constructed by low-rank approximating all the off-diagonal blocks in a non-nested manner.

"testHODLR.cpp" performs HODLR matrix vector product. The matrix entries are to be defined in the function "getMatrixEntry(i,j)" of the kernel.hpp file.

The vector to be multiplied to the matrix is to be defined as vector "b" of "testHODLR.cpp" file. Currently for test purposes, it is defined as a vector containing ones and zeros at random locations.

The particles are currently distributed uniformly. To change the location of particles, define the vector "particles" of object "mykernel" of file "testHODLR.cpp" accordingly.

It takes these inputs at run time:

  1. cubeRootN: It determines the system size, that is N=pow(cubeRootN,3)

  2. nParticlesInLeafAlong1D: it determines the maximum number of particles in a leaf node, which is pow(nParticlesInLeafAlong1D,3).

  3. L: half side length of the cube which is the computational domain

  4. TOL_POW: tolerance set for ACA routine

  5. Qchoice: Use this to select the kernel you want to use. For various choices look at the kernel.hpp file.

  For instance, below is a sample execution of the code and its output.

  It is always good to clean using make clean before running the code, i.e.,

  	make -f Makefile3D.mk clean

  Then make the file

  	make -f Makefile3D.mk

  Run the generated executable as, for instance,

   ./testHODLR 20 6 1 7 7

   Number of particles is: 8000

   Time taken to create the tree is: 0.003958

   Time taken to construct HODLR3D representation is: 0.831133

   Time taken to do Mat-Vec product is: 0.007521

   Memory in GB: 0.00284596

   CR: 0.355745

   max rank: 46

   relative forward error: 4.96199e-08

------------------------------------------------------------------------------------------------------------------------

"testHODLRsolve.cpp" solves an integral equation. The integral equation is discretized and the resulting linear system is solved using GMRES. The matrix-vector product encountered in each of GMRES' iterations is accelerated using the HODLR matrix-vector product.

For test purposes, the true solution is considered to be a vector of zeros and ones at random locations. The true right-hand side is computed exactly (upto roundoff) and the solution is found using GMRES.

The 'Qchoice' is to be set to 16.

It takes these inputs at run time:

  1. cubeRootN: It determines the system size, that is N=pow(cubeRootN,3)

  2. nParticlesInLeafAlong1D: it determines the maximum number of particles in a leaf node, which is pow(nParticlesInLeafAlong1D,3).

  3. L: half side length of the cube which is the computational domain

  4. TOL_POW: tolerance set for ACA routine

  5. Qchoice: 16.

  For instance, below is a sample execution of the code and its output.

  It is always good to clean using make clean before running the code, i.e.,

  	make -f Makefile3D.mk clean

  Then make the file

  	make -f Makefile3D.mk

  Run the generated executable as, for instance,

  ./testHODLRsolve 20 6 1 7 16

  Number of particles is: 8000

  Time taken to create the tree is: 0.003719

  time GMRES: 1.20493

  GMRES residual err: 1.50231e-11

  GMRES no. of iterations: 8

  relative forward error in solution: 1.69799e-09
