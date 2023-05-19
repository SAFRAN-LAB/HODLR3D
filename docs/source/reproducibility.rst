***************
Reproducibility
***************

This part of the documentation helps in reproducing the results shown in the article titled "HODLR3D: Hierarchical matrices for :math:`N`-body problems in three dimensions".

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


Various benchmarks of HODLR3D matrix-vector product in comparison with those of HODLR and :math:`\mathcal{H}` matrix-vector products for the kernel :math:`\frac{1}{r}`
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
To reproduce the results illustrated in Figure 9 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, H, and HODLR.

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

H
^

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


Various benchmarks of HODLR3D matrix-vector product in comparison with those of HODLR and :math:`\mathcal{H}` matrix-vector products for the kernel :math:`\frac{1}{r^4}`
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the results illustrated in Figure 10 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, H, and HODLR.

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

H
^

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


Various benchmarks of HODLR3D matrix-vector product in comparison with those of HODLR and :math:`\mathcal{H}` matrix-vector products for the kernel :math:`\frac{cos(r)}{r}`
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the results illustrated in Figure 11 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, H, and HODLR.

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

H
^

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


Various benchmarks of the HODLR3D accelerated iterative solver for the integral equation in comparison with those of HODLR and :math:`\mathcal{H}`
--------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the results illustrated in Figure 12 of the article, follow the instructions given here.

The following values are inputed at run-time for the three codes HODLR3D, H, and HODLR.

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

H
^

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
