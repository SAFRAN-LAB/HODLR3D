********
Tutorial
********

For the sake of this tutorial, we are going to be using the ``testHODLR3D.cpp`` file that is listed under ``examples/`` since it demonstrates the features of this library. For the most part, comments in the file demonstrate intended functionality. However, we go over the main functions that may be of interest to a user on this page.

**NOTE**: It is assumed that you have already completed the installation process of getting the dependencies.

Setting Parameters in CMakeLists.txt
------------------------------------

There are some variables that need to be set by the user at the top of the ``CMakeLists.txt`` file:

- ``INPUT_FILE``: This is the input ``.cpp`` file that needs to be compiled. For this tutorial, it's going to be set to ``examples/testHODLR3D.cpp``.
- ``OUTPUT_EXECUTABLE``: This is the name that the final build executable is given. Here we are just going to set is as ``testHODLR3D``.

Running the program::
---------------------

For this particular tutorial, the problem parameters are passed to the executable during runtime. We have the lines::

  cubeRootN		            =	atoi(argv[1]);
	nParticlesInLeafAlong1D	=	atoi(argv[2]);
	L		                   	=	atof(argv[3]);
	TOL_POW                 = atoi(argv[4]);
	Qchoice                 = atoi(argv[5]);

Each of the above are described below

1. cubeRootN = pow(N,1/3) where N is the system size. In this tutorial, the code takes N to be perfect cubes

2. nParticlesInLeafAlong1D is a parameter that controls the maximum leaf size, where leaf size=pow(nParticlesInLeafAlong1D,3)

3. L is the half side length of square cube, which is the computational domain

4. TOL_POW = -log(tolerance for the ACA routine), which is approximately the number of digits of accuracy we want.

5. Qchoice is the choice of kernel function to be used from among those defined in function getMatrixEntry of userkernel class.

For instance, running ``./testHODLR3D 40 8 1.0 6 7`` would correspond to solving the problem with parameters :cubeRootN=40, nParticlesInLeafAlong1D=8, L=1.0, TOL_POW=6, Qchoice=7.

Creating a Derived Class of ``kernel``:
---------------------------------------

The matrix that needs to be solved for is abstracted through this derived class of ``kernel``. The main method that needs to be set for this class is ``getMatrixEntry`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix. Many kernel functions are defined, choose the that is needed using the Qchoice parameter at run-time. If a kernel function that is not already defined, one can define it at this function ``getMatrixEntry``. The ``getMatrixEntry`` is defined as follows in the ``testHODLR3D.cpp`` file::

  double userkernel::getMatrixEntry(const unsigned i, const unsigned j) {
	 if (i==j) {
	 	return 0.0;
	 }
	 else {
		if (Qchoice == 0)
			return RBF_Logarithm(i, j);
		else if (Qchoice == 1)
			return RBF_Exponential(i, j);
		else if (Qchoice == 2)
			return RBF_Inverse_Quadric(i, j);
		else if (Qchoice == 3)
			return RBF_Quadric(i, j);
		else if (Qchoice == 4)
			return RBF_Inverse_Multi_Quadric(i, j);
		else if (Qchoice == 5)
			return RBF_Gaussian(i, j);
		else if (Qchoice == 6)
			return RBF_Multi_quadric(i, j);
		else if (Qchoice == 7)
			return Laplacian_3D(i, j);
		else if (Qchoice == 8)
			return oneOverR4(i, j);
		else if (Qchoice == 9)
			return Laplacian_2D(i, j);
		else if (Qchoice == 10)
			return RBF_spline(i, j);
		else if (Qchoice == 11)
			return kernel_besselJ(i, j);
		else if (Qchoice == 12)
			return kernel_besselY(i, j);
		else if (Qchoice == 13)
			return Helmholtz_cos(i, j);
		else if (Qchoice == 14)
			return Feynman(i, j);
		else if (Qchoice == 15)
			return Yukawa(i, j);
		}
   }

Defining location of particles:
-------------------------------

In this tutorial, we have initialized ``particles``, at the tensor product grid of uniformly distributed :math:`pow(cubeRootN,3)` number of particles in :math:`(-L, L)^3` which are then passed to the constructor of object ``HODLR3D``::

  void set_Uniform_Nodes(int cubeRootN, double L, std::vector<pts3D>& particles) {
    std::vector<double> Nodes1D;
    for (int k=0; k<cubeRootN; ++k) {
      Nodes1D.push_back(-L+2.0*L*(k+1.0)/(cubeRootN+1.0));
    }
    pts3D temp1;
    for (int j=0; j<cubeRootN; ++j) {
      for (int k=0; k<cubeRootN; ++k) {
        for (int i=0; i<cubeRootN; ++i) {
          temp1.x	=	Nodes1D[k];
          temp1.y	=	Nodes1D[j];
          temp1.z	=	Nodes1D[i];
          particles.push_back(temp1);
        }
      }
    }
  }

One can also define particles at a distribution of choice and pass that to the constructor of object ``HODLR3D``.

Defining vector ``b``:
----------------------

In this tutorial, we have defined ``b``, the vector that is to be multiplied to the matrix, as a vector consisting of ones and zeros at random locations. This choice of b is considered to test the code in a fast way::

  Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
  int n = N/500; //randomly choosing n different indices where b is set to 1, b at the rest of the indices is set to 0
  srand(time(NULL));
  std::set<int> s;
  while(s.size() < n) {
    int index	=	rand()%N;
    s.insert(index);
  }
  std::set<int>::iterator it;
  for (it = s.begin(); it != s.end(); it++) {
    b(*it) = 1.0;
  }

Creating the Instance of ``HODLR3D``:
-------------------------------------

The main operations of this library are carried out through the ``HODLR3D`` class. The parameters that are taken for the constructor are N, MinParticlesInLeaf, TOL_POW, loc::

  HODLR3D *H = new HODLR3D(nParticlesInLeafAlong1D, L, TOL_POW, Qchoice, particles);


We will now proceed to demonstrate the individual methods available under this class.

``assemble``
^^^^^^^^^^^^

Assemble the matrix in HODRL3D structure; i.e. it finds the low rank representation of the low-rank matrix sub-blocks::

  H->assemble();

``MatVecProduct``
^^^^^^^^^^^^^^^^^^^

Multiplies the matrix that is defined through object ``userkernel`` with the vector ``b``::

  H->MatVecProduct(b, AFMM_Ab);
