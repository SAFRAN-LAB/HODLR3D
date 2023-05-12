********
Tutorial
********

For the sake of this tutorial, we are going to be using the ``tutorial.cpp`` file that is listed under ``examples/`` since it demonstrates the features of this library. For the most part, comments in the file demonstrate intended functionality. However, we go over the main functions that may be of interest to a user on this page.

**NOTE**: It is assumed that you have already completed the installation process of getting the dependencies.

Setting Parameters in CMakeLists.txt
------------------------------------

There are some variables that need to be set by the user at the top of the ``CMakeLists.txt`` file:

- ``INPUT_FILE``: This is the input ``.cpp`` file that needs to be compiled. For this tutorial, it's going to be set to ``examples/tutorial.cpp``.
- ``OUTPUT_EXECUTABLE``: This is the name that the final build executable is given. Here we are just going to set is as ``tutorial``.

Running the program::
---------------------

For this particular tutorial, the problem parameters are passed to the executable during runtime. We have the lines::

    N                  = atoi(argv[1]); // Size of the Matrix in consideration
    MinParticlesInLeaf = atoi(argv[2]); // Minimum particles that can be present in a leaf of KD tree
    TOL_POW            = atoi(argv[3]); // Tolerance of problem

This means that the first argument would be the matrix size considered, the second one would be the Minimum particles that can be present in a leaf of KD tree, and the final argument is approximately the number of digits of accuracy we want. For instance, running ``./tutorial 10000 32 12`` would correspond to solving the problem with parameters :math:`N=1000, MinParticlesInLeaf=32, \epsilon=10^{-12}`.

Creating a Derived Class of ``kernel``:
---------------------------------------

The matrix that needs to be solved for is abstracted through this derived class of ``kernel``. The main method that needs to be set for this class is ``getMatrixEntry`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix. For instance, for the :math:`1 / R` kernel, this would be set as::

  class userkernel: public kernel {

  public:

  getMatrixEntry(const unsigned i, const unsigned j) {

  	pts2D r1 = particles_X[i];

  	pts2D r2 = particles_X[j];

  	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);

  	double R	=	sqrt(R2);

  	if (R < a) {

  		return R/a;
  	}
  	else {
  		return a/R;
  	}
  }
  	~userkernel() {};
  };


In this tutorial, we have initialized ``loc``, a random set of points in :math:`(-1, 1)` which are then passed to the constructor of object ``AFMM``::

  Eigen::MatrixXd loc(N,Dimension);

  for (size_t j = 0; j < N; j++) {

    for (size_t k = 0; k < Dimension; k++) {

      loc(j,k) = 2.0*double(rand())/double(RAND_MAX)-1.0; // domain: [-1,1]x[-1,1]
    }
  }

Defining vector ``b``:
----------------------

In this tutorial, we have defined ``b``, the vector that is to be multiplied to the matrix, as a random set of points::

  Eigen::VectorXd b = Eigen::VectorXd::Random(N);

Creating the Instance of ``AFMM``:
----------------------------------

The main operations of this library are carried out through the ``AFMM`` class. The parameters that are taken for the constructor are N, MinParticlesInLeaf, TOL_POW, loc::

  AFMM* afmm = new AFMM(N, MinParticlesInLeaf, TOL_POW, loc);

We will now proceed to demonstrate the individual methods available under this class.

``assemble``
^^^^^^^^^^^^

Assemble the matrix in AFMM structure; i.e. it finds the low rank representation of the appropriate matrix sub-blocks::

  afmm->assemble();

``Mat-Vec product``
^^^^^^^^^^^^^^^^^^^

Multiplies the matrix that is defined through object ``userkernel`` with the vector ``b``::

  outputVec = afmm->computeMatVecProduct(b);
