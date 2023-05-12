*************************
Installation and Building
*************************

Downloading the Source
-----------------------

:math:`\texttt{HODLR3Dlib}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/SAFRAN-LAB/HODLR3D.git --recursive

The ``--recursive`` flag is argument to ensure that all submodules are also checked out.

Dependencies
-------------

- Eigen Linear Algebra Library (get it `here <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.

**NOTE**: On MacOS, the default compiler is `clang` which doesn't have OpenMP support. You will have to use g++ to make use of the speedups from OpenMP::

    user@computer HODLR3D$ brew install g++-10
    user@computer HODLR3D$ export CXX=g++

Installation
-------------

Manually Installing
^^^^^^^^^^^^^^^^^^^

Then, set the environment variable ``EIGEN_PATH`` to the location of your Eigen installation. This is needed by the CMake script.::

    user@computer HODLR3D$ export EIGEN_PATH=path/to/eigen/

Testing
-------

Now, we need to ensure that all the functions of the libraries function as intended. For this purpose, we will be running the script ``examples/testHODLR3D.cpp``.
Key in the file ```examples/testHODLR3D.cpp`` as input under ``INPUT_FILE`` in ``HODLR3Dlib/CMakeLists.txt``. Here you also set the name of the output executable, say ``testHODLR3D``, under ``OUTPUT_EXECUTABLE_NAME``.
Then create a directory called ``build`` and navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable.
To check this on your computer, run the following lines::

    user@computer HODLR3D$ mkdir build && cd build
    user@computer build$ cmake ..
    user@computer build$ make
    user@computer build$ ./testHODLR3D

Building and Executing
----------------------

Key in the required ``.cpp`` to be used as input under ``INPUT_FILE`` in ``HODLR3Dlib/CMakeLists.txt``. Here you also set the name of the output executable under ``OUTPUT_EXECUTABLE_NAME``. Then navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable::

    user@computer build$ cmake path/to/HODLR3D/
    user@computer build$ make -j n_threads
    user@computer build$ ./executable
