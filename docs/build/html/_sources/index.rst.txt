.. role:: underline
    :class: underline

Welcome to HODLR3Dlib's Documentation!
**************************************

About :math:`\texttt{HODLR3Dlib}`:
==================================

HODLR3D (Hierarchical Off Diagonal Low Rank 3D) is a new class of Hierarchical matrix representation for matrices arising out of `N`-body problems in 3D.

:math:`\texttt{HODLR3Dlib}` is a library to compute fast matrix-vector products using HODLR3D.

It relies on the fact that certain off-diagonal blocks in the matrix can be well-approximated by low-rank matrices.

Low-rank approximation of the appropriate blocks is obtained using ACA.

The domain has been hierachically partitioned using a uniform oct-tree.

The code is written in C++ and features an easy-to-use interface, where the user provides the following inputs:

- a ``kernel`` object which abstracts data of the matrix through a member function ``getMatrixEntry(int i, int j)`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

- locations of nodes in the domain through an Eigen matrix ``loc``

- the vector ``b`` to be multiplied to the matrix

The algorithm has been parallelized using OpenMP.

Obtains :math:`A x` at a cost of :math:`\mathcal{O}\left(N\right)`

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   benchmarks

Other Links
===========

Learn more about :math:`\texttt{AFMM2Dlib}` by visiting the

* Code Repository: http://github.com/sivaramambikasaran/AFMM2D
* Documentation: http://afmm2d.rtfd.io
