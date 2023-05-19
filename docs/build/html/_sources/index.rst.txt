.. role:: underline
    :class: underline

Welcome to HODLR3Dlib's Documentation!
**************************************

About :math:`\texttt{HODLR3Dlib}`:
==================================

HODLR3D (Hierarchical Off Diagonal Low Rank 3D) is a new class of Hierarchical matrix representation for matrices arising out of `N`-body problems in 3D.

:math:`\texttt{HODLR3Dlib}` is a library to compute fast matrix-vector products using HODLR3D.

It relies on the fact that certain off-diagonal blocks in the matrix can be well-approximated by low-rank matrices.

The domain has been hierarchically partitioned using a uniform oct-tree.

The low-rank approximation of the appropriate blocks is obtained using ACA.

The code is written in C++ and features an easy-to-use interface, where the user provides the following inputs:

- a ``kernel`` object which abstracts data of the matrix through a member function ``getMatrixEntry(int i, int j)`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

- locations of nodes in the domain through a std::vector ``particles``

- the vector ``b`` to be multiplied to the matrix

It obtains :math:`A x` in almost linear complexity.

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   benchmarks
   reproducibility

Other Links
===========

Learn more about :math:`\texttt{HODLR3Dlib}` by visiting the

* Code Repository: https://github.com/SAFRAN-LAB/HODLR3D
* Documentation: https://hodlr3d.readthedocs.io/en/latest/
