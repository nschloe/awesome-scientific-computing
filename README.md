# Awesome Numerical Software

A curated list of awesome software for numerical analysis.


## Contents

- [Basic linear algebra](#basic-linear-algebra)
- [Frameworks](#frameworks)
- [Finite Elements](#finite-elements)
- [Meshing](#meshing)
- [Sparse linear solvers](#sparse-linear-solvers)
- [Miscellaneous](#miscellaneous)


## Basic linear algebra

 - [BLAS](http://www.netlib.org/blas/)
   (Fortran, public domain).
   The BLAS (Basic Linear Algebra Subprograms) are routines that provide
   standard building blocks for performing basic vector and matrix operations.

 - [OpenBLAS](https://www.openblas.net/)
   (Fortran, BSD, [GitHub](https://github.com/xianyi/OpenBLAS)).
   Optimized BLAS library based on GotoBLAS2.

 - [BLIS](https://github.com/flame/blis)
   (C++, BSD, [GitHub](https://github.com/flame/blis)).
   Portable software framework for instantiating high-performance BLAS-like
   dense linear algebra libraries.

 - [LAPACK](http://www.netlib.org/lapack/)
   (Fortran, BSD, [GitHub](https://github.com/Reference-LAPACK/lapack)).
   Routines for solving systems of simultaneous linear equations, least-squares
   solutions of linear systems of equations, eigenvalue problems, and singular
   value problems.

 - [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
   (C++, MPL2, [BitBucket](https://bitbucket.org/eigen/eigen) (Mercurial)).
   C++ template library for linear algebra: matrices, vectors,
   numerical solvers, and related algorithms.


## Frameworks

 - [PETSc](https://www.mcs.anl.gov/petsc/)
   (C, 2-clause BSD license, [BitBucket](https://bitbucket.org/petsc/petsc/src)).
   Suite of data structures and routines for the scalable (parallel) solution
   of scientific applications modeled by partial differential equations. It
   supports MPI, and GPUs through CUDA or OpenCL, as well as hybrid MPI-GPU
   parallelism.

 - [DUNE Numerics](https://www.dune-project.org/)
   (C++, GPL2, [GitLab](https://gitlab.dune-project.org/core/)).
   DUNE, the Distributed and Unified Numerics Environment is a modular toolbox
   for solving partial differential equations (PDEs) with grid-based methods.
   It supports the easy implementation of methods like Finite Elements (FE),
   Finite Volumes (FV), and also Finite Differences (FD).

 - [Trilinos](https://trilinos.org/)
   (mostly C++, mostly BSD, [GitHub](https://github.com/trilinos/)).
   Algorithms and enabling technologies within an object-oriented software
   framework for the solution of large-scale, complex multi-physics engineering
   and scientific problems. A unique design feature of Trilinos is its focus on
   packages.

 - [SciPy](https://www.scipy.org/)
   (Python, mostly BSD, [GitHub](https://github.com/scipy/scipy/)).
   SciPy (pronounced "Sigh Pie") is open-source software for mathematics,
   science, and engineering. It includes modules for statistics, optimization,
   integration, linear algebra, Fourier transforms, signal and image
   processing, ODE solvers, and more.

 - [NumPy](http://www.numpy.org/)
   (Python, BSD, [GitHub](https://github.com/numpy/numpy)).
   NumPy is the fundamental package needed for scientific computing with Python.


## Finite Elements

 - [FEniCS](https://fenicsproject.org/)
   (C++/Python, LGPLv3, [BitBucket](https://bitbucket.org/fenics-project/)).
   Open-source computing platform for solving partial differential equations
   (PDEs). FEniCS enables users to quickly translate scientific models into
   efficient finite element code.

 - [libMesh](https://libmesh.github.io/)
   (C++, LGPL 2.1, [GitHub](https://github.com/libMesh/libmesh)).
   A framework for the numerical simulation of partial differential equations
   using arbitrary unstructured discretizations on serial and parallel
   platforms. A major goal of the library is to provide support for adaptive
   mesh refinement (AMR) computations in parallel.

 - [deal.II](http://dealii.org/)
   (C++, LGPL 2.1, [GitHub](https://github.com/dealii/dealii)).
   A C++ software library supporting the creation of finite element codes and
   an open community of users and developers.

 - [FreeFem++](http://www.freefem.org/)
   (C++, LGPL 2.1, [GitHub](https://github.com/FreeFem/FreeFem-sources)).
   Partial differential equation solver. It has its own language. freefem
   scripts can solve multiphysics non linear systems in 2D and 3D.


## Meshing

 - [gmsh](http://gmsh.info/)
   (C++, GPL, [GitLab](https://gitlab.onelab.info/gmsh/gmsh)).
   Three-dimensional finite element mesh generator with built-in pre- and
   post-processing facilities.

   [A separate Python interface exists.](https://github.com/nschloe/pygmsh)

 - [MeshPy](https://mathema.tician.de/software/meshpy/)
   (Python, MIT, [GitHub](https://github.com/inducer/meshpy)).
   Quality triangular and tetrahedral mesh generation.

 - [Netgen](https://github.com/NGSolve/netgen)
   (C++ LGPL 2.1, [GitHub](https://github.com/NGSolve/netgen)).
   Automatic 3D tetrahedral mesh generator.

 - [meshio](https://github.com/nschloe/meshio)
   (Python, MIT, [GitHub](https://github.com/nschloe/meshio)).
   I/O for various mesh formats, file conversion.

 - [CGAL](https://www.cgal.org/)
   (C++, mixed LGPL/GPL, [GitHub](https://github.com/CGAL/cgal)).
   The Computational Geometry Algorithms Library (CGAL) is a C++ library that
   aims to provide easy access to efficient and reliable algorithms in
   computational geometry.

   Separate Python interfaces exist
   ([pygalmesh](https://github.com/nschloe/pygalmesh),
   [mshr](https://bitbucket.org/fenics-project/mshr/)).

 - [MOAB](http://sigma.mcs.anl.gov/moab-library/)
   (C++, mostly LGPL3, [BitBucket](https://bitbucket.org/fathomteam/moab/)).
   A component for representing and evaluating mesh data.

 - [NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
   (C/C++/Fortran/Java/Python, [custom open source](https://www.unidata.ucar.edu/software/netcdf/copyright.html), [GitHub](https://github.com/Unidata/netcdf-c/)).
   NetCDF is a set of software libraries and self-describing,
   machine-independent data formats that support the creation, access, and
   sharing of array-oriented scientific data.

 - [HDF5](https://support.hdfgroup.org/HDF5/)
   (C/Fortran, BSD).
   A data model, library, and file format for storing and managing data. It
   supports an unlimited variety of data types, and is designed for flexible and
   efficient I/O and for high volume and complex data.

 - [XDMF](http://www.xdmf.org)
   (C++, [GitLab](https://gitlab.kitware.com/xdmf/xdmf)).
   eXtensible Data Model and Format to exchange scientific data between High
   Performance Computing codes and tools.


## Sparse linear solvers

 - [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)
   (C, mostly BSD, [GitHub](https://github.com/xiaoyeli/superlu)).
   General purpose library for the direct solution of large, sparse,
   nonsymmetric systems of linear equations.

 - [MUMPS](http://mumps.enseeiht.fr/)
   (Fortran, CeCILL-C).
   MUltifrontal Massively Parallel Solver. A package for solving systems of
   linear equations of the form Ax = b, where A is a square sparse matrix that
   can be either unsymmetric, symmetric positive definite, or general
   symmetric, on distributed memory computers.

 - [KryPy](https://github.com/andrenarchy/krypy)
   (Python, MIT, [GitHub](https://github.com/andrenarchy/krypy)).
   Module for Krylov subspace methods for the solution of linear algebraic
   systems. This includes enhanced versions of CG, MINRES and GMRES as well as
   methods for the efficient solution of sequences of linear systems.


## Miscellaneous

 - [FFTW](http://www.fftw.org/)
   (C, GPL2, [GitHub](https://github.com/FFTW/fftw3)).
   Library for computing the discrete Fourier transform (DFT) in one or more
   dimensions, of arbitrary input size, and of both real and complex data (as
   well as of even/odd data, i.e. the discrete cosine/sine transforms or
   DCT/DST).

 - [Qhull](http://www.qhull.org/)
   (C/C++, [custom open source license](http://www.qhull.org/COPYING.txt), [GitHub](https://github.com/qhull/qhull/)).
   Computes the convex hull, Delaunay triangulation, Voronoi diagram, halfspace
   intersection about a point, furthest-site Delaunay triangulation, and
   furthest-site Voronoi diagram.

 - [GNU Octave](https://www.gnu.org/software/octave/)
   (GPL).
   Powerful mathematics-oriented syntax with built-in plotting and
   visualization tools. Free and open-source MATLAB clone.

 - [GSL](https://www.gnu.org/software/gsl/)
   (C/C++, GPL, [Savannah](https://savannah.gnu.org/projects/gsl)).
   Provides a wide range of mathematical routines such as random number
   generators, special functions and least-squares fitting. There are over 1000
   functions in total.

 - [OpenFOAM](https://www.openfoam.com/)
   (C++, GPL3, [GitHub](https://github.com/OpenFOAM/OpenFOAM-dev)).
   Free, open source CFD (computational fluid dynamics) software.

 - [ParaView](https://www.paraview.org/)
   (C++, BSD, [GitLab](https://gitlab.kitware.com/paraview/paraview)).
   Open-source, multi-platform data analysis and visualization application
   based on Visualization Toolkit (VTK).

 - [PyAMG](https://pyamg.github.io/)
   (Python, MIT, [GitHub](https://github.com/pyamg/pyamg)).
   Algebraic Multigrid Solvers in Python.

 - [quadpy](https://github.com/nschloe/quadpy)
   (Pyhon, MIT, [GitHub](https://github.com/nschloe/quadpy)).
   Numerical integration (quadrature, cubature) in Python.


## License

[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

To the extent possible under law, [Nico Schlömer](https://github.com/nschloe)
has waived all copyright and related or neighboring rights to this work.
