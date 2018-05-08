# Awesome Numerical Software
[![Awesome](https://awesome.re/badge-flat.svg)](https://github.com/sindresorhus/awesome)

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
   (Fortran, public domain) -
   Standard building blocks for performing basic vector and matrix operations.

 - [OpenBLAS](https://www.openblas.net/)
   (Fortran, BSD, [GitHub](https://github.com/xianyi/OpenBLAS)) -
   Optimized BLAS library based on GotoBLAS2.

 - BLIS
   (C++, BSD, [GitHub](https://github.com/flame/blis)) -
   Portable software framework for instantiating high-performance BLAS-like
   dense linear algebra libraries.

 - [LAPACK](http://www.netlib.org/lapack/)
   (Fortran, BSD, [GitHub](https://github.com/Reference-LAPACK/lapack)) -
   Routines for solving systems of simultaneous linear equations, least-squares
   solutions of linear systems of equations, eigenvalue problems, and singular
   value problems.

 - [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
   (C++, MPL2, [BitBucket](https://bitbucket.org/eigen/eigen)) -
   C++ template library for linear algebra: matrices, vectors, numerical
   solvers, and related algorithms.


## Frameworks

 - [PETSc](http://www.mcs.anl.gov/petsc/)
   (C, 2-clause BSD license, [BitBucket](https://bitbucket.org/petsc/petsc/src)) -
   Suite of data structures and routines for the scalable (parallel) solution
   of scientific applications modeled by partial differential equations.

 - [DUNE Numerics](https://www.dune-project.org/)
   (C++, GPL2, [GitLab](https://gitlab.dune-project.org/core/)) -
   Modular toolbox for solving partial differential equations (PDEs) with
   grid-based methods.

 - Trilinos
   (mostly C++, mostly BSD, [GitHub](https://github.com/trilinos/)) -
   Algorithms and enabling technologies within an object-oriented software
   framework for the solution of large-scale, complex multi-physics engineering
   and scientific problems.

 - [SciPy](https://www.scipy.org/)
   (Python, mostly BSD, [GitHub](https://github.com/scipy/scipy/)) -
   Python modules for statistics, optimization, integration, linear algebra,
   Fourier transforms, signal and image processing, ODE solvers, and more.

 - [NumPy](http://www.numpy.org/)
   (Python, BSD, [GitHub](https://github.com/numpy/numpy)) -
   Fundamental package needed for scientific computing with Python.


## Finite Elements

 - [FEniCS](https://fenicsproject.org/)
   (C++/Python, LGPLv3, [BitBucket](https://bitbucket.org/fenics-project/)) -
   Open-source computing platform for solving partial differential equations
   (PDEs).

 - [libMesh](https://libmesh.github.io/)
   (C++, LGPL 2.1, [GitHub](https://github.com/libMesh/libmesh)) -
   Framework for the numerical simulation of partial differential equations
   using arbitrary unstructured discretizations on serial and parallel
   platforms.

 - [deal.II](http://dealii.org/)
   (C++, LGPL 2.1, [GitHub](https://github.com/dealii/dealii)) -
   Software library supporting the creation of finite element codes and an open
   community of users and developers.


## Meshing

 - [gmsh](http://gmsh.info/)
   (C++, GPL, [GitLab](https://gitlab.onelab.info/gmsh/gmsh)) -
   Three-dimensional finite element mesh generator with built-in pre- and
   post-processing facilities.
   [A separate Python interface exists.](https://github.com/nschloe/pygmsh)

 - [MeshPy](https://mathema.tician.de/software/meshpy/)
   (Python, MIT, [GitHub](https://github.com/inducer/meshpy)) -
   Quality triangular and tetrahedral mesh generation.

 - [Netgen/NGSolve](https://ngsolve.org/)
   (C++, LGPL 2.1, [GitHub](https://github.com/NGSolve/netgen)) -
   High performance multiphysics finite element software.

 - meshio
   (Python, MIT, [GitHub](https://github.com/nschloe/meshio)) -
   I/O for various mesh formats, file conversion.

 - [CGAL](https://www.cgal.org/)
   (C++, mixed LGPL/GPL, [GitHub](https://github.com/CGAL/cgal)) -
   Efficient and reliable algorithms in computational geometry.
   Separate Python interfaces for meshing exist
   ([pygalmesh](https://github.com/nschloe/pygalmesh),
   [mshr](https://bitbucket.org/fenics-project/mshr/)).

 - [MOAB](http://sigma.mcs.anl.gov/moab-library/)
   (C++, mostly LGPL3, [BitBucket](https://bitbucket.org/fathomteam/moab/)) -
   Representing and evaluating mesh data.

 - [NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
   (C/C++/Fortran/Java/Python, [custom open source](https://www.unidata.ucar.edu/software/netcdf/copyright.html), [GitHub](https://github.com/Unidata/netcdf-c/)) -
   Software libraries and self-describing, machine-independent data formats
   that support the creation, access, and sharing of array-oriented scientific
   data.

 - [HDF5](https://support.hdfgroup.org/HDF5/)
   (C/Fortran, BSD) -
   Data model, library, and file format for storing and managing data.

 - [XDMF](http://www.xdmf.org/index.php/Main_Page)
   (C++, [GitLab](https://gitlab.kitware.com/xdmf/xdmf)) -
   eXtensible Data Model and Format to exchange scientific data between High
   Performance Computing codes and tools.


## Sparse linear solvers

 - [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)
   (C, mostly BSD, [GitHub](https://github.com/xiaoyeli/superlu)) -
   General purpose library for the direct solution of large, sparse,
   nonsymmetric systems of linear equations.

 - [MUMPS](http://mumps.enseeiht.fr/)
   (Fortran, CeCILL-C) -
   Solving systems of linear equations of the form Ax = b, where A is a square
   sparse matrix that can be either unsymmetric, symmetric positive definite,
   or general symmetric, on distributed memory computers.

 - KryPy
   (Python, MIT, [GitHub](https://github.com/andrenarchy/krypy)) -
   Krylov subspace methods for the solution of linear algebraic systems.


## Miscellaneous

 - [FFTW](http://www.fftw.org/)
   (C, GPL2, [GitHub](https://github.com/FFTW/fftw3)) -
   Computes the discrete Fourier transform (DFT) in one or more dimensions, of
   arbitrary input size, and of both real and complex data (as well as of
   even/odd data, i.e. the discrete cosine/sine transforms or DCT/DST).

 - [Qhull](http://www.qhull.org/)
   (C/C++, [custom open source license](http://www.qhull.org/COPYING.txt), [GitHub](https://github.com/qhull/qhull/)) -
   Computes the convex hull, Delaunay triangulation, Voronoi diagram, halfspace
   intersection about a point, furthest-site Delaunay triangulation, and
   furthest-site Voronoi diagram.

 - [GSL](https://www.gnu.org/software/gsl/)
   (C/C++, GPL, [Savannah](https://savannah.gnu.org/projects/gsl)) -
   Wide range of mathematical routines such as random number generators,
   special functions and least-squares fitting. There are over 1000 functions
   in total.

 - [OpenFOAM](https://www.openfoam.com/)
   (C++, GPL3, [GitHub](https://github.com/OpenFOAM/OpenFOAM-dev)) -
   Free, open source CFD (computational fluid dynamics) software.

 - [ParaView](https://www.paraview.org/)
   (C++, BSD, [GitLab](https://gitlab.kitware.com/paraview/paraview)) -
   Open-source, multi-platform data analysis and visualization application
   based on Visualization Toolkit (VTK).

 - [PyAMG](https://pyamg.github.io/)
   (Python, MIT, [GitHub](https://github.com/pyamg/pyamg)) -
   Algebraic Multigrid Solvers in Python.

 - quadpy
   (Pyhon, MIT, [GitHub](https://github.com/nschloe/quadpy)) -
   Numerical integration (quadrature, cubature) in Python.


## License

[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

To the extent possible under law, [Nico Schlömer](https://github.com/nschloe)
has waived all copyright and related or neighboring rights to this work.
