# Awesome Numerical Software [![Awesome](https://awesome.re/badge-flat.svg)](https://github.com/sindresorhus/awesome)

A curated list of awesome software for numerical analysis.


## Contents

- [Basic linear algebra](#basic-linear-algebra)
- [Frameworks](#frameworks)
- [Finite Elements](#finite-elements)
- [Meshing](#meshing)
- [Sparse linear solvers](#sparse-linear-solvers)
- [Miscellaneous](#miscellaneous)


## Basic linear algebra

 - [BLAS](http://www.netlib.org/blas/) -
   Standard building blocks for performing basic vector and matrix operations.
   [Fortran, public domain, [GitHub](https://github.com/Reference-LAPACK/lapack/tree/master/BLAS)]

 - [OpenBLAS](https://www.openblas.net/) -
   Optimized BLAS library based on GotoBLAS2.
   [Fortran, BSD, [GitHub](https://github.com/xianyi/OpenBLAS)]

 - [BLIS](https://github.com/flame/blis) -
   Portable software framework for instantiating high-performance BLAS-like
   dense linear algebra libraries.
   [C++, BSD, GitHub]

 - [LAPACK](http://www.netlib.org/lapack/) -
   Routines for solving systems of simultaneous linear equations, least-squares
   solutions of linear systems of equations, eigenvalue problems, and singular
   value problems.
   [Fortran, BSD, [GitHub](https://github.com/Reference-LAPACK/lapack)]

 - [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) -
   C++ template library for linear algebra: matrices, vectors, numerical
   solvers, and related algorithms.
   [C++, MPL 2, [BitBucket](https://bitbucket.org/eigen/eigen)]


## Frameworks

 - [PETSc](http://www.mcs.anl.gov/petsc/) -
   Suite of data structures and routines for the scalable (parallel) solution
   of scientific applications modeled by partial differential equations.
   [C, 2-clause BSD, [BitBucket](https://bitbucket.org/petsc/petsc/src)]

 - [DUNE Numerics](https://www.dune-project.org/) -
   Modular toolbox for solving partial differential equations with grid-based
   methods.
   [C++, GPL 2, [GitLab](https://gitlab.dune-project.org/core/)]

 - [Trilinos](https://trilinos.org/) -
   Algorithms and enabling technologies for the solution of large-scale,
   complex multi-physics engineering and scientific problems.
   [mostly C++, mostly BSD, [GitHub](https://github.com/trilinos/)]

 - [SciPy](https://www.scipy.org/) -
   Python modules for statistics, optimization, integration, linear algebra,
   Fourier transforms, signal and image processing, ODE solvers, and more.
   [Python, mostly BSD, [GitHub](https://github.com/scipy/scipy/)]

 - [NumPy](http://www.numpy.org/) -
   Fundamental package needed for scientific computing with Python.
   [Python, BSD, [GitHub](https://github.com/numpy/numpy)]


## Finite Elements

 - [FEniCS](https://fenicsproject.org/) -
   Open-source computing platform for solving partial differential equations
   with Python and C++ interfaces.
   [C++/Python, LGPL 3, [BitBucket](https://bitbucket.org/fenics-project/)]

 - [libMesh](https://libmesh.github.io/) -
   Framework for the numerical simulation of partial differential equations
   using arbitrary unstructured discretizations on serial and parallel
   platforms.
   [C++, LGPL 2.1, [GitHub](https://github.com/libMesh/libmesh)]

 - [deal.II](http://dealii.org/) -
   Software library supporting the creation of finite element codes and an open
   community of users and developers.
   [C++, LGPL 2.1, [GitHub](https://github.com/dealii/dealii)]

 - [Netgen/NGSolve](https://ngsolve.org/) -
   High performance multiphysics finite element software.
   [C++, LGPL 2.1, [GitHub](https://github.com/NGSolve/netgen)]

 - [Firedrake](https://www.firedrakeproject.org/) -
   Automated system for the portable solution of partial differential equations
   using the finite element method.
   [Python, LGPL 3, [GitHub](https://github.com/firedrakeproject/firedrake)]

 - [MOOSE](http://www.mooseframework.org) -
   Multiphysics Object Oriented Simulation Environment.
   [C/Python, LGPL 2.1, [GitHub](https://github.com/idaholab/moose)]

 - [MFEM](http://mfem.org/) -
   Free, lightweight, scalable C++ library for finite element methods.
   [C++, LGPL 2.1, [GitHub](https://github.com/mfem/mfem)]


## Meshing

 - [Gmsh](http://gmsh.info/) -
   Three-dimensional finite element mesh generator with built-in pre- and
   post-processing facilities.
   [C++, GPL, [GitLab](https://gitlab.onelab.info/gmsh/gmsh)]

 - [pygmsh](https://github.com/nschloe/pygmsh) -
   Python interface for Gmsh.
   [Python, MIT, GitHub]

 - [MeshPy](https://mathema.tician.de/software/meshpy/) -
   Quality triangular and tetrahedral mesh generation.
   [Python, MIT, [GitHub](https://github.com/inducer/meshpy)]

 - [meshio](https://github.com/nschloe/meshio) -
   I/O for various mesh formats, file conversion.
   [Python, MIT, GitHub]

 - [CGAL](https://www.cgal.org/) -
   Efficient and reliable algorithms in computational geometry.
   [C++, mixed LGPL/GPL, [GitHub](https://github.com/CGAL/cgal)]

 - [pygalmesh](https://github.com/nschloe/pygalmesh) -
   Python interface for CGAL's 3D meshing capabilities.
   [Python, MIT, GitHub]

 - [mshr](https://bitbucket.org/fenics-project/mshr/) -
   The mesh generation component of FEniCS.
   [Python, GPL 3, BitBucket]

 - [MOAB](http://sigma.mcs.anl.gov/moab-library/) -
   Representing and evaluating mesh data.
   [C++, mostly LGPL3, [BitBucket](https://bitbucket.org/fathomteam/moab/)]

 - [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) -
   Software libraries and data formats for the creation, access, and sharing of
   array-oriented scientific data.
   [C/C++/Fortran/Java/Python, [custom open-source
   license](https://www.unidata.ucar.edu/software/netcdf/copyright.html),
   [GitHub](https://github.com/Unidata/netcdf-c/)]

 - [HDF5](https://support.hdfgroup.org/HDF5/) -
   Data model, library, and file format for storing and managing data.
   [C/Fortran, BSD]

 - [XDMF](http://www.xdmf.org/index.php/Main_Page) -
   eXtensible Data Model and Format to exchange scientific data between High
   Performance Computing codes and tools.
   [C++, [GitLab](https://gitlab.kitware.com/xdmf/xdmf)]

 - [TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen) -
   Quality Tetrahedral Mesh Generator and a 3D Delaunay Triangulator.
   [C++, AGPLv3]

 - [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) -
   Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.
   [C, *nonfree software*]


## Sparse linear solvers

 - [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) -
   General purpose library for the direct solution of large, sparse,
   nonsymmetric systems of linear equations.
   [C, mostly BSD, [GitHub](https://github.com/xiaoyeli/superlu)]

 - [KryPy](https://github.com/andrenarchy/krypy) -
   Krylov subspace methods for the solution of linear algebraic systems.
   [Python, MIT, GitHub]

 - [PyAMG](https://pyamg.github.io/) -
   Algebraic Multigrid Solvers in Python.
   [Python, MIT, [GitHub](https://github.com/pyamg/pyamg)]


## Miscellaneous

 - [FFTW](http://www.fftw.org/) -
   Computes the discrete Fourier transform in one or more dimensions, of
   arbitrary input size, and of both real and complex data.
   [C, GPL2, [GitHub](https://github.com/FFTW/fftw3)]

 - [Qhull](http://www.qhull.org/) -
   Computes the convex hull, Delaunay triangulation, Voronoi diagram, halfspace
   intersection about a point, furthest-site Delaunay triangulation, and
   furthest-site Voronoi diagram.
   [C/C++, [custom open source license](http://www.qhull.org/COPYING.txt),
   [GitHub](https://github.com/qhull/qhull/)]

 - [GSL](https://www.gnu.org/software/gsl/) -
   Wide range of mathematical routines such as random number generators,
   special functions, and least-squares fitting.
   [C/C++, GPL 3, [Savannah](https://savannah.gnu.org/projects/gsl)]

 - [OpenFOAM](https://www.openfoam.com/) -
   Free, open source CFD (computational fluid dynamics) software.
   [C++, GPL 3, [GitHub](https://github.com/OpenFOAM/OpenFOAM-dev)]

 - [ParaView](https://www.paraview.org/) -
   Open-source, multi-platform data analysis and visualization application
   based on VTK.
   [C++, BSD, [GitLab](https://gitlab.kitware.com/paraview/paraview)]

 - [quadpy](https://github.com/nschloe/quadpy) -
   Numerical integration (quadrature, cubature) in Python.
   [Python, MIT, GitHub]

 - [FiPy](https://www.ctcms.nist.gov/fipy/) -
   Finite-volume PDF solver.
   [Python, [custom open-source
   license](https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software),
   [GitHub](https://github.com/usnistgov/fipy)]

 - [accupy](https://github.com/nschloe/accupy) -
   Accurate sums and dot products for Python.
   [Python, MIT, GitHub]


## License

[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

To the extent possible under law, [Nico Schlömer](https://github.com/nschloe)
has waived all copyright and related or neighboring rights to this work.
