# -*- coding: utf-8 -*-
#
import os.path

import dmsh
import meshio
import optimesh
import numpy

# from dolfin import (
#     Mesh,
#     XDMFFile,
#     mpi_comm_world,
#     Constant,
#     DirichletBC,
#     inner,
#     grad,
#     dx,
#     FunctionSpace,
#     TestFunction,
#     TrialFunction,
#     Function,
#     solve
# )


def _create_mesh(filename):
    poly = dmsh.Polygon(
        [
            [-295.0, 0.0],
            [-160.0, 0.0],
            [-50.0, 110.0],
            [-50.0, 190.0],
            [+50.0, 190.0],
            [+50.0, 110.0],
            [+160.0, 0.0],
            [+295.0, 0.0],
            [+405.0, 110.0],
            [+405.0, 235.0],
            [+200.0, 430.0],
            [+170.0, 400.0],
            [+355.0, 235.0],
            [-355.0, 235.0],
            [-170.0, 400.0],
            [-200.0, 430.0],
            [-405.0, 235.0],
            [-405.0, 110.0],
        ]
    )

    geo = dmsh.Union(
        [
            poly,
            dmsh.Circle([-295.0, 110.0], 110.0),
            dmsh.Circle([+295.0, 110.0], 110.0),
            dmsh.Circle([-160.0, 110.0], 110.0),
            dmsh.Circle([+160.0, 110.0], 110.0),
        ]
    )

    X, cells = dmsh.generate(
        geo,
        35.0,
        # show=True
    )

    X, cells = optimesh.cvt.quasi_newton_uniform_full(
        X, cells, 1.0e-3, 100, verbose=True
    )

    X = numpy.column_stack([X[:, 0], X[:, 1], numpy.zeros(X.shape[0])])

    meshio.write_points_cells(filename, X, {"triangle": cells})
    return


def _main():
    this_dir = os.path.dirname(os.path.abspath(__file__))

    filename = os.path.join(this_dir, "sunglasses.xdmf")
    _create_mesh(filename)

    # mesh = Mesh()
    # with XDMFFile(mpi_comm_world(), filename) as f:
    #     f.read(mesh)

    # # Create mesh and define function space
    # V = FunctionSpace(mesh, "Lagrange", 1)

    # def boundary(x):
    #     return x[1] < 1.0e-10

    # # Define boundary condition
    # u0 = Constant(0.0)
    # bc = DirichletBC(V, u0, "on_boundary")

    # # Define variational problem
    # u = TrialFunction(V)
    # v = TestFunction(V)
    # f = Constant(1.0)
    # a = inner(grad(u), grad(v)) * dx
    # L = f * v * dx

    # # Compute solution
    # u = Function(V)
    # solve(a == L, u, bc)

    # filename = os.path.join(this_dir, "solution.xdmf")
    # xf = XDMFFile(filename)
    # xf.write(u)
    return


if __name__ == "__main__":
    _main()
