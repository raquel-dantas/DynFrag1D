import numpy as np
import itertools

import DFMesh


# Local stiffness and mass matrices and local force vector
k_elem = DFMesh.young_modulus * DFMesh.area * np.array([[1.0, -1.0], [-1.0, 1.0]])
m_elem = DFMesh.rho * DFMesh.area / 2 * np.diag([1, 1])
f_elem = np.array([0.0, 0.0])


def getGlobalIndex(elem_index, local_dof):
    """Returns the global index of a local dof."""

    return DFMesh.connect[elem_index][local_dof]


def applyBC(M, F, elem_index):
    """Apply boundary conditions.\n
    Arguments:\n
    M -- Global stiffness matrix;\n
    F -- Global load vector.\n"""

    phi = 1.0
    # Penalty number
    bignumber = 10.0**30
    dof = DFMesh.connect[elem_index][0]
    value, bc_type = DFMesh.bc_dict[DFMesh.materials[elem_index]]
    if bc_type == "dirichlet":
        F[dof] += bignumber * float(value) * phi
        M[dof, dof] += float("inf")
    if bc_type == "velocity":
        M[dof, dof] += float("inf")


def contributeEl(M, F, elem_index):
    """Computes the contribution of element in the global mass matrix, and load vector.\n
    Arguments:\n
    M -- Global mass matrix;\n
    F -- Global load vector."""

    # Contribution to M and F
    i_loc = 0
    i_gl = getGlobalIndex(elem_index, i_loc)
    F[i_gl] += f_elem[i_loc]
    M[i_gl, i_gl] += m_elem[i_loc, i_loc] * DFMesh.getElemLength(elem_index)

    i_loc = 1
    i_gl = getGlobalIndex(elem_index, i_loc)
    F[i_gl] += f_elem[i_loc]
    M[i_gl, i_gl] += m_elem[i_loc, i_loc] * DFMesh.getElemLength(elem_index)


def contribute(M, F, elem_index):
    material_id = DFMesh.materials[elem_index]
    if material_id == 0:
        contributeEl(M, F, elem_index)
    elif material_id in DFMesh.bc_dict:
        applyBC(M, F, elem_index)


def globalSystem():
    """Returns global stiffness and mass matrices, and global load vector."""

    n_dofs = max(list(itertools.chain.from_iterable(DFMesh.connect))) + 1
    F = np.zeros(n_dofs)
    M = np.diag(np.zeros(DFMesh.n_elements * 2))
    M.fill(0)

    # Assembly of elements
    n_elements = len(DFMesh.connect)
    [contribute(M, F, i_el) for i_el in range(n_elements)]

    return M, F
