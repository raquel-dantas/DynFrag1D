import numpy as np
import itertools
import DFMesh

# class DFMesh:
#     DFMesh.L
#     connect = list()




# Global index of local DOF:
def Gl_index(elem_index, local_dof):
    return DFMesh.connect[elem_index][local_dof]

# Boundary conditions
def Apply_bc(K, F, elem_index):
    phi = 1.0
    bignumber = 10.0**30
    dof = DFMesh.connect[elem_index][0]
    value, bc_type = DFMesh.bc_dict[DFMesh.materials[elem_index]]
    if bc_type == "dirichlet":
        K[dof, dof] += bignumber*phi*phi
        F[dof] += bignumber*float(value)*phi


def Contribute_el(K, M, F, elem_index):
    """Computes the contribution of element in the global stiffness and mass matrices, and load vector.\n
    Arguments:\n
    K -- Global stiffness matrix which will be alterated to take account the contribution;\n
    M -- Global stiffness matrix;\n
    F -- Global load vector;\n
    elem_index -- element index"""
    
    # Element stiffness matrix and element load vector
    h = DFMesh.L/DFMesh.n_el
    k_elem = DFMesh.E*DFMesh.A/h * np.array([[1.0, -1.0], [-1.0, 1.0]])
    m_elem = DFMesh.rho*DFMesh.A*h/2 * np.array([[1.0, 0.0], [0.0, 1.0]])
    f_elem = np.array([0.0, 0.0])
    # Contribution to K, M and F
    for i_loc in range(2):
        i_gl = Gl_index(elem_index, i_loc)
        F[i_gl] += f_elem[i_loc]
        for j_loc in range(2):
            j_gl = Gl_index(elem_index, j_loc)
            K[i_gl, j_gl] += k_elem[i_loc, j_loc]
            M[i_gl, j_gl] += m_elem[i_loc, j_loc]

# Global stiffness matrix and mass matrix
def GlobalSystem():
    """ Returns global stiffness and mass matrices, and global load vector.    
    """
    n_dofs = max(list(itertools.chain.from_iterable(DFMesh.connect))) + 1
    K = np.zeros((n_dofs, n_dofs))
    F = np.zeros((n_dofs))
    M = np.zeros((n_dofs, n_dofs))
    # Assembly of elements
    n_el = len(DFMesh.connect)
    for i_el in range(n_el):
        mat_id = DFMesh.materials[i_el]
        if mat_id == 0:
            Contribute_el(K, M, F, i_el)
        elif mat_id == 2:
            Apply_bc(K, F, i_el)

    return K, M, F