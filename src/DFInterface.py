import numpy as np
import DFMesh
import DFFem


def InsertInterface(el_left, el_right, u, v, acel):
    """ Insert a interface element. Updates: number of DOFs, connect list, materials ID, node ID, and vectors u, v e acel vector. \n
    Returns u, v and acel \n
    Arguments:\n
    el_left -- element immediatly at left side from the interface element;\n
    el_right -- element immediatly at left side from the interface element;\n
    u -- displacement vector;\n
    v -- velocity vector;\n
    acel -- aceleration vector."""

    dof_broken = DFMesh.connect[el_left][1]
    # Update the connect list
    n_dofs = u.shape[0]
    DFMesh.connect[el_right][0] = n_dofs
    dof_right = n_dofs
    dof_left = DFMesh.connect[el_left][1]
    DFMesh.connect.append([dof_left, dof_right])
    # Update materials type
    DFMesh.materials.append(1)
    # Update node ID
    DFMesh.node_id.append(DFMesh.node_id[el_left][1])
    # Update u, v and acel arrays
    u = np.append(u, u[dof_broken])
    v = np.append(v, v[dof_broken])
    acel = np.append(acel, acel[dof_broken])

    return u, v, acel


def CohesiveLaw(jump_u,el_index):
    """ Returns the stress for interface element through a linear cohesive law. \n
    Arguments:\n
    jump_u -- jump in the displacement between the DOFs at right and left sizes of the interface element;\n
    el_index -- cohesive element index."""

    if DFMesh.delta_max[el_index] > jump_u:
        Tmax = DFMesh.stress_c * (1.0 - DFMesh.delta_max[el_index]/DFMesh.delta_c)
        stress = Tmax/DFMesh.delta_max[el_index] * jump_u
    else:
        stress = DFMesh.stress_c * (1.0 - jump_u/DFMesh.delta_c)
        DFMesh.delta_max[el_index] = min(jump_u,DFMesh.delta_c)
    
    # After contact law implementation update this point in the code
    stress = max(stress,0.0)

    return stress


def DamageParameter(el_index):
    """Returns the damage parameter for an interface element.\n"""

    if DFMesh.materials[el_index] == 1:
        return min(1.0,DFMesh.delta_max[el_index]/DFMesh.delta_c)
    else:
        return 0.0


def ForceInterface(u):
    """ Returns the force of interfaces (flambda).\n
    Arguments:\n
    u -- displacemnt vector."""

    n_dofs = u.shape[0]
    flambda = np.zeros(n_dofs)
    for el in range(len(DFMesh.materials)):
        if DFMesh.materials[el] == 1:
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            flambda[DFMesh.connect[el][0]] = CohesiveLaw(jump_u,el) * DFMesh.A
            flambda[DFMesh.connect[el][1]] = -flambda[DFMesh.connect[el][0]]

    return flambda


def ForceInt(u):
    """ Returns the internal force vector (ku + flambda)\n
    Arguments:\n
    u -- displacemnt vector."""
    
    # Element stiffness matrix
    h = DFMesh.L/DFMesh.n_el
    k_elem = DFMesh.E*DFMesh.A/h * np.array([[1.0, -1.0], [-1.0, 1.0]])
    n_dofs = u.shape[0]
    fint = np.zeros(n_dofs)
    for el in range(DFMesh.n_el):
        u_loc = np.array([u[DFFem.Gl_index(el, 0)], u[DFFem.Gl_index(el, 1)]])
        fint_loc = np.matmul(k_elem, u_loc)

        for i_loc in range(2):
            i_gl = DFFem.Gl_index(el, i_loc)
            fint[i_gl] += fint_loc[i_loc]

    fint += ForceInterface(u)

    return fint