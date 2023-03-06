import numpy as np

import DFMesh
import DFFem


def insertInterface(el_left, el_right, u, v, acel):
    """Insert a interface element. Updates: number of DOFs, connect list, materials ID, node ID, and vectors u, v e acel vector. \n
    Returns u, v and acel \n
    Arguments:\n
    el_left -- element immediatly at left side from the interface element;\n
    el_right -- element immediatly at left side from the interface element;\n
    u -- displacement vector;\n
    v -- velocity vector;\n
    acel -- aceleration vector."""

    # Identify the crack
    dof_broken = DFMesh.connect[el_left][1]

    # Update the connectivity list
    n_dofs = u.shape[0]
    DFMesh.connect[el_right][0] = n_dofs
    dof_right = n_dofs
    dof_left = DFMesh.connect[el_left][1]
    DFMesh.connect.append([dof_left, dof_right])

    # Assign interface material
    DFMesh.materials.append(1)

    # Update node id
    DFMesh.node_id.append(DFMesh.node_id[el_left][1])

    # Update u, v and acel arrays
    u = np.append(u, u[dof_broken])
    v = np.append(v, v[dof_broken])
    acel = np.append(acel, acel[dof_broken])

    return u, v, acel


def stressCohesiveLaw(jump_u, el_index):
    """Returns the stress for interface element through a linear cohesive law. \n
    Arguments:\n
    jump_u -- jump in the displacement between the DOFs at right and left sizes of the interface element;\n
    el_index -- cohesive element index."""

    # i is the crack position to relate with the critical values array
    i = DFMesh.connect[el_index][0] - 1

    if jump_u >= 0:
        if DFMesh.jump_max[el_index] > jump_u:
            stress_max = DFMesh.stress_critical[i] * (
                1.0 - DFMesh.jump_max[el_index] / DFMesh.crack_critical[i]
            )
            stress_interface = stress_max / DFMesh.jump_max[el_index] * jump_u
        else:
            stress_interface = DFMesh.stress_critical[i] * (
                1.0 - min(jump_u / DFMesh.crack_critical[i], 1.0)
            )
            DFMesh.jump_max[el_index] = min(jump_u, DFMesh.crack_critical[i])
    else:
        stress_interface = DFMesh.contact_penalty[i] * jump_u

    return stress_interface


def getDamageParameter():
    """Returns the damage parameter for an interface element.\n"""

    for el_index in range(len(DFMesh.materials)):
        if DFMesh.materials[el_index] == 1:
            i = DFMesh.connect[el_index][0] - 1
            return min(1.0, DFMesh.jump_max[el_index] / DFMesh.crack_critical[i])
        else:
            return 0.0


def forceInterface(u):
    """Returns the force on interfaces (force_lambda).\n
    Arguments:\n
    u -- displacemnt vector."""

    n_dofs = u.shape[0]
    force_lambda = np.zeros(n_dofs)

    for el in range(len(DFMesh.materials)):
        if DFMesh.materials[el] == 1:
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            force_lambda[DFMesh.connect[el][0]] = -stressCohesiveLaw(jump_u, el) * DFMesh.area
            force_lambda[DFMesh.connect[el][1]] = -force_lambda[DFMesh.connect[el][0]]

    return force_lambda


def internalForce(u):
    """Returns the internal force vector (ku + force_lambda)\n
    Arguments:\n
    u -- displacemnt vector for all dofs."""

    n_dofs = u.shape[0]
    fint = np.zeros(n_dofs)

    for el in range(DFMesh.n_elements):
        # u_loc returns a vector contained u the current element
        u_loc = np.array([u[DFFem.getGlobalIndex(el, 0)], u[DFFem.getGlobalIndex(el, 1)]])
        fint_loc = np.matmul(DFFem.k_elem, u_loc) / DFMesh.getElemLength(el)
        # Contribution of each dof in the internal force vector
        for i_loc in range(2):
            i_gl = DFFem.getGlobalIndex(el, i_loc)
            fint[i_gl] += fint_loc[i_loc]

    # The internal force is the sum of the force from the linear and interface elements
    fint += forceInterface(u)

    return fint
