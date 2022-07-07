import numpy as np
import DFMesh
import DFInterface
import DFFem


def Newmark_exp(M, u, v, acel, p_next, dt):
    """Apply Newmarks explicity integration scheme. Returns vectors with displacement, velocity and aceleration in all dofs for the next time step.\n
    Arguments: \n
    K -- Global stiffness matrix; \n
    M -- Global mass matrix; \n
    u -- displacement in the current time step;\n
    v -- velocity in the current time step; \n
    acel -- velocity in the current time step;\n
    p_next -- applied forces of next time step; \n
    dt -- time increment."""

    # Newmark explicity constants
    gamma = 0.5
    beta = 0.0

    # Number of degrees of freedom
    dofs = M.shape[0]
    
    # Initiation of variables
    u_next = np.zeros((dofs))
    acel_next = np.zeros((dofs))
    v_next = np.zeros((dofs))
    
    # u_next returns a vector with the displacement in all dofs for the next time step 
    u_next = u + dt*v + ((1.0/2.0)*dt**2)*acel

    # Solution of the linear problem: acel_next returns a vector with the acceleration in all dofs for the next time step 
    Minv = np.linalg.inv(M)
    f_int = DFInterface.InternalForce(u_next)
    inertia = p_next - f_int 
    acel_next = np.matmul(Minv, inertia)
    
    # v_next returns a vector with the velocity in all dofs for the next time step 
    v_next = v + (1 - gamma)*dt*acel + gamma*dt*acel_next

    # If there is a dirichlet conditions apply:
    for i_el in range(len(DFMesh.connect)):
        if DFMesh.materials[i_el] == 2:
            value, bc_type = DFMesh.bc_dict[DFMesh.materials[i_el]]
            if bc_type == "dirichlet":
                i_dof = DFFem.Gl_index(i_el,0)
                u_next[i_dof] = value
                acel_next[i_dof] = 0.0
                v_next[i_dof] = 0.0

    return u_next, v_next, acel_next