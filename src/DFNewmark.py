import numpy as np
import copy
import DFMesh
import DFFem
import DFInterface
import DFDiffuseDamage


def explicitScheme(M, u, v, acel, d, p_next, dt):
    """Apply Newmarks explicit integration scheme. Returns vectors with displacement, velocity and aceleration in all dofs for the next time step.\n
    Arguments: \n
    M -- Global mass matrix; \n
    u -- displacement in the current time-step;\n
    v -- velocity in the current time-step; \n
    acel -- velocity in the current time-step;\n
    d -- damage at the current time-step; \n
    p_next -- applied forces of next time step; \n
    dt -- time increment."""

    # Newmark explicity constants
    gamma = 0.5
    beta = 0.0

    # Number of degrees of freedom
    dofs = u.shape[0]

    # Initiation of variables
    u_next = np.zeros(dofs)
    vp = np.zeros(dofs)
    v_next = np.zeros(dofs)
    d_next = np.zeros(DFMesh.n_elements)
    acel_next = np.zeros(dofs)

    # Compute displacement at next time-step 
    u_next = u + dt * v + (0.5 * dt**2) * acel
    
    # Compute velocity predictor (vp)
    vp = v + (1.0 - gamma) * dt * acel

    if DFMesh.use_cohesive_elements == True:
        d_next = d
        f_int = DFInterface.internalForce(u_next)
    

    if DFMesh.use_lipfield == True:

        # u_next_copy = copy.deepcopy(u_next)
        # d_copy = copy.deepcopy(d)
        # Compute damage at the next time-step (d_next)
        
        d_next = DFDiffuseDamage.computeDamageNextStep_useProjection(u_next, d, predictor_method='Newton', projection_method='FM')
        # d_next = DFDiffuseDamage.computeDamageNextStep_noProjection(u_next, d)
        
        # tolerance = 10e-8
        # assert(np.linalg.norm(u_next_copy - u_next) < tolerance)
        # assert(np.linalg.norm(d_copy - d) < tolerance)
        
        # Solution of the linear problem: acel_next returns a vector with the acceleration in all dofs for the next time step
        f_int = DFDiffuseDamage.internalForce(u_next, d_next)
    
    
    inertia = p_next - f_int
    # Compute aceleration and velocity at next time-step 
    acel_next = np.linalg.solve(M[:dofs, :dofs], inertia)
    v_next = vp + gamma * dt * acel_next

    # If there is a dirichlet conditions apply:
    for i_el in range(len(DFMesh.connect)):
        if DFMesh.materials[i_el] == 2:
            value, bc_type = DFMesh.bc_dict[DFMesh.materials[i_el]]
            if bc_type == "dirichlet":
                i_dof = DFFem.getGlobalIndex(i_el, 0)
                u_next[i_dof] = value
                acel_next[i_dof] = 0.0
                v_next[i_dof] = 0.0

    return u_next, v_next, acel_next, d_next
