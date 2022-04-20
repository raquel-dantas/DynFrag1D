import numpy as np
import DFMesh
import DFInterface
import DFFem

# Newmark's method explicit
def Newmark_exp(K, M, C, u, v, acel, p_next, dt, gamma):

    # Degrees of freedom
    dofs = K.shape[0]
    
    u_next = np.zeros((dofs))
    acel_next = np.zeros((dofs))
    v_next = np.zeros((dofs))
    
    u_next = u + dt*v + ((1.0/2.0)*dt**2)*acel

    for i in range(dofs):
        coord = DFMesh.u0[i] - u_next[i]
        if coord < 0:
            u_next[i] = 0

    # Solution of the linear problem:
    Minv = np.linalg.inv(M)
    f_int = DFInterface.ForceInt(u_next)
    inertia = p_next - f_int 
    acel_next = np.matmul(Minv, inertia)
    
    v_next = v + (1 - gamma)*dt*acel + gamma*dt*acel_next
    

    for i_el in range(len(DFMesh.connect)):
        if DFMesh.materials[i_el] == 2:
            value, bc_type = DFMesh.bc_dict[DFMesh.materials[i_el]]
            if bc_type == "dirichlet":
                i_dof = DFFem.Gl_index(i_el,0)
                u_next[i_dof] = value
                acel_next[i_dof] = 0.0
                v_next[i_dof] = 0.0
                

    return u_next, v_next, acel_next