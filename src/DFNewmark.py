import numpy as np
import DFMesh
import DFInterface
import DFFem

# Newmark's method explicit
def Newmark_exp(K, M, C, u, v, acel, p_next, dt, gamma):

    # Degrees of freedom
    dofs = K.shape[0]
    
    up = np.zeros((dofs))
    vp = np.zeros((dofs))
    u_next = np.zeros((dofs))
    acel_next = np.zeros((dofs))
    v_next = np.zeros((dofs))

    u_next = u + dt*v + ((1.0/2.0)*dt**2)*acel
    # Predictor vector vp
    vp = v + (1 - gamma)*dt*acel

    # Solution of the linear problem:
    term1 = M + gamma*dt/2*C 
    Minv = np.linalg.inv(term1)
    f_int = DFInterface.ForceInt(u_next)
    term2 = p_next - f_int - np.matmul(C, vp)
    acel_next = np.matmul(Minv, term2)
    

    for i_el in range(len(DFMesh.connect)):
        if DFMesh.materials[i_el] == 2:
            value, bc_type = DFMesh.bc_dict[DFMesh.materials[i_el]]
            if bc_type == "dirichlet":
                i_dof = DFFem.Gl_index(i_el,0)
                u_next[i_dof] = value
                acel_next[i_dof] = 0.0
                v_next[i_dof] = 0.0
                
    # Correctors
    v_next = vp + gamma*dt*acel_next

    return u_next, v_next, acel_next