import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
import DFMesh
import DFFem
import DFDamage


def Newmark_exp(M, u, v, acel, d, p_next, dt):
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
    dofs = u.shape[0]

    # Initiation of variables
    u_next = np.zeros((dofs))
    vp = np.zeros((dofs))
    # dp = np.zeros((DFMesh.n_el))
    d_next = np.zeros((DFMesh.n_el))
    acel_next = np.zeros((dofs))
    v_next = np.zeros((dofs))


    # u_next returns a vector with the displacement in all dofs for the next time step
    u_next = u + dt*v + ((1.0/2.0)*dt**2)*acel

    # Velocity predictor
    vp = v + (1. - gamma)*dt*acel

    # Strain for u_next
    strain = [(u_next[DFMesh.connect[el][1]] - u_next[DFMesh.connect[el][0]]) / DFMesh.ElemLength(el) for el in range(DFMesh.n_el)]

    # Compute damage next time-step
    def func(d): return DFDamage.w*sum([
            (0.5*(1. - d[el])**2 * 
            DFMesh.E*strain[el]**2 + 
            DFDamage.Yc[el] * 
            DFDamage.h(DFDamage.lamb[el], d[el])) * 
            DFMesh.hun*0.5
            for el in range(DFMesh.n_el)
            ])
    

    A = scipy.sparse.eye(DFMesh.n_el-1, DFMesh.n_el) - scipy.sparse.eye(DFMesh.n_el-1, DFMesh.n_el, 1)
    b = DFMesh.hun/DFDamage.l
    const = LinearConstraint(A, -b * np.ones(DFMesh.n_el-1), b * np.ones(DFMesh.n_el-1))

    dlip_opt = minimize(
        fun=func,
        x0=d,
        method='SLSQP',
        bounds=zip(d, [1.]*DFMesh.n_el),
        tol=1e-6,
        constraints=const,
    )
    d_next = dlip_opt.x
    
    if dlip_opt.success == False:
        raise Exception('optimization failed')

    if np.linalg.norm(d_next)>0:
        pass

    # Solution of the linear problem: acel_next returns a vector with the acceleration in all dofs for the next time step
    f_int = DFFem.InternalForce(u_next, d_next)
    inertia = p_next - f_int
    acel_next = np.linalg.solve(M[:dofs, :dofs], inertia)
    # print(acel_next)
    # v_next returns a vector with the velocity in all dofs for the next time step
    v_next = vp + gamma*dt*acel_next
    
    # If there is a dirichlet conditions apply:
    for i_el in range(len(DFMesh.connect)):
        if DFMesh.materials[i_el] == 2:
            value, bc_type = DFMesh.bc_dict[DFMesh.materials[i_el]]
            if bc_type == "dirichlet":
                i_dof = DFFem.Gl_index(i_el, 0)
                u_next[i_dof] = value
                acel_next[i_dof] = 0.0
                v_next[i_dof] = 0.0

    return u_next, v_next, acel_next, d_next
