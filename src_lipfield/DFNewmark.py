import numpy as np
import scipy
from scipy.optimize import minimize
import DFMesh
import DFInterface
import DFFem


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
    dp = np.zeros((DFMesh.n_el))
    d_next = np.zeros((DFMesh.n_el))
    acel_next = np.zeros((dofs))
    v_next = np.zeros((dofs))
    
    # u_next returns a vector with the displacement in all dofs for the next time step 
    u_next = u + dt*v + ((1.0/2.0)*dt**2)*acel

    # Velocity predictor
    vp = v + (1 - gamma)*dt*acel

    epsilon = [(u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]) / DFMesh.ElemLength(el) for el in range(DFMesh.n_el)]

    # Damage prediction

    Yc = DFMesh.sigma_c**2/(2*DFMesh.E)
    lip_constraint = 2.21*10**-6
    lamb = 2*Yc*lip_constraint/DFMesh.Gc
    h = lambda d: (2*d-d**2)/(1-d+lamb*d**2)**2
    # TODO: numerical integration on func
    func = lambda d: (0.5*(1-d)**2*DFMesh.E*epsilon**2 + Yc*h(d) )
    result = minimize(
        fun=func, 
        x0=d, 
        method='SLSQP', 
        bounds=zip(d,[1]*DFMesh.n_el),
        tol=1e-5
    )
    ddash = result.x

    # Coords integration point 
    x = [(DFMesh.node_coord[el+1]+ DFMesh.node_coord[el])*0.5 for el in range(DFMesh.n_el)]
    # y = x
    # minimize(lambda j: ddash[el]+abs(y[j]-x[el])/lip_constraint, x0= ) 


    # TODO: Constraints of the minimizations may have to be put on the arguments but also on the minimized functional

    # Projections
    pi_lower =  [minimize(
            lambda y: ddash[el]+abs(y-x[el])/lip_constraint, 
            x0=0.,
            method='SLSQP', 
            bounds=(DFMesh.x0, DFMesh.xf), 
            tol=1e-5
        ).fun
        for el in range(DFMesh.n_el)
    ]
    pi_upper =  [minimize(
            lambda y: ddash[el]-abs(y-x[el])/lip_constraint, 
            x0=0.,
            method='SLSQP', 
            bounds=(DFMesh.x0, DFMesh.xf), 
            tol=1e-5
        ).fun
        for el in range(DFMesh.n_el)
    ]

    # Compute d_next
    # Lipchitz constraints
    # Constraint type ineq means that is a non negative result
    # TODO: DFMesh.ElemLength(i) is assuming a uniform mesh
    cons = ({'type': 'ineq', 
            'fun': lambda i: -(d[i] - d[i+1] -DFMesh.ElemLength(i)/lip_constraint)},
            {'type': 'ineq', 
            'fun': lambda i: -(d[i] - d[i-1] -DFMesh.ElemLength(i)/lip_constraint)})



    # Solution of the linear problem: acel_next returns a vector with the acceleration in all dofs for the next time step 
    # Minv = np.linalg.inv(M)
    f_int = DFInterface.InternalForce(u_next)
    inertia = p_next - f_int 
    acel_next = np.linalg.solve(M[:dofs,:dofs], inertia)
    
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