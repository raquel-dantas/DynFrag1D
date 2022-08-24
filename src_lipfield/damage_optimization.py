from scipy.optimize import minimize
from matplotlib import pyplot as plt
import numpy as np
import itertools
import copy
import inspect


# Inputs 

E = 275.0*10**9     # Young's module (Pa)
rho = 2750.0        # Density (kg/m3)
Gc = 100.0          # Fracture energy (N/m)
stress_critical = 300.0*10**4   # Limit stress / critical stress (Pa)

A = 1*10**-3        # Cross sectional area (m2)
L = 1.05*10**-3     # Lenght of the bar (m)
x0 = 0              # Left extremitiy x coordinate / 0-initial
xf = L              # Rigth extremitiy x coordinate / f-final
n_el = 21           # Number of linear elements (n_el)
hun = L/n_el        # Size of the elemenets (h) for a uniform mesh (un)

strain_rate = 10.0**5   # Applied strain rate (s-1)
vel = strain_rate*L/2   # Applied velocity

# Set time increment
dt_crit = 0.2*hun/((E/rho)**0.5)    # Critical time step
dt = dt_crit*0.5                    # Adopted time step (s)
time_simulation = 4.0*10**-7        # Total time of simulation (s)
n_steps = int(time_simulation/dt)   # Number of time-steps 



# Set BC's

# Material id convention:
# 0 : line element
# 1 : interface element
# 2 : Support left node
# 3 : Support right node
# 4 : Velocity applied left node
# 5 : Velocity applied right node
materials = [0] * n_el
materials.append(4)
materials.append(5)

# BC dictionary
bc_dict = {
    2: (0, "dirichlet"),
    4: (-vel, "velocity"),
    5: (vel, "velocity")
}



# Mesh 

# node_id[elem_index][local_node], returns the global node id
node_id = [[i,i+1] for i in range(n_el)]

# Identify each node has apllied BCs
node_id.append([0])         # Applied velocity at left boundary
node_id.append([n_el])      # Applied velocity at right boundary

# Connect[el][j] returns the global index of local dof 'j' of the element 'el'
connect = copy.deepcopy(node_id)

n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1  # Number of degree of freedom
n_points = n_dofs                           # Number of points in the mesh
lp = np.linspace(x0, xf, n_points)           # Points coordinates for a uniform mesh
node_coord = lp                              # Node coordinates
x = [hun*i + hun*0.5 for i in range(n_el)]   # Coordenates of integration points



# A critical stress field to fabricate concentrated damage at L/2
a = L/2
k = 4
b = (L/2)**0.5/(stress_critical*(k-1))
sigmac = [(abs(node_coord[i]-a)**0.5)/b + stress_critical for i in range(n_el)]



# Set initial values

acel0 = np.zeros((n_dofs))          # Initial acceleration
p = np.zeros((n_steps+1, n_dofs))   # External forces
C = np.zeros((n_dofs, n_dofs))      # Damping

# Initial velocity (v0): velocity profile (vel) is a function v(x)
v0 = np.array([strain_rate*x for x in node_coord])
# v0 = np.zeros((n_dofs))

# Initial displacement
if strain_rate < 5.0 * 10.0**3:
    # Apply initial displacement neq zero to save computational time during pre-crack phase
    u0 = np.array([0.98*stress_critical*x / E for x in node_coord])
else:
    u0 = np.zeros((n_dofs))

# Initial damage 
# d0 = [0., 0.04525, 0.086, 0.154, 0.18125, 0.204, 0.22225, 0.236, 0.24525, 0.25, 0.25025, 0.246, 0.23725, 
# 0.224, 0.20625, 0.184, 0.15725, 0.126, 0.09025, 0.05, 0.00525]

d0 = [0.,0.0473944, 0.0898275, 0.127299, 0.15981, 0.187359, 0.209947, 0.227574,  0.24024,  0.247944, 0.250687,0.248469, 0.24129,  0.229149, 0.212047, 0.189984, 0.16296, 0.130974, 0.0940275, 0.0521194, 0.00525]


def ElemLength(elem_index):
    """Returns the element length (hel)."""

    hel = node_coord[elem_index+1] - node_coord[elem_index]

    return hel



# FEM

k_elem = E*A* np.array([[1.0, -1.0], [-1.0, 1.0]])  # Local stiffness matirx
m_elem = rho*A/2 * np.diag([1,1])                   # Local mass matrix (lumped)
f_elem = np.array([0.0, 0.0])                       # Local force vector
M = np.diag(np.zeros(n_el*2))                       # Initialize M 

def LocalSystem():
    """Returns local stiffness and mass matrices for linear 1D element."""

    return k_elem, m_elem



def Gl_index(elem_index, local_dof):
    """Returns the global index of a local dof."""

    return connect[elem_index][local_dof]



def Apply_bc(M, F, elem_index):
    """Apply boundary conditions.\n
    Arguments:\n
    M -- Global stiffness matrix;\n
    F -- Global load vector.\n"""

    phi = 1.0
    bignumber = 10.0**30    # Penalty number
    dof = connect[elem_index][0]
    value, bc_type = bc_dict[materials[elem_index]]
    if bc_type == "dirichlet":
        F[dof] += bignumber*float(value)*phi
        M[dof, dof] += float("inf")
    if bc_type == "velocity":
        M[dof, dof] += float("inf")



def Contribute_el(M, F, elem_index):
    """Computes the contribution of element in the global stiffness and mass matrices, and load vector.\n
    Arguments:\n
    M -- Global stiffness matrix;\n
    F -- Global load vector."""
    
    i_loc = 0
    i_gl = Gl_index(elem_index, i_loc)
    F[i_gl] += f_elem[i_loc]
    M[i_gl, i_gl] += m_elem[i_loc, i_loc]*ElemLength(elem_index)
    
    i_loc = 1
    i_gl = Gl_index(elem_index, i_loc)
    F[i_gl] += f_elem[i_loc]
    M[i_gl, i_gl] += m_elem[i_loc, i_loc]*ElemLength(elem_index)



def Contribute(M, F, elem_index):

        mat_id = materials[elem_index]
        if mat_id == 0:
            Contribute_el(M, F, elem_index)
        elif mat_id in bc_dict:
            Apply_bc(M, F, elem_index)



def GlobalSystem():
    """ Returns global stiffness and mass matrices, and global load vector."""

    n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1

    # Initiation of variables 
    F = np.zeros((n_dofs))
    M.fill(0)

    # Assembly of elements
    n_el = len(connect)
    [Contribute(M,F,i_el) for i_el in range(n_el)]

    return M, F


def DecreasingFun(damage, elem_index):

    return (1-damage[elem_index])**2


def InternalForce(u):
    """ Returns the internal force vector (ku)\n
    Arguments:\n
    u -- displacemnt vector for all dofs."""
    
    n_dofs = u.shape[0]
    fint = np.zeros(n_dofs)

    for el in range(n_el):
        # u_loc returns a vector contained u for a local dof
        u_loc = np.array([u[Gl_index(el, 0)], u[Gl_index(el, 1)]])
        # Decreasing function
        g = DecreasingFun(d,el)
        # Internal local force vector  
        fint_loc = np.matmul(k_elem*g, u_loc) / ElemLength(el)
        # Contribution of each dof in the internal force vector
        for i_loc in range(2):
            i_gl = Gl_index(el, i_loc)
            fint[i_gl] += fint_loc[i_loc]

    return fint



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
    dofs = u.shape[0]

    # Initiation of variables
    u_next = np.zeros((dofs))       # Displacement next time-step
    vp = np.zeros((dofs))           # Velocity predictor
    acel_next = np.zeros((dofs))    # Aceleration next time-step
    v_next = np.zeros((dofs))       # Velocity next time-step

    # u_next returns a vector with the displacement in all dofs for the next time step
    u_next = u + dt*v + ((1.0/2.0)*dt**2)*acel

    # Velocity predictor
    vp = v + (1 - gamma)*dt*acel

    # Solution of the linear problem: acel_next returns a vector with the acceleration in all dofs for the next time step
    f_int = InternalForce(u_next)
    inertia = p_next - f_int
    acel_next = np.linalg.solve(M[:dofs, :dofs], inertia)

    # v_next returns a vector with the velocity in all dofs for the next time step
    v_next = v + (1 - gamma)*dt*acel + gamma*dt*acel_next
    
    return u_next, v_next, acel_next


def GetStrain(u):
    
    # Total number of elements (linear + cohesive)
    numel = len(materials)
    # Initiation of variables
    strain = np.zeros(numel)

    for el in range(numel):
        if materials[el] == 0:
            # Strain[u,L] returns the strain value at each linear element 'el' 
            strain[el] = (u[connect[el][1]] - u[connect[el][0]]) / ElemLength(el)

    return strain


def PlotDamage(x, dn, ddash, upper, lower, dlip):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')

    plt.title(str("Damage field"))
    plt.xlabel(str("x"))
    plt.ylabel(str("D"))
    plt.plot(x, dn, label='dn')
    plt.plot(x, ddash, label='ddash')
    plt.plot(x, upper, label='upper')
    plt.plot(x, lower, label='lower')
    plt.plot(x, dlip, label='dlip')
    plt.legend()
    plt.show()

def PlotDamage_ContEl(dn, ddash, upper, lower, dlip):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')

    plt.title(str("Damage field"))
    plt.xlabel(str("x"))
    plt.ylabel(str("D"))

    x = np.array([[node_coord[el] , node_coord[el+1]] for el in range(n_el)])
    x = x.flatten()

    dn_elem = np.array([[dn[el], dn[el]] for el in range(n_el)])
    dn_elem = dn_elem.flatten()

    ddash_elem = np.array([[ddash[el], ddash[el]] for el in range(n_el)])
    ddash_elem = ddash_elem.flatten()

    upper_elem = np.array([[upper[el], upper[el]] for el in range(n_el)])
    upper_elem = upper_elem.flatten()

    lower_elem = np.array([[lower[el], lower[el]] for el in range(n_el)])
    lower_elem = lower_elem.flatten()

    dlip_elem = np.array([[dlip[el], dlip[el]] for el in range(n_el)])
    dlip_elem = dlip_elem.flatten()

    plt.plot(x, dn_elem, label='dn')
    plt.plot(x, ddash_elem, label='ddash')
    plt.plot(x, upper_elem, label='upper')
    plt.plot(x, lower_elem, label='lower')
    plt.plot(x, dlip_elem, label='dlip')
    plt.legend()
    plt.show()

# Main

# Initiation of variables
u = u0
v = v0
acel = acel0
d = d0

for n in range(n_steps):

    
    # Get K, M and F
    M, F = GlobalSystem()

    # u,v,acel returns a vector for u,v and acel at every dof at the n step
    u, v, acel = Newmark_exp(M, u, v, acel, F, dt)
    dn = d

    # Post process (stress, strain, energies)
    strain = GetStrain(u)


    # Damage regularization
    # l = 0.0001/2                        # Lip-field regularization lenght
    l = 0.00067
    Yc = [sigmac[el]**2/(2*E) for el in range(n_el)]       # Critical energy release rate
    lamb = [2*Yc[el]*l/Gc for el in range(n_el)]       
    def h(lamb,d): return (2*d-d**2)/(1-d+lamb*d**2)**2
    # Functional 
    w = 2.  # Weight quadrature
    def func(d): return w*sum([
            (0.5*(1-d[el])**2*E*strain[el]**2 + Yc[el]*h(lamb[el],d[el]))*ElemLength(el)/2.
            for el in range(n_el)
        ])


    # Damage prediction
    ddash = minimize(
        fun=func,
        x0=dn,
        method='SLSQP',
        bounds=zip(dn, [1]*n_el),
        tol=1e-6,
    ).x

    # Upper projection
    upper = [-minimize(
        lambda y: -ddash[np.searchsorted(node_coord, y[0])-1] + abs(x[el]-y[0])/l,
        x0=0.5*L,
        method='SLSQP',
        bounds=[(x0, xf)],
        tol=1e-6
    ).fun
        for el in range(n_el)
    ]

    # Lower projection
    lower = [minimize(
        lambda y: ddash[np.searchsorted(node_coord, y[0])-1] + abs(x[el]-y[0])/l,
        x0=0.5*L,
        method='SLSQP',
        bounds=[(x0, xf)],
        tol=1e-6
    ).fun
        for el in range(n_el)
    ]

    # Damage lip-field regularization
    # Constraints for regularized damage for time step n+1 (dlip)
    const22a = [{'type': 'ineq', 'fun': lambda d:
            -(d[i] - d[i+1] - ElemLength(i)/l)} for i in range(0,n_el-2)]
    const22b = [{'type': 'ineq', 'fun': lambda d:
            -(d[i] - d[i-1] - ElemLength(i)/l)} for i in range(1,n_el-1)]
    const = tuple([*const22a,*const22b])

    # Regularized damage for time step n+1 (dlip) (Optimization 21b)
    dlip = minimize(
        fun=func,
        x0=dn,
        method='SLSQP',
        bounds=zip(dn, [1]*n_el),
        tol=1e-6,
        constraints=const
    ).x

    # if n>200:
    #     d = dlip
    d = dlip
    
    PlotDamage(x, dn, ddash, upper, lower, dlip)
    PlotDamage_ContEl(dn, ddash, upper, lower, dlip)