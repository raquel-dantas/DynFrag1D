import numpy as np
import itertools
import copy



E = 275.0*10**9     # Young's module (Pa)
rho = 2750.0        # Density (kg/m3)
Gc = 100.0          # Fracture energy (N/m)
stress_critical = 300.0*10**6   # Limit stress / critical stress (Pa)

A = 1*10**-3        # Cross sectional area (m2)
L = 1.05*10**-3     # Lenght of the bar (m)
x0 = 0              # Left extremitiy x coordinate / 0-initial
xf = L              # Rigth extremitiy x coordinate / f-final
n_el = 10           # Number of linear elements (n_el)
hun = L/n_el        # Size of the elemenets (h) for a uniform mesh (un) 





# Set applied velocity 

# Applied strain rate (s-1)
# strain_rate = 10.0**2
# strain_rate = 10.0**3
# strain_rate = 10.0**4
strain_rate = 10.0**5

# Applied velocity
vel = strain_rate*L/2


# Set BC's

# Material id convention:
# 0 : line elemnt
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
# node_id.append([0])         # Support at left boundary
node_id.append([0])         # Applied velocity at left boundary
node_id.append([n_el])      # Applied velocity at right boundary


# Connect[el][j] returns the global index of local dof 'j' of the element 'el'
connect = copy.deepcopy(node_id)

# Number of degree of freedom
n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1
n_points = n_dofs       # Number of points in the mesh

# Points coordinates for a uniform mesh
l = np.linspace(x0, xf, n_points)
node_coord = l
# Points coordinates for non-uniform mesh 
# node_coord = np.array([x + np.random.uniform(low=-0.4, high=0.4) * hun for x in l])
# node_coord[0] = x0
# node_coord[n_el] = xf




# Interface parameters

# A critical stress field to fabricate concentrated damage at L/2
# a = L/2
# k = 3.5
# b = (L/2)**0.5/(stress_critical*(k-1))
# sigmac = [(abs(node_coord[i]-a)**0.5)/b + stress_critical for i in range(n_el-1)]
delta_critical = (2.0*Gc)/stress_critical       # Limit crack oppening
alpha_critical = (stress_critical**2 + 4.5 * strain_rate**(2/3) * E * Gc**(2/3) * rho**(1/3)) / (4.5 * Gc)
# sigmac returns a array with an random distribution of critical stress at the interface elements in order to consider heterogeneities 
sigmac = np.random.uniform(low=stress_critical-1, high=stress_critical+1, size=(n_el-1))
deltac = [(2.0*Gc)/sigmac[i] for i in range(n_el-1)]
# Contact penalty
alpha = [(sigmac[i]**2 + 4.5 * strain_rate**(2/3) * E * Gc**(2/3) * rho**(1/3)) / (4.5 * Gc) for i in range(n_el-1)]
# Initialization of maximum jump u between two linear elements (delta_max)
delta_max = np.zeros((len(materials)*2))

# Set time increment
dt_crit = 0.2*hun/((E/rho)**0.5)    # Critical time step
dt = dt_crit*0.4                    # Adopted time step (s)
time_simulation = 4.0*10**-7        # Total time of simulation (s)
n_steps = int(time_simulation/dt)   # Number of time-steps 
# Time related to peak stress
time_peakstress = stress_critical / (E * strain_rate)   # Time peak stress 
nstep_peak = int(time_peakstress/dt)                    # Time-step to peak stress 







# Set initial values

# Initial velocity (v0): velocity profile (vel) is a function v(x)
v0 = np.array([strain_rate*x for x in node_coord])

# Initial displacement
if strain_rate < 5.0 * 10.0**3:
    # Apply initial displacement neq zero to save computational time during pre-crack phase
    u0 = np.array([0.98*stress_critical*x / E for x in node_coord])
else:
    u0 = np.zeros((n_dofs))

 
acel0 = np.zeros((n_dofs))          # Initial acceleration
p = np.zeros((n_steps+1, n_dofs))   # External forces
C = np.zeros((n_dofs, n_dofs))      # Damping





def ElemLength(elem_index):
    """Returns the element length (hel)."""

    hel = node_coord[elem_index+1] - node_coord[elem_index]

    return hel



def GetEl(connect, dof_id):
    elid = 0
    for el in connect:
        locdof = 0
        for dof in el:
            if dof == dof_id:
                return elid, locdof
            locdof = locdof + 1
        elid = elid + 1


def ListDofCoord():
    """Returns a list of coordinates by dof.\n
    DofCoord[i] returns [x,y,z] of dof i.
    """

    ndofs = max(list(itertools.chain.from_iterable(connect))) + 1
    DofCoord = np.zeros((ndofs,3))
    for el in range(n_el):
        dof = connect[el][0]
        DofCoord[dof,0] = node_coord[node_id[el][0]]
        dof = connect[el][1]
        DofCoord[dof,0] = node_coord[node_id[el][1]]
    return DofCoord


