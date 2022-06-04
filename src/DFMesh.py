import numpy as np
import itertools
import copy 


# Young's module (Pa)
E = 275.0*10**9  
# Density (kg/m3)
rho = 2750.0  

# Lenght of the bar (m)
L = 50*10**-3  
x0 = -L/2
xf = L/2

# Number of linear elements (n_el)
n_el = 1500
# Lenght of each linear element (h)
h = L/n_el

# Limit stress / critical stress (stress_c) (Pa)
stress_c = 300.0*10**6 
# Assuming a random distribution of critical stress in the linear elements
diststress_c = np.random.uniform(low=299*10**6, high=301*10**6, size=(n_el))

# Applied strain rate (s-1)
# strain_rate = 10.0**2 
# strain_rate = 10.0**3  
strain_rate = 10.0**4  
# strain_rate = 10.0**5 

# Critical time step 
dt_crit = h/((E/rho)**0.5)
# Adopted time step (s)
dt = dt_crit*0.1  

# Time peak stress (stress_c)
time_peakstress = stress_c / (E * strain_rate)
print(time_peakstress)
nstep_peak = int(time_peakstress/dt)
print(nstep_peak)
# Total time of simulation (s)
time_simulation = 6.0*10**-7
# Number of time steps (n_steps)
n_steps = int(time_simulation/dt)
print(dt)
print(n_steps)

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

# node_id[elem_index][local_node], returns the global node id
node_id = [[i,i+1] for i in range(n_el)]

# Identify each node has apllied BCs
# node_id.append([0]) # Support at left boundary
node_id.append([0]) # Applied velocity at left boundary
node_id.append([n_el]) # Applied velocity at right boundary 

# Connect[el][j] returns the global index of local dof 'j' of the element 'el'
connect = copy.deepcopy(node_id) 

# Applied velocity
vel = strain_rate*L/2 

# BC dictionary
bc_dict = {
    2: (0, "dirichlet"),
    4: (-vel, "velocity"),
    5: (vel, "velocity")
}
# Number of degree of freedom 
n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1

# Cross sectional area (m2)
A = 1*10**-3  
# A = 1.0  
# Fracture energy (N/m)
Gc = 100.0 

# Limit fracture oppening
delta_c = (2.0*Gc)/stress_c
# Assuming the random distribution of critical stress
distdelta_c = np.zeros(n_el)
for el in range(n_el):
    distdelta_c[el] = (2.0*Gc)/diststress_c[el]



# Initial values

# Initial velocity (v0): velocity profile (vel) is a function v(x)
n_points = n_dofs
l = np.linspace(-L/2, L/2, n_points)
v0 = np.array([strain_rate*x for x in l])
v0 = np.round(v0, 8)

# Initial displacement (u0)
if strain_rate < 5.0 * 10.0**3:
    u0 = np.array([0.98*stress_c*x / E for x in l])
else:
    u0 = np.zeros((n_dofs))

# Initial acceleration (acel0)
acel0 = np.zeros((n_dofs))
# External forces (p)
p = np.zeros((n_steps+1, n_dofs))
# Damping
C = np.zeros((n_dofs, n_dofs))
# Initialization of maximum jump u between two linear elements (delta_max)
delta_max = np.zeros((len(materials)*2))
# Contact penalty
alpha = (stress_c**2 + 4.5 * strain_rate**(2/3) * E * Gc**(2/3) * rho**(1/3)) / (4.5 * Gc)
print(alpha)
distalpha = np.zeros(n_el)
for el in range(n_el):
    distalpha[el] = (diststress_c[el]**2 + 4.5 * strain_rate**(2/3) * E * Gc**(2/3) * rho**(1/3)) / (4.5 * Gc)
    




def NodeCoord(node_id):
    return x0 + node_id*h


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
        DofCoord[dof,0] = NodeCoord(node_id[el][0])
        dof = connect[el][1]
        DofCoord[dof,0] = NodeCoord(node_id[el][1])
    return DofCoord
