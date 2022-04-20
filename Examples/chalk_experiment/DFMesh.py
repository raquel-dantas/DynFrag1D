import numpy as np
import itertools
import copy 

# Lenght of the bar (L)
L = 80*10**-3  # (m)
x0 = 0
xf = L
# Number of linear elements (n_el)
n_el = 5
# Lenght of each linear element (h)
h = L/n_el

# Material id convention: 
# 0 : line elemnt
# 1 : interface element
# 2 : Support left node
# 3 : Support right node 
# 4 : Velocity applied left node
# 5 : Velocity applied right node
materials = [0] * n_el
# materials.append(2) 
# materials.append(4)
# materials.append(5)

# node_id[elem_index][local_node], returns the global node id
node_id = []
for i in range(n_el):
    node_id.append([i, i+1])

# Identify each node has apllied BCs
# node_id.append([0]) # Support
# node_id.append([0]) # Applied velocity at left boundary
# node_id.append([n_el]) # Applied velocity at right boundary 

# Connect[el][j] returns the global index of local dof 'j' of the element 'el'
connect = copy.deepcopy(node_id) 

# Applied strain rate and veloctities
# strain_rate = 10.0**3  # (s-1)
# vel = strain_rate*L/2 
# vel = strain_rate*L

# BC dictionary
bc_dict = {
    2: (0, "dirichlet"),
    4: (0, "velocity"),
    5: (0, "velocity")
}
# Number of degree of freedom 
n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1

# Young's module 
E = 275.0*10**3  # (Pa)
# Cross sectional area
A = 1*10**-4  # (m2)
# Fracture energy
Gc = 100.0  # (N/m)
# Limit stress
stress_c = 300.0*10**6  # (Pa)
# Limit fracture oppening
delta_c = (2.0*Gc)/stress_c
# Density
rho = 2750.0*10  # (kg/m3)

# n_steps = 5000
dt_crit = h/(2*3**0.5*np.pi*(E/rho)**0.5)
dt = dt_crit*0.1  # (s)
# time_simulation = n_steps*dt # (s)
# print(time_simulation)

# Time integration
time_simulation = 0.2 # (s)
# Critical time step
# dt_crit = h/((E/rho)**0.5)
# dt_crit = h/(2*3**0.5*np.pi*(E/rho)**0.5)
# Adopted time step
# dt = dt_crit*0.1  # (s)
# Number of time steps (n_steps)
n_steps = int(time_simulation/dt)
# print(dt)
print(n_steps)
# Newmark explicity constants
gamma = 0.5
beta = 0.0

# Initial values

# Initial high
h0 = 20.0 * 10**-2 #(m)
# h0 = 0 #(m)
# Initial displacement (u0)
u0 = np.zeros((n_dofs))
n_points = n_dofs
l = np.linspace(0, L, n_points)
u0 = np.array([h0 + x for x in l])
u0 = np.round(u0, 8)
print(u0)
# Initial velocity (v0): 
v0 = np.zeros((n_dofs))

# Initial acceleration (acel0)
acel0 = np.zeros((n_dofs))

# Inital load (p):
# Gravity
g = 9.81 #(m/s)
weight = -rho*L*A*g/n_points
# weight = -2000
p = np.full((n_steps+1, n_dofs),weight)
# print(p)
# Damping
C = np.zeros((n_dofs, n_dofs))
# Initialization of maximum jump u between two linear elements (delta_max)
delta_max = np.zeros((len(materials)*2))

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
    ndofs = max(list(itertools.chain.from_iterable(connect))) + 1
    DofCoord = np.zeros((ndofs,3))
    for el in range(n_el):
        dof = connect[el][0]
        DofCoord[dof,0] = NodeCoord(node_id[el][0])
        dof = connect[el][1]
        DofCoord[dof,0] = NodeCoord(node_id[el][1])
    return DofCoord
