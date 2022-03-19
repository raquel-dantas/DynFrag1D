import numpy as np
import itertools
import copy 

# Mesh
L = 50*10**-3  # (m)
n_el = 3
h = L/n_el

# Material id convention: 
# 0(line elemnt)/1(interface element)/2(BC left node)/3(BC right node)/ 4(Velocity right node)/5(Velocity right node)
materials = [0] * n_el
materials.append(2)
materials.append(4)
# node_id[elem_index][local_node], returns the global node id
node_id = []
for i in range(n_el):
    node_id.append([i, i+1])
# Node BC
node_id.append([0]) # Support
node_id.append([n_el]) # Applied velocity at right boundary 

# Connect[el][j] returns the global index of local dof 'j' of the element 'el'
connect = copy.deepcopy(node_id) 

# Applied strain rate and veloctities
strain_rate = 10.0**1  # (s-1)
vel = strain_rate*L

# BC
bc_dict = {
    2: (0, "dirichlet"),
    4: (vel, "velocity"),
    5: (-vel, "velocity")
}
n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1



# Input parameters

# Material parameters
E = 275.0*10**9  # (Pa)

A = 1*10**-3  # (unit area)
Gc = 100.0  # (N/m)
stress_c = 300.0*10**5  # (Pa)
rho = 75.0*10**3  # (kg/m3)
delta_c = (2.0*Gc)/stress_c

# Time steps
n_steps = 85



# Initial displacement
u0 = np.zeros((n_dofs))

# Initial velocity (v0): velocity profile (vel) is a function v(x)
n_points = n_dofs
l = np.linspace(0, L, n_points)
v0 = np.array([vel/L*x for x in l])
v0 = np.round(v0, 8)


# Initial acceleration (acel0)
acel0 = np.zeros((n_dofs))
# Load (p)
p = np.zeros((n_steps+1, n_dofs))

# Initialization of maximum jump u between two linear elements
delta_max = np.zeros((len(materials)*2))

# Time integration
dt_crit = h/((E/rho)**0.5)
dt = dt_crit*0.1  # (s)
gamma = 0.5
beta = 0.0
C = np.zeros((n_dofs, n_dofs))
