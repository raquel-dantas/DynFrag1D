import numpy as np
import itertools
import copy
import pickle

import input_files.input_data as inputdata


# Import or set inputs

# Assign material properties
young_modulus = inputdata.young_modulus
rho = inputdata.density
fracture_energy = inputdata.fracture_energy
stress_limit = inputdata.stress_limit
generate_limit_stress_variation = inputdata.generate_limit_stress_variation

# Assign geometry
area = inputdata.area
bar_length = inputdata.bar_length
x0 = inputdata.x0
xf = inputdata.xf

# Assign load
strain_rate = inputdata.strain_rate
applied_vel = strain_rate * bar_length * 0.5
# applied_vel = strain_rate * bar_length 

# Assign method
use_cohesive_elements = inputdata.use_cohesive_elements
use_lipfield = inputdata.use_lipfield

# Assign mesh inputs
uniform_mesh = inputdata.uniform_mesh
create_mesh = inputdata.create_mesh
n_elements = inputdata.number_elements

# Assign time of simulation
time_simulation = inputdata.time_simulation


# Assign material to elements according to the following material convention
# Material id convention:
# 0 : line element
# 1 : interface element
# 2 : Support left node
# 3 : Support right node
# 4 : Velocity applied left node
# 5 : Velocity applied right node
materials = [0] * n_elements
materials.append(4)
# materials.append(2)
materials.append(5)


# BC's dictionary
bc_dict = {
    2: (0, "dirichlet"),
    4: (-applied_vel, "velocity"),
    5: (+applied_vel, "velocity"),
}


# Mesh

# Compute the size of elements (h) for uniform mesh
h_uniform = bar_length / n_elements

node_id = [[i, i + 1] for i in range(n_elements)] # node_id : returns the global node id for all elements
# Append the nodes with apllied BCs to match material convention
node_id.append([0])  # Applied velocity at id = 0
node_id.append([n_elements])  # Applied velocity at id = nb_elements

connect = copy.deepcopy(node_id) # connect: returns the connectivity of all elements
n_dofs = max(list(itertools.chain.from_iterable(connect))) + 1
n_points = n_dofs # n_points: number of geometric points in the mesh

if create_mesh:
    # uniform_coord: returns the points coodinate for a uniform mesh
    uniform_coord = np.linspace(x0, xf, n_points)
    node_coord = uniform_coord
    if uniform_mesh == False:
        node_coord = np.array(
            [
                x + np.random.uniform(low=-0.4, high=0.4) * h_uniform
                for x in uniform_coord
            ]
        )
        node_coord[0] = x0
        node_coord[n_elements] = xf

else:
    # import the coordinates from a picke file
    with open("src/input_files/mesh.pickle", "rb") as handle:
        node_coord = pickle.load(handle)

# intpoit_coord: returns the coordinates of integration points
intpoint_coord = [node_coord[i + 1] - 0.5*node_coord[i] for i in range(n_elements)]


# Interface parameters
# crack_limit: limit crack oppening
crack_limit = (2.0 * fracture_energy) / stress_limit
# contact_penalty_limit: contact penalty computed based in the stress_limit
contact_penalty_limit = (
    stress_limit**2
    + 4.5
    * strain_rate ** (2 / 3)
    * young_modulus
    * fracture_energy ** (2 / 3)
    * rho ** (1 / 3)
) / (4.5 * fracture_energy)

# stress_critical: returns a random distribution of limit stress at the interface elements 
if generate_limit_stress_variation == True:
    stress_critical = np.random.uniform(
        low=stress_limit - 1.0 * 10**6,
        high=stress_limit + 1.0 * 10**6,
        size=(n_elements),
    )
else:
    with open("src/input_files/random_stress_critical.pickle", "rb") as handle:
        stress_critical = pickle.load(handle)

# crack_critical: returns the critical crack opening for the values in stress_critical
crack_critical = [(2.0 * fracture_energy) / stress_critical[i] for i in range(n_elements)]

# contact_penalty: returns the penalty based on the values in stress_critical
contact_penalty = [
    (stress_critical[i] ** 2 + 4.5 * strain_rate ** (2 / 3) * young_modulus * fracture_energy ** (2 / 3) * rho ** (1 / 3))
    / (4.5 * fracture_energy)
    for i in range(n_elements - 1)
]

# jump_max: returns the maximum displacement jump historically archieve between two consecutive linear elements
jump_max = np.zeros(len(materials) * 2)


# Set time 

# Compute the critical time step
if uniform_mesh:
    smallest_element = h_uniform
else:
    smallest_element = (0.2 * h_uniform)  
dt_crit = smallest_element / ((young_modulus / rho) ** 0.5)

dt = dt_crit * 0.5  # Adopted time step (s)

n_steps = int(time_simulation / dt)  # Number of time-steps
time_peakstress = stress_limit / (young_modulus * strain_rate)  # Time peak stress
nstep_peak = int(time_peakstress / dt)  # Time-step to peak stress


# Set initial values

v0 = np.array([strain_rate * i for i in node_coord])  # Initial velocity 
acel0 = np.zeros(n_dofs) # Initial acceleration 
p = np.zeros((n_steps + 1, n_dofs))  # External forces 
C = np.zeros((n_dofs, n_dofs))  # Damping 
if use_lipfield == True:
    d0 = np.zeros(n_elements)  # Initial damage

# Initial displacement
# Apply initial displacement neq zero to save computational time during pre-crack phase for low strain-rates
if strain_rate < 5.0 * 10.0**3:
    u0 = np.array([0.98 * stress_limit * i / young_modulus for i in node_coord])
else:
    u0 = np.zeros(n_dofs)




def getElemLength(elem_index):
    """Returns the element length (hel)."""

    hel = node_coord[elem_index + 1] - node_coord[elem_index]

    return hel


def getEl(connect, dof_id):
    el_id = 0
    for el in connect:
        local_dof = 0
        for dof in el:
            if dof == dof_id:
                return el_id, local_dof
            local_dof = local_dof + 1
        el_id = el_id + 1



def listDofCoord():
    """Returns a list of coordinates by dof.\n
    dof_coord[i] returns [x,y,z] of dof i."""

    ndofs = max(list(itertools.chain.from_iterable(connect))) + 1
    dof_coord = np.zeros((ndofs, 3))
    for el in range(n_elements):
        dof = connect[el][0]
        dof_coord[dof, 0] = node_coord[node_id[el][0]]
        dof = connect[el][1]
        dof_coord[dof, 0] = node_coord[node_id[el][1]]
    return dof_coord
