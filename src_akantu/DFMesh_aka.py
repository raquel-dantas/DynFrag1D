import subprocess
import pickle
import akantu as aka
import numpy as np

import input_files.input_data_aka as inputdata
filepath = inputdata.filepath

# Assign material properties
young_modulus = inputdata.young_modulus
rho = inputdata.density
fracture_energy = inputdata.fracture_energy
stress_limit = inputdata.stress_limit
generate_limit_stress_variation = inputdata.generate_limit_stress_variation
if generate_limit_stress_variation == False:
    stress_limit_file_name = inputdata.stress_limit_file_name

# Assign geometry
area = inputdata.area
bar_length = inputdata.bar_length
x0 = inputdata.x0
xf = inputdata.xf

# Assign load
strain_rate = inputdata.strain_rate
if inputdata.half_bar == True:
    applied_vel = strain_rate * bar_length 
else:
    applied_vel = strain_rate * bar_length * 0.5


# Assign mesh inputs
uniform_mesh = inputdata.uniform_mesh
create_mesh = inputdata.create_mesh
n_elements = inputdata.number_elements

# Assign time of simulation
time_simulation = inputdata.time_simulation



# Mesh

# Compute the size of elements (h) for uniform mesh
h_uniform = bar_length / (n_elements * 0.5) 

geometry_file = f"""
Point(1) = {{ {x0}, 0, 0, {h_uniform} }};
Point(2) = {{ {xf}, 0, 0, {h_uniform} }};
Point(3) = {{ {xf}, {h_uniform}, 0, {h_uniform} }};
Point(4) = {{ {x0}, {h_uniform}, 0, {h_uniform} }};

Line(1) = {{1,2}};
Line(2) = {{2,3}};
Line(3) = {{3,4}};
Line(4) = {{4,1}};

Line Loop(1) = {{1,2,3,4}};
Plane Surface(1) = {{1}};
Physical Surface(1) = {{1}};
Physical Line("Yblocked") = {{1,3}};
Physical Line("right") = {{2}};
Physical Line("left") = {{4}};

Transfinite Surface {1} Left;
"""

subprocess.run("mkdir src_akantu/LOG", shell=True)
subprocess.run("rm -r paraview_backup", shell=True)
subprocess.run("mv paraview paraview_backup", shell=True)

with open('src_akantu/LOG/bar.geo', 'w') as f:
    f.write(geometry_file)

ret = subprocess.run("gmsh -2 -order 1 -o src_akantu/LOG/bar.msh src_akantu/LOG/bar.geo", shell=True)
if ret.returncode:
    print("FATAL    : Gmsh error")
else:
    print("Info    : Mesh generated")

# Contact penalty
alpha = (stress_limit**2 + 4.5 * strain_rate**(2/3) * young_modulus * fracture_energy**(2/3) * rho**(1/3)) / (4.5 * fracture_energy)

material_file = f"""
seed = 1.0
model solid_mechanics_model_cohesive [

    material elastic [
        name = linear
        rho = {rho} # Density (kg/m3)
        E = {young_modulus}  # Young's module (Pa)
        nu = 0.  
        finite_deformation = true
    ]

    material cohesive_linear [
        name = insertion
        sigma_c = {stress_limit} uniform [-1e6, 1e6] # critical stress (Pa)
        G_c = {fracture_energy} # Fracture energy (N/m)
        beta = 0.
        penalty = {alpha}
    ]
]
"""

subprocess.run("mkdir src_akantu/LOG", shell=True)

with open('src_akantu/LOG/material.dat', 'w') as f:
    f.write(material_file)


# Read material file to akantu
aka.parseInput('src_akantu/LOG/material.dat')

# Read mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read('src_akantu/LOG/bar.msh')

# Get number of nodes
n_nodes = mesh.getNbNodes()
# Get connectivity list
connect = mesh.getConnectivity(aka._triangle_3)
# Get coordinates 
node_coord = mesh.getNodes()



if create_mesh == False:
    # import the coordinates from a picke file
    with open(inputdata.mesh_file_name, "rb") as handle:
        node_coord_x = pickle.load(handle)
    for i in range(len(node_coord_x)-2):
        node_coord[i + 4, 0] = node_coord_x[i + 1]
        node_coord[n_elements + 1 - i , 0] = node_coord_x[i + 1]

