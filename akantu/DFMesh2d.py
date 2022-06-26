import akantu as aka
import numpy as np
import subprocess


# Geometry parameters

# Lenght of the bar (m)
L = 50*10**-3  
x0 = -L/2
xf = L/2
# Number of triangular elements (n_el)
n_el = 2
# Lenght of each linear element (h)
h = L/(n_el/2)
# Cross sectional area (m2)
A = 1*10**-3 

# Material parameters

# Young's module (Pa)
E = 275.0*10**9  
# Density (kg/m3)
rho = 2750.0  
# Limit stress / critical stress (stress_c) (Pa)
stress_c = 300.0*10**6



# Mesh (Triangles elements)

geometry_file = f"""
Point(1) = {{ {x0}, 0, 0, {h} }};
Point(2) = {{ {xf}, 0, 0, {h} }};
Point(3) = {{ {xf}, {h}, 0, {h} }};
Point(4) = {{ {x0}, {h}, 0, {h} }};

Line(1) = {{1,2}};
Line(2) = {{2,3}};
Line(3) = {{3,4}};
Line(4) = {{4,1}};

Line Loop(1) = {{1,2,3,4}};
Plane Surface(1) = {{1}};
Physical Surface(1) = {{1}};
Physical Line("YBlocked") = {{1,3}};
Physical Line("right") = {{2}};
Physical Line("left") = {{4}};

Transfinite Surface {1} Left;
"""

subprocess.run("mkdir LOG", shell=True)
subprocess.run("rm -r paraview_backup", shell=True)
subprocess.run("mv paraview paraview_backup", shell=True)

with open('LOG/bar.geo', 'w') as f:
    f.write(geometry_file)

ret = subprocess.run("gmsh -2 -order 1 -o LOG/bar.msh LOG/bar.geo", shell=True)
if ret.returncode:
    print("FATAL    : Gmsh error")
else:
    print("Info    : Mesh generated")



material_file = f"""
seed = 1.0
model solid_mechanics_model_cohesive [

    material elastic [
        name = linear
        rho = {rho} # Density (kg/m3)
        E = {E}  # Young's module (Pa)
        nu = 0.3  
        finite_deformation = true
    ]

    material cohesive_linear [
        name = insertion
        sigma_c = {stress_c} uniform [-1e6, 1e6] # critical stress (Pa)
        G_c = 100.0 # Fracture energy (N/m)
        beta = 0.
        penalty = 1e11
    ]
]
"""

subprocess.run("mkdir LOG", shell=True)

with open('LOG/material.dat', 'w') as f:
    f.write(material_file)






