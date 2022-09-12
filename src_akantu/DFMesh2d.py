import akantu as aka
import numpy as np
import subprocess


# Material
E = 275.0*10**9     # Young's module (Pa)
rho = 2750.0        # Density (kg/m3)
Gc = 100.0          # Fracture energy (N/m)
stress_critical = 300.0*10**6   # Limit stress / critical stress (Pa)
# Geometry
A = 1*10**-3        # Cross sectional area (m2)
L = 1.05*10**-3     # Lenght of the bar (m)
x0 = 0              # Left extremitiy x coordinate / 0-initial
xf = L              # Rigth extremitiy x coordinate / f-final
n_el = 10           # Number of linear elements (n_el)
hun = L/n_el        # Size of the elemenets (h) for a uniform mesh (un) 




# Mesh (Triangles elements)

geometry_file = f"""
Point(1) = {{ {x0}, 0, 0, {hun} }};
Point(2) = {{ {xf}, 0, 0, {hun} }};
Point(3) = {{ {xf}, {hun}, 0, {hun} }};
Point(4) = {{ {x0}, {hun}, 0, {hun} }};

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
        sigma_c = {stress_critical} uniform [-1e6, 1e6] # critical stress (Pa)
        G_c = 100.0 # Fracture energy (N/m)
        beta = 0.
        penalty = 1e17
    ]
]
"""

subprocess.run("mkdir LOG", shell=True)

with open('LOG/material.dat', 'w') as f:
    f.write(material_file)






