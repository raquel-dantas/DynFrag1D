import akantu as aka
import numpy as np
import subprocess


# Lenght of the bar (m)
L = 50*10**-3  
x0 = -L/2
xf = L/2
# Number of linear elements (n_el)
n_el = 5
# Lenght of each linear element (h)
h = L/n_el


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

subprocess.run("mkdir LOG")

with open('LOG/bar.geo', 'w') as f:
    f.write(geometry_file)

ret = subprocess.run("gmsh -2 -order 1 -o LOG/bar.msh LOG/bar.geo", shell=True)
if ret.returncode:
    print("FATAL    : Gmsh error")
else:
    print("Info    : Mesh generated")
