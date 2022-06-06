from operator import index
from textwrap import indent
import akantu as aka
import numpy as np


# Read material file
aka.parseInput('akantu/material.dat')

# Read mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read('LOG/bar.msh')

# Create model
model = aka.SolidMechanicsModel(mesh)
model.

# Initial values
model.initFull(_analysis_method=aka._explicit_lumped_mass)


# Applied strain rate (s-1)
# strain_rate = 10.0**2
# strain_rate = 10.0**3
strain_rate = 10.0**4
# strain_rate = 10.0**5


# Configure solver

critical_time_step = model.getStableTimeStep()

n_nodes = mesh.getNbNodes()

u0 = model.getDisplacement()
v0 = model.getVelocity()
# model.setDisplacement(u0)

# for n in range():node = 0; node < mesh.getNbNodes(); ++node) {
#     disp(node, 0) = 0.1;
#     velo(node, 1) = 1.;
# }

v0[:,0] = np.array([strain_rate * x for x,y in mesh.getNodes()])

# Young's module (Pa)
E = 275.0*10**9 

material = model.getMaterial("linear")
stress_c = 300.0*10**6 
index = model.
# Initial displacement (u0)
if strain_rate < 5.0 * 10.0**3:
    u0[:,0] = np.array([0.98*stress_c*x / E for x,y in mesh.getNodes()])



