import akantu as aka
import numpy as np
import subprocess

# Material parameters
# Young's module (Pa)
E = 275.0*10**9  
# Density (kg/m3)
rho = 2750.0  
# Limit stress / critical stress (stress_c) (Pa)
stress_c = 300.0*10**6

# Material file
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

# Read material
aka.parseInput('LOG/material.dat')




# Geometry parameters
# Lenght of the bar (m)
L = 50*10**-3  
x0 = -L/2
xf = L/2
# Number of triangular elements (n_el)
n_el = 2
# Lenght of each linear element (h)
h = L/(n_el/2)

# Mesh file (Triangles elements)
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

# Read mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read('LOG/bar.msh')
# Get connectivity list
connect = mesh.getConnectivity(aka._triangle_3)
# Get coordinates 
coords = mesh.getNodes()




# Create model
model = aka.SolidMechanicsModelCohesive(mesh)
model.initFull(_analysis_method=aka._static, _is_extrinsic=True)

# Configure solver
solver = model.getNonLinearSolver('static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)
model.initNewSolver(aka._explicit_lumped_mass)

# Dynamic insertion of cohesive elements
model.updateAutomaticInsertion()

# Critical time step
dt_crit = model.getStableTimeStep()
# Adopted time step
dt = dt_crit*0.1
# Total time of simulation (s)
time_simulation = 6.0*10**-6
# Number of time steps (n_steps)
n_steps = int(time_simulation/dt)


# Apply Dirichlet BC to block displacements at y direction on top and bottom
model.applyBC(aka.FixedValue(0., aka._y), 'YBlocked')

# Applied strain rate (s-1)
strain_rate = 10.0**5
# Applied velocity at the boundary
vel = strain_rate*L/2 


# Apply Velocity at the extremities
class FixedVelocity (aka.DirichletFunctor):
    """Fixed velocity at the boundaries."""

    def __init__(self, axis, vel):
        super().__init__(axis)
        self.axis = axis
        self.time = 0
        self.vel = vel
    
    def set_time(self, t):
        self.time = t
    
    def __call__(self, node, flags, disp, coord):
        flags[int(self.axis)] = True
        disp[int(self.axis)] = self.vel*self.time

functor_left = FixedVelocity(aka._x, -vel)
functor_right = FixedVelocity(aka._x, vel)
model.applyBC(functor_left, 'left')
model.applyBC(functor_right, 'right')


# VTK plot
model.setBaseName('bar')
model.addDumpFieldVector('displacement')
# VTK plot for Cohesive model
model.setBaseNameToDumper('cohesive elements', 'cohesive')
model.addDumpFieldVectorToDumper('cohesive elements', 'displacement')


for n in range(n_steps):

    # Apply velocity at the extremities
    functor_left.set_time(dt*n)
    functor_right.set_time(dt*n)
    model.applyBC(functor_left, 'left')
    model.applyBC(functor_right, 'right')

    model.dump()
    model.dump('cohesive elements')

    # Run simulation
    model.checkCohesiveStress()
    model.solveStep('explicit_lumped')

    # Outputs
    u = model.getDisplacement()[:,0]

