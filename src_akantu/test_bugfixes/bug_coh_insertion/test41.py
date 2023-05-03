import akantu as aka
import numpy as np
import subprocess

# Material parameters
E = 275.0*10**9         # Young's module (Pa)
rho = 27.500            # Density (kg/m3)
stress_c = 300.0*10**8  # Limit stress / critical stress (stress_c) (Pa)

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
aka.parseInput('LOG/material.dat')


# Geometry parameters
L = 50*10**-3   # Lenght of the bar (m)
x0 = -L/2
xf = L/2
n_el = 10        # Number of triangular elements (n_el)
h = L/(n_el/2)  # Lenght of each linear element (h)

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

dt_crit = model.getStableTimeStep() # Critical time step
dt = dt_crit*0.1                    # Adopted time step
model.setTimeStep(dt)
time_simulation = 6.0*10**-6        # Total time of simulation (s)
n_steps = int(time_simulation/dt)   # Number of time steps

# Apply Dirichlet BC to block displacements at y direction on top and bottom
model.applyBC(aka.FixedValue(0., aka._y), 'YBlocked')

# Constant velocity boundary condition
# Applied strain rate (s-1)
strain_rate = 10.0**5
# Applied velocity at the boundary
vel = strain_rate*L/2 
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

# VTK plot setup
model.setBaseName('bar')
model.addDumpFieldVector('displacement')
# VTK plot setup for Cohesive model
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
