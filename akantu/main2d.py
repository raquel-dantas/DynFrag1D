from operator import index
from textwrap import indent
import akantu as aka
import numpy as np
import DFMesh2d
import DFModel
import DFPosprocess2d
import DFPlot2d


# Read material file
aka.parseInput('LOG/material.dat')

# Read mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read('LOG/bar.msh')

# Create model
model = aka.SolidMechanicsModelCohesive(mesh)
model.initFull(_analysis_method=aka._static, _is_extrinsic=True)

# Configure static solver
solver = model.getNonLinearSolver('static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

# Solver (explicit Newmark with lumped mass)
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
print(n_steps)

# Apply Dirichlet BC to block dispacements at y direction on top and botton of the elements
model.applyBC(aka.FixedValue(0., aka._y), 'YBlocked')


# Applied strain rate (s-1)

# strain_rate = 10.0**2
# strain_rate = 10.0**3
# strain_rate = 10.0**4
strain_rate = 10.0**5

# Applied velocity
vel = strain_rate*DFMesh2d.L/2 



# Apply constant velocity at the boundaries
functor_left = DFModel.FixedVelocity(aka._x, -vel)
functor_right = DFModel.FixedVelocity(aka._x, vel)
model.applyBC(functor_left, 'left')
model.applyBC(functor_right, 'right')


# Initial value 

n_nodes = mesh.getNbNodes()
u0 = model.getDisplacement()
v0 = model.getVelocity()

# Initial velocity profile
v0[:,0] = np.array([strain_rate * x for x,y in mesh.getNodes()])

# DFPlot2d.PlotByCoord(v0[:,0],mesh,'initial vel')

# Apply initial conditions:
model.getVelocity()[:] = v0



# VTK plot
model.setBaseName('bar')
model.addDumpFieldVector('displacement')
model.addDumpFieldVector('velocity')
model.addDumpField('strain')
model.addDumpField('stress')
model.addDumpField('blocked_dofs')
model.addDumpField('material_index')

# VTK plot for Cohesive model
model.setBaseNameToDumper('cohesive elements', 'cohesive')
model.addDumpFieldVectorToDumper('cohesive elements', 'displacement')
model.addDumpFieldToDumper('cohesive elements', 'damage')
model.addDumpFieldVectorToDumper('cohesive elements', 'tractions')
model.addDumpFieldVectorToDumper('cohesive elements', 'opening')







# Initiation of variables
Epot = np.zeros(n_steps)
Ekin = np.zeros(n_steps)
Edis = np.zeros(n_steps)
Erev = np.zeros(n_steps)
Econ = np.zeros(n_steps)
Wext = np.zeros(n_steps)
work = 0.0
fp_left = 0.0
fp_right = 0.0
avg_stress = np.zeros(n_steps)

# vp = np.zeros()




for n in range(n_steps):


    # Apply velocity at the boundaries 
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
    v = model.getVelocity()[:,0]
    # coord = mesh.getNodes()
    # DFPlot2d.PlotByCoord(v0[:,0],mesh,'v')
    a = model.getAcceleration()[:,0]
    fint = model.getInternalForce()[:,0]
    stress = model.getMaterial(0).getStress(aka._triangle_3)
    stress_xx = stress[:,0]
    avg_stress[n] = np.mean(stress_xx)

 
    # Energy balance
    Epot[n] = model.getEnergy('potential')
    Ekin[n] = model.getEnergy('kinetic')
    Edis[n] = model.getEnergy('dissipated')
    Erev[n] = model.getEnergy('reversible')
    Econ[n] = model.getEnergy('cohesive contact')
    Wext[n], fp_left, fp_right = DFPosprocess2d.ExternalWork(mesh, fint, fp_left, fp_right, work, vel, dt)
    work = Wext[n]






    # Fragmentation
    
# print(Epot)
# print(Ekin)
# print(Wext)

# Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPosprocess2d.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps)

# Power [Energy, time] returns the energy difference between consecutive time steps
PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot = DFPosprocess2d.Power(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps)

# DFPlot2d.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot, time_simulation, n_steps)
# DFPlot2d.PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot, time_simulation, n_steps)

# DFPlot2d.PlotAverageStressBar(avg_stress, time_simulation, n_steps)

