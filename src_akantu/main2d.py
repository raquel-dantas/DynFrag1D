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
# Get connectivity list
connect = mesh.getConnectivity(aka._triangle_3)
# Get coordinates 
coords = mesh.getNodes()



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



# Set time increment
dt_crit = model.getStableTimeStep()     # Critical time step (s)
dt = dt_crit*0.1                        # Adopted time step
model.setTimeStep(dt)
time_simulation = 1.0*10**-7            # Total time of simulation (s)
n_steps = int(time_simulation/dt)       # Number of time steps



# Apply Dirichlet BC to block dispacements at y direction on top and botton of the elements
model.applyBC(aka.FixedValue(0., aka._y), 'YBlocked')



# Applied strain rate (s-1)
# strain_rate = 10.0**2
# strain_rate = 10.0**3
strain_rate = 10.0**5
# strain_rate = 10.0**5



# Set velocity applied at the boundaries
# Value of the velocity at the bounary
vel = strain_rate*DFMesh2d.L/2 
# Apply constant velocity at the boundaries
functor_left = DFModel.FixedVelocity(aka._x, -vel)
functor_right = DFModel.FixedVelocity(aka._x, vel)
model.applyBC(functor_left, 'left')
model.applyBC(functor_right, 'right')



# Initial values
n_nodes = mesh.getNbNodes()     # Number of nodes of the mesh
u0 = model.getDisplacement()    # Initial displacement
v0 = model.getVelocity()        # Initial velocity
# Set initial velocity profile
v0[:,0] = np.array([strain_rate * x for x,y in mesh.getNodes()])
# DFPlot2d.PlotByCoord(v0[:,0],mesh,'initial vel')
# Apply initial velocity values
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
Epot = np.zeros(n_steps)            # Potential energy 
Ekin = np.zeros(n_steps)            # Kinetic energy 
Edis = np.zeros(n_steps)            # Dissipated energy 
Erev = np.zeros(n_steps)            # Reversible energy 
Econ = np.zeros(n_steps)            # Contact energy 
Wext = np.zeros(n_steps)            # External work 
work = 0.0                          # External work from the previous time-step 
fp_left = 0.0                       # Equivalent force applied at the left boundary in the previous time-step
fp_right = 0.0                      # Equivalent force applied at the right boundary in the previous time-step
avg_stress = np.zeros(n_steps)      # Average stress in the whole bar


# Initiation of fragmentatio postprocess variables
nfrag = np.zeros(n_steps)                           # Number of fragments  
mean_nelfrag = np.zeros(n_steps)                    # Mean number of elements per fragment

mean_sfrag = np.zeros(n_steps)                      # Mean fragment size
sfrag_all_steps = [[] for _ in range(n_steps)]      # Fragment size for all fragments at all time-steps

mean_velfrag = np.zeros(n_steps)                    # Mean fragments velocity
velfrag_all_steps = [[] for _ in range(n_steps)]    # Velocitis of all fragments at all time-steps

mfrag_all_steps = [[] for _ in range(n_steps)]      # Mass of all fragments ata all time-steps 


# Create files to save outputs
f = str(nfrag)
number_fragments = f
with open('LOG/number_fragments_dynfrag_akantu.txt','w') as f: 
    f.write("--\n")


for n in range(n_steps):

    # Apply velocity at the boundaries 
    functor_left.set_time(dt*n)
    functor_right.set_time(dt*n)
    model.applyBC(functor_left, 'left')
    model.applyBC(functor_right, 'right')

    # VTK files
    model.dump()
    model.dump('cohesive elements')

    # Run simulation
    model.checkCohesiveStress()
    model.solveStep('explicit_lumped')


    u = model.getDisplacement()[:,0]
    v = model.getVelocity()[:,0]
    acel = model.getAcceleration()[:,0]
    # DFPlot2d.PlotByCoord(v0[:,0],mesh,'v')
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


    # Fragmentation data
    fragmentdata = aka.FragmentManager(model)
    fragmentdata.computeAllData()

    # Number of fragments
    nfrag[n] = fragmentdata.getNbFragment()                             # Number of fragments at time-step n   
    mean_nelfrag[n] = np.mean(fragmentdata.getNbElementsPerFragment())  # Mean number of elements per fragment at time-step n 
    
    # Fragments size (assuming uniform mesh)
    sfrag = np.zeros(fragmentdata.getNbFragment())                      
    sfrag = fragmentdata.getNbElementsPerFragment() * DFMesh2d.hun      # Sizes of all fragments 
    sfrag_all_steps[n] = sfrag                                          # Store the fragments sizes of time-step n
    mean_sfrag[n] = mean_nelfrag[n] * DFMesh2d.hun                      # Mean size of fragments 
    
    # Fragments velocities
    velfrag = np.zeros(fragmentdata.getNbFragment())
    velfrag = fragmentdata.getVelocity()                                # Velocities of all fragments 
    velfrag_all_steps[n] = velfrag                                      # Store the velocities of time-step n
    mean_velfrag[n] = np.mean(velfrag)                                  # Mean velocity of fragments
    
    # Fragments mass
    mfrag = np.zeros(fragmentdata.getNbFragment())
    mfrag = fragmentdata.getMass()                                      # Fragments mass of all fragments
    mfrag_all_steps[n] = mfrag                                          # Store the mass of time-step n

    # write outputs in .txt
    f = str(nfrag)
    number_fragments = f
    with open('LOG/number_fragments_dynfrag_akantu.txt','w') as f: 
        f.write("--\n")



# Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPosprocess2d.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps)

# Power [Energy, time] returns the energy difference between consecutive time steps
PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot = DFPosprocess2d.Power(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps)



# Plots

DFPlot2d.PlotAverageStressBar(avg_stress, time_simulation, n_steps)

DFPlot2d.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, time_simulation, n_steps)

DFPlot2d.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot, time_simulation, n_steps)

DFPlot2d.PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot, time_simulation, n_steps)



# Plots fragmentation data

DFPlot2d.PlotNumberFragments(nfrag, time_simulation, n_steps)

DFPlot2d.PlotAvgFragmentSize(mean_sfrag, time_simulation, n_steps)

DFPlot2d.PlotFragmentSizeHistogram(sfrag)



# Save average results for convergence study
f = str(mean_sfrag[n])
average_fragment_size = f
with open('LOG/average_fraglen_dynfrag_akantu.txt','w') as f: 
    f.write(average_fragment_size)
    
f = str(nfrag[n])
number_fragments = f
with open('LOG/number_fragments_dynfrag_akantu.txt','w') as f: 
    f.write(number_fragments)

f = str(Edis[n])
Edis = f
with open('LOG/final_diss_energy_dynfrag_akantu.txt','w') as f: 
    f.write(Edis)