import akantu as aka
import numpy as np
import pickle
from matplotlib import pyplot as plt
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
# Add solver 
model.initNewSolver(aka._explicit_lumped_mass)


# Insertion of cohesive elements
# Use automatic insertion
# model.updateAutomaticInsertion()

# Configure the insertion

# Facets mesh and connectivity
mesh_facets = mesh.getMeshFacets()
connect_facets = mesh_facets.getConnectivities()

# Get the cohesive inserter and check facets
cohesive_inserter = model.getElementInserter()
check_facets = cohesive_inserter.getCheckFacets()
up = np.array([0.0, 1.0])
# Loop by all facet types used in the simulation 
for facet_type in connect_facets.elementTypes(dim=(spatial_dimension - 1)):
    conn_facettype = connect_facets(facet_type)
    check_facettype = check_facets(facet_type)
    for el, conn_facettype in enumerate(conn_facettype):
        # Check the direction of the vector 
        dir_vec = coords[conn_facettype[1], :] - coords[conn_facettype[0], :]
        direction = (dir_vec / np.linalg.norm(dir_vec)).dot(up)
        # If the direction is not 1 it means that is a diagonal facet, then assign False
        if abs(direction) < 0.9:
            check_facettype[el] = False



# Set time increment
dt_crit = model.getStableTimeStep()     # Critical time step (s)
dt = dt_crit*0.1                        # Adopted time step
model.setTimeStep(dt)
n_steps = int(DFMesh2d.time_simulation/dt)       # Number of time steps


# Apply Dirichlet BC to block dispacements at y direction on top and botton of the elements
model.applyBC(aka.FixedValue(0., aka._y), 'YBlocked')



# Set velocity applied at the boundaries
# Value of the velocity at the bounary
vel = DFMesh2d.strain_rate * DFMesh2d.L/2 
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
v0[:,0] = np.array([DFMesh2d.strain_rate * x for x,y in mesh.getNodes()])
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
nfrag = np.zeros(n_steps)           # Number of fragments  
mean_nelfrag = np.zeros(n_steps)    # Mean number of elements per fragment
sfrag = DFMesh2d.L                  # Mean size of fragment
avg_sfrag = np.zeros(n_steps)      # Mean fragment size
datahist = []                       # Data of histogram of fragment size
mean_velfrag = np.zeros(n_steps)    # Mean fragments velocity





for n in range(n_steps):

    # Apply velocity at the boundaries 
    functor_left.set_time(dt*n)
    functor_right.set_time(dt*n)
    model.applyBC(functor_left, 'left')
    model.applyBC(functor_right, 'right')

    # VTK files
    if n%100 == 0:
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

    # coh = model.getMaterial(1)
    # d = model.getMaterial(1).getInternalReal('damage')
    # d = d(aka._cohesive_2d_4)
    # coh_id = model.getMaterial('insertion').getElementFilter()(aka._cohesive_2d_4)

    # Fragmentation data
    fragmentdata = aka.FragmentManager(model)
    fragmentdata.computeAllData()

    # Number of fragments
    nfrag[n] = fragmentdata.getNbFragment()                             # Number of fragments at time-step n   
    mean_nelfrag[n] = np.mean(fragmentdata.getNbElementsPerFragment())  # Mean number of elements per fragment at time-step n 
    
    # Fragments size (assuming uniform mesh)
    sfrag = np.zeros(fragmentdata.getNbFragment())                      
    sfrag = fragmentdata.getNbElementsPerFragment() * DFMesh2d.hun      # Sizes of all fragments 
    avg_sfrag[n] = (mean_nelfrag[n]%2 + (mean_nelfrag[n] - mean_nelfrag[n]%2)/2 ) * DFMesh2d.hun                      # Mean size of fragments 
    datahist = plt.hist(sfrag,10)


    # Fragments velocities
    velfrag = np.zeros(fragmentdata.getNbFragment())
    velfrag = fragmentdata.getVelocity()                # Velocities of all fragments 
    mean_velfrag[n] = np.mean(velfrag)                  # Mean velocity of fragments
    
    # Fragments mass
    mfrag = np.zeros(fragmentdata.getNbFragment())
    mfrag = fragmentdata.getMass()                      # Fragments mass of all fragments




# Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPosprocess2d.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps)

# Power [Energy, time] returns the energy difference between consecutive time steps
PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot = DFPosprocess2d.Power(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps)



# Plots 
DFPlot2d.PlotAverageStressBar(avg_stress, DFMesh2d.time_simulation, n_steps)
DFPlot2d.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, DFMesh2d.time_simulation, n_steps)
DFPlot2d.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot, DFMesh2d.time_simulation, n_steps)
DFPlot2d.PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot, DFMesh2d.time_simulation, n_steps)
DFPlot2d.PlotNumberFragments(nfrag, DFMesh2d.time_simulation, n_steps)
DFPlot2d.PlotFragmentSizeHistogram(sfrag)


# Save results 
# Number of fragments
with open('LOG/src_akantu_number_fragments.pickle', 'wb') as handle:
    pickle.dump(nfrag, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Average fragment size
with open('LOG/src_akantu_average_fragment_size.pickle', 'wb') as handle:
    pickle.dump(avg_sfrag, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Histogram fragment size
with open('LOG/src_akantu_datahist_fragment_size.pickle', 'wb') as handle:
    pickle.dump(datahist, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Average stress for the bar 
with open('LOG/src_akantu_avg_stress.pickle', 'wb') as handle:
    pickle.dump(avg_stress, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Energy
with open('LOG/src_akantu_epot.pickle', 'wb') as handle:
    pickle.dump(Epot, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_ekin.pickle', 'wb') as handle:
    pickle.dump(Ekin, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_edis.pickle', 'wb') as handle:
    pickle.dump(Edis, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_erev.pickle', 'wb') as handle:
    pickle.dump(Erev, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_econ.pickle', 'wb') as handle:
    pickle.dump(Econ, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_wext.pickle', 'wb') as handle:
    pickle.dump(Wext, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Variation of energy
with open('LOG/src_akantu_var_epot.pickle', 'wb') as handle:
    pickle.dump(varEpot, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_var_ekin.pickle', 'wb') as handle:
    pickle.dump(varEkin, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_var_edis.pickle', 'wb') as handle:
    pickle.dump(varEdis, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_var_erev.pickle', 'wb') as handle:
    pickle.dump(varErev, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_var_econ.pickle', 'wb') as handle:
    pickle.dump(varEcon, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('LOG/src_akantu_var_wext.pickle', 'wb') as handle:
    pickle.dump(varWext, handle, protocol=pickle.HIGHEST_PROTOCOL)









    
