import akantu as aka
import numpy as np
import progressbar

import DFMesh_aka as DFMesh
import DFApplyVel_aka as DFApplyVel


# Create model a model
dynfrag = aka.SolidMechanicsModelCohesive(DFMesh.mesh)
dynfrag.initFull(_analysis_method=aka._static, _is_extrinsic=True)
# Configure static solver
solver = dynfrag.getNonLinearSolver("static")
solver.set("max_iterations", 100)
solver.set("threshold", 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)
# Add solver
dynfrag.initNewSolver(aka._explicit_lumped_mass)


# Configure cohesive element insertion
use_automatic_insertion = False
if use_automatic_insertion == True:
    dynfrag.updateAutomaticInsertion()
else:
    mesh_facets = DFMesh.mesh.getMeshFacets()
    connect_facets = mesh_facets.getConnectivities()

    # Get the cohesive inserter and check facets
    cohesive_inserter = dynfrag.getElementInserter()
    check_facets = cohesive_inserter.getCheckFacets()
    up = np.array([0.0, 1.0])
    # Loop by all facet types used in the simulation 
    for facet_type in connect_facets.elementTypes(dim=(DFMesh.spatial_dimension - 1)):
        conn_facettype = connect_facets(facet_type)
        check_facettype = check_facets(facet_type)
        for el, conn_facettype in enumerate(conn_facettype):
            # Check the direction of the vector 
            dir_vec = DFMesh.node_coord[conn_facettype[1], :] - DFMesh.node_coord[conn_facettype[0], :]
            direction = (dir_vec / np.linalg.norm(dir_vec)).dot(up)
            # If the direction is not 1 it means that is a diagonal facet, then assign False
            if abs(direction) < 0.9:
                check_facettype[el] = False



# Set time increment
dt_crit = dynfrag.getStableTimeStep()     # Critical time step (s)
dt = dt_crit*0.1                        # Adopted time step
dynfrag.setTimeStep(dt)
n_steps = int(DFMesh.time_simulation/dt)       # Number of time steps


# Apply Dirichlet BC to block dispacements at y direction on top and botton of the elements
dynfrag.applyBC(aka.FixedValue(0., aka._y), 'Yblocked')

# Apply constant velocity at the boundaries
left_extremity = DFApplyVel.FixedVelocity(aka._x, - DFMesh.applied_vel)
right_extremity = DFApplyVel.FixedVelocity(aka._x, DFMesh.applied_vel)
dynfrag.applyBC(left_extremity, 'left')
dynfrag.applyBC(right_extremity, 'right')


# Initial values
u0 = dynfrag.getDisplacement()
v0 = dynfrag.getVelocity() 
v0[:,0] = np.array([DFMesh.strain_rate * x for x,y in DFMesh.mesh.getNodes()])
dynfrag.getVelocity()[:] = v0

# Initiation of variables
energy_potential = np.zeros(n_steps)
energy_kinetic = np.zeros(n_steps)
energy_dissipated = np.zeros(n_steps)
energy_reversible = np.zeros(n_steps)
energy_contact = np.zeros(n_steps)
external_work = np.zeros(n_steps)
work_previous_step = 0.0
energies = [
        ["energy potential", energy_potential],
        ["energy kinetic", energy_kinetic],
        ["energy dissipated", energy_dissipated],
        ["energy reversible", energy_reversible],
        ["energy contact", energy_contact],
        ["external work", external_work],
    ]

avg_stress_bar = np.zeros(n_steps)

u_all_steps = [u0]
damage_all_steps = []
fraglen_all_steps = []
stress_all_steps = []


n_fragments = np.zeros(n_steps)
elements_per_frag = []
avg_frag_size = np.zeros(n_steps)
data_histogram_frag_size = []
data_bc = [0,0]


def initProgressBar():
    bar = progressbar.ProgressBar(
        maxval=50,
        widgets=[progressbar.Bar("=", "[", "]"), " ", progressbar.Percentage()],
    )
    bar.start()
    return bar


def endProgressBar(bar):
    bar.finish()
    print("\n")


def updateProgressBar(n, bar):
    progress = int(bar.maxval * float(n / n_steps))
    bar.update(progress)


def applyVel(current_time_step):
    left_extremity.set_time(dt*current_time_step)
    right_extremity.set_time(dt*current_time_step)
    dynfrag.applyBC(left_extremity, 'left')
    dynfrag.applyBC(right_extremity, 'right')