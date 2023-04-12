import akantu as aka
import numpy as np
import progressbar
import pickle

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
    facets_coords = mesh_facets.getNodes()

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
            dir_vec = (
                DFMesh.node_coord[conn_facettype[1], :]
                - DFMesh.node_coord[conn_facettype[0], :]
            )
            direction = (dir_vec / np.linalg.norm(dir_vec)).dot(up)
            # If the direction is not 1 it means that is a diagonal facet, then assign False
            if abs(direction) < 0.99:
                check_facettype[el] = False

if DFMesh.generate_limit_stress_variation == False:
    stress_critical = dynfrag.getMaterial(1).getInternalReal("sigma_c")
    stress_critical = stress_critical(aka._segment_2)
    with open(DFMesh.stress_limit_file_name, "rb") as handle:
        random_stress = pickle.load(handle)

    update_stress = []
    stress_value = []
    for i in check_facettype:
        update_stress.extend([i, i])
    for i in random_stress:
        stress_value.extend([i, i])

    i = 0
    deb = 0
    for facet in range(len(stress_critical)):
        if update_stress[facet] == True:
            deb += 1
            stress_critical[facet] = stress_value[i]
            i += 1


# Set time increment
dt_crit = dynfrag.getStableTimeStep()
dt = dt_crit * 0.1
dynfrag.setTimeStep(dt)
n_steps = int(DFMesh.time_simulation / dt)


# Apply Dirichlet BC to block dispacements at y direction on top and botton of the elements
dynfrag.applyBC(aka.FixedValue(0.0, aka._y), "Yblocked")

# Apply constant velocity at the boundaries
left_extremity = DFApplyVel.FixedVelocity(aka._x, -DFMesh.applied_vel)
right_extremity = DFApplyVel.FixedVelocity(aka._x, DFMesh.applied_vel)
dynfrag.applyBC(left_extremity, "left")
dynfrag.applyBC(right_extremity, "right")


# Initial values
u0 = dynfrag.getDisplacement()
v0 = dynfrag.getVelocity()
v0[:, 0] = np.array([DFMesh.strain_rate * x for x, y in DFMesh.mesh.getNodes()])
dynfrag.getVelocity()[:] = v0
data_bc = [0, 0]

# if DFMesh.strain_rate < 5.0 * 10.0**3:
#     u0[:,0] = np.array([0.98 * DFMesh.stress_limit * x / DFMesh.young_modulus for x, y in DFMesh.mesh.getNodes()])
#     dynfrag.getDisplacement()[:] = u0



# Initiation of variables
energy_potential = 0.0
energy_kinetic = 0.0
energy_dissipated = 0.0
energy_reversible = 0.0
energy_contact = 0.0
external_work = 0.0
work_previous_step = 0.0
energies = [
    ["energy potential", energy_potential],
    ["energy kinetic", energy_kinetic],
    ["energy dissipated", energy_dissipated],
    ["energy reversible", energy_reversible],
    ["energy contact", energy_contact],
    ["external work", external_work],
]



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
    left_extremity.set_time(dt * current_time_step)
    right_extremity.set_time(dt * current_time_step)
    dynfrag.applyBC(left_extremity, "left")
    dynfrag.applyBC(right_extremity, "right")
