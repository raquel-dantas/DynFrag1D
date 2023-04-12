import numpy as np
import akantu as aka

import DFModel_aka as DFModel
import DFMesh_aka as DFMesh


def getDamageParameter():

    d = DFModel.dynfrag.getMaterial(1).getInternalReal("damage")
    d = d(aka._cohesive_2d_4)
    mat = DFModel.dynfrag.getMaterial(1)
    coh_id = DFModel.dynfrag.getMaterial("insertion").getElementFilter()(
        aka._cohesive_2d_4
    )
    # model.getFEEngine().computeIntegrationPointsCoordinate(quad_coords, model.getMaterial("cohesive material").getElementFilter())
    return d


# d
# (facet, localNode) -> damage

# def getDamagePerFacet(cohID, localNodeID):
#     return d[cohID + localNodeID]

# cohID -> facetID

# coord -> ... -> damage

# nodeindex -> damage

# facetIndex -> globalNode (feito)
# globalNode -> coord (feito)

# facetIndex -> localNode


# localNode -> coord     (feito)


def computeEnergies(work_previous_step, fint_current_step):

    energy_potential = DFModel.dynfrag.getEnergy("potential")
    energy_kinetic = DFModel.dynfrag.getEnergy("kinetic")
    energy_dissipated = DFModel.dynfrag.getEnergy("dissipated")
    energy_reversible = DFModel.dynfrag.getEnergy("reversible")
    energy_contact = DFModel.dynfrag.getEnergy("cohesive contact")

    # External work
    nodes_bc_left = DFMesh.mesh.getElementGroup("left").getNodeGroup().getNodes()
    nodes_bc_right = DFMesh.mesh.getElementGroup("right").getNodeGroup().getNodes()
    # Internal for at the current time step
    fint_current_step_bcleft = -np.sum(fint_current_step[nodes_bc_left])
    fint_current_step_bcright = -np.sum(fint_current_step[nodes_bc_right])
    # The reaction force is taken as an average between the internal force in the current and previous time-set
    freact_previous_step_bcleft = DFModel.data_bc[0]
    freact_previous_step_bcright = DFModel.data_bc[1]

    freact_bcleft = (fint_current_step_bcleft + freact_previous_step_bcleft) * 0.5
    freact_bcright = (fint_current_step_bcright + freact_previous_step_bcright) * 0.5
    external_work = (
        work_previous_step
        + (freact_bcleft * -DFMesh.applied_vel + freact_bcright * DFMesh.applied_vel)
        * DFModel.dt
    )

    DFModel.data_bc = [freact_bcleft, freact_bcright]

    energies_step = [
        ["energy potential", energy_potential],
        ["energy kinetic", energy_kinetic],
        ["energy dissipated", energy_dissipated],
        ["energy reversible", energy_reversible],
        ["energy contact", energy_contact],
        ["external work", external_work],
    ]

    return energies_step


def updateEnergies(energies, n, work_previous_step, fint_current_step):
    energies_step = computeEnergies(work_previous_step, fint_current_step)
    for i in range(len(energies)):
        energies[i][1][n] = energies_step[i][1]
    return energies


def getEnergy(energies, energy_name):
    for i in range(len(energies)):
        if energies[i][0] == energy_name:
            energy = energies[i][1]
            return energy


def computeVariationEnergy(energies):
    """Returns the variation of energies between the current time step and the time step 0."""

    var_energy_potential = np.zeros(DFModel.n_steps)
    var_energy_kinetic = np.zeros(DFModel.n_steps)
    var_energy_dissipated = np.zeros(DFModel.n_steps)
    var_energy_contact = np.zeros(DFModel.n_steps)
    var_energy_reversible = np.zeros(DFModel.n_steps)
    var_energy_total = np.zeros(DFModel.n_steps)
    var_external_work = np.zeros(DFModel.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    energy_reversible = getEnergy(energies, "energy reversible")
    energy_contact = getEnergy(energies, "energy contact")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFModel.n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_energy_reversible[n] = energy_reversible[n] - energy_reversible[0]
        var_energy_contact[n] = energy_contact[n] - energy_contact[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (
            var_energy_potential[n]
            + var_energy_kinetic[n]
            + var_energy_dissipated[n]
            + var_energy_reversible[n]
            + var_energy_contact[n]
        )

    var_energies = [
        ["var energy potential", var_energy_potential],
        ["var energy kinetic", var_energy_kinetic],
        ["var energy dissipated", var_energy_dissipated],
        ["var energy reversible", var_energy_reversible],
        ["var energy contact", var_energy_contact],
        ["var external work", var_external_work],
        ["var energy total", var_energy_total],
    ]

    return var_energies


def computePower(energies):
    """Returns the variation of energies between two consecutives time steps."""

    power_potential = np.zeros(DFModel.n_steps)
    power_kinetic = np.zeros(DFModel.n_steps)
    power_dissipated = np.zeros(DFModel.n_steps)
    power_reversible = np.zeros(DFModel.n_steps)
    power_contact = np.zeros(DFModel.n_steps)
    power_external_work = np.zeros(DFModel.n_steps)
    power_total = np.zeros(DFModel.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    energy_reversible = getEnergy(energies, "energy reversible")
    energy_contact = getEnergy(energies, "energy contact")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFModel.n_steps):
        power_potential[n] = energy_potential[n] - energy_potential[n - 1]
        power_kinetic[n] = energy_kinetic[n] - energy_kinetic[n - 1]
        power_dissipated[n] = energy_dissipated[n] - energy_dissipated[n - 1]
        power_reversible[n] = energy_reversible[n] - energy_reversible[n - 1]
        power_contact[n] = energy_contact[n] - energy_contact[n - 1]
        power_external_work[n] = external_work[n] - external_work[n - 1]
        power_total[n] = power_external_work[n] - (
            power_potential[n]
            + power_kinetic[n]
            + power_dissipated[n]
            + power_reversible[n]
            + power_contact[n]
        )

    power = [
        ["power potential", power_potential],
        ["power kinetic", power_kinetic],
        ["power dissipated", power_dissipated],
        ["power reversible", power_reversible],
        ["power contact", power_contact],
        ["power external work", power_external_work],
        ["power total", power_total],
    ]

    return power
