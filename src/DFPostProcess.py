from matplotlib.pyplot import connect
import numpy as np

import DFMesh
import DFInterface
import DFFem
import DFDiffuseDamage


def computeStress(u, d):
    """Returns the stress vector for all elements, and the average stress vector between two consecutives line elements.\n
    Arguments:\n
    u -- displacemnt vector; \n
    d -- damage field."""

    # Total number of elements (linear + cohesive)
    n_elements = len(DFMesh.materials)
    # Initiation of variables
    strain = np.zeros(n_elements)
    stress = np.zeros(n_elements)

    for el in range(n_elements):

        if DFMesh.use_cohesive_elements == True:
            if DFMesh.materials[el] == 0:
                strain[el] = (
                    u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
                ) / DFMesh.getElemLength(el)
                stress[el] = DFMesh.young_modulus * strain[el]

            elif DFMesh.materials[el] == 1:
                # jump_u returns the jump in the displacement between two consecutive linear elements
                jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
                # Stress retuns the stress value at each cohesive el
                stress[el] = DFInterface.stressCohesiveLaw(jump_u, el)

        if DFMesh.use_lipfield == True:
            if DFMesh.materials[el] == 0:
                g = (1.0 - d[el]) ** 2
                strain[el] = (
                    u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
                ) / DFMesh.getElemLength(el)
                stress[el] = DFMesh.young_modulus * g * strain[el]

    # average_stress returns the average stress between two consecutive linear elements
    average_stress = (
        lambda el: (stress[el] + stress[el + 1]) / 2.0
        if DFMesh.connect[el][1] == DFMesh.connect[el + 1][0]
        else 0
    )
    average_stress_neighbors = [
        average_stress(el) for el in range(DFMesh.n_elements - 1)
    ]

    return stress, average_stress_neighbors


def logStress(time_step, stress_evolution, current_stress):
    """Returns a matrix that contains the stress for all elements at all time steps (cols).\n
    Arguments:\n
    time_step: current time step of the analysis;\n
    stress_evolution: stress evolution array;\n
    current_stress: the stress vector of the current time step."""

    n_elements = len(DFMesh.materials)
    for el in range(n_elements):
        stress_evolution[el, time_step] = current_stress[el]

    return stress_evolution


def stressBar(current_stress):
    """Returns the average stress at the whole bar at each time step.\n
    Arguments:\n
    current_stress: the stress vector of the current time step."""

    return sum(current_stress) / len(DFMesh.materials)


def saveResultsAtBC(u, d):

    for bc in range(len(DFMesh.materials)):
        if DFMesh.materials[bc] == 4 or DFMesh.materials[bc] == 5:
            if DFMesh.materials[bc] == 4:
                el_bc = 0
                uprevious_bc_left = np.array(
                    [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
                )
                if DFMesh.use_lipfield == True:
                    dprevious_bc_left = d[el_bc]
            else:
                el_bc = DFMesh.n_elements - 1
                uprevious_bc_right = np.array(
                    [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
                )
                if DFMesh.use_lipfield == True:
                    dprevious_bc_right = d[el_bc]

    uprevious_bc = [uprevious_bc_left, uprevious_bc_right]
    if DFMesh.use_cohesive_elements == True:
        dprevious_bc = None
        data_bc = [uprevious_bc, dprevious_bc]
    if DFMesh.use_lipfield == True:
        dprevious_bc = [dprevious_bc_left, dprevious_bc_right]
        data_bc = [uprevious_bc, dprevious_bc]

    return data_bc


def computeEnergiesCZM(u, v, stress, data_bc, work_previous_step):
    """Returns potential, kinetic, dissipated, reversible, contact and external energies.\n
    Arguments:\n
    uprevious_bc_left -- displacement from previous time step at left boundary element;\n
    uprevious_bc_right -- displacement from previous time step at right boundary element;\n
    u -- displacements vector;\n
    v -- velocities vector;\n
    work_previous_step -- external energy from the previous time step."""

    energy_potential = 0.0
    energy_kinetic = 0.0
    energy_dissipated = 0.0
    energy_reversible = 0.0
    energy_contact = 0.0
    external_work = work_previous_step

    for el in range(len(DFMesh.materials)):

        if DFMesh.materials[el] == 0:

            # Get local displacement and velocity
            u_local = np.array([u[DFMesh.connect[el][0]], u[DFMesh.connect[el][1]]])
            v_local = np.array([v[DFMesh.connect[el][0]], v[DFMesh.connect[el][1]]])

            energy_potential += (
                0.5 * np.dot(np.matmul(DFFem.k_elem, u_local), u_local)
            ) / DFMesh.getElemLength(el)

            energy_kinetic += (
                0.5 * np.dot(np.matmul(DFFem.m_elem, v_local), v_local)
            ) * DFMesh.getElemLength(el)

        if DFMesh.materials[el] == 1:

            energy_dissipated += (
                0.5 * DFMesh.stress_limit * DFMesh.jump_max[el] * DFMesh.area
            )

            # Get jump in displacment and current stress
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            stress_cohesive = stress[el]

            if jump_u >= 0:
                energy_reversible += 0.5 * stress_cohesive * jump_u * DFMesh.area
            else:
                i = DFMesh.connect[el][0] - 1
                energy_contact += (
                    0.5 * DFMesh.contact_penalty[i] * jump_u**2 * DFMesh.area
                )

        if DFMesh.materials[el] == 4 or DFMesh.materials[el] == 5:
            # vel_extremity is the velocity applied on the extremity
            # e_lbc is the element index of the applied velocity
            # uprvious_local is the local displacement of el_bc from the previous time step

            uprevious_bc = data_bc[0]
            if DFMesh.materials[el] == 4:
                vel_extremity = np.array([-DFMesh.applied_vel, 0])
                el_bc = 0
                uprevious_local = uprevious_bc[0]

            else:
                vel_extremity = np.array([0, DFMesh.applied_vel])
                el_bc = DFMesh.n_elements - 1
                uprevious_local = uprevious_bc[1]

            hel_bc = DFMesh.getElemLength(el_bc)
            u_local = np.array(
                [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
            )

            fint_current_step = np.matmul(DFFem.k_elem, u_local) / hel_bc
            fint_previous_step = np.matmul(DFFem.k_elem, uprevious_local) / hel_bc
            # The reaction force is taken as an average between the internal force in the current and previous time-set
            freact = (fint_current_step + fint_previous_step) * 0.5
            # Stress on the boundary
            stress_boundary = freact / DFMesh.area

            # Work is power integrated in time
            work = np.dot(stress_boundary, vel_extremity) * DFMesh.dt * DFMesh.area
            external_work += work

    energies_step = [
        ["energy potential", energy_potential],
        ["energy kinetic", energy_kinetic],
        ["energy dissipated", energy_dissipated],
        ["energy reversible", energy_reversible],
        ["energy contact", energy_contact],
        ["external work", external_work],
    ]

    return energies_step


def computeEnergiesLipfield(u, v, d, data_bc, work_previous_step):
    """Returns potential, kinetic, dissipated, reversible, contact and external energies.\n
    Arguments:\n
    uprevious_bc_left -- displacement from previous time step at left boundary element;\n
    uprevious_bc_right -- displacement from previous time step at right boundary element;\n
    u -- displacements vector;\n
    v -- velocities vector;\n
    d -- damage field;
    dprevious_bc_left -- damage from previous time step at left boundary element;\n
    dprevious_bc_right -- damage from previous time step at right boundary element;\n
    work_previous_step -- external energy from the previous time step."""

    energy_potential = 0.0
    energy_kinetic = 0.0
    energy_dissipated = 0.0
    external_work = work_previous_step

    for el in range(len(DFMesh.materials)):

        if DFMesh.materials[el] == 0:

            # Get local displacement and velocity
            g = (1.0 - d[el]) ** 2
            u_local = np.array([u[DFMesh.connect[el][0]], u[DFMesh.connect[el][1]]])
            v_local = np.array([v[DFMesh.connect[el][0]], v[DFMesh.connect[el][1]]])

            energy_potential += (
                0.5 * g * np.dot(np.matmul(DFFem.k_elem, u_local), u_local)
            ) / DFMesh.getElemLength(el)

            energy_kinetic += (
                0.5 * np.dot(np.matmul(DFFem.m_elem, v_local), v_local)
            ) * DFMesh.getElemLength(el)

            energy_dissipated += (
                DFDiffuseDamage.Yc[el]
                * DFDiffuseDamage.h(DFDiffuseDamage.lamb_const[el], d[el])
                * DFMesh.area
                * DFMesh.getElemLength(el)
            )

        if DFMesh.materials[el] == 4 or DFMesh.materials[el] == 5:
            # vel_extremity is the velocity applied on the extremity
            # e_lbc is the element index of the applied velocity
            # uprevious_local is the local displacement of el_bc from the previous time step

            uprevious_bc = data_bc[0]
            dprevious_bc = data_bc[1]
            if DFMesh.materials[el] == 4:
                vel_extremity = np.array([-DFMesh.applied_vel, 0])
                el_bc = 0
                g = (1.0 - d[el_bc]) ** 2
                uprevious_local = uprevious_bc[0]
                gprevious = (1.0 - dprevious_bc[0]) ** 2

            else:
                vel_extremity = np.array([0, DFMesh.applied_vel])
                el_bc = DFMesh.n_elements - 1
                g = (1.0 - d[el_bc]) ** 2
                uprevious_local = uprevious_bc[1]
                gprevious = (1.0 - dprevious_bc[1]) ** 2

            hel_bc = DFMesh.getElemLength(el_bc)
            u_local = np.array(
                [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
            )
            fint_current_step = np.matmul(g * DFFem.k_elem, u_local) / hel_bc
            fint_previous_step = (
                np.matmul(gprevious * DFFem.k_elem, uprevious_local) / hel_bc
            )
            # The reaction force is taken as an average between the internal force in the current and previous time-set
            freact = (fint_current_step + fint_previous_step) * 0.5
            # Stress on the boundary
            stress_boundary = freact / DFMesh.area

            # Work is power integrated in time
            work = np.dot(stress_boundary, vel_extremity) * DFMesh.dt * DFMesh.area
            external_work += work

    energies_step = [
        ["energy potential", energy_potential],
        ["energy kinetic", energy_kinetic],
        ["energy dissipated", energy_dissipated],
        ["external work", external_work],
    ]

    return energies_step


def updateEnergies(energies, n, u, v, d, stress, data_bc, work_previous_step):

    if DFMesh.use_cohesive_elements == True:
        energies_step = computeEnergiesCZM(u, v, stress, data_bc, work_previous_step)
    if DFMesh.use_lipfield == True:
        energies_step = computeEnergiesLipfield(u, v, d, data_bc, work_previous_step)

    for i in range(len(energies)):
        energies[i][1][n] = energies_step[i][1]
    return energies


def getEnergy(energies, energy_name):

    for i in range(len(energies)):
        if energies[i][0] == energy_name:
            energy = energies[i][1]
            return energy
    # else:
    #     raise Exception("energy name don't match!")


def computeVarEnergiesCZM(energies):
    """Returns the variation of energies between the current time step and the time step 0."""

    var_energy_potential = np.zeros(DFMesh.n_steps)
    var_energy_kinetic = np.zeros(DFMesh.n_steps)
    var_energy_dissipated = np.zeros(DFMesh.n_steps)
    var_energy_contact = np.zeros(DFMesh.n_steps)
    var_energy_reversible = np.zeros(DFMesh.n_steps)
    var_energy_total = np.zeros(DFMesh.n_steps)
    var_external_work = np.zeros(DFMesh.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    energy_reversible = getEnergy(energies, "energy reversible")
    energy_contact = getEnergy(energies, "energy contact")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFMesh.n_steps):
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


def computeVarEnergiesLipfield(energies):
    """Returns the variation of energies between the current time step and the time step 0."""

    var_energy_potential = np.zeros(DFMesh.n_steps)
    var_energy_kinetic = np.zeros(DFMesh.n_steps)
    var_energy_dissipated = np.zeros(DFMesh.n_steps)
    var_external_work = np.zeros(DFMesh.n_steps)
    var_energy_total = np.zeros(DFMesh.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFMesh.n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (
            var_energy_potential[n] + var_energy_kinetic[n] + var_energy_dissipated[n]
        )

    var_energies = [
        ["var energy potential", var_energy_potential],
        ["var energy kinetic", var_energy_kinetic],
        ["var energy dissipated", var_energy_dissipated],
        ["var external work", var_external_work],
        ["var energy total", var_energy_total],
    ]

    return var_energies


def computeVariationEnergy(energies):

    if DFMesh.use_cohesive_elements == True:
        var_energies = computeVarEnergiesCZM(energies)
    if DFMesh.use_lipfield == True:
        var_energies = computeVarEnergiesLipfield(energies)
    return var_energies


def computePowerCZM(energies):
    """Returns the variation of energies between two consecutives time steps."""

    power_potential = np.zeros(DFMesh.n_steps)
    power_kinetic = np.zeros(DFMesh.n_steps)
    power_dissipated = np.zeros(DFMesh.n_steps)
    power_reversible = np.zeros(DFMesh.n_steps)
    power_contact = np.zeros(DFMesh.n_steps)
    power_external_work = np.zeros(DFMesh.n_steps)
    power_total = np.zeros(DFMesh.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    energy_reversible = getEnergy(energies, "energy reversible")
    energy_contact = getEnergy(energies, "energy contact")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFMesh.n_steps):
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


def computePowerLipfield(energies):
    """Returns the variation of energies between two consecutives time steps."""

    power_potential = np.zeros(DFMesh.n_steps)
    power_kinetic = np.zeros(DFMesh.n_steps)
    power_dissipated = np.zeros(DFMesh.n_steps)
    power_external_work = np.zeros(DFMesh.n_steps)
    power_total = np.zeros(DFMesh.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFMesh.n_steps):
        power_potential[n] = energy_potential[n] - energy_potential[n - 1]
        power_kinetic[n] = energy_kinetic[n] - energy_kinetic[n - 1]
        power_dissipated[n] = energy_dissipated[n] - energy_dissipated[n - 1]
        power_external_work[n] = external_work[n] - external_work[n - 1]
        power_total[n] = power_external_work[n] - (
            power_potential[n] + power_kinetic[n] + power_dissipated[n]
        )

    power = [
        ["power potential", power_potential],
        ["power kinetic", power_kinetic],
        ["power dissipated", power_dissipated],
        ["power external work", power_external_work],
        ["power total", power_total],
    ]

    return power


def computePower(energies):

    if DFMesh.use_cohesive_elements == True:
        power = computePowerCZM(energies)
    if DFMesh.use_lipfield == True:
        power = computePowerLipfield(energies)
    return power
