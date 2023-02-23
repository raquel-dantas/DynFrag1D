from matplotlib.pyplot import connect
import numpy as np

import DFMesh
import DFInterface
import DFFem
import DFDiffuseDamage


def postProcess(u, d):
    """Returns the strain for linear elements, the stress vector for all elements, and the average stress vector between two consecutives line elements.\n
    Arguments:\n
    u -- displacemnt vector."""

    # Total number of elements (linear + cohesive)
    numel = len(DFMesh.materials)
    # Initiation of variables
    strain = np.zeros(numel)
    stress = np.zeros(numel)

    for el in range(numel):

        if DFMesh.use_lipfield == True:
            if DFMesh.materials[el] == 0:
                g = (1.0 - d[el]) ** 2
                strain[el] = (
                    u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
                ) / DFMesh.ElemLength(el)
                stress[el] = DFMesh.E * g * strain[el]

        else:
            if DFMesh.materials[el] == 0:
                strain[el] = (
                    u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
                ) / DFMesh.ElemLength(el)
                stress[el] = DFMesh.E * strain[el]

            elif DFMesh.materials[el] == 1:
                # jump_u returns the jump in the displacement between two consecutive linear elements
                jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
                # Stress retuns the stress value at each cohesive el
                stress[el] = DFInterface.CohesiveLaw(jump_u, el)

    # average_stress returns the average stress between two consecutive linear elements
    average = (
        lambda el: (stress[el] + stress[el + 1]) / 2.0
        if DFMesh.connect[el][1] == DFMesh.connect[el + 1][0]
        else 0
    )
    average_stress = [average(el) for el in range(DFMesh.n_el - 1)]

    return strain, stress, average_stress


def logStress(time_step, stress_evl, current_stress):
    """Returns a matrix that contains the stress for all elements (lin) at all time steps (cols).\n
    Arguments:\n
    time_step: current time step of the analysis;\n
    stress_evl: current stress evolution matrix;\n
    current_stress: the stress vector of the current time step."""

    numel = len(DFMesh.materials)
    for el in range(numel):
        stress_evl[el, time_step] = current_stress[el]

    return stress_evl


def stressBar(current_stress, els_step):
    """Returns the average stress at the whole bar at each time step.\n
    Arguments:\n
    current_stress: the stress vector of the current time step.\n
    els_step: number of elements (linear + cohesive) at the current time step."""

    sum_stress = 0.0
    numel = len(DFMesh.materials)
    for el in range(numel):
        sum_stress += current_stress[el]
    avg_stress_bar = sum_stress / els_step

    return avg_stress_bar


def computeEnergiesCZM(uprevious_bc_left, uprevious_bc_right, u, v, stress, work_previous_step):
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
            ) / DFMesh.ElemLength(el)

            energy_kinetic += (
                0.5 * np.dot(np.matmul(DFFem.m_elem, v_local), v_local)
            ) * DFMesh.ElemLength(el)


        if DFMesh.materials[el] == 1:

            energy_dissipated += 0.5 * DFMesh.stress_critical * DFMesh.delta_max[el] * DFMesh.area

            # Get jump in displacment and current stress
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            stress_cohesive = stress[el]

            if jump_u >= 0:
                energy_reversible += 0.5 * stress_cohesive * jump_u * DFMesh.area
            else:
                i = DFMesh.connect[el][0] - 1
                energy_contact += 0.5 * DFMesh.alpha[i] * jump_u**2 * DFMesh.area


        if DFMesh.materials[el] == 4 or DFMesh.materials[el] == 5:
            # vel_extremity is the velocity applied on the extremity
            # e_lbc is the element index of the applied velocity
            # uprvious_local is the local displacement of el_bc from the previous time step

            if DFMesh.materials[el] == 4:
                vel_extremity = np.array([-DFMesh.vel, 0])
                el_bc = 0
                uprevious_local = uprevious_bc_left

            else:
                vel_extremity = np.array([0, DFMesh.vel])
                el_bc = DFMesh.n_elements - 1
                uprevious_local = uprevious_bc_right

            hel_bc = DFMesh.getElemLength(el_bc)
            u_local = np.array([u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]])

            fint_current_step = np.matmul(DFFem.k_elem, u_local) / hel_bc
            fint_previous_step = np.matmul(DFFem.k_elem, uprevious_local) / hel_bc
            # The reaction force is taken as an average between the internal force in the current and previous time-set
            freact = (fint_current_step + fint_previous_step) * 0.5
            # Stress on the boundary
            stress_boundary = freact / DFMesh.area

            # Work is power integrated in time
            work = np.dot(stress_boundary, vel_extremity) * DFMesh.dt * DFMesh.area
            external_work += work

    return energy_potential, energy_kinetic, energy_dissipated, energy_reversible, energy_contact, external_work


def computeEnergiesLipfield(uprevious_bc_left, uprevious_bc_right, u, v, d, dprevious_bc_left, dprevious_bc_right, work_previous_step):
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
    external_work = work_previous_step


    for el in range(len(DFMesh.materials)):

        if DFMesh.materials[el] == 0:

            # Get local displacement and velocity 
            g = (1.0 - d[el]) ** 2
            u_local = np.array([u[DFMesh.connect[el][0]], u[DFMesh.connect[el][1]]])
            v_local = np.array([v[DFMesh.connect[el][0]], v[DFMesh.connect[el][1]]])

            energy_potential += (
                0.5 * g * np.dot(np.matmul(DFFem.k_elem, u_local), u_local)
            ) / DFMesh.ElemLength(el)

            energy_kinetic += (
                0.5 * np.dot(np.matmul(DFFem.m_elem, v_local), v_local)
            ) * DFMesh.ElemLength(el)

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

            if DFMesh.materials[el] == 4:
                vel_extremity = np.array([-DFMesh.vel, 0])
                el_bc = 0
                g = (1.0 - d[el_bc]) ** 2
                uprevious_local = uprevious_bc_left
                gprevious = (1.0 - dprevious_bc_left) ** 2

            else:
                vel_extremity = np.array([0, DFMesh.vel])
                el_bc = DFMesh.n_elements - 1
                g = (1.0 - d[el_bc]) ** 2
                uprevious_local = uprevious_bc_right
                gprevious = (1.0 - dprevious_bc_right) ** 2

            hel_bc = DFMesh.getElemLength(el_bc)
            u_local = np.array([u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]])
            fint_current_step = np.matmul(g * DFFem.k_elem, u_local) / hel_bc
            fint_previous_step = np.matmul(gprevious * DFFem.k_elem, uprevious_local) / hel_bc
            # The reaction force is taken as an average between the internal force in the current and previous time-set
            freact = (fint_current_step + fint_previous_step) * 0.5
            # Stress on the boundary
            stress_boundary = freact / DFMesh.area

            # Work is power integrated in time
            work = np.dot(stress_boundary, vel_extremity) * DFMesh.dt * DFMesh.are
            external_work += work

    return energy_potential, energy_kinetic, energy_dissipated, external_work




def computeVarEnergyCZM(energy_potential, energy_kinetic, energy_dissipated, Erev, energy_contact, external_work):
    """Returns the variation of energies between the current time step and the time step 0."""

    var_energy_potential = np.zeros(DFMesh.n_steps)
    var_energy_kinetic = np.zeros(DFMesh.n_steps)
    var_energy_dissipated = np.zeros(DFMesh.n_steps)
    var_energy_contact = np.zeros(DFMesh.n_steps)
    varErev = np.zeros(DFMesh.n_steps)
    varEtot = np.zeros(DFMesh.n_steps)
    var_external_work = np.zeros(DFMesh.n_steps)
    
    for n in range(1, DFMesh.n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        varErev[n] = Erev[n] - Erev[0]
        var_energy_contact[n] = energy_contact[n] - energy_contact[0]
        var_external_work[n] = external_work[n] - external_work[0]
        varEtot[n] = var_external_work[n] - (
            var_energy_potential[n] + var_energy_kinetic[n] + var_energy_dissipated[n] + varErev[n] + var_energy_contact[n]
        )

    return var_energy_potential, var_energy_kinetic,  var_energy_dissipated, varErev, var_energy_contact, var_external_work, varEtot


def computeVarEnergyLipfield(energy_potential, energy_kinetic, energy_dissipated, external_work):
    """Returns the variation of energies between the current time step and the time step 0."""

    var_energy_potential = np.zeros(DFMesh.n_steps)
    var_energy_kinetic = np.zeros(DFMesh.n_steps)
    var_energy_dissipated = np.zeros(DFMesh.n_steps)
    var_external_work = np.zeros(DFMesh.n_steps)
    var_energy_total= np.zeros(DFMesh.n_steps)
    
    for n in range(1, DFMesh.n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (var_energy_potential[n] + var_energy_kinetic[n] + var_energy_dissipated[n])

    return var_energy_potential, var_energy_kinetic, var_energy_dissipated, var_external_work, var_energy_total


def PowerCZM(Epot, Ekin, Edis, Erev, Econ, Wext):
    """Returns the variation of energies between two consecutives time steps."""

    PEkin = np.zeros((DFMesh.n_steps))
    PEpot = np.zeros((DFMesh.n_steps))
    PEdis = np.zeros((DFMesh.n_steps))
    PWext = np.zeros((DFMesh.n_steps))
    PErev = np.zeros((DFMesh.n_steps))
    PEcon = np.zeros((DFMesh.n_steps))
    PEtot = np.zeros((DFMesh.n_steps))
    for n in range(1, DFMesh.n_steps):
        PEpot[n] = Epot[n] - Epot[n - 1]
        PEkin[n] = Ekin[n] - Ekin[n - 1]
        PEdis[n] = Edis[n] - Edis[n - 1]
        PErev[n] = Erev[n] - Erev[n - 1]
        PEcon[n] = Econ[n] - Econ[n - 1]
        PWext[n] = Wext[n] - Wext[n - 1]
        PEtot[n] = PWext[n] - (PEpot[n] + PEkin[n] + PEdis[n] + PErev[n] + PEcon[n])

    return PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot


def PowerLipfield(Epot, Ekin, Edis, Wext):
    """Returns the variation of energies between two consecutives time steps."""

    PEkin = np.zeros((DFMesh.n_steps))
    PEpot = np.zeros((DFMesh.n_steps))
    PEdis = np.zeros((DFMesh.n_steps))
    PWext = np.zeros((DFMesh.n_steps))
    PEtot = np.zeros((DFMesh.n_steps))
    for n in range(1, DFMesh.n_steps):
        PEpot[n] = Epot[n] - Epot[n - 1]
        PEkin[n] = Ekin[n] - Ekin[n - 1]
        PEdis[n] = Edis[n] - Edis[n - 1]
        PWext[n] = Wext[n] - Wext[n - 1]
        PEtot[n] = PWext[n] - (PEpot[n] + PEkin[n] + PEdis[n])

    return PEkin, PEpot, PEdis, PWext, PEtot