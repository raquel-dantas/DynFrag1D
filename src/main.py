import numpy as np
import pickle
import progressbar
import time

import DFMesh
import DFFem
import src.DFPostProcess as DFPostProcess
import DFNewmark
import DFInterface
import DFPlot
import DFFragmentation


def runSimulation(strain_rate):

    bar = progressbar.ProgressBar(
        maxval=50,
        widgets=[progressbar.Bar("=", "[", "]"), " ", progressbar.Percentage()],
    )
    bar.start()

    # Initiation of variables
    u = DFMesh.u0
    v = DFMesh.v0
    acel = DFMesh.acel0
    energy_potential = np.zeros(DFMesh.n_steps)
    energy_kinetic = np.zeros(DFMesh.n_steps)
    energy_dissipated = np.zeros(DFMesh.n_steps)
    external_work = np.zeros(DFMesh.n_steps)
    work_previous_step = 0.0
    els_step = DFMesh.n_el
    stress_evolution = np.zeros((2 * len(DFMesh.materials), DFMesh.n_steps))
    avg_stress_bar = np.zeros(DFMesh.n_steps)
    uprevious_bc_left = np.array([0, 0])
    uprevious_bc_right = np.array([0, 0])
    n_fragments = np.zeros(DFMesh.n_steps)
    avg_frag_size = np.zeros(DFMesh.n_steps)
    data_histogram_frag_size = []

    if DFMesh.use_lipfield == True:
        d = DFMesh.d0
        dprevious_bc_left = 0.0
        dprevious_bc_right = 0.0

    if DFMesh.use_1d_cohesive_elements == True:
        d = None
        energy_reversible = np.zeros(DFMesh.n_steps)
        energy_contact = np.zeros(DFMesh.n_steps)

    for n in range(DFMesh.n_steps):

        progress = int(bar.maxval * float(n / DFMesh.n_steps))
        bar.update(progress)

        strain, stress, average_stress_neighbors = DFPostProcess.postProcess(u, d)
        stress_evolution = DFPostProcess.logStress(n, stress_evolution, stress)
        avg_stress_bar[n] = DFPostProcess.stressBar(stress)
        M, F = DFFem.globalSystem()

        # Plots at each time step
        # DFPlot.PlotByDOF(v)
        # DFPlot.PlotByElement(stress)
        # DFPlot.PlotVTK('animation/test',n,u,stress)

        if DFMesh.use_1d_cohesive_elements == True:
            d = [
                DFInterface.getDamageParameter(el)
                for el in range(len(DFMesh.materials))
            ]
            # DFPlot.PlotByInterface(D)
            (
                energy_potential[n],
                energy_kinetic[n],
                energy_dissipated[n],
                energy_reversible[n],
                energy_contact[n],
                external_work[n],
            ) = DFPostProcess.computeEnergiesCZM(
                uprevious_bc_left, uprevious_bc_right, u, v, stress, work_previous_step
            )

        if DFMesh.use_lipfield == True:
            (
                energy_potential[n],
                energy_kinetic[n],
                energy_dissipated[n],
                external_work[n],
            ) = DFPostProcess.computeEnergiesLipfield(
                uprevious_bc_left,
                uprevious_bc_right,
                u,
                v,
                d,
                dprevious_bc_left,
                dprevious_bc_right,
                work_previous_step,
            )

        work_previous_step = external_work[n]

        # Fragmentation data
        n_fragments[n] = DFFragmentation.getNumberFragments(d)
        frag_lengths, avg_frag_size[n] = DFFragmentation.getFragmentSizes(d)
        data_histogram_frag_size = DFFragmentation.getFragmentSizeHistogramData(
            frag_lengths, 10
        )

        # up_bc is the previous displacement vector for the local dofs in the boundary elements (left and right)
        uprevious_bc_left = np.array([0, 0])
        uprevious_bc_right = np.array([0, 0])
        for bc in range(len(DFMesh.materials)):
            if DFMesh.materials[bc] == 4 or DFMesh.materials[bc] == 5:
                if DFMesh.materials[bc] == 4:
                    el_bc = 0
                    uprevious_bc_left = np.array(
                        [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
                    )
                else:
                    el_bc = DFMesh.n_el - 1
                    uprevious_bc_right = np.array(
                        [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
                    )

        # Time integration
        u, v, acel = DFNewmark.explicitScheme(M, u, v, acel, d, F, DFMesh.dt)

        if DFMesh.use_1d_cohesive_elements == True:
            for el in range(DFMesh.n_el - 1):
                if average_stress_neighbors[el] > DFMesh.stress_critical[el]:
                    # Fracture happens: creat new interface element
                    u, v, acel = DFInterface.insertInterface(el, el + 1, u, v, acel)
                    els_step = els_step + 1

    bar.finish()
    print("\n")

    DFPlot.plotAverageStressBar(avg_stress_bar)
    DFPlot.plotNumberFragments(n_fragments)
    DFPlot.plotFragmentSizeHistogram(frag_lengths, 10)

    # Energy balance
    if DFMesh.use_1d_cohesive_elements == True:
        (
            var_energy_potential,
            var_energy_kinetic,
            var_energy_dissipated,
            var_energy_reversible,
            var_energy_contact,
            var_external_work,
            var_energy_total,
        ) = DFPostProcess.computeVarEnergiesCZM(
            energy_potential,
            energy_kinetic,
            energy_dissipated,
            energy_reversible,
            energy_contact,
            external_work,
        )

        (
            power_potential,
            power_kinetic,
            power_dissipated,
            power_reversible,
            power_contact,
            power_external_work,
            power_total,
        ) = DFPostProcess.computePowerCZM(
            energy_potential,
            energy_kinetic,
            energy_dissipated,
            energy_reversible,
            energy_contact,
            external_work,
        )

        DFPlot.plotEnergiesCZM(
            energy_potential,
            energy_kinetic,
            energy_dissipated,
            energy_reversible,
            energy_contact,
            external_work,
        )
        DFPlot.plotVarEnergiesCZM(
            var_energy_potential,
            var_energy_kinetic,
            var_energy_dissipated,
            var_energy_reversible,
            var_energy_contact,
            var_external_work,
            var_energy_total,
        )
        DFPlot.plotPowerCZM(
            power_potential,
            power_kinetic,
            power_dissipated,
            power_reversible,
            power_contact,
            power_external_work,
            power_total,
        )

    if DFMesh.use_lipfield == True:

        var_energy_potential,
        var_energy_kinetic,
        var_energy_dissipated,
        var_external_work,
        var_energy_total = DFPostProcess.computeVarEnergiesLipfield(
            energy_potential, energy_kinetic, energy_dissipated, external_work
        )

        power_potential,
        power_kinetic,
        power_dissipated,
        power_external_work,
        power_total = DFPostProcess.computePowerLipfield(energy_potential, energy_kinetic, energy_dissipated, external_work)

        DFPlot.plotEnergiesLipfield(
            energy_potential, energy_kinetic, energy_dissipated, external_work
        )

        DFPlot.plotVarEnergiesLipfield(
            var_energy_potential,
            var_energy_kinetic,
            var_energy_dissipated,
            var_external_work,
            var_energy_total,
        )
        DFPlot.plotPowerLipfield(
            power_potential,
            power_kinetic,
            power_dissipated,
            power_external_work,
            power_total,
        )

    # Save results
    if DFMesh.use_1d_cohesive_elements == True:

        DFPlot.saveResultsCZM(n_fragments)
        DFPlot.saveResultsCZM(avg_frag_size)
        DFPlot.saveResultsCZM(data_histogram_frag_size)

        DFPlot.saveResultsCZM(avg_stress_bar)

        DFPlot.saveResultsCZM(energy_potential)
        DFPlot.saveResultsCZM(energy_kinetic)
        DFPlot.saveResultsCZM(energy_dissipated)
        DFPlot.saveResultsCZM(energy_reversible)
        DFPlot.saveResultsCZM(energy_contact)
        DFPlot.saveResultsCZM(external_work)

        DFPlot.saveResultsCZM(var_energy_potential)
        DFPlot.saveResultsCZM(var_energy_kinetic)
        DFPlot.saveResultsCZM(var_energy_dissipated)
        DFPlot.saveResultsCZM(var_energy_reversible)
        DFPlot.saveResultsCZM(var_energy_contact)
        DFPlot.saveResultsCZM(var_external_work)

        DFPlot.saveResultsCZM(power_potential)
        DFPlot.saveResultsCZM(power_kinetic)
        DFPlot.saveResultsCZM(power_dissipated)
        DFPlot.saveResultsCZM(power_reversible)
        DFPlot.saveResultsCZM(power_contact)
        DFPlot.saveResultsCZM(power_external_work)


    if DFMesh.use_lipfield == True:

        DFPlot.saveResultsLipfield(n_fragments)
        DFPlot.saveResultsLipfield(avg_frag_size)
        DFPlot.saveResultsLipfield(data_histogram_frag_size)

        DFPlot.saveResultsLipfield(avg_stress_bar)

        DFPlot.saveResultsLipfield(energy_potential)
        DFPlot.saveResultsLipfield(energy_kinetic)
        DFPlot.saveResultsLipfield(energy_dissipated)
        DFPlot.saveResultsLipfield(external_work)

        DFPlot.saveResultsLipfield(var_energy_potential)
        DFPlot.saveResultsLipfield(var_energy_kinetic)
        DFPlot.saveResultsLipfield(var_energy_dissipated)
        DFPlot.saveResultsLipfield(var_external_work)

        DFPlot.saveResultsLipfield(power_potential)
        DFPlot.saveResultsLipfield(power_kinetic)
        DFPlot.saveResultsLipfield(power_dissipated)
        DFPlot.saveResultsLipfield(power_external_work)

    


if __name__ == "__main__":
    start_time = time.time()
    runSimulation(DFMesh.strain_rate)
    total_time = time.time() - start_time
    print("--- %s seconds ---" % (time.time() - start_time))
    if DFMesh.use_1d_cohesive_elements == True:
        DFPlot.saveResultsCZM(total_time)
    if DFMesh.use_lipfield == True:
        DFPlot.saveResultsLipfield(total_time)
