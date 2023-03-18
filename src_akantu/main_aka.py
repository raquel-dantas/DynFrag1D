import akantu as aka
import numpy as np
import time

import DFMesh_aka as DFMesh
import DFModel_aka as DFModel
import DFPostProcess_aka as DFPostProcess
import DFPlot_aka as DFPlot


def runSimulation(strain_rate):

    bar = DFModel.initProgressBar()

    # Initiation of variables

    avg_stress_bar = DFModel.avg_stress_bar
    work_previous_step = DFModel.work_previous_step
    energies = DFModel.energies

    n_fragments = DFModel.n_fragments
    elements_per_fragment = DFModel.elements_per_frag
    avg_frag_size = DFModel.avg_frag_size
    fraglen_all_steps = DFModel.fraglen_all_steps
    damage_all_steps = DFModel.damage_all_steps

    # VTK plot
    DFPlot.addPlotVtk()

    for n in range(DFModel.n_steps):

        DFModel.updateProgressBar(n, bar)

        DFPlot.addVtkFiles(n)

        # Apply velocity at the boundaries
        DFModel.applyVel(n)

        # Run simulation
        DFModel.dynfrag.checkCohesiveStress()
        DFModel.dynfrag.solveStep("explicit_lumped")

        u = DFModel.dynfrag.getDisplacement()[:, 0]
        v = DFModel.dynfrag.getVelocity()[:, 0]
        acel = DFModel.dynfrag.getAcceleration()[:, 0]

        fint = DFModel.dynfrag.getInternalForce()[:, 0]
        stress = DFModel.dynfrag.getMaterial(0).getStress(aka._triangle_3)
        stress_xx = stress[:, 0]
        avg_stress_bar[n] = np.mean(stress_xx)

        # Energy balance
        energies = DFPostProcess.updateEnergies(energies, n, work_previous_step, fint)
        work_previous_step = DFPostProcess.getEnergy(energies, "external work")[n]

        d = DFPostProcess.getDamageParameter()
        damage_all_steps.append(d)

        # Fragmentation data
        fragment_data = aka.FragmentManager(DFModel.dynfrag)
        fragment_data.computeAllData()

        # Number of fragments
        n_fragments[n] = fragment_data.getNbFragment()
        elements_per_fragment.append([fragment_data.getNbElementsPerFragment()])
        mean_els_per_fragment = np.mean(fragment_data.getNbElementsPerFragment())

        # Fragments size (assuming uniform mesh)
        frag_lengths = np.zeros(fragment_data.getNbFragment())
        frag_lengths = fragment_data.getNbElementsPerFragment() * DFMesh.h_uniform
        fraglen_all_steps.append(frag_lengths)
        # avg_sfrag = (mean_els_per_fragment%2 + (mean_els_per_fragment - mean_els_per_fragment%2)*0.5 ) * DFMesh.h_uniform

    DFModel.endProgressBar(bar)

    # Energy balance
    var_energies = DFPostProcess.computeVariationEnergy(energies)
    power = DFPostProcess.computePower(energies)

    # Plots
    # DFPlot.plotDamage(d)
    DFPlot.plotAverageStressBar(avg_stress_bar)
    DFPlot.plotNumberFragments(n_fragments)
    DFPlot.plotFragmentSizeHistogram(frag_lengths, 10)
    DFPlot.plotEnergies(energies)
    DFPlot.plotVarEnergies(var_energies)
    DFPlot.plotPower(power)
    # Save results
    DFPlot.saveResults(u)
    DFPlot.saveResults(v)
    DFPlot.saveResults(acel)
    DFPlot.saveResults(d)
    DFPlot.saveResults(damage_all_steps)
    DFPlot.saveResults(fraglen_all_steps)
    DFPlot.saveResults(n_fragments)
    DFPlot.saveResults(avg_frag_size)
    DFPlot.saveResults(avg_stress_bar)
    DFPlot.saveResults(energies)
    DFPlot.saveResults(var_energies)
    DFPlot.saveResults(power)
    time_data = [DFMesh.time_simulation, DFModel.dt, DFModel.n_steps]
    DFPlot.saveResults(time_data)


if __name__ == "__main__":
    start_time = time.time()
    runSimulation(DFMesh.strain_rate)
    total_time = time.time() - start_time
    print("--- %s seconds ---" % (time.time() - start_time))
    DFPlot.saveResults(total_time)
