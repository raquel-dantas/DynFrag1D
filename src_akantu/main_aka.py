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

    work_previous_step = DFModel.work_previous_step
    energies = DFModel.energies

    # VTK plot
    # DFPlot.addPlotVtk()

    for n in range(DFModel.n_steps):

        DFModel.updateProgressBar(n, bar)

        # DFPlot.addVtkFiles(n)

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
        avg_stress_bar = np.mean(stress_xx)

        # Energy balance
        energies = DFPostProcess.computeEnergies(work_previous_step, fint)
        work_previous_step = DFPostProcess.getEnergy(energies, "external work")

        d = DFPostProcess.getDamageParameter()

        # Fragmentation data
        fragment_data = aka.FragmentManager(DFModel.dynfrag)
        fragment_data.computeAllData()

        # Number of fragments
        n_fragments = fragment_data.getNbFragment()

        # Fragments size (assuming uniform mesh)
        frag_lengths = np.zeros(fragment_data.getNbFragment())
        frag_lengths = fragment_data.getNbElementsPerFragment() 

        # if n%10 == 0:
        results = [
                ["displacement", u],
                ["velocity", v],
                ["acceleration", acel],
                ["damage", d],
                ["stress", stress],
                ["avg_stress_bar", avg_stress_bar],
                ["energies", energies],
                ["n_fragments", n_fragments],
                ["frag_lengths", frag_lengths],
            ]
        DFPlot.saveResultsCurrentStep(results, n)

    DFModel.endProgressBar(bar)

    # Energy balance
    # var_energies = DFPostProcess.computeVariationEnergy(energies)
    # power = DFPostProcess.computePower(energies)

    # Plots
    # DFPlot.plotDamage(d)
    # DFPlot.plotAverageStressBar(avg_stress_bar)
    # DFPlot.plotNumberFragments(n_fragments)
    # DFPlot.plotFragmentSizeHistogram(frag_lengths, 10)
    # DFPlot.plotEnergies(energies)
    # DFPlot.plotVarEnergies(var_energies)
    # DFPlot.plotPower(power)
    # # Save results
    # DFPlot.saveResults(u)
    # DFPlot.saveResults(v)
    # DFPlot.saveResults(acel)
    # DFPlot.saveResults(d)
    # DFPlot.saveResults(damage_all_steps)
    # DFPlot.saveResults(fraglen_all_steps)
    # DFPlot.saveResults(n_fragments)
    # DFPlot.saveResults(avg_frag_size)
    # DFPlot.saveResults(avg_stress_bar)
    # DFPlot.saveResults(energies)
    # DFPlot.saveResults(var_energies)
    # DFPlot.saveResults(power)
    time_data = [DFMesh.time_simulation, DFModel.dt, DFModel.n_steps]
    DFPlot.saveResults(time_data)


if __name__ == "__main__":
    start_time = time.time()
    runSimulation(DFMesh.strain_rate)
    total_time = time.time() - start_time
    print("--- %s seconds ---" % (time.time() - start_time))
    DFPlot.saveResults(total_time)
