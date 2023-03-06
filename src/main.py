import numpy as np
import time

import DFMesh
import DFFem
import DFPostProcess
import DFNewmark
import DFInterface
import DFPlot
import DFFragmentation
import DFModel


def runSimulation(strain_rate):

    bar = DFModel.initProgressBar()

    # Initiation of variables

    u = DFModel.u
    v = DFModel.v
    acel = DFModel.acel
    d = DFModel.d
    avg_stress_bar = DFModel.avg_stress_bar

    data_bc = DFModel.data_bc
    work_previous_step = DFModel.work_previous_step
    energies = DFModel.energies

    n_fragments = DFModel.n_fragments
    avg_frag_size = DFModel.avg_frag_size
    data_histogram_frag_size = []

    for n in range(DFModel.n_init, DFModel.n_final):

        DFModel.updateProgressBar(n, bar)

        stress, average_stress_neighbors = DFPostProcess.computeStress(u, d)
        avg_stress_bar[n] = DFPostProcess.stressBar(stress)
        M, F = DFFem.globalSystem()

        # DFPlot.PlotVTK('animation/test',n,u,stress)

        if DFMesh.use_cohesive_elements == True:
            d = DFInterface.getDamageParameter()

        energies = DFPostProcess.updateEnergies(
            energies, n, u, v, d, stress, data_bc, work_previous_step
        )
        work_previous_step = DFPostProcess.getEnergy(energies, "external work")[n]

        # Fragmentation data
        n_fragments[n] = DFFragmentation.getNumberFragments(d)
        frag_lengths, avg_frag_size[n] = DFFragmentation.getFragmentSizes(d)
        data_histogram_frag_size = DFFragmentation.getFragSizeHistogramData(
            frag_lengths
        )

        # Save results previous step at BC to compute external work
        data_bc = DFPostProcess.saveResultsAtBC(u, d)

        # Time integration
        u, v, acel, d = DFNewmark.explicitScheme(M, u, v, acel, d, F, DFMesh.dt)

        if DFMesh.use_cohesive_elements == True:
            u, v, acel = DFInterface.verifyStress(average_stress_neighbors, u, v, acel)

    DFModel.endProgressBar(bar)

    DFPlot.plotAverageStressBar(avg_stress_bar)
    DFPlot.plotNumberFragments(n_fragments)
    DFPlot.plotFragmentSizeHistogram(frag_lengths, 10)

    # Energy balance
    var_energies = DFPostProcess.computeVariationEnergy(energies)
    power = DFPostProcess.computePower(energies)

    DFPlot.plotEnergies(energies)
    DFPlot.plotVarEnergies(var_energies)
    DFPlot.plotPower(power)

    if DFMesh.use_lipfield == True:
        DFPlot.plotByIntPoint(d)

    # Save results
    DFPlot.saveResults(u)
    DFPlot.saveResults(v)
    DFPlot.saveResults(acel)
    DFPlot.saveResults(d)

    DFPlot.saveResults(n_fragments)
    DFPlot.saveResults(avg_frag_size)
    DFPlot.saveResults(data_histogram_frag_size)

    DFPlot.saveResults(avg_stress_bar)
    DFPlot.saveResults(energies)
    DFPlot.saveResults(var_energies)
    DFPlot.saveResults(power)


if __name__ == "__main__":
    start_time = time.time()
    runSimulation(DFMesh.strain_rate)
    total_time = time.time() - start_time
    print("--- %s seconds ---" % (time.time() - start_time))
    if DFMesh.use_cohesive_elements == True:
        DFPlot.saveResultsCZM(total_time)
    if DFMesh.use_lipfield == True:
        DFPlot.saveResultsLipfield(total_time)
