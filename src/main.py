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
    data_bc = DFModel.data_bc
    work_previous_step = DFModel.work_previous_step


    for n in range(DFModel.n_init, DFModel.n_final):

        DFModel.updateProgressBar(n, bar)

        stress, average_stress_neighbors = DFPostProcess.computeStress(u, d)
        avg_stress_bar = DFPostProcess.stressBar(stress)
        M, F = DFFem.globalSystem()

        if DFMesh.use_cohesive_elements == True:
            d = DFInterface.getDamageParameter()

        energies = DFPostProcess.updateEnergies(u, v, d, stress, data_bc, work_previous_step )
        work_previous_step = DFPostProcess.getEnergy(energies, "external work")

        # Fragmentation data
        n_fragments = DFFragmentation.getNumberFragments(d)
        frag_lengths, avg_frag_size = DFFragmentation.getFragmentSizes(d)

        # Save results previous step at BC to compute external work
        data_bc = DFPostProcess.saveResultsAtBC(u, d)

        # Time integration
        u, v, acel, d = DFNewmark.explicitScheme(M, u, v, acel, d, F, DFMesh.dt)

        if DFMesh.use_cohesive_elements == True:
            u, v, acel = DFInterface.verifyStress(average_stress_neighbors, u, v, acel)

        if n%10 == 0:
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

    time_data = [DFMesh.time_simulation, DFMesh.dt, DFModel.n_final]
    DFPlot.saveResults(time_data)
    if DFMesh.use_cohesive_elements == True:
        materials = DFMesh.materials
        connect = DFMesh.connect
        DFPlot.saveResults(materials)
        DFPlot.saveResults(connect)


if __name__ == "__main__":
    start_time = time.time()
    runSimulation(DFMesh.strain_rate)
    computational_time = time.time() - start_time
    print("--- %s seconds ---" % (time.time() - start_time))
    DFPlot.saveResults(computational_time)
