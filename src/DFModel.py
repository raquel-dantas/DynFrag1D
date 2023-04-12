import numpy as np
import pickle
import progressbar

import DFMesh
import DFPostProcess


# Initiation of variables
n_init = 0
n_final = DFMesh.n_steps

u = DFMesh.u0
v = DFMesh.v0
acel = DFMesh.acel0
d = DFMesh.d0
data_bc = DFPostProcess.saveResultsAtBC(u, d)

work_previous_step = 0.0





def readPreviousResults(previous_simulation_file):
    
    with open(previous_simulation_file, "rb") as handle:
        previous_results = pickle.load(handle)

    return previous_results



def getResults(results, variable_name):

    for i in range(len(results)):
        if results[i][0] == variable_name:
            variable = results[i][1]
            return variable







if DFMesh.continue_simulation_from_step == True:
    n_init = DFMesh.n_init
    n_final = DFMesh.n_steps

    previous_results = readPreviousResults(DFMesh.previous_simulation)

    # Load previous results
    u = getResults(previous_results, "displacement")
    v = getResults(previous_results, "velocity")
    acel = getResults(previous_results, "acceleration")
    d = getResults(previous_results, "damage")
    energies = getResults(previous_results, "energies")
    work_previous_step = DFPostProcess.getEnergy(energies, "external work")
    data_bc = DFPostProcess.saveResultsAtBC(u, d)


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
    progress = int(bar.maxval * float(n / DFMesh.n_steps))
    bar.update(progress)
