import numpy as np
import pickle
import progressbar

import DFMesh
import DFPostProcess
import input_files.input_data as inputdata


# Initiation of variables
n_init = 0
n_final = DFMesh.n_steps

u = DFMesh.u0
v = DFMesh.v0
acel = DFMesh.acel0
d = DFMesh.d0
data_bc = DFPostProcess.saveResultsAtBC(u, d)

# energy_potential = np.zeros(DFMesh.n_steps)
# energy_kinetic = np.zeros(DFMesh.n_steps)
# energy_dissipated = np.zeros(DFMesh.n_steps)
# external_work = np.zeros(DFMesh.n_steps)
energy_potential = 0.0
energy_kinetic = 0.0
energy_dissipated = 0.0
external_work = 0.0
work_previous_step = 0.0

# stress_evolution = np.zeros((2 * len(DFMesh.materials), DFMesh.n_steps))
# avg_stress_bar = np.zeros(DFMesh.n_steps)

u_all_steps = [DFMesh.u0]
damage_all_steps = [DFMesh.d0]
# fraglen_all_steps = []
# stress_all_steps = []


# n_fragments = np.zeros(DFMesh.n_steps)
# avg_frag_size = np.zeros(DFMesh.n_steps)
# data_histogram_frag_size = []

if DFMesh.use_cohesive_elements == True:
    energy_reversible = np.zeros(DFMesh.n_steps)
    energy_contact = np.zeros(DFMesh.n_steps)
    energies = [
        ["energy potential", energy_potential],
        ["energy kinetic", energy_kinetic],
        ["energy dissipated", energy_dissipated],
        ["energy reversible", energy_reversible],
        ["energy contact", energy_contact],
        ["external work", external_work],
    ]
if DFMesh.use_lipfield == True:
    energies = [
        ["energy potential", energy_potential],
        ["energy kinetic", energy_kinetic],
        ["energy dissipated", energy_dissipated],
        ["external work", external_work],
    ]


def readPreviousResultsLipfield(variable_name):

    with open("src/input_files/lipfield_" + variable_name + ".pickle", "rb") as handle:
        variable = pickle.load(handle)
    return variable

def readPreviousResultsCZM(variable_name):

    with open("src/input_files/czm_" + variable_name + ".pickle", "rb") as handle:
        variable = pickle.load(handle)
    return variable

def readPreviousResults(variable_name):

    if DFMesh.use_cohesive_elements == True:
        variable = readPreviousResultsCZM(variable_name)
    if DFMesh.use_lipfield == True:
        variable = readPreviousResultsLipfield(variable_name)

    return variable




if inputdata.continue_simulation_from_step == True:
    n_init = inputdata.initial_step
    n_final = DFMesh.n_steps

    # Load previous results
    u = readPreviousResults('u')
    v = readPreviousResults('v')
    acel = readPreviousResults('acel')
    d = readPreviousResults('d')
    n_fragments = readPreviousResults('n_fragments')
    avg_frag_size = readPreviousResults('avg_frag_size')
    avg_stress_bar = readPreviousResults('avg_stress_bar')
    energies = readPreviousResults('energies')
    work_previous_step = DFPostProcess.getEnergy(energies, "external work")[n_init]
    u_all_steps = readPreviousResults('u_all_steps')
    damage_all_steps = readPreviousResults('damage_all_steps')
    stress_all_steps = readPreviousResults('stress_all_steps')
    fraglen_all_steps = readPreviousResults('fraglen_all_steps')
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

