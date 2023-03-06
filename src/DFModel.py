import numpy as np
import pickle
import progressbar

import DFMesh
import input_files.input_data as inputdata


# Initiation of variables
n_init = 0
n_final = DFMesh.n_steps

u = DFMesh.u0
v = DFMesh.v0
acel = DFMesh.acel0
uprevious_bc_left = np.array([0, 0])
uprevious_bc_right = np.array([0, 0])

energy_potential = np.zeros(DFMesh.n_steps)
energy_kinetic = np.zeros(DFMesh.n_steps)
energy_dissipated = np.zeros(DFMesh.n_steps)
external_work = np.zeros(DFMesh.n_steps)
work_previous_step = 0.0

# stress_evolution = np.zeros((2 * len(DFMesh.materials), DFMesh.n_steps))
avg_stress_bar = np.zeros(DFMesh.n_steps)


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



if inputdata.continue_simulation_from_step == True:
    n_init = inputdata.initial_step
    n_final = DFMesh.n_steps

    # Load previous results
    if DFMesh.use_lipfield == True:

        with open("src/input_files/lipfield_u.pickle", "rb") as handle:
            u = pickle.load(handle)

        with open("src/input_files/lipfield_v.pickle", "rb") as handle:
            v = pickle.load(handle)

        with open("src/input_files/lipfield_acel.pickle", "rb") as handle:
            acel = pickle.load(handle)

        with open("src/input_files/lipfield_d.pickle", "rb") as handle:
            d = pickle.load(handle)

        with open("src/input_files/lipfield_energy_potential.pickle", "rb") as handle:
            energy_potential_previous = pickle.load(handle)

        with open("src/input_files/lipfield_energy_kinetic.pickle", "rb") as handle:
            energy_kinetic_previous = pickle.load(handle)

        with open("src/input_files/lipfield_energy_dissipated.pickle", "rb") as handle:
            energy_dissipated_previous = pickle.load(handle)

        with open("src/input_files/lipfield_external_work.pickle", "rb") as handle:
            external_work_previous = pickle.load(handle)

        with open("src/input_files/lipfield_avg_stress_bar.pickle", "rb") as handle:
            avg_stress_bar_previous = pickle.load(handle)
        
        with open("src/input_files/lipfield_n_fragments.pickle", "rb") as handle:
            n_fragments_previous = pickle.load(handle)

        with open("src/input_files/lipfield_avg_frag_size.pickle", "rb") as handle:
            avg_frag_size_previous = pickle.load(handle)


        for n in range(n_init):
            energy_potential[n] = energy_potential_previous[n]
            energy_kinetic[n] = energy_kinetic_previous[n]
            energy_dissipated[n] = energy_dissipated_previous[n]
            external_work[n] = external_work_previous[n]
            avg_stress_bar[n] = avg_stress_bar_previous[n]
            n_fragments[n] = n_fragments_previous[n]
            avg_frag_size[n] = avg_frag_size_previous[n]


        work_previous_step = external_work[n_init]



        for bc in range(len(DFMesh.materials)):
            if DFMesh.materials[bc] == 4 or DFMesh.materials[bc] == 5:
                if DFMesh.materials[bc] == 4:
                    el_bc = 0
                    uprevious_bc_left = np.array(
                        [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
                    )
                    dprevious_bc_left = d[el_bc]
                else:
                    el_bc = DFMesh.n_elements - 1
                    uprevious_bc_right = np.array(
                        [u[DFMesh.connect[el_bc][0]], u[DFMesh.connect[el_bc][1]]]
                    )
                    dprevious_bc_right = d[el_bc]

    






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