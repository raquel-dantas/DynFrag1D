import numpy as np

import DFProcessResults
import DFPlotResults



# This script is for a mesh study: compare several meshes
# Initialize what you want to compare:
avg_stress_10to5_uniform_mesh = []
n_fragments_10to5_uniform_mesh  = []
avg_stress_10to5_non_uniform_mesh  = []
n_fragments_10to5_non_uniform_mesh = []

avg_stress_10to4_uniform_mesh = []
n_fragments_10to4_uniform_mesh  = []
avg_stress_10to4_non_uniform_mesh  = []
n_fragments_10to4_non_uniform_mesh = []

avg_stress_10to3_uniform_mesh = []
n_fragments_10to3_uniform_mesh  = []
avg_stress_10to3_non_uniform_mesh  = []
n_fragments_10to3_non_uniform_mesh = []


# Path to the output files: strain-rate of 10to5
filepath = "output_mesh_study/10to5/"

# 10to5: Lip-field 625 elements Uniform mesh
simulation_name = "LIP 625 elements"
file_address = filepath + "10to5_lipfield_uniform_625/lipfield_"
time_data = DFProcessResults.getTimeData(file_address)
time_simulation = time_data[0]
dt = time_data[1]
n_steps = time_data[2]
n_files = int(n_steps / 10 + 1)
avg_stress_bar, energies, n_fragments = DFProcessResults.getResultsAllSteps(file_address, n_files, n_steps)
time = np.linspace(0, time_simulation, n_files)
avg_stress_10to5_uniform_mesh.append([simulation_name, time, avg_stress_bar])
n_fragments_10to5_uniform_mesh.append([simulation_name, time, n_fragments])

# 10to5: Lip-field 1250 elements Uniform mesh
simulation_name = "LIP 1250 elements"
file_address = filepath + "10to5_lipfield_uniform_1250/lipfield_"
time_data = DFProcessResults.getTimeData(file_address)
time_simulation = time_data[0]
dt = time_data[1]
n_steps = time_data[2]
n_files = int(n_steps / 10 + 1)
avg_stress_bar, energies, n_fragments = DFProcessResults.getResultsAllSteps(file_address, n_files, n_steps)
time = np.linspace(0, time_simulation, n_files)
avg_stress_10to5_uniform_mesh.append([simulation_name, time, avg_stress_bar])
n_fragments_10to5_uniform_mesh.append([simulation_name, time, n_fragments])

# 10to5: Lip-field 2500 elements Uniform mesh
simulation_name = "LIP 2500 elements"
file_address = filepath + "10to5_lipfield_uniform_2500/lipfield_"
time_data = DFProcessResults.getTimeData(file_address)
time_simulation = time_data[0]
dt = time_data[1]
n_steps = time_data[2]
n_files = int(n_steps / 10 + 1)
avg_stress_bar, energies, n_fragments = DFProcessResults.getResultsAllSteps(file_address, n_files, n_steps)
time = np.linspace(0, time_simulation, n_files)
avg_stress_10to5_uniform_mesh.append([simulation_name, time, avg_stress_bar])
n_fragments_10to5_uniform_mesh.append([simulation_name, time, n_fragments])

# 10to5: Lip-field 7500 elements Uniform mesh
simulation_name = "LIP 7500 elements"
file_address = filepath + "10to5_lipfield_uniform_7500/lipfield_"
time_data = DFProcessResults.getTimeData(file_address)
time_simulation = time_data[0]
dt = time_data[1]
n_steps = time_data[2]
n_files = int(n_steps / 10 + 1)
avg_stress_bar, energies, n_fragments = DFProcessResults.getResultsAllSteps(file_address, n_files, n_steps)
time = np.linspace(0, time_simulation, n_files)
avg_stress_10to5_uniform_mesh.append([simulation_name, time, avg_stress_bar])
n_fragments_10to5_uniform_mesh.append([simulation_name, time, n_fragments])













DFPlotResults.plotResults(
    avg_stress_10to5_uniform_mesh,
    label_x="time (s)",
    label_y="Average stress at the bar (MPa)",
)
DFPlotResults.plotResults(
    n_fragments_10to5_uniform_mesh,
    label_x="time (s)",
    label_y="N",
)

var_energies = DFProcessResults.computeVarEnergiesCZM(
    energies, n_files, 1250)
DFPlotResults.plotVarEnergiesCZM(
    var_energies,
    time_simulation,
    title="CZM: Variation of energy"
)

# DFPlotResults.plotFragmentSizeHistogram(sfrag[n_files -1])