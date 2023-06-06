import numpy as np

import DFProcessResults
import DFPlotResults



# This script is to quick check results
stress = []
n_frag = []

# src_akantu
# Path to the output files
simulation_name = "CZM 1250 elements"
filepath = "output/" 
file_address = filepath + "akantu_"

time_data = DFProcessResults.getTimeData(file_address)
time_simulation = time_data[0]
dt = time_data[1]
n_steps = time_data[2]

n_files = int(n_steps / 10 + 1)

avg_stress_bar, energies, n_fragments = DFProcessResults.getResultsAllSteps(file_address, n_files, n_steps)
time = np.linspace(0, time_simulation, n_files)
stress.append([simulation_name, time, avg_stress_bar])
n_frag.append([simulation_name, time, n_fragments])

DFPlotResults.plotResults(
    stress,
    label_x="time (s)",
    label_y="Average stress at the bar (MPa)",
)
DFPlotResults.plotResults(
    n_frag,
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