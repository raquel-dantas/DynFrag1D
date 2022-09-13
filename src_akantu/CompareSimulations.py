import numpy as np
import DFPlotCompare




# Inputs
# sim1 -- src_czm_interface
# sim2 -- src_akantu

time_simulation = 6.0*10**-6
n_steps = 500

avg_stress_bar_sim1 = []
avg_stress_bar_sim2 = []

Epot_sim1 = []
Ekin_sim1 = []
Edis_sim1 = []
Erev_sim1 = []
Econ_sim1 = []
Wext_sim1 = []

Epot_sim2 = []
Ekin_sim2 = []
Edis_sim2 = []
Erev_sim2 = []
Econ_sim2 = []
Wext_sim2 = []

varEpot_sim1 = []
varEkin_sim1 = [] 
varEdis_sim1 = [] 
varErev_sim1 = [] 
varEcon_sim1 = []
varWext_sim1 = []

varEpot_sim2 = []
varEkin_sim2 = [] 
varEdis_sim2 = [] 
varErev_sim2 = [] 
varEcon_sim2 = []
varWext_sim2 = []

nfrag_sim1 = []
nfrag_sim2 = []

sfrag_sim1 = []
sfrag_sim2 = []

frag_sizes_sim1 = []
frag_sizes_sim2 = []


# Plots

DFPlotCompare.PlotCompareAverageStressBar(avg_stress_bar_sim1, avg_stress_bar_sim2, time_simulation, n_steps)

DFPlotCompare.PlotCompareEnergies(Epot_sim1, Ekin_sim1, Edis_sim1, Erev_sim1, Econ_sim1, Wext_sim1,Epot_sim2, Ekin_sim2, Edis_sim2, Erev_sim2, Econ_sim2, Wext_sim2, time_simulation, n_steps)

DFPlotCompare.PlotCompareVarEnergies(varEpot_sim1, varEkin_sim1, varEdis_sim1, varErev_sim1, varEcon_sim1, varWext_sim1, varEpot_sim2, varEkin_sim2, varEdis_sim2, varErev_sim2, varEcon_sim2, varWext_sim2, time_simulation, n_steps)

DFPlotCompare.PlotCompareNumberFragments(nfrag_sim1, nfrag_sim2, time_simulation, n_steps)

DFPlotCompare.PlotCompareAvgFragmentSize(sfrag_sim1, sfrag_sim2, time_simulation, n_steps)

DFPlotCompare.PlotCompareFragmentSizeHistogram(frag_sizes_sim1,frag_sizes_sim2)
