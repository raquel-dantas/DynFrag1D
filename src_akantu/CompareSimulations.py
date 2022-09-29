import numpy as np
import DFPlotCompare
import pickle



# Inputs
# src_czm_interface
# src_akantu




# sim1 -- src_akantu_1000el
# sim2 -- src_akantu_5000el
# sim3 -- src_akantu_7000el
# sim4 -- src_akantu_10000el
# sim5 -- scr_akantu_14000el

time_simulation = 7.0*10**-7
n_steps_sim1 = 1386
n_steps_sim2 = 6932
n_steps_sim3 = 9705
n_steps_sim4 = 13864
n_steps_sim5 = 19410

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

with open('LOG/number_fragments_1000.pickle', 'rb') as handle:
    nfrag_sim1 = pickle.load(handle)
with open('LOG/number_fragments_5000.pickle', 'rb') as handle:
    nfrag_sim2 = pickle.load(handle)
with open('LOG/number_fragments_7000.pickle', 'rb') as handle:
    nfrag_sim3 = pickle.load(handle)
with open('LOG/number_fragments_10000.pickle', 'rb') as handle:
    nfrag_sim4 = pickle.load(handle)
with open('LOG/number_fragments_14000.pickle', 'rb') as handle:
    nfrag_sim5 = pickle.load(handle)



sfrag_sim1 = []
sfrag_sim2 = []

frag_sizes_sim1 = []
frag_sizes_sim2 = []


# Plots

# DFPlotCompare.PlotCompareAverageStressBar(avg_stress_bar_sim1, avg_stress_bar_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareEnergies(Epot_sim1, Ekin_sim1, Edis_sim1, Erev_sim1, Econ_sim1, Wext_sim1,Epot_sim2, Ekin_sim2, Edis_sim2, Erev_sim2, Econ_sim2, Wext_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareVarEnergies(varEpot_sim1, varEkin_sim1, varEdis_sim1, varErev_sim1, varEcon_sim1, varWext_sim1, varEpot_sim2, varEkin_sim2, varEdis_sim2, varErev_sim2, varEcon_sim2, varWext_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

DFPlotCompare.PlotCompareNumberFragments(nfrag_sim1, nfrag_sim2, nfrag_sim3, nfrag_sim4, nfrag_sim5, time_simulation, n_steps_sim1, n_steps_sim2, n_steps_sim3, n_steps_sim4, n_steps_sim5)

# DFPlotCompare.PlotCompareAvgFragmentSize(sfrag_sim1, sfrag_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareFragmentSizeHistogram(frag_sizes_sim1,frag_sizes_sim2)
