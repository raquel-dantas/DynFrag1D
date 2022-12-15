import numpy as np
import DFPlotCompare
import DFPlot
import pickle
from matplotlib import pyplot as plt


# Compare number of fragments for 10to5 using czm_interface for uniform and non-uniform mesh

# Uniform mesh
with open('LOG/10to5_src_czm_250el_number_fragments.pickle', 'rb') as handle:
    nfrag_250el_un = pickle.load(handle)
    nfrag_250el_un = max(nfrag_250el_un)

with open('LOG/10to5_src_czm_500el_number_fragments.pickle', 'rb') as handle:
    nfrag_500el_un = pickle.load(handle)
    nfrag_500el_un = max(nfrag_500el_un)
 
with open('LOG/10to5_src_czm_750el_number_fragments.pickle', 'rb') as handle:
    nfrag_750el_un = pickle.load(handle)
    nfrag_750el_un = max(nfrag_750el_un)

with open('LOG/10to5_src_czm_1000el_number_fragments.pickle', 'rb') as handle:
    nfrag_1000el_un = pickle.load(handle)
    nfrag_1000el_un = max(nfrag_1000el_un)

with open('LOG/10to5_src_czm_1250el_number_fragments.pickle', 'rb') as handle:
    nfrag_1250el_un = pickle.load(handle)
    nfrag_1250el_un = max(nfrag_1250el_un)

with open('LOG/10to5_src_czm_1500el_number_fragments.pickle', 'rb') as handle:
    nfrag_1500el_un = pickle.load(handle)
    nfrag_1500el_un = max(nfrag_1500el_un)

with open('LOG/10to5_src_czm_2000el_number_fragments.pickle', 'rb') as handle:
    nfrag_2000el_un = pickle.load(handle)
    nfrag_2000el_un = max(nfrag_2000el_un)

with open('LOG/10to5_src_czm_3000el_number_fragments.pickle', 'rb') as handle:
    nfrag_3000el_un = pickle.load(handle)
    nfrag_3000el_un = max(nfrag_3000el_un)

with open('LOG/10to5_src_czm_4000el_number_fragments.pickle', 'rb') as handle:
    nfrag_4000el_un = pickle.load(handle)
    nfrag_4000el_un = max(nfrag_4000el_un)

# Non-uniform mesh
with open('LOG/10to5_src_czm_250elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_250el_nun = pickle.load(handle)
    nfrag_250el_nun = max(nfrag_250el_nun)

with open('LOG/10to5_src_czm_500elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_500el_nun = pickle.load(handle)
    nfrag_500el_nun = max(nfrag_500el_nun)
 
with open('LOG/10to5_src_czm_750elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_750el_nun = pickle.load(handle)
    nfrag_750el_nun = max(nfrag_750el_nun)

with open('LOG/10to5_src_czm_1000elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_1000el_nun = pickle.load(handle)
    nfrag_1000el_nun = max(nfrag_1000el_nun)

with open('LOG/10to5_src_czm_1250elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_1250el_nun = pickle.load(handle)
    nfrag_1250el_nun = max(nfrag_1250el_nun)

with open('LOG/10to5_src_czm_1500elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_1500el_nun = pickle.load(handle)
    nfrag_1500el_nun = max(nfrag_1500el_nun)

with open('LOG/10to5_src_czm_2000elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_2000el_nun = pickle.load(handle)
    nfrag_2000el_nun = max(nfrag_2000el_nun)

with open('LOG/10to5_src_czm_3000elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_3000el_nun = pickle.load(handle)
    nfrag_3000el_nun = max(nfrag_3000el_nun)

with open('LOG/10to5_src_czm_4000elnun_number_fragments.pickle', 'rb') as handle:
    nfrag_4000el_nun = pickle.load(handle)
    nfrag_4000el_nun = max(nfrag_4000el_nun)


meshes = [250, 500, 750, 1000, 1250,1500, 2000, 3000, 4000]
nfrags_un = [nfrag_250el_un, nfrag_500el_un, nfrag_750el_un, nfrag_1000el_un, nfrag_1250el_un, nfrag_1500el_un, nfrag_2000el_un, nfrag_3000el_un, nfrag_4000el_un]
nfrags_nun = [nfrag_250el_nun, nfrag_500el_nun, nfrag_750el_nun, nfrag_1000el_nun, nfrag_1250el_nun, nfrag_1500el_nun, nfrag_2000el_nun, nfrag_3000el_nun, nfrag_4000el_nun]




def PlotConvergenceNumfragMeshes(meshes, nfrags_un, nfrags_nun):
   
    # plt.title(str("Number of fragments convergence"))
    hfont = {'fontname':'sans-serif'}
    # fig, axes = plt.subplots()
    # axes.grid(True, which='both')
    # axes.axhline(y=0, color='k')
    plt.rcParams.update({'font.size': 15})
    plt.xlabel('Number of elements', **hfont)
    plt.ylabel('Number of fragments', **hfont)
    nnodes = [meshes[i]+1 for i in range(len(meshes))]
    
    plt.plot(meshes, nfrags_un, label='Uniform mesh', marker='o')
    plt.plot(meshes, nfrags_nun, label='Non-uniform mesh',marker='s')
    plt.legend()
    plt.savefig("LOG/convergence_nfrags_akantu.pdf")
    plt.show()


PlotConvergenceNumfragMeshes(meshes, nfrags_un, nfrags_nun)




































# print(nfrag_1000el[n_steps_1000el-1])
# print(nfrag_sim2[n_steps_sim2-1])
# print(nfrag_sim3[n_steps_sim3-1])
# print(nfrag_sim4[n_steps_sim4-1])
# print(nfrag_sim5[n_steps_sim5-1])

# avg_stress_bar_sim1 = []
# avg_stress_bar_sim2 = []

# Epot_sim1 = []
# Ekin_sim1 = []
# Edis_sim1 = []
# Erev_sim1 = []
# Econ_sim1 = []
# Wext_sim1 = []

# Epot_sim2 = []
# Ekin_sim2 = []
# Edis_sim2 = []
# Erev_sim2 = []
# Econ_sim2 = []
# Wext_sim2 = []

# varEpot_sim1 = []
# varEkin_sim1 = [] 
# varEdis_sim1 = [] 
# varErev_sim1 = [] 
# varEcon_sim1 = []
# varWext_sim1 = []



# with open('LOG/number_fragments_500_czmint.pickle', 'rb') as handle:
#     nfrag_sim6 = pickle.load(handle)
# with open('LOG/number_fragments_1000_czmint.pickle', 'rb') as handle:
#     nfrag_sim7 = pickle.load(handle)
# with open('LOG/number_fragments_2000_czmint.pickle', 'rb') as handle:
#     nfrag_sim8 = pickle.load(handle)
# with open('LOG/number_fragments_2500_czmint.pickle', 'rb') as handle:
#     nfrag_sim9 = pickle.load(handle)
# with open('LOG/number_fragments_3000_czmint.pickle', 'rb') as handle:
#     nfrag_sim10 = pickle.load(handle)





# print(nfrag_sim6[n_steps_sim6-1])
# print(nfrag_sim7[n_steps_sim7-1])
# print(nfrag_sim8[n_steps_sim8-1])
# print(nfrag_sim9[n_steps_sim9-1])
# print(nfrag_sim10[n_steps_sim10-1])


# sfrag_sim1 = []
# sfrag_sim2 = []

# frag_sizes_sim1 = []
# frag_sizes_sim2 = []


# Plots

# DFPlotCompare.PlotCompareAverageStressBar(avg_stress_bar_sim1, avg_stress_bar_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareEnergies(Epot_sim1, Ekin_sim1, Edis_sim1, Erev_sim1, Econ_sim1, Wext_sim1,Epot_sim2, Ekin_sim2, Edis_sim2, Erev_sim2, Econ_sim2, Wext_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareVarEnergies(varEpot_sim1, varEkin_sim1, varEdis_sim1, varErev_sim1, varEcon_sim1, varWext_sim1, varEpot_sim2, varEkin_sim2, varEdis_sim2, varErev_sim2, varEcon_sim2, varWext_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareNumberFragments(nfrag_sim1, nfrag_sim2, nfrag_sim3, nfrag_sim4, nfrag_sim5, time_simulation, n_steps_sim1, n_steps_sim2, n_steps_sim3, n_steps_sim4, n_steps_sim5)

# DFPlotCompare.PlotCompareAvgFragmentSize(sfrag_sim1, sfrag_sim2, time_simulation, n_steps_sim1, n_steps_sim2)

# DFPlotCompare.PlotCompareFragmentSizeHistogram(frag_sizes_sim1,frag_sizes_sim2)




