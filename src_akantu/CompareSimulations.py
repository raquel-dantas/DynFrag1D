import numpy as np
import DFPlotCompare
import DFPlot2d
import pickle
from matplotlib import pyplot as plt


# n_steps_1000el = 3961
# n_steps_2000el = 7922
# n_steps_3000el = 11883
# n_steps_4000el = 15845
# n_steps_5000el = 19806
# n_steps_6000el = 23767
# n_steps_7000el = 27729
# n_steps_8000el = 31690
# n_steps_9000el = 35651
# n_steps_10000el = 39613

# 10to3
with open('LOG/10to3_1000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_1000el = pickle.load(handle)
    nfrag_10to3_1000el = max(nfrag_10to3_1000el)

with open('LOG/10to3_2000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_2000el = pickle.load(handle)
    nfrag_10to3_2000el = max(nfrag_10to3_2000el)
 
with open('LOG/10to3_3000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_3000el = pickle.load(handle)
    nfrag_10to3_3000el = max(nfrag_10to3_3000el)

with open('LOG/10to3_4000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_4000el = pickle.load(handle)
    nfrag_10to3_4000el = max(nfrag_10to3_4000el)

with open('LOG/10to3_5000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_5000el = pickle.load(handle)
    nfrag_10to3_5000el = max(nfrag_10to3_5000el)

with open('LOG/10to3_6000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_6000el = pickle.load(handle)
    nfrag_10to3_6000el = max(nfrag_10to3_6000el)

with open('LOG/10to3_7000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_7000el = pickle.load(handle)
    nfrag_10to3_7000el = max(nfrag_10to3_7000el)

with open('LOG/10to3_8000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_8000el = pickle.load(handle)
    nfrag_10to3_8000el = max(nfrag_10to3_8000el)

with open('LOG/10to3_9000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_9000el = pickle.load(handle)
    nfrag_10to3_9000el = max(nfrag_10to3_9000el)

with open('LOG/10to3_10000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to3_10000el = pickle.load(handle)
    nfrag_10to3_10000el = max(nfrag_10to3_10000el)

meshes_10to3 = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
nfrags_10to3 = [nfrag_10to3_1000el, nfrag_10to3_2000el, nfrag_10to3_3000el, nfrag_10to3_4000el, nfrag_10to3_5000el, nfrag_10to3_6000el, nfrag_10to3_7000el, nfrag_10to3_8000el, nfrag_10to3_9000el, nfrag_10to3_10000el]



# 10to4
with open('LOG/10to4_1000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_1000el = pickle.load(handle)
    nfrag_10to4_1000el = max(nfrag_10to4_1000el)

with open('LOG/10to4_2000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_2000el = pickle.load(handle)
    nfrag_10to4_2000el = max(nfrag_10to4_2000el)

with open('LOG/10to4_3000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_3000el = pickle.load(handle)
    nfrag_10to4_3000el = max(nfrag_10to4_3000el)

with open('LOG/10to4_4000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_4000el = pickle.load(handle)
    nfrag_10to4_4000el = max(nfrag_10to4_4000el)

with open('LOG/10to4_5000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_5000el = pickle.load(handle)
    nfrag_10to4_5000el = max(nfrag_10to4_5000el)

with open('LOG/10to4_6000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_6000el = pickle.load(handle)
    nfrag_10to4_6000el = max(nfrag_10to4_6000el)

with open('LOG/10to4_7000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_7000el = pickle.load(handle)
    nfrag_10to4_7000el = max(nfrag_10to4_7000el)

with open('LOG/10to4_8000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_8000el = pickle.load(handle)
    nfrag_10to4_8000el = max(nfrag_10to4_8000el)

with open('LOG/10to4_9000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_9000el = pickle.load(handle)
    nfrag_10to4_9000el = max(nfrag_10to4_9000el)

with open('LOG/10to4_10000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to4_10000el = pickle.load(handle)
    nfrag_10to4_10000el = max(nfrag_10to4_10000el)


meshes_10to4 = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
nfrags_10to4 = [nfrag_10to4_1000el, nfrag_10to4_2000el, nfrag_10to4_3000el, nfrag_10to4_4000el, nfrag_10to4_5000el, nfrag_10to4_6000el, nfrag_10to4_7000el, nfrag_10to4_8000el, nfrag_10to4_9000el, nfrag_10to4_10000el]



# 10to5
with open('LOG/10to5_1000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_1000el = pickle.load(handle)
    nfrag_10to5_1000el = max(nfrag_10to5_1000el)

with open('LOG/10to5_2000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_2000el = pickle.load(handle)
    nfrag_10to5_2000el = max(nfrag_10to5_2000el)

with open('LOG/10to5_3000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_3000el = pickle.load(handle)
    nfrag_10to5_3000el = max(nfrag_10to5_3000el)

with open('LOG/10to5_4000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_4000el = pickle.load(handle)
    nfrag_10to5_4000el = max(nfrag_10to5_4000el)

with open('LOG/10to5_5000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_5000el = pickle.load(handle)
    nfrag_10to5_5000el = max(nfrag_10to5_5000el)

with open('LOG/10to5_6000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_6000el = pickle.load(handle)
    nfrag_10to5_6000el = max(nfrag_10to5_6000el)

with open('LOG/10to5_7000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_7000el = pickle.load(handle)
    nfrag_10to5_7000el = max(nfrag_10to5_7000el)

with open('LOG/10to5_8000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_8000el = pickle.load(handle)
    nfrag_10to5_8000el = max(nfrag_10to5_8000el)

with open('LOG/10to5_9000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_9000el = pickle.load(handle)
    nfrag_10to5_9000el = max(nfrag_10to5_9000el)

with open('LOG/10to5_10000el_src_akantu_number_fragments.pickle', 'rb') as handle:
    nfrag_10to5_10000el = pickle.load(handle)
    nfrag_10to5_10000el = max(nfrag_10to5_10000el)


meshes_10to5 = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
nfrags_10to5 = [nfrag_10to5_1000el, nfrag_10to5_2000el, nfrag_10to5_3000el, nfrag_10to5_4000el, nfrag_10to5_5000el, nfrag_10to5_6000el, nfrag_10to5_7000el, nfrag_10to5_8000el, nfrag_10to5_9000el, nfrag_10to5_10000el]

meshes = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]


def PlotConvergenceNumfragMeshesAkantu(meshes, nfrags_10to3, nfrags_10to4,  nfrags_10to5 ):
   
    plt.title(str("Number of fragments convergence"))
    plt.xlabel(str("Number of elements"))
    plt.ylabel(str("Number of fragments"))
    nnodes = [meshes[i]+1 for i in range(len(meshes))]
    
    plt.plot(meshes, nfrags_10to3, label=10**3)
    plt.plot(meshes, nfrags_10to4, label=10**4)
    plt.plot(meshes, nfrags_10to5, label=10**5)
    plt.legend()
    # plt.savefig("LOG/convergence_nfrags_akantu.svg")
    plt.show()


PlotConvergenceNumfragMeshesAkantu(meshes, nfrags_10to3, nfrags_10to4, nfrags_10to5)




































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




