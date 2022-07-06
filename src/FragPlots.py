import numpy as np
import DFPlot
import DFFragmentation


# Analytical estimations of fragment size

points_interval = 50
min_strainrate = 10**1
max_strainrate= 10**10
values_strainrate = np.linspace(min_strainrate, max_strainrate, num=points_interval)
s_grady = DFFragmentation.GradyFragSize(values_strainrate)
s_gc = DFFragmentation.GlenChudnoviskFragSize(values_strainrate)
s_zmr = DFFragmentation.ZhouMolinariRameshFragSize(values_strainrate)

DFPlot.PlotAnalyticals(s_grady, s_gc, s_zmr, values_strainrate)


# Convergence study

meshes = np.array(500,1000,1500,2000,2500)
energies = np.array()
num_frags = np.array()

DFPlot.PlotConvergenceEnergy(energies,meshes)
DFPlot.PlotConvergenceNumfrag(num_frags,meshes)
