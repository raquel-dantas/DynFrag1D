import matplotlib
import DFMesh
import DFPlot
import numpy as np
from matplotlib import pyplot as plt




def GradyAnalyticalFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulaion proposed by Grady (1982)for a given interval of strain rate.    
    """
    s = np.zeros(len(values_strainrate))
    s = np.array([((24.* DFMesh.Gc)/(DFMesh.rho * epsilon**2))**(1/3) for epsilon in values_strainrate])

    return s



points_interval = 10
min_strainrate = 10**1
max_strainrate= 10**10
values_strainrate = np.linspace(min_strainrate, max_strainrate, num=points_interval)
s = GradyAnalyticalFragSize(values_strainrate)

DFPlot.Plotlog(values_strainrate, s, "strain rate (s-1)", "Frag. size (m)", "Grady (1982)")