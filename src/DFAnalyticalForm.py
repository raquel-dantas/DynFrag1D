import matplotlib
import DFMesh
import DFPlot
import numpy as np
from matplotlib import pyplot as plt




def GradyFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Grady (1982) for a given interval of strain rate."""

    s = np.zeros(len(values_strainrate))
    s = np.array([((24.* DFMesh.Gc)/(DFMesh.rho * epsilon**2))**(1/3) for epsilon in values_strainrate])

    return s


def GlenChudnoviskFragSize(values_strainrate):

    for epsilon in range(len(values_strainrate)):
        alpha = 3.*DFMesh.stress_c / (DFMesh.rho * DFMesh.E * epsilon**2)
        beta = 3.*DFMesh.G_c / (2.* DFMesh.rho * DFMesh.E * epsilon**2)
        phi = np.arcsinh(beta*(3./alpha)**(3/2))
        s = 4.*(alpha/3.)*np.sinh(phi/3.)

    return s



points_interval = 10
min_strainrate = 10**1
max_strainrate= 10**10
values_strainrate = np.linspace(min_strainrate, max_strainrate, num=points_interval)
s_grady = GradyFragSize(values_strainrate)
s_gc = GlenChudnoviskFragSize(values_strainrate)


def PlotMultLog():

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Stress evolution")
    plt.xlabel("Time (s)")
    plt.ylabel("Stress (Pa)")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    for el in range(len(DFMesh.materials)):
        plt.plot(x, stress_evl[el], label=el)
    plt.legend()
    plt.show()




DFPlot.Plotlog(values_strainrate, s, "strain rate (s-1)", "Frag. size (m)", "Grady (1982)")