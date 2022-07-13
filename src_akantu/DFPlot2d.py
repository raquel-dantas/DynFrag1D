import numpy as np
import inspect
from matplotlib import pyplot as plt

def Plot(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.plot(x, y)
    plt.show()


def PlotByCoord(func, mesh, labely):
    """Plot a vector of values that corresponds to each DOF of the mesh"""

    title = labely

    x = np.array([x for x,y in mesh.getNodes()])
    y = func
    x = x.flatten()
    y = y.flatten()
    Plot(x,y,"x",labely,title)



def PlotAverageStressBar(average_stress_bar, time_simulation, n_steps):
    """Plot a vector of values that corresponds to the average stress between all elements at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Average stress bar")
    plt.xlabel("Time (s)")
    plt.ylabel("Average Stress (Pa)")

    x = np.linspace(0, time_simulation, n_steps)
    y = average_stress_bar
    plt.plot(x, y)
    plt.show()



def PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, time_simulation, n_steps):
    """Plot energies values per time"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Energies")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy  (N/m)")
        
    x = np.linspace(0, time_simulation, n_steps)
    plt.plot(x, Epot, label='Epot')
    plt.plot(x, Ekin, label='Ekin')
    plt.plot(x, Edis, label='Edis')
    plt.plot(x, Erev, label='Erev')
    plt.plot(x, Econ, label='Econ')
    plt.plot(x, Wext, label='Wext')
    plt.legend()
    plt.show()

def PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot, time_simulation, n_steps):
    """Plot variation of energy from time t to t0."""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str("Variation of energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Variation of energy"))

    x = np.linspace(0, time_simulation, n_steps)
    plt.plot(x, varEpot, label='varEpot')
    plt.plot(x, varEkin, label='varEkin')
    plt.plot(x, varEdis, label='varEdis')
    plt.plot(x, varErev, label='varErev')
    plt.plot(x, varEcon, label='varEcon')
    plt.plot(x, -varWext, label='-varWext')
    plt.plot(x, varEtot, label='varEtot')
    plt.legend()
    plt.show()


def PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot, time_simulation, n_steps):
    """Plot variation of energy between time steps"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str("Power"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Power"))

    x = np.linspace(0, time_simulation, n_steps)
    plt.plot(x, PEpot, label='varEpot')
    plt.plot(x, PEkin, label='varEkin')
    plt.plot(x, PEdis, label='varEdis')
    plt.plot(x, PErev, label='varErev')
    plt.plot(x, PEcon, label='varEcon')
    plt.plot(x, -PWext, label='-varWext')
    plt.plot(x, PEtot, label='varEtot')
    plt.legend()
    plt.show()