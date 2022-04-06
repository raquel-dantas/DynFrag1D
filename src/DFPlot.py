import numpy as np
from matplotlib import pyplot as plt
import inspect
import DFMesh
import DFFem

# Plot function
def Plot(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.plot(x, y)
    plt.show()


def retrieve_name(var):
    """Gets the name of the argument passed to it, as you coded it in your python script"""

    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]


def PlotByDOF(func):
    """Plot a vector of values that corresponds to each DOF of the mesh"""

    labely = retrieve_name(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    h = DFMesh.L/n_oneD_elements
    x0 = 0
    x = np.array([[x0+h*el, x0+h*(el+1)] for el in range(n_oneD_elements)])
    x = x.flatten()
    y = np.array([[func[DFFem.Gl_index(el,0)],func[DFFem.Gl_index(el,1)]] for el in range(n_oneD_elements)])
    y = y.flatten()
    Plot(x,y,"x",labely,title)


def PlotByElement(func):
    """Plot a vector of values that corresponds to each element of the mesh """

    labely = retrieve_name(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    h = DFMesh.L/n_oneD_elements
    x0 = 0
    x = np.array([[x0+h*el, x0+h*(el+1)] for el in range(n_oneD_elements)])
    x = x.flatten()
    y = np.array([[func[el], func[el]] for el in range(n_oneD_elements)])
    y = y.flatten()
    Plot(x,y,"x",labely,title)


# Plot points
def PlotScatter(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.scatter(x, y)
    plt.show()


def PlotByInterface(func):
    """Plot a vector of values that corresponds to each interface element of the mesh """

    labely = retrieve_name(func)[0]
    title = labely
    n_int_elements = sum(1 for el in DFMesh.materials if el == 1)
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    h = DFMesh.L/n_oneD_elements
    x = np.array([i*h for i in range(0, n_oneD_elements + 1)])
    y = np.zeros(x.shape[0])
    for el in range(n_oneD_elements, len(DFMesh.materials)):
        if DFMesh.materials[el] == 1:
            j = DFMesh.connect[el][0]
            y[j] = func[el]
    PlotScatter(x,y,"x",labely,title)

def PlotAverageStressBar(average_stress_bar):
    """Plot a vector of values that corresponds to the average stress between all elements at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Average stress bar")
    plt.xlabel("Time (s)")
    plt.ylabel("Average Stress (Pa)")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    y = average_stress_bar
    plt.plot(x, y)
    plt.show()


# Plot more than one function

def PlotStressByTime(stress_evl):
    """Plot the stress on the elements at each time step"""

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


def PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext):
    """Plot energies values per time"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Energies")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy  (N/m)")
        
    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, Epot, label='Epot')
    plt.plot(x, Ekin, label='Ekin')
    plt.plot(x, Edis, label='Edis')
    # plt.plot(x, Erev, label='Erev')
    # plt.plot(x, Econ, label='Econ')
    plt.plot(x, Wext, label='Wext')
    plt.legend()
    plt.show()


def PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot):
# def PlotVarEnergy( varEtot):
    """Plot variation of energy per time"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str("Variation of energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Variation of energy"))

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, varEpot, label='varEpot')
    plt.plot(x, varEkin, label='varEkin')
    plt.plot(x, varEdis, label='varEdis')
    # plt.plot(x, varErev, label='varErev')
    # plt.plot(x, varEcon, label='varEcon')
    plt.plot(x, -varWext, label='-varWext')
    plt.plot(x, varEtot, label='varEtot')
    plt.legend()
    plt.show()