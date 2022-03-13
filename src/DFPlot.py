import numpy as np
from matplotlib import pyplot as plt
import inspect
import DFMesh
import DFFem

# Plot functions

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

def PlotStressByTime(n_steps, stress_evl):
    """Plot the stress on the elements at each time step"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')

    x = np.linspace(0, n_steps, n_steps)

    for el in range(len(DFMesh.materials)):
        plt.plot(x, stress_evl[el])

    plt.show()

# Plot energy

def PlotEnergy(n_steps, E_kin, E_pot, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    x = np.linspace(0, n_steps, n_steps)
    plt.plot(x, E_kin, label='E_kin')
    plt.plot(x, E_pot, label='E_pot')
    plt.legend()
    plt.show()


def PlotEnergyTotal(n_steps, E_kin, E_pot, E_tot, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    x = np.linspace(0, n_steps, n_steps)
    plt.plot(x, E_kin, label='E_kin')
    plt.plot(x, E_pot, label='E_pot')
    plt.plot(x, E_tot, label='E_tot')
    plt.legend()
    plt.show()