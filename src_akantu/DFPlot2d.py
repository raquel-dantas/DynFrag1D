import numpy as np
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



def Plotlog(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.xscale("log")
    plt.yscale("log")
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
    plt.savefig("LOG/average_stress_bar_dynfrag_akantu.svg")
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
    plt.savefig("LOG/energies_dynfrag_akantu.svg")
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
    plt.savefig("LOG/var_energies_dynfrag_akantu.svg")
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
    plt.savefig("LOG/power_dynfrag_akantu.svg")
    plt.show()



def PlotNumberFragments(nfrag, time_simulation, n_steps):
    """Plot a vector of values that corresponds to the number of fragments at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Number of fragments")
    plt.xlabel("Time (s)")
    plt.ylabel("N")

    x = np.linspace(0, time_simulation, n_steps)
    y = nfrag
    plt.plot(x, y)
    plt.savefig("LOG/number_fragments_dynfrag_akantu.svg")
    plt.show()


def PlotAvgFragmentSize(avg_frag_sizes, time_simulation, n_steps):
    """Plot a vector of values that corresponds to the fragments length at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Fragments sizes")
    plt.xlabel("Time (s)")
    plt.ylabel("m")

    x = np.linspace(0, time_simulation, n_steps)
    y = avg_frag_sizes
    plt.plot(x, y)
    plt.savefig("LOG/size_fragments_dynfrag_akantu.svg")
    plt.show()



def PlotFragmentSizeHistogram(frag_sizes):

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Fragment size distribution")
    plt.xlabel("Fragment size (mm)")
    plt.ylabel("Number of fragments")
    plt.hist(frag_sizes*10**3,5) # Usinf mm
    plt.savefig("LOG/fragment_size_distribution_dynfrag_akantu.svg")
    plt.show()



def PlotConvergenceEnergy(energies_un, energies_nun, meshes):
    """Plot total number of nodes of the mesh x energy dissipated.\n
    Arguments: \n
    energies -- dissipated energies for different mesh sizes;\n
    meshes -- number of linear elements."""

    plt.title(str("Energy convergence"))
    plt.xlabel(str("Number of nodes"))
    plt.ylabel(str("Energy dissipated (N/m)"))
    nnodes = [meshes[i]+1 for i in range(len(meshes))]

    plt.plot(nnodes, energies_un, label='Uniform mesh')
    plt.plot(nnodes, energies_nun, label='Non-uniform mesh')
    plt.legend()
    plt.savefig("LOG/convergence_energy_dynfrag_akantu.svg")
    plt.show()



def PlotLogConvergenceEnergy(energies_un, energies_nun, meshes):
    """Plot total number of nodes of the mesh x energy dissipated.\n
    Arguments: \n
    energies -- dissipated energies for different mesh sizes;\n
    meshes -- number of linear elements."""

    plt.title(str("Energy convergence"))
    plt.xlabel(str("Number of nodes"))
    plt.ylabel(str("Energy dissipated (N/m)"))
    plt.xscale("log")
    plt.yscale("log")

    plt.xlim(10**2,5*10**5)

    nnodes = [meshes[i]+1 for i in range(len(meshes))]

    plt.plot(nnodes, energies_un, label='Uniform mesh')
    plt.plot(nnodes, energies_nun, label='Non-uniform mesh')
    plt.legend()
    plt.savefig("LOG/logconvergence_energy_dynfrag_akantu.svg")
    plt.show()



def PlotConvergenceNumfrag(nfrags_un, nfrags_nun, meshes):
    """Plot total number of nodes of the mesh x final number of fragments.\n
    Arguments: \n
    energies -- dissipated energies for different mesh sizes;\n
    meshes -- number of linear elements."""

    plt.title(str("Nfrag convergence"))
    plt.xlabel(str("Number of nodes"))
    plt.ylabel(str("Number of fragments"))
    nnodes = [meshes[i]+1 for i in range(len(meshes))]

    plt.plot(nnodes, nfrags_un, label='Uniform mesh')
    plt.plot(nnodes, nfrags_nun, label='Non-uniform mesh')
    plt.legend()
    plt.savefig("LOG/convergence_nfrags_dynfrag_akantu.svg")
    plt.show()