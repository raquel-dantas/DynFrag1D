import numpy as np
from matplotlib import pyplot as plt
import inspect
import DFMesh
import DFFem



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

def PlotScatter(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.scatter(x, y)
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
    x = np.array([[DFMesh.node_coord[el] , DFMesh.node_coord[el+1]] for el in range(n_oneD_elements)])
    x = x.flatten()
    y = np.array([[func[DFFem.Gl_index(el,0)],func[DFFem.Gl_index(el,1)]] for el in range(n_oneD_elements)])
    y = y.flatten()
    Plot(x,y,"x",labely,title)



def PlotByElement(func):
    """Plot a vector of values that corresponds to each element of the mesh """

    labely = retrieve_name(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    x = np.array([[DFMesh.node_coord[el] , DFMesh.node_coord[el+1]] for el in range(n_oneD_elements)])
    x = x.flatten()
    y = np.array([[func[el], func[el]] for el in range(n_oneD_elements)])
    y = y.flatten()
    Plot(x,y,"x",labely,title)



def PlotByInterface(func):
    """Plot a vector of values that corresponds to each interface element of the mesh """

    labely = retrieve_name(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    x = [DFMesh.node_coord[el] for el in range(n_oneD_elements+1)]
    y = np.zeros(len(x))
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
    plt.savefig("LOG/average_stress_bar.svg")
    plt.show()



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
    plt.plot(x, Erev, label='Erev')
    plt.plot(x, Econ, label='Econ')
    plt.plot(x, Wext, label='Wext')
    plt.legend()
    plt.savefig("LOG/energies.svg")
    plt.show()



def PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot):
    """Plot variation of energy from time t to t0."""

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
    plt.plot(x, varErev, label='varErev')
    plt.plot(x, varEcon, label='varEcon')
    plt.plot(x, -varWext, label='-varWext')
    plt.plot(x, varEtot, label='varEtot')
    plt.legend()
    plt.savefig("LOG/var_energies.svg")
    plt.show()



def PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot):
    """Plot variation of energy between time steps"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title(str("Power"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Power"))

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, PEpot, label='varEpot')
    plt.plot(x, PEkin, label='varEkin')
    plt.plot(x, PEdis, label='varEdis')
    plt.plot(x, PErev, label='varErev')
    plt.plot(x, PEcon, label='varEcon')
    plt.plot(x, -PWext, label='-varWext')
    plt.plot(x, PEtot, label='varEtot')
    plt.legend()
    plt.savefig("LOG/power.svg")
    plt.show()



def PlotNumberFragments(nfrag):
    """Plot a vector of values that corresponds to the number of fragments at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Number of fragments")
    plt.xlabel("Time (s)")
    plt.ylabel("N")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    y = nfrag
    plt.plot(x, y)
    plt.savefig("LOG/number_fragments.svg")
    plt.show()



def PlotAvgFragmentSize(avg_frag_sizes):
    """Plot a vector of values that corresponds to the fragments length at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Fragments sizes")
    plt.xlabel("Time (s)")
    plt.ylabel("m")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    y = avg_frag_sizes
    plt.plot(x, y)
    plt.savefig("LOG/size_fragments.svg")
    plt.show()



def PlotFragmentSizeHistogram(frag_sizes):

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Fragment size distribution")
    plt.xlabel("Fragment size (m)")
    plt.ylabel("Number of fragments")

    plt.hist(frag_sizes,10)
    plt.savefig("LOG/fragment_size_distribution.svg")
    plt.show()



def PlotVTK(prefix, timestep, u, stress):
    ndofs = len(u)
    filename = prefix + '.' + str(timestep) + '.vtk'
    header = '''# vtk DataFile Version 3.0
Dynamic fragmentation 
ASCII

DATASET UNSTRUCTURED_GRID
'''
    coord = DFMesh.ListDofCoord()
    points = 'POINTS ' + str(len(coord)) + ' float\n' + np.array2string(coord).replace('[','').replace(']','') + '\n'

    cellist = np.zeros((DFMesh.n_el,3),dtype=int)
    for i in range(DFMesh.n_el):
        cellist[i,0] = len(DFMesh.connect[i])
        cellist[i,1] = DFMesh.connect[i][0]
        cellist[i,2] = DFMesh.connect[i][1]

    cells = '\nCELLS ' + str(DFMesh.n_el) + ' ' + str(DFMesh.n_el*3) + '\n' + np.array2string(cellist).replace('[','').replace(']','')+ '\n'

    celltypes = '\nCELL_TYPES ' + str(DFMesh.n_el) + '\n' + '\n'.join(map(str,[3]*DFMesh.n_el))+ '\n'

    displacement = np.array([[u[i],0.,0.] for i in range(len(u))])
    displacement = '\nVECTORS displacement float\n' +  np.array2string(displacement).replace('[','').replace(']','')+ '\n'

    el_avgdisp = [(u[DFMesh.connect[el][0]] + u[DFMesh.connect[el][1]])*0.5 for el in range(DFMesh.n_el)]
    
    avgdisp = np.zeros((ndofs,3))
    for el in range(DFMesh.n_el):
        avgdisp[DFMesh.connect[el][0],0] = el_avgdisp[el]
        avgdisp[DFMesh.connect[el][1],0] = el_avgdisp[el]

    avgdisp = '\nVECTORS AvgDisplacement float\n' +  np.array2string(avgdisp).replace('[','').replace(']','')+ '\n'

    stressplot = 'StressX 1 '+str(DFMesh.n_el)+' float\n' +'\n'.join(map(str,stress[:DFMesh.n_el]))+ '\n'

    output = open(filename, 'w')
    output.write(header)
    output.write(points)
    output.write(cells)
    output.write(celltypes)
    output.write('\nCELL_DATA ' + str(DFMesh.n_el) + '\n')
    output.write('FIELD FieldData 1\n')
    output.write(stressplot)
    output.write('\nPOINT_DATA ' + str(len(u)) + '\n')
    output.write(displacement)
    output.write(avgdisp)
    output.close()



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
    plt.savefig("LOG/convergence_energy.svg")
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
    nnodes = [meshes[i]+1 for i in range(len(meshes))]

    plt.plot(nnodes, energies_un, label='Uniform mesh')
    plt.plot(nnodes, energies_nun, label='Non-uniform mesh')
    plt.legend()
    plt.savefig("LOG/logconvergence_energy.svg")
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
    plt.savefig("LOG/convergence_nfrags.svg")
    plt.show()



def PlotAnalyticals(grady, gc, zmr, values_strainrate):
    """Plot analytical estimation of fragment size given by analytical models.\n
    Arguments:\n
    grady -- estmations by Grady(1982);\n
    gc - estmations by Glen and Chudnovisk(1986);
    zmr -- estmations by Zhou, Molinari and Ramesh (2006)."""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Fragment size")
    plt.xlabel("Strain rate (s-1)")
    plt.ylabel("s")
    plt.xscale("log")
    plt.yscale("log")

    plt.plot(values_strainrate, grady, label="Grady(1982)")
    plt.plot(values_strainrate, gc, label="Glen and Chudnovisk (1986)")
    plt.plot(values_strainrate, zmr, label="Zhou, Molinari and Ramesh (2006)")

    plt.legend()
    plt.show()

        
