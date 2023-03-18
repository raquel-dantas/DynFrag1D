import numpy as np
import pickle
import inspect
from matplotlib import pyplot as plt

import DFMesh_aka as DFMesh
import DFPostProcess_aka as DFPostProcess
import DFModel_aka as DFModel


def plot(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.plot(x, y)
    # plt.show()


def plotlog(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(x, y)
    # plt.show()


def retrieveName(var):
    """Gets the name of the argument passed to it, as you coded it in your python script"""

    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]


def plotByCoord(func, labely):
    """Plot a vector of values that corresponds to each DOF of the mesh"""

    title = labely
    x = np.array([x for x, y in DFMesh.getNodes()])
    y = func
    x = x.flatten()
    y = y.flatten()
    plot(x, y, "x", labely, title)


def plotAverageStressBar(average_stress_bar):
    """Plot a vector of values that corresponds to the average stress between all elements at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Average stress bar")
    plt.xlabel("Time (s)")
    plt.ylabel("Average Stress (Pa)")

    x = np.linspace(0, DFMesh.time_simulation, DFModel.n_steps)
    y = average_stress_bar
    plt.plot(x, y)
    plt.savefig(DFMesh.filepath + "aka_average_stress_bar.svg")


def plotEnergies(energies):
    """Plot energies values per time"""

    energy_potential = DFPostProcess.getEnergy(energies, "energy potential")
    energy_kinetic = DFPostProcess.getEnergy(energies, "energy kinetic")
    energy_dissipated = DFPostProcess.getEnergy(energies, "energy dissipated")
    energy_reversible = DFPostProcess.getEnergy(energies, "energy reversible")
    energy_contact = DFPostProcess.getEnergy(energies, "energy contact")
    external_work = DFPostProcess.getEnergy(energies, "external work")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Energies")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy  (N/m)")

    x = np.linspace(0, DFMesh.time_simulation, DFModel.n_steps)
    plt.plot(x, energy_potential, label="Epot")
    plt.plot(x, energy_kinetic, label="Ekin")
    plt.plot(x, energy_dissipated, label="Edis")
    plt.plot(x, energy_reversible, label="Erev")
    plt.plot(x, energy_contact, label="Econ")
    plt.plot(x, external_work, label="Wext")
    plt.legend()

    plt.savefig(DFMesh.filepath + "aka_energies_dynfrag_akantu.svg")


def plotVarEnergies(var_energies):
    """Plot variation of energy from time t to t0."""

    var_energy_potential = DFPostProcess.getEnergy(var_energies, "var energy potential")
    var_energy_kinetic = DFPostProcess.getEnergy(var_energies, "var energy kinetic")
    var_energy_dissipated = DFPostProcess.getEnergy(
        var_energies, "var energy dissipated"
    )
    var_energy_reversible = DFPostProcess.getEnergy(
        var_energies, "var energy reversible"
    )
    var_energy_contact = DFPostProcess.getEnergy(var_energies, "var energy contact")
    var_external_work = DFPostProcess.getEnergy(var_energies, "var external work")
    var_energy_total = DFPostProcess.getEnergy(var_energies, "var energy total")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Variation of energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Variation of energy"))

    x = np.linspace(0, DFMesh.time_simulation, DFModel.n_steps)
    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, var_energy_reversible, label="varErev")
    plt.plot(x, var_energy_contact, label="varEcon")
    plt.plot(x, -var_external_work, label="-varWext")
    plt.plot(x, var_energy_total, label="varEtot")
    plt.legend()
    plt.savefig(DFMesh.filepath + "aka_var_energies.svg")


def plotPower(power):
    """Plot variation of energy between time steps"""

    power_potential = DFPostProcess.getEnergy(power, "power potential")
    power_kinetic = DFPostProcess.getEnergy(power, "power kinetic")
    power_dissipated = DFPostProcess.getEnergy(power, "power dissipated")
    power_reversible = DFPostProcess.getEnergy(power, "power reversible")
    power_contact = DFPostProcess.getEnergy(power, "power contact")
    power_external_work = DFPostProcess.getEnergy(power, "power external work")
    power_total = DFPostProcess.getEnergy(power, "power total")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Power"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Power"))

    x = np.linspace(0, DFMesh.time_simulation, DFModel.n_steps)
    plt.plot(x, power_potential, label="pEpot")
    plt.plot(x, power_kinetic, label="pEkin")
    plt.plot(x, power_dissipated, label="pEdis")
    plt.plot(x, power_reversible, label="pErev")
    plt.plot(x, power_contact, label="pEcon")
    plt.plot(x, -power_external_work, label="-pWext")
    plt.plot(x, power_total, label="pEtot")
    plt.legend()
    plt.savefig(DFMesh.filepath + "aka_power.svg")


def plotDamage(d):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Damage")
    plt.xlabel("x")
    plt.ylabel("d")

    x = np.array([x for x, y in DFModel.facets_coords])
    plt.plot(x,d)
    plt.savefig(DFMesh.filepath + "aka_damage.svg")



def plotNumberFragments(nfrag):
    """Plot a vector of values that corresponds to the number of fragments at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Number of fragments")
    plt.xlabel("Time (s)")
    plt.ylabel("N")

    x = np.linspace(0, DFMesh.time_simulation, DFModel.n_steps)
    y = nfrag
    plt.plot(x, y)
    plt.savefig(DFMesh.filepath + "aka_number_fragments.svg")


def plotAvgFragmentSize(avg_frag_sizes):
    """Plot a vector of values that corresponds to the fragments length at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Fragments sizes")
    plt.xlabel("Time (s)")
    plt.ylabel("m")

    x = np.linspace(0, DFModel.time_simulation, DFModel.n_steps)
    y = avg_frag_sizes
    plt.plot(x, y)
    plt.savefig(DFMesh.filepath + "aka_avg_size_fragments.svg.svg")
    plt.show()


def plotFragmentSizeHistogram(frag_sizes, n_columns):

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Fragment size distribution")
    plt.xlabel("Fragment size (m)")
    plt.ylabel("Number of fragments")

    plt.hist(frag_sizes, n_columns)

    plt.savefig(DFMesh.filepath + "aka_fragment_size_histogram.svg")


def plotConvergenceEnergy(energies_un, energies_nun, meshes):
    """Plot total number of nodes of the mesh x energy dissipated.\n
    Arguments: \n
    energies -- dissipated energies for different mesh sizes;\n
    meshes -- number of linear elements."""

    plt.title(str("Energy convergence"))
    plt.xlabel(str("Number of nodes"))
    plt.ylabel(str("Energy dissipated (N/m)"))
    nnodes = [meshes[i] + 1 for i in range(len(meshes))]

    plt.plot(nnodes, energies_un, label="Uniform mesh")
    plt.plot(nnodes, energies_nun, label="Non-uniform mesh")
    plt.legend()
    plt.savefig("LOG/convergence_energy_dynfrag_akantu.svg")
    plt.show()


def plotLogConvergenceEnergy(energies_un, energies_nun, meshes):
    """Plot total number of nodes of the mesh x energy dissipated.\n
    Arguments: \n
    energies -- dissipated energies for different mesh sizes;\n
    meshes -- number of linear elements."""

    plt.title(str("Energy convergence"))
    plt.xlabel(str("Number of nodes"))
    plt.ylabel(str("Energy dissipated (N/m)"))
    plt.xscale("log")
    plt.yscale("log")

    plt.xlim(10**2, 5 * 10**5)

    nnodes = [meshes[i] + 1 for i in range(len(meshes))]

    plt.plot(nnodes, energies_un, label="Uniform mesh")
    plt.plot(nnodes, energies_nun, label="Non-uniform mesh")
    plt.legend()
    plt.savefig("LOG/logconvergence_energy_dynfrag_akantu.svg")
    plt.show()


def plotConvergenceNumfrag(nfrags_un, nfrags_nun, meshes):
    """Plot total number of nodes of the mesh x final number of fragments.\n
    Arguments: \n
    energies -- dissipated energies for different mesh sizes;\n
    meshes -- number of linear elements."""

    plt.title(str("Nfrag convergence"))
    plt.xlabel(str("Number of nodes"))
    plt.ylabel(str("Number of fragments"))
    nnodes = [meshes[i] + 1 for i in range(len(meshes))]

    plt.plot(nnodes, nfrags_un, label="Uniform mesh")
    plt.plot(nnodes, nfrags_nun, label="Non-uniform mesh")
    plt.legend()
    plt.savefig("LOG/convergence_nfrags_dynfrag_akantu.svg")
    plt.show()


def addPlotVtk():
    DFModel.dynfrag.setBaseName('bar')
    DFModel.dynfrag.addDumpFieldVector('displacement')
    DFModel.dynfrag.addDumpFieldVector('velocity')
    DFModel.dynfrag.addDumpField('strain')
    DFModel.dynfrag.addDumpField('stress')
    DFModel.dynfrag.addDumpField('blocked_dofs')
    DFModel.dynfrag.addDumpField('material_index')

    # VTK plot for Cohesive model
    DFModel.dynfrag.setBaseNameToDumper('cohesive elements', 'cohesive')
    DFModel.dynfrag.addDumpFieldVectorToDumper('cohesive elements', 'displacement')
    DFModel.dynfrag.addDumpFieldToDumper('cohesive elements', 'damage')
    DFModel.dynfrag.addDumpFieldVectorToDumper('cohesive elements', 'tractions')
    DFModel.dynfrag.addDumpFieldVectorToDumper('cohesive elements', 'opening')

def addVtkFiles(n_steps_to_save):
    if n_steps_to_save % 100 == 0:
        DFModel.dynfrag.dump()
        DFModel.dynfrag.dump("cohesive elements")


def saveResults(variable):
    variable_name = retrieveName(variable)[0]
    with open(DFMesh.filepath + "aka_" + variable_name + ".pickle", "wb") as handle:
        pickle.dump(variable, handle, protocol=pickle.HIGHEST_PROTOCOL)


def saveResultsCurrentStep(results, current_step):

    step_id = "step_" + str(current_step) + "_"
    with open(DFMesh.filepath + "lipfield_" + step_id + ".pickle", "wb") as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)