import numpy as np
from matplotlib import pyplot as plt
import inspect
import pickle

import DFMesh
import DFFem
import DFPostProcess


def plot(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    # for col in range(len(y)):
    #     plt.plot(x, y[col])
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


def plotScatter(x, y, labelx, labely, title):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str(title))
    plt.xlabel(str(labelx))
    plt.ylabel(str(labely))
    plt.scatter(x, y)
    # plt.show()


def retrieveName(var):
    """Gets the name of the argument passed to it, as you coded it in your python script"""

    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]


def plotByDOF(func):
    """Plot a vector of values that corresponds to each DOF of the mesh"""

    labely = retrieveName(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    x = np.array(
        [
            [DFMesh.node_coord[el], DFMesh.node_coord[el + 1]]
            for el in range(n_oneD_elements)
        ]
    )
    x = x.flatten()
    y = np.array(
        [
            [func[DFFem.Gl_index(el, 0)], func[DFFem.Gl_index(el, 1)]]
            for el in range(n_oneD_elements)
        ]
    )
    y = y.flatten()

    plot(x, y, "x", labely, title)


def plotByElement(func):
    """Plot a vector of values that corresponds to each element of the mesh"""

    labely = retrieveName(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    x = np.array(
        [
            [DFMesh.node_coord[el], DFMesh.node_coord[el + 1]]
            for el in range(n_oneD_elements)
        ]
    )
    x = x.flatten()
    y = np.array([[func[el], func[el]] for el in range(n_oneD_elements)])
    y = y.flatten()

    plot(x, y, "x", labely, title)


def plotByInterface(func):
    """Plot a vector of values that corresponds to each interface element of the mesh"""

    labely = retrieveName(func)[0]
    title = labely
    n_oneD_elements = sum(1 for el in DFMesh.materials if el == 0)
    x = [DFMesh.node_coord[el] for el in range(n_oneD_elements + 1)]
    y = np.zeros(len(x))
    for el in range(n_oneD_elements, len(DFMesh.materials)):
        if DFMesh.materials[el] == 1:
            j = DFMesh.connect[el][0]
            y[j] = func[el]

    plotScatter(x, y, "x", labely, title)


def plotByIntPoint(func):
    """Plot a vector of values that corresponds to each integration point"""

    labely = retrieveName(func)[0]
    plot(DFMesh.intpoint_coord, func, "x", labely, labely)


def plotDamage(damage):


    if DFMesh.use_cohesive_elements == True:
        plotByInterface(damage)
        plt.savefig("LOG/czm_damage.svg")

    if DFMesh.use_lipfield == True:
        plotByElement(damage)
        plt.savefig("LOG/lipfield_damage.svg")


def plotAverageStressBar(average_stress_bar):
    """Plot a vector of values that corresponds to the average stress between all elements at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Average stress bar")
    plt.xlabel("Time (s)")
    plt.ylabel("Average Stress (Pa)")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    y = average_stress_bar

    plt.plot(x, y)
    if DFMesh.use_cohesive_elements == True:
        plt.savefig("LOG/czm_average_stress_bar.svg")
    if DFMesh.use_lipfield == True:
        plt.savefig("LOG/lipfield_average_stress_bar.svg")


def plotStressByTime(stress_evolution):
    """Plot the stress on the elements at each time step"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Stress evolution")
    plt.xlabel("Time (s)")
    plt.ylabel("Stress (Pa)")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)

    for el in range(len(DFMesh.materials)):
        plt.plot(x, stress_evolution[el], label=el)
    plt.legend()
    plt.show()


def plotEnergiesCZM(energies):
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

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, energy_potential, label="Epot")
    plt.plot(x, energy_kinetic, label="Ekin")
    plt.plot(x, energy_dissipated, label="Edis")
    plt.plot(x, energy_reversible, label="Erev")
    plt.plot(x, energy_contact, label="Econ")
    plt.plot(x, external_work, label="Wext")
    plt.legend()

    plt.savefig("LOG/czm_energies.svg")
        


def plotEnergiesLipfield(energies):
    """Plot energies values x time"""

    energy_potential = DFPostProcess.getEnergy(energies, "energy potential")
    energy_kinetic = DFPostProcess.getEnergy(energies, "energy kinetic")
    energy_dissipated = DFPostProcess.getEnergy(energies, "energy dissipated")
    external_work = DFPostProcess.getEnergy(energies, "external work")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Energies")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy  (N/m)")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)

    plt.plot(x, energy_potential, label="Potential")
    plt.plot(x, energy_kinetic, label="Ekin")
    plt.plot(x, energy_dissipated, label="Edis")
    plt.plot(x, external_work, label="Wext")
    plt.legend()
    plt.savefig("LOG/lipfield_energies.svg")


def plotEnergies(energies):

    if DFMesh.use_cohesive_elements == True:
        plotEnergiesCZM(energies)

    if DFMesh.use_lipfield == True:
        plotEnergiesLipfield(energies)


def plotVarEnergiesCZM(var_energies):
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

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, var_energy_reversible, label="varErev")
    plt.plot(x, var_energy_contact, label="varEcon")
    plt.plot(x, -var_external_work, label="-varWext")
    plt.plot(x, var_energy_total, label="varEtot")
    plt.legend()
    plt.savefig("LOG/czm_var_energies.svg")
    # plt.show()


def plotVarEnergiesLipfield(var_energies):
    """Plot variation of energy from time t to t0."""

    var_energy_potential = DFPostProcess.getEnergy(var_energies, "var energy potential")
    var_energy_kinetic = DFPostProcess.getEnergy(var_energies, "var energy kinetic")
    var_energy_dissipated = DFPostProcess.getEnergy(
        var_energies, "var energy dissipated"
    )
    var_external_work = DFPostProcess.getEnergy(var_energies, "var external work")
    var_energy_total = DFPostProcess.getEnergy(var_energies, "var energy total")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Variation of energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Variation of energy"))

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)

    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, -var_external_work, label="-varWext")
    plt.plot(x, var_energy_total, label="varEtot")
    plt.legend()
    plt.savefig("LOG/lipfield_var_energies.svg")
    # plt.show()


def plotVarEnergies(var_energies):

    if DFMesh.use_cohesive_elements == True:
        plotVarEnergiesCZM(var_energies)

    if DFMesh.use_lipfield == True:
        plotVarEnergiesLipfield(var_energies)


def plotPowerCZM(power):
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

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, power_potential, label="pEpot")
    plt.plot(x, power_kinetic, label="pEkin")
    plt.plot(x, power_dissipated, label="pEdis")
    plt.plot(x, power_reversible, label="pErev")
    plt.plot(x, power_contact, label="pEcon")
    plt.plot(x, -power_external_work, label="-pWext")
    plt.plot(x, power_total, label="pEtot")
    plt.legend()
    plt.savefig("LOG/czm_power.svg")


def plotPowerLipfield(power):
    """Plot variation of energy between time steps"""

    power_potential = DFPostProcess.getEnergy(power, "power potential")
    power_kinetic = DFPostProcess.getEnergy(power, "power kinetic")
    power_dissipated = DFPostProcess.getEnergy(power, "power dissipated")
    power_external_work = DFPostProcess.getEnergy(power, "power external work")
    power_total = DFPostProcess.getEnergy(power, "power total")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Power"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Power"))

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    plt.plot(x, power_potential, label="pEpot")
    plt.plot(x, power_kinetic, label="pEkin")
    plt.plot(x, power_dissipated, label="pEdis")
    plt.plot(x, -power_external_work, label="-pWext")
    plt.plot(x, power_total, label="pEtot")
    plt.legend()
    plt.savefig("LOG/lipfield_power.svg")


def plotPower(power):

    if DFMesh.use_cohesive_elements == True:
        plotPowerCZM(power)

    if DFMesh.use_lipfield == True:
        plotPowerLipfield(power)


def plotNumberFragments(nfrag):
    """Plot a vector of values that corresponds to the number of fragments at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Number of fragments")
    plt.xlabel("Time (s)")
    plt.ylabel("N")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    y = nfrag

    plt.plot(x, y)
    if DFMesh.use_cohesive_elements == True:
        plt.savefig("LOG/czm_number_fragments.svg")
    if DFMesh.use_lipfield == True:
        plt.savefig("LOG/lipfield_number_fragments.svg")


def plotAvgFragmentSize(avg_frag_sizes):
    """Plot a vector of values that corresponds to the fragments length at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Fragments sizes")
    plt.xlabel("Time (s)")
    plt.ylabel("m")

    x = np.linspace(0, DFMesh.time_simulation, DFMesh.n_steps)
    y = avg_frag_sizes

    plt.plot(x, y)

    if DFMesh.use_cohesive_elements == True:
        plt.savefig("LOG/czm_avg_size_fragments.svg")
    if DFMesh.use_lipfield == True:
        plt.savefig("LOG/lipfield_avg_size_fragments.svg")

def plotFragmentSizeHistogram(frag_sizes, n_columns):

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Fragment size distribution")
    plt.xlabel("Fragment size (m)")
    plt.ylabel("Number of fragments")

    plt.hist(frag_sizes, n_columns)

    if DFMesh.use_cohesive_elements == True:
        plt.savefig("LOG/czm_fragment_size_histogram.svg")

    if DFMesh.use_lipfield == True:
        plt.savefig("LOG/lipfield_fragment_size_histogram.svg")


def plotVTK(prefix, timestep, u, stress):
    ndofs = len(u)
    filename = prefix + "." + str(timestep) + ".vtk"
    header = """# vtk DataFile Version 3.0
Dynamic fragmentation 
ASCII

DATASET UNSTRUCTURED_GRID
"""
    coord = DFMesh.listDofCoord()
    points = (
        "POINTS "
        + str(len(coord))
        + " float\n"
        + np.array2string(coord).replace("[", "").replace("]", "")
        + "\n"
    )

    cellist = np.zeros((DFMesh.n_elements, 3), dtype=int)
    for i in range(DFMesh.n_elements):
        cellist[i, 0] = len(DFMesh.connect[i])
        cellist[i, 1] = DFMesh.connect[i][0]
        cellist[i, 2] = DFMesh.connect[i][1]

    cells = (
        "\nCELLS "
        + str(DFMesh.n_elements)
        + " "
        + str(DFMesh.n_elements * 3)
        + "\n"
        + np.array2string(cellist).replace("[", "").replace("]", "")
        + "\n"
    )

    celltypes = (
        "\nCELL_TYPES "
        + str(DFMesh.n_elements)
        + "\n"
        + "\n".join(map(str, [3] * DFMesh.n_elements))
        + "\n"
    )

    displacement = np.array([[u[i], 0.0, 0.0] for i in range(len(u))])
    displacement = (
        "\nVECTORS displacement float\n"
        + np.array2string(displacement).replace("[", "").replace("]", "")
        + "\n"
    )

    el_avgdisp = [
        (u[DFMesh.connect[el][0]] + u[DFMesh.connect[el][1]]) * 0.5
        for el in range(DFMesh.n_elements)
    ]

    avgdisp = np.zeros((ndofs, 3))
    for el in range(DFMesh.n_elements):
        avgdisp[DFMesh.connect[el][0], 0] = el_avgdisp[el]
        avgdisp[DFMesh.connect[el][1], 0] = el_avgdisp[el]

    avgdisp = (
        "\nVECTORS AvgDisplacement float\n"
        + np.array2string(avgdisp).replace("[", "").replace("]", "")
        + "\n"
    )

    stressplot = (
        "StressX 1 "
        + str(DFMesh.n_elements)
        + " float\n"
        + "\n".join(map(str, stress[: DFMesh.n_elements]))
        + "\n"
    )

    output = open(filename, "w")
    output.write(header)
    output.write(points)
    output.write(cells)
    output.write(celltypes)
    output.write("\nCELL_DATA " + str(DFMesh.n_elements) + "\n")
    output.write("FIELD FieldData 1\n")
    output.write(stressplot)
    output.write("\nPOINT_DATA " + str(len(u)) + "\n")
    output.write(displacement)
    output.write(avgdisp)
    output.close()


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
    plt.savefig("LOG/convergence_energy.svg")
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
    plt.savefig("LOG/logconvergence_energy.svg")
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
    plt.savefig("LOG/convergence_nfrags.svg")
    plt.show()


def plotLogAnalyticals(grady, gc, zmr, values_strainrate):
    """Plot analytical estimation of fragment size given by analytical models.\n
    Arguments:\n
    grady -- estmations by Grady(1982);\n
    gc - estmations by Glen and Chudnovisk(1986);\n
    zmr -- estmations by Zhou, Molinari and Ramesh (2006)."""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
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


def saveResultsCZM(variable_name, variable):
    with open("LOG/czm_" + variable_name + ".pickle", "wb") as handle:
        pickle.dump(variable, handle, protocol=pickle.HIGHEST_PROTOCOL)


def saveResultsLipfield(variable_name, variable):
    with open("LOG/lipfield_" + variable_name + ".pickle", "wb") as handle:
        pickle.dump(variable, handle, protocol=pickle.HIGHEST_PROTOCOL)


def saveResults(variable):

    variable_name = retrieveName(variable)[0]
    if DFMesh.use_cohesive_elements == True:
        saveResultsCZM(variable_name, variable)

    if DFMesh.use_lipfield == True:
        saveResultsLipfield(variable_name, variable)
