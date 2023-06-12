import numpy as np
from matplotlib import pyplot as plt

import DFProcessResults


# Plot functions


def plotResults(
    results_variable,
    x_values,
    label_x: str,
    label_y: str,
    plot_title: str,
    save_plot: bool,
    save_filename: str,
):
    """Plot results values for a given variable
    """

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(plot_title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    plt.plot( x_values, results_variable)
    plt.legend()
    if save_plot == True:
        plt.savefig(save_filename + ".svg")
    plt.show()


def plotResultsComparison(
    results: list,
    label_x: str,
    label_y: str,
    plot_title: str,
    save_plot: bool,
    save_filename: str,
):
    """Plot results values in the following form: \n
    results = [label_simulation , x_values, y_values]
    """

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(plot_title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    nb_simulations = len(results)

    for i in range(nb_simulations):
        name_simulation = results[i][0]
        x_values = results[i][1]
        y_values = results[i][2] 
        plt.plot(x_values, y_values, label=name_simulation)
    plt.legend()
    if save_plot == True:
        plt.savefig(save_filename + ".svg")
    plt.show()


def plotVarEnergiesCZM(
    var_energies, time_simulation, plot_title: str, save_plot: bool, save_filename: str
):
    """Plot variation of energy from time t to t0."""

    unit_convertion = 1e-3  # From J to kJ

    var_energy_potential = (
        DFProcessResults.getEnergy(var_energies, "var energy potential")
        * unit_convertion
    )
    var_energy_kinetic = (
        DFProcessResults.getEnergy(var_energies, "var energy kinetic") * unit_convertion
    )
    var_energy_dissipated = (
        DFProcessResults.getEnergy(var_energies, "var energy dissipated")
        * unit_convertion
    )
    var_energy_reversible = (
        DFProcessResults.getEnergy(var_energies, "var energy reversible")
        * unit_convertion
    )
    var_energy_contact = (
        DFProcessResults.getEnergy(var_energies, "var energy contact") * unit_convertion
    )
    var_external_work = (
        DFProcessResults.getEnergy(var_energies, "var external work") * unit_convertion
    )
    var_energy_total = (
        DFProcessResults.getEnergy(var_energies, "var energy total") * unit_convertion
    )

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(plot_title)
    plt.xlabel(str("Time (s)"))
    plt.ylabel("Variation of energy ($ kJ/ {m^2} $)")

    x = time_simulation

    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, var_energy_reversible, label="varErev")
    plt.plot(x, var_energy_contact, label="varEcon")
    plt.plot(x, -var_external_work, label="-varWext")
    plt.plot(x, var_energy_total, label="varEtot")
    plt.legend()
    if save_plot == True:
        plt.savefig(save_filename + ".svg")
    plt.show()


def plotVarEnergiesLipfield(
    var_energies, time_simulation, plot_title: str, save_plot: bool, save_filename: str
):
    """Plot variation of energy from time t to t0."""

    unit_convertion = 1e-3  # From J to kJ

    var_energy_potential = (
        DFProcessResults.getEnergy(var_energies, "var energy potential")
        * unit_convertion
    )
    var_energy_kinetic = (
        DFProcessResults.getEnergy(var_energies, "var energy kinetic") * unit_convertion
    )
    var_energy_dissipated = (
        DFProcessResults.getEnergy(var_energies, "var energy dissipated")
        * unit_convertion
    )
    var_external_work = (
        DFProcessResults.getEnergy(var_energies, "var external work") * unit_convertion
    )
    var_energy_total = (
        DFProcessResults.getEnergy(var_energies, "var energy total") * unit_convertion
    )

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str(plot_title))
    plt.xlabel(str("Time (s)"))
    plt.ylabel("Variation of energy ($ kJ/ {m^2} $)")

    x = time_simulation

    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, -var_external_work, label="-varWext", color="saddlebrown")
    plt.plot(x, var_energy_total, label="varEtot", color="orchid")
    plt.legend()
    if save_plot == True:
        plt.savefig(save_filename + ".svg")
    plt.show()


def plotFragmentSizeHistogram(
    fragment_sizes, n_cols, plot_title: str, save_plot: bool, save_filename: str
):
    unit_convertion = 1e3  # From m to mm

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(plot_title)
    plt.xlabel("Fragment size (mm)")
    plt.ylabel("Number of fragments")
    plt.xlim(xmin=0, xmax=1.4)
    plt.ylim(ymin=0, ymax=28)

    plt.hist(
        fragment_sizes * unit_convertion,
        n_cols,
        histtype="stepfilled",
        alpha=0.4,
    )
    plt.legend()
    if save_plot == True:
        plt.savefig(save_filename + ".svg")
    plt.show()


def plotConvergence(
    results: list,
    meshes: list,
    label_y: str,
    plot_title: str,
    save_plot: bool,
    save_filename: str,
):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(plot_title)
    plt.xlabel("Number of elements")
    plt.ylabel(label_y)

    nb_simulations = len(results)
    x_values = np.zeros(nb_simulations)
    final_values = np.zeros(nb_simulations)

    for i in range(nb_simulations):
        x_values[i] = meshes[i]
        final_values[i] = max(results[i][2])
    plt.plot(x_values, final_values, marker=".")

    plt.legend()
    if save_plot == True:
        plt.savefig(save_filename + ".svg")
    plt.show()


def plotConvergenceComparison(
    results1, results2, meshes, title, label_x, label_y, save_filename
):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    nb_simulations = len(meshes)
    x_values = np.zeros(nb_simulations)
    final_values1 = np.zeros(nb_simulations)
    final_values2 = np.zeros(nb_simulations)

    for i in range(nb_simulations):
        x_values[i] = meshes[i]
        final_values1[i] = max(results1[i][2])
        final_values2[i] = max(results2[i][2])
    plt.plot(x_values, final_values1, marker=".", label="Uniform mesh")
    plt.plot(x_values, final_values2, marker=".", label="Non-uniform mesh")

    plt.legend()
    plt.savefig(save_filename)
    plt.show()
