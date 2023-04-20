import numpy as np
from matplotlib import pyplot as plt
import pickle


# Read results from previous simulation


def readPreviousResults(previous_simulation_file):
    """Reads previous simulation results saved in the results variable for a given time-step"""

    with open(previous_simulation_file, "rb") as handle:
        previous_results = pickle.load(handle)

    return previous_results


def getResults(results, variable_name):
    """Get a vriable results from a given time-step"""
    for i in range(len(results)):
        if results[i][0] == variable_name:
            variable = results[i][1]
            return variable


def readResultsVariable(file_address, variable_name):
    with open(file_address + variable_name + ".pickle", "rb") as handle:
        variable = pickle.load(handle)
    return variable


def plotResults(results, label_x, label_y):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    nb_simulations = len(results)

    for i in range(nb_simulations):
        name_simulation = results[i][0]
        x_values = results[i][1]
        y_values = results[i][2]
        plt.plot(x_values, y_values, label=name_simulation)
    plt.legend()

    plt.show()
    


def plotCompareSimulations(results, title, label_x, label_y, save_filename):
    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    nb_simulations = len(results)

    for i in range(nb_simulations):
        name_simulation = results[i][0]
        x_values = results[i][1]
        y_values = results[i][2]
        plt.plot(x_values, y_values, label=name_simulation)
    plt.legend()
    # plt.savefig(save_filename)
    plt.show()

def getEnergy(energies, energy_name):

    for i in range(len(energies)):
        if energies[i][0] == energy_name:
            energy = energies[i][1]
            return energy

def computeVarEnergiesCZM(energies, n_steps, n_elements):
    """Returns the variation of energies between the current time step and the time step 0."""

    energy_potential = np.zeros(n_steps)
    energy_kinetic = np.zeros(n_steps)
    energy_dissipated = np.zeros(n_steps)
    energy_contact = np.zeros(n_steps)
    energy_reversible = np.zeros(n_steps)
    energy_total = np.zeros(n_steps)
    external_work = np.zeros(n_steps)
    
    var_energy_potential = np.zeros(n_steps)
    var_energy_kinetic = np.zeros(n_steps)
    var_energy_dissipated = np.zeros(n_steps)
    var_energy_contact = np.zeros(n_steps)
    var_energy_reversible = np.zeros(n_steps)
    var_energy_total = np.zeros(n_steps)
    var_external_work = np.zeros(n_steps)

    h = 50*10**-3 / n_elements
    for i in range(n_steps):

        energy_potential[i] = getEnergy(energies[i], "energy potential") / h
        energy_kinetic[i] = getEnergy(energies[i], "energy kinetic") / h
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated") / h
        energy_reversible[i] = getEnergy(energies[i], "energy reversible") / h
        energy_contact[i] = getEnergy(energies[i], "energy contact") / h
        external_work[i] = getEnergy(energies[i], "external work") / h

    for n in range(1, n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_energy_reversible[n] = energy_reversible[n] - energy_reversible[0]
        var_energy_contact[n] = energy_contact[n] - energy_contact[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (
            var_energy_potential[n]
            + var_energy_kinetic[n]
            + var_energy_dissipated[n]
            + var_energy_reversible[n]
            + var_energy_contact[n]
        )

    var_energies = [
        ["var energy potential", var_energy_potential],
        ["var energy kinetic", var_energy_kinetic],
        ["var energy dissipated", var_energy_dissipated],
        ["var energy reversible", var_energy_reversible],
        ["var energy contact", var_energy_contact],
        ["var external work", var_external_work],
        ["var energy total", var_energy_total],
    ]

    return var_energies


def computeVarEnergiesLipfield(energies, n_steps):
    """Returns the variation of energies between the current time step and the time step 0."""

    energy_potential = np.zeros(n_steps)
    energy_kinetic = np.zeros(n_steps)
    energy_dissipated = np.zeros(n_steps)
    external_work = np.zeros(n_steps)
    
    var_energy_potential = np.zeros(n_steps)
    var_energy_kinetic = np.zeros(n_steps)
    var_energy_dissipated = np.zeros(n_steps)
    var_external_work = np.zeros(n_steps)
    var_energy_total = np.zeros(n_steps)

    for i in range(n_steps):
        energy_potential[i] = getEnergy(energies[i], "energy potential")
        energy_kinetic[i] = getEnergy(energies[i], "energy kinetic")
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated")
        external_work[i] = getEnergy(energies[i], "external work")

    for n in range(1, n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (
            var_energy_potential[n] + var_energy_kinetic[n] + var_energy_dissipated[n]
        )

    var_energies = [
        ["var energy potential", var_energy_potential],
        ["var energy kinetic", var_energy_kinetic],
        ["var energy dissipated", var_energy_dissipated],
        ["var external work", var_external_work],
        ["var energy total", var_energy_total],
    ]

    return var_energies


        
def plotVarEnergiesCZM(var_energies, time, save_filename, title):
    """Plot variation of energy from time t to t0."""

    var_energy_potential = getEnergy(var_energies, "var energy potential") / 10**3
    var_energy_kinetic = getEnergy(var_energies, "var energy kinetic") / 10**3
    var_energy_dissipated = getEnergy(
        var_energies, "var energy dissipated"
    ) / 10**3
    var_energy_reversible = getEnergy(
        var_energies, "var energy reversible"
    ) / 10**3
    var_energy_contact = getEnergy(var_energies, "var energy contact") / 10**3
    var_external_work = getEnergy(var_energies, "var external work") / 10**3
    var_energy_total = getEnergy(var_energies, "var energy total") / 10**3

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(title)
    plt.xlabel(str("Time (s)"))
    plt.ylabel("Variation of energy ($ kJ/ {m^2} $)")
    plt.ylim(-30,20)
    plt.xlim(0, 3.5*10**-7)

    x = time 

    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, var_energy_reversible, label="varErev")
    plt.plot(x, var_energy_contact, label="varEcon")
    plt.plot(x, -var_external_work, label="-varWext")
    plt.plot(x, var_energy_total, label="varEtot")
    plt.legend()
    # plt.savefig(save_filename)
    plt.show()


def plotVarEnergiesLipfield(var_energies, time, save_filename, title):
    """Plot variation of energy from time t to t0."""

    var_energy_potential = getEnergy(var_energies, "var energy potential") / 10**3
    var_energy_kinetic = getEnergy(var_energies, "var energy kinetic") / 10**3
    var_energy_dissipated = getEnergy(
        var_energies, "var energy dissipated"
    ) / 10**3
    var_external_work = getEnergy(var_energies, "var external work")/ 10**3
    var_energy_total = getEnergy(var_energies, "var energy total")/ 10**3

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str(title))
    plt.xlabel(str("Time (s)"))
    plt.ylabel("Variation of energy ($ kJ/ {m^2} $)")
    plt.ylim(-30,20)
    plt.xlim(0, 3.5*10**-7)

    x = time

    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, -var_external_work, label="-varWext", color='saddlebrown')
    plt.plot(x, var_energy_total, label="varEtot", color='orchid' )
    plt.legend()
    # plt.savefig(save_filename)
    plt.show()