import numpy as np
import pickle
from matplotlib import pyplot as plt
from itertools import groupby
from operator import itemgetter


# Read results functions


def readResultsTimeStep(previous_simulation_file: str):
    """Returns the a list with previous simulation results saved in a .pickle for a time-step"""

    with open(previous_simulation_file, "rb") as handle:
        previous_results = pickle.load(handle)

    return previous_results


def getResults(results: list, variable_name: str):
    """Get a variable results from a given time-step"""

    for i in range(len(results)):
        if results[i][0] == variable_name:
            variable = results[i][1]
            return variable


def getResultsAllSteps(file_path: str, n_files, n_steps):
    """Returns the results for all time-steps for average stress at the bar, energies values and number of fragments."""

    damage = []
    avg_stress_bar = np.zeros(n_files)
    energies = []
    n_fragments = np.zeros(n_files)

    i = 0
    for n in range(n_steps):
        if n%10==0: 
            results = readResultsTimeStep(file_path + "step_" + str(n) + "_.pickle")
            avg_stress_bar[i] = getResults(results, "avg_stress_bar")
            energies.append(getResults(results, "energies"))
            damage.append(getResults(results, "damage"))
            n_fragments[i] = getNumberFragments(damage[i])
            i += 1

    return avg_stress_bar, energies, n_fragments


def getResultsAllStepsCZM_dumptype2(file_path: str, n_files, n_steps):
    """Returns the results for all time-steps for average stress at the bar, energies values and number of fragments."""

    damage = []
    avg_stress_bar = np.zeros(n_files)
    energies = []
    n_fragments = np.zeros(n_files)

    i = 0
    for n in range(n_steps):
        if n%20==0: 
            results = readResultsTimeStep(file_path + "step_" + str(n) + "_.pickle")
            avg_stress_bar[i] = getResults(results, "avg_stress_bar")
            energies.append(getResults(results, "energies"))
            n_fragments[i] = getResults(results, "n_fragments")
            i += 1

    return avg_stress_bar, energies, n_fragments

def getResultsAllStepsCZM(file_path: str, n_files, n_steps):
    """Returns the results for all time-steps for average stress at the bar, energies values and number of fragments."""

    damage = []
    avg_stress_bar = np.zeros(n_files)
    energies = []
    n_fragments = np.zeros(n_files)

    i = 0
    for n in range(n_steps):
        if n%10==0: 
            results = readResultsTimeStep(file_path + "step_" + str(n) + "_.pickle")
            avg_stress_bar[i] = getResults(results, "avg_stress_bar")
            energies.append(getResults(results, "energies"))
            n_fragments[i] = getResults(results, "n_fragments")
            i += 1

    return avg_stress_bar, energies, n_fragments


def getResultFinalStep(variable_name: str, file_path: str, n_steps):
    
    results = readResultsTimeStep(file_path + "step_" + str(n_steps) + "_.pickle")
    variable = getResults(results, variable_name)

    return variable



def getTimeData(file_address: str):
    """Return the time data from a simulation in the followin form: \n
    time_data = [time_simulation, time_increment, number_time_steps]"""

    with open(file_address + "time_data" + ".pickle", "rb") as handle:
        time_data = pickle.load(handle)
    return time_data


#  Post process functions


def getEnergy(energies, energy_name):
    for i in range(len(energies)):
        if energies[i][0] == energy_name:
            energy = energies[i][1]
            return energy


def getDissipatedEnergy(energies, n_files):
    
    energy_dissipated = np.zeros(n_files)
    for i in range(n_files):
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated")
    return energy_dissipated





def getNumberFragments(damage):
    """Compute the number of fragments based on the damage field."""

    elements_full_damage = []

    for el in range(len(damage)):
        if damage[el] > 0.99:
            elements_full_damage.append(el)

        # Divide groups of consecutives elements with damage > 0.99
        cracks_full_damage = []
        for i, subgroup in groupby(
            enumerate(elements_full_damage), lambda index: index[0] - index[1]
        ):
            cracks_full_damage.append(list(map(itemgetter(1), subgroup)))

        n_fragments = 1 + len(cracks_full_damage)

    return n_fragments


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

    h = (
        50 * 10**-3 / n_elements
    )  # Correction of unities to compare with src_akantu (J/m2)

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


def computeVarEnergiesLipfield(energies, n_files):
    """Returns the variation of energies between the current time step and the time step 0."""

    energy_potential = np.zeros(n_files)
    energy_kinetic = np.zeros(n_files)
    energy_dissipated = np.zeros(n_files)
    external_work = np.zeros(n_files)

    var_energy_potential = np.zeros(n_files)
    var_energy_kinetic = np.zeros(n_files)
    var_energy_dissipated = np.zeros(n_files)
    var_external_work = np.zeros(n_files)
    var_energy_total = np.zeros(n_files)

    for i in range(n_files):
        energy_potential[i] = getEnergy(energies[i], "energy potential")
        energy_kinetic[i] = getEnergy(energies[i], "energy kinetic")
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated")
        external_work[i] = getEnergy(energies[i], "external work")

    for n in range(1, n_files):
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
