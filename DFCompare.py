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
    """Get a vriable results from a given time-step """
    for i in range(len(results)):
        if results[i][0] == variable_name:
            variable = results[i][1]
            return variable


def readResultsVariable(file_address, variable_name):
    with open(file_address + variable_name + ".pickle", "rb") as handle:
        variable = pickle.load(handle)
    return variable


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




def plotConvergence(results, meshes, title, label_x, label_y, save_filename):

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    nb_simulations = len(meshes)
    x_values = np.zeros(nb_simulations)
    final_values = np.zeros(nb_simulations)

    for i in range(nb_simulations):
        x_values[i] = meshes[i]
        final_values[i] = max(results[i][2])
    plt.plot(x_values, final_values, marker='.')

    plt.legend()
    # plt.savefig(save_filename)
    plt.show()



def plotConvergenceComparison(results1, results2, meshes, title, label_x, label_y, save_filename):

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
    plt.plot(x_values, final_values1, marker='.')
    plt.plot(x_values, final_values2, marker='.')

    plt.legend()
    # plt.savefig(save_filename)
    plt.show()