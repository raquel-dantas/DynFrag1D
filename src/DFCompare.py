import numpy as np
from matplotlib import pyplot as plt
import pickle


# Read results from previous simulation


def readResults(file_address, variable_name):
    file = file_address + variable_name + ".pickle"
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
    plt.savefig(save_filename)
    plt.show()




def plotConvergence(results,meshes, title, label_x, label_y, save_filename):

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
    plt.savefig(save_filename)
    plt.show()