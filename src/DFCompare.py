import numpy as np
from matplotlib import pyplot as plt
import pickle


# Read results from previous simulation


def readResults(file_address, variable_name):
    file = file_address + variable_name + ".pickle"
    with open(file_address + variable_name + ".pickle", "rb") as handle:
        variable = pickle.load(handle)
    return variable


def plotCompareSimulations(results, x_values, title, label_x, label_y, save_filename):

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    nb_simulations = len(results)

    for i in range(nb_simulations):
        name_simulation = results[i][0]
        y_values = results[i][1]
        plt.plot(x_values, y_values, label=name_simulation)
    plt.legend()
    plt.savefig(save_filename)
    plt.show()


#  Compare symmetry

# file_address = 'LOG/test_symmetry/uniform_mesh/lipfield/half_bar_300el/lipfield_'

# u_lip_half = readResults(file_address, "u")
# v_lip_half = readResults(file_address, "v")
# acel_lip_half = readResults(file_address, "acel")
# d_lip_half = readResults(file_address, "d")
# n_fragments_lip_half = readResults(file_address, "n_fragments")
# avg_frag_size_lip_half = readResults(file_address, "avg_frag_size")
# avg_stress_bar_lip_half = readResults(file_address, "avg_stress_bar")
# energies_lip_half = readResults(file_address, "energies")
# u_all_steps_lip_half = readResults(file_address, "u_all_steps")
# damage_all_steps_lip_half = readResults(file_address, "damage_all_steps")
# stress_all_steps_lip_half = readResults(file_address, "stress_all_steps")
# fraglen_all_steps_lipfield = readResults(file_address, "fraglen_all_steps")
# time_data = readResults(file_address, "time_data")
# time_simulation = time_data[0]
# n_steps = time_data[2]


# file_address = 'LOG/test_symmetry/uniform_mesh/lipfield/whole_bar_600el/lipfield_'
# u_lip_whole = readResults(file_address, "u")
# v_lip_whole = readResults(file_address, "v")
# acel_lip_whole =readResults(file_address, "acel")
# d_lip_whole = readResults(file_address, "d")
# n_fragments_lip_whole = readResults(file_address, "n_fragments")
# avg_frag_size_lip_whole = readResults(file_address, "avg_frag_size")
# avg_stress_bar_lip_whole = readResults(file_address, "avg_stress_bar")
# energies_lip_whole = readResults(file_address, "energies")
# u_all_steps_lip_whole = readResults(file_address, "u_all_steps")
# damage_all_steps_lip_whole = readResults(file_address, "damage_all_steps")
# # stress_all_steps_lip_whole = readResults(file_address, "stress_all_steps")

# time = np.linspace(0, time_simulation, n_steps)

# stress_bar = [['whole bar', avg_stress_bar_lip_whole],
#               ['half bar', avg_stress_bar_lip_half]]





# plotCompareSimulations(
#     stress_bar,
#     time,
#     title="Test use symmetry",
#     label_x="time (s)",
#     label_y="Average stress at the bar",
#     save_filename="LOG/test_symmetry_avg_stress.svg",
# )

# Compare symmetry
