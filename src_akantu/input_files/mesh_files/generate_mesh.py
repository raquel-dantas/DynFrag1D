import numpy as np
import pickle

# This file generate a non-uniform mesh and save the nodes coordinates in a pickle file to be read in DFMesh.py if create_mesh == False

n_elements = 5
n_points = n_elements + 1

bar_length = 50 * 10**-3  # (m)
x0 = -0.5 * bar_length 
xf = 0.5 * bar_length  

use_uniform_mesh = False

h_uniform = bar_length / n_elements

uniform_node_coord = np.linspace(x0, xf, n_points)
non_uniform_node_coord = uniform_node_coord
non_uniform_node_coord = np.array(
    [x + np.random.uniform(low=-0.4, high=0.4) * h_uniform for x in uniform_node_coord]
)
non_uniform_node_coord[0] = x0
non_uniform_node_coord[n_elements] = xf

if use_uniform_mesh == True:
    with open('src/input_files/mesh_uniform_625.pickle', 'wb') as handle:
        pickle.dump(uniform_node_coord, handle, protocol=pickle.HIGHEST_PROTOCOL)

else:
    with open('src/input_files/mesh_non_uniform_test.pickle', 'wb') as handle:
        pickle.dump(non_uniform_node_coord, handle, protocol=pickle.HIGHEST_PROTOCOL)


