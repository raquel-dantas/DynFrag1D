import numpy as np
import pickle

# This file generate a non-uniform mesh and save the nodes coordinates in a pickle file to be read in DFMesh.py if create_mesh == False

n_elements = 250
n_points = n_elements + 1

bar_length = 5 * 10**-3  # (m)
x0 = -0.5 * bar_length 
xf = 0.5 * bar_length  

h_uniform = bar_length / n_elements

uniform_coord = np.linspace(x0, xf, n_points)
node_coord = uniform_coord
node_coord = np.array(
    [x + np.random.uniform(low=-0.4, high=0.4) * h_uniform for x in uniform_coord]
)
node_coord[0] = x0
node_coord[n_elements] = xf


with open('src/input_files/non_uniform_mesh.pickle', 'wb') as handle:
    pickle.dump(node_coord, handle, protocol=pickle.HIGHEST_PROTOCOL)


