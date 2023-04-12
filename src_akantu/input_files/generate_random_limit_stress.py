import numpy as np
import pickle

# This file generate a random variation around stress limit value pickle file to be read in DFMesh.py if generate_limit_stress_variation == False

n_elements = 5

stress_limit = 300.0 * 10**6  # (Pa)

stress_critical = np.random.uniform(
    low=stress_limit - 1.0 * 10**6,
    high=stress_limit + 1.0 * 10**6,
    size=(n_elements),
)

with open('src/input_files/random_stress_test.pickle', 'wb') as handle:
    pickle.dump(stress_critical, handle, protocol=pickle.HIGHEST_PROTOCOL)

# for analyse symmmetry :
# stress_critical_half = stress_critical[300:]
# with open('src/input_files/random_stress_critical_300.pickle', 'wb') as handle:
#     pickle.dump(stress_critical_half, handle, protocol=pickle.HIGHEST_PROTOCOL)
