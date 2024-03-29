# Input model 1D dynamic fragmentation of an expanding ring
import subprocess

# Method
use_cohesive_elements = False
use_lipfield = True


# Type of mesh
uniform_mesh = True
# if there is an mesh file for input set create_mesh = False
create_mesh = True
if create_mesh == False:
    mesh_file_name = "input_files/mesh_non_uniform.pickle"


# Material
young_modulus = 275e9  # (Pa)
density = 2.75e3  # (kg/m3)
fracture_energy = 100.0  # (N/m)
stress_limit = 300e6  # (Pa)

# if there is already a file with the random values for limit stress set generate_limit_stress_variation = False
generate_limit_stress_variation = True
if generate_limit_stress_variation == False:
    stress_limite_file_name = "input_files/random_stress_critical.pickle"


# Geometry
bar_length = 50e-3  # (m)
x0 = -0.5 * bar_length  # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 5000
area = 1.0  # Cross sectional area
regularization_length = bar_length/2500 # Considering at least 5 elements in a 2500 elements mesh


# Load
strain_rate = 1e5  # (s-1)


# Time
time_simulation = 2.0e-7  # Total time of simulation (s)


# if there is previous data to continue the simulation set continue_simulation_from_step = True and give the time to start the simulation and the files path
initial_step = 0
continue_simulation_from_step = False
if continue_simulation_from_step == True:
    previous_simulation = "output/akantu_step__.pickle"


# if use symmetry we have to add the bc proper
half_bar = False

subprocess.Popen("mkdir output", shell=True)
filepath_save_results = "output/"
