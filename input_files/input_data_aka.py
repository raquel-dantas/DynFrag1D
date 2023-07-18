# Input model for src_akantu
import subprocess

# Type of mesh
uniform_mesh = False
# if there is an mesh file for input set create_mesh = False
create_mesh = False
if create_mesh == False:
    mesh_file_name = "input_files/mesh_files/mesh_non_uniform_7500.pickle"

# Material
young_modulus = 275e9  # (Pa)
density = 2.75e3  # (kg/m3)
fracture_energy = 100.0  # (N/m)
stress_limit = 300e6  # (Pa)

# if there is already a file with the random values for limit stress set generate_limit_stress_variation = False
generate_limit_stress_variation = False
if generate_limit_stress_variation == False:
    stress_limit_file_name = "input_files/random_stress_files/random_stress_critical_7500.pickle"


# Geometry
bar_length = 50e-3  # (m)
x0 = -0.5 * bar_length  # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 7500 * 2 # Total number of triangular elements
area = 1.0  # Cross sectional area (m2) (Equal to element size )


# Load
strain_rate = 1e3  # (s-1)


# Time
time_simulation = 3.5e-6   # Total time of simulation (s)

# if there is previous data to continue the simulation set continue_simulation_from_step = True and give the time to start the simulation and the files path
initial_step = 0
continue_simulation_from_step = False
if continue_simulation_from_step == True:
    previous_simulation = "output_mesh_study/10to3/10to3_akantu_uniform_7500/akantu_step_0_.pickle"


# if use symmetry we have to add the bc proper
half_bar = False

filepath_save_results = "output_mesh_study/10to3/10to3_akantu_non_uniform_7500/"

