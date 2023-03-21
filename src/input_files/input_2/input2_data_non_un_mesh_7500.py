import subprocess
# Input model 1D dynamic fragmentation of an expanding ring


# Method
use_cohesive_elements = False
use_lipfield = True


# Type of mesh
uniform_mesh = False
# if there is an mesh file for input set create_mesh = False
create_mesh = False  
if create_mesh == False:
    mesh_file_name = "src/input_files/mesh_files/mesh_non_uniform_7500.pickle"


# Material
young_modulus = 275.0 * 10**9  # (Pa)
density = 2750.0  # (kg/m3)
fracture_energy = 100.0  # (N/m)
stress_limit = 300.0 * 10**6  # (Pa)

# if there is already a file with the random values for limit stress set generate_limit_stress_variation = False
generate_limit_stress_variation = False  
if generate_limit_stress_variation == False:
    stress_limite_file_name = (
        "src/input_files/random_stress_files/random_stress_critical_7500.pickle"
    )


# Geometry
bar_length = 50.0 * 10**-3  # (m)
x0 = -0.5 * bar_length  # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 7500
area = 1.0  # Cross sectional area (m2) (Equal to element size )


# Load
strain_rate = 10.0**4  # (s-1)


# Time
time_simulation = 2.5 * 10**-7  # Total time of simulation (s)

# if there is previous data to continue the simulation set continue_simulation_from_step = True and give the time to start the simulation and the files path
initial_step = 11740
continue_simulation_from_step = True
if continue_simulation_from_step == True:
    previous_simulation = "output/lipfield_step_11740_.pickle"

half_bar = False
# if use symmetry we have to add the bc properly

# filepath_save_results = "LOG/mesh_study/lipfield_non_uniform_mesh/7500el/"

subprocess.Popen("mkdir output", shell=True)
filepath_save_results = "output/"
