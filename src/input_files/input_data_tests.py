# Input model 1D dynamic fragmentation of an expanding ring


# Method
use_cohesive_elements = False
use_lipfield = True


# Type of mesh
uniform_mesh = True
# if there is an mesh file for input set create_mesh = False
create_mesh = True 
if create_mesh == False:
    mesh_file_name = 'filename.pickle'


# Material
young_modulus = 275.0*10**9   # (Pa)
density = 2750.0 # (kg/m3)
fracture_energy = 300.  # (N/m)
stress_limit = 300.0 * 10**6  # (Pa)

# if there is already a file with the random values for limit stress set generate_limit_stress_variation = False
generate_limit_stress_variation = True 
if generate_limit_stress_variation == False:
    stress_limite_file_name = 'filename.pickle'


# Geometry
bar_length = 50 * 10** -3  # (m)
x0 = -0.5 * bar_length  # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 1250
area = 1.  # Cross sectional area 


# Load
strain_rate = 10.0**4  # (s-1)


# Time
time_simulation = 3.5 * 10**-7  # Total time of simulation (s)


# if there is previous data to continue the simulation set continue_simulation_from_step = True and give the time to start the simulation and the files path
initial_step = 498
continue_simulation_from_step = True
if continue_simulation_from_step == True:
    previous_simulation = "LOG/lipfield_test/lipfield_step_498_.pickle"


half_bar = False
# if use symmetry we have to add the bc proper

filepath_save_results = "LOG/lipfield_test/"