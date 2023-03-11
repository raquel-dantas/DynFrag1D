# Input model 1D dynamic fragmentation of an expanding ring


# Method
use_cohesive_elements = False
use_lipfield = True


# Type of mesh
uniform_mesh = True
create_mesh = True # if there is an mesh file for input set create_mesh = False
if create_mesh == False:
    mesh_file_name = 'src/input_files/filename.pickle'


# Material
young_modulus = 275.0*10**9   # (Pa)
density = 2750.0 # (kg/m3)
fracture_energy = 100.  # (N/m)
stress_limit = 300.0 * 10**6  # (Pa)

generate_limit_stress_variation = False # if there is already a file with the random values for limit stress set generate_limit_stress_variation = False
if generate_limit_stress_variation == False:
    stress_limite_file_name = 'src/input_files/random_stress_critical_5000.pickle'


# Geometry
bar_length = 50. * 10** -3  # (m)
x0 = - 0.5 * bar_length # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 5000
area = 1.  # Cross sectional area (m2) (Equal to element size )


# Load
strain_rate = 10.0**4  # (s-1)


# Time
time_simulation = 3.5 * 10**-7  # Total time of simulation (s)


initial_step = 0
continue_simulation_from_step = False
# if there is previous data to continue the simulation from a previous simulation set 
# continue_simulation_from_step = True and give the time to start the simulation

half_bar = False
# if use symmetry we have to add the bc properly