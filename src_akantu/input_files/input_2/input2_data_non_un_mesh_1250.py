# Input model for src_akantu

# Type of mesh
uniform_mesh = False
# if there is an mesh file for input set create_mesh = False
create_mesh = False 
if create_mesh == False:
    mesh_file_name = 'src_akantu/input_files/mesh_files/mesh_non_uniform_1250.pickle'

# Material
young_modulus = 275.0*10**9   # (Pa)
density = 2750.0 # (kg/m3)
fracture_energy = 100.  # (N/m)
stress_limit = 300.0 * 10**6  # (Pa)

# if there is already a file with the random values for limit stress set generate_limit_stress_variation = False
generate_limit_stress_variation = False 
if generate_limit_stress_variation == False:
    stress_limit_file_name = 'src_akantu/input_files/random_stress_files/random_stress_critical_1250.pickle'


# Geometry
bar_length = 50 * 10** -3  # (m)
x0 = -0.5 * bar_length  # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 1250 * 2 # Total number of triangular elements
area = 1.  # Cross sectional area (m2) (Equal to element size )


# Load
strain_rate = 10.0**4  # (s-1)


# Time
time_simulation = 4.5 * 10**-7   # Total time of simulation (s)

# if there is previous data to continue the simulation from a previous simulation set 
# continue_simulation_from_step = True and give the time to start the simulation
initial_step = 0
continue_simulation_from_step = False


half_bar = False
# if use symmetry we have to add the bc proper

filepath = "src_akantu/LOG/mesh_study/czm_non_uniform_mesh/1250el/"

