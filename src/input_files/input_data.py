# Input model 1D dynamic fragmentation of an expanding ring


# Method
use_1d_cohesive_elements = True
use_akantu = False
use_lipfield = False


# Type of mesh
uniform_mesh = False
create_mesh = True
# if there is an mesh file for input set create_mesh = False


# Material
young_modulus = 275.0 * 10**9  # (Pa)
density = 2750.0  # (kg/m3)
fracture_energy = 100.0  # (N/m)
stress_limit = 300.0 * 10**6  # (Pa)
generate_limit_stress_variation = True
# if there is already a file with the random values for limit stress set generate_limit_stress_variation = False


# Geometry
bar_length = 50 * 10**-3  # (m)
x0 = -0.5 * bar_length  # Left extremitiy x coordinate / 0-initial
xf = 0.5 * bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 100  
area = bar_length / number_elements  # Cross sectional area (m2) (Equal to element size )


# Load
strain_rate = 10.0**4  # (s-1)


# Time
time_simulation = 3.0 * 10**-7  # Total time of simulation (s)
