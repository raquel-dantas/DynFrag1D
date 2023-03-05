# Input model 1D dynamic fragmentation of an expanding ring


# Method
use_1d_cohesive_elements = False
use_lipfield = True


# Type of mesh
uniform_mesh = True
create_mesh = True # if there is an mesh file for input set create_mesh = False


# Material
young_modulus = 380 * 10**9  # (Pa)
density = 3.9 * 10 **3  # (kg/m3)
fracture_energy = 83.13  # (N/m)
stress_limit = 100.0 * 10**6  # (Pa)
generate_limit_stress_variation = True # if there is already a file with the random values for limit stress set generate_limit_stress_variation = False


# Geometry
bar_length = 2 * 10** -3  # (m)
x0 = 0  # Left extremitiy x coordinate / 0-initial
xf = bar_length  # Rigth extremitiy x coordinate / f-final
number_elements = 500  
area = 3 * 10 ** -7  # Cross sectional area (m2) (Equal to element size )


# Load
strain_rate = 10.0**4  # (s-1)


# Time
time_simulation = 3.0 * 10**-7  # Total time of simulation (s)
