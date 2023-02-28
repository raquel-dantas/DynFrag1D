# Input model for src_lipfield

# Material
E = 275.0 * 10**9               # Young's module (Pa)
rho = 2750.0                    # Density (kg/m3)
Gc = 100.0                      # Fracture energy (N/m)
stress_critical = 300.0 * 10**6 # Limit stress / critical stress (Pa)

# Geometry
L = 5 * 10**-3  # Lenght of the bar (m)
x0 = -0.5 * L   # Left extremitiy x coordinate / 0-initial
xf = 0.5 * L    # Rigth extremitiy x coordinate / f-final
n_el = 500      # Number of linear elements (n_el)
hun = L / n_el  # Size of the elements (h) for a uniform mesh (un)
# A = hun         # Cross sectional area (m2)
A = 1*10**-3         # Cross sectional area (m2)

# Load
strain_rate = 10.0**4  # (s-1)

# Time
dt_crit = hun / ((E / rho) ** 0.5)  # Critical time step
dt = dt_crit * 0.1                  # Adopted time step (s)
time_simulation = 3.0 * 10**-7      # Total time of simulation (s)
n_steps = int(time_simulation / dt) # Number of time-steps
