# Input model for src_czm_interface

# Material
E = 380.0 * 10**9              # Young's module (Pa)
rho = 3.9 * 10**3              # Density (kg/m3)
Gc = 83.13                     # Fracture energy (N/m)
stress_critical = 1.0 * 10**9  # Limit stress / critical stress (Pa)

# Geometry
A = 2.0 * 10**-7  # Cross sectional area (m2)
L = 2.0 * 10**-3  # Lenght of the bar (m)
x0 = -0.5 * L     # Left extremitiy x coordinate / 0-initial
xf = 0.5 * L      # Rigth extremitiy x coordinate / f-final
n_el = 100        # Number of linear elements (n_el)
hun = L / n_el    # Size of the elemenets (h) for a uniform mesh (un)

# Load
strain_rate = 10.0**5  # (s-1)

# Time
dt_crit = hun / ((E / rho) ** 0.5)   # Critical time step
dt = dt_crit * 0.9                   # Adopted time step (s)
time_simulation = 3.0 * 10**-7       # Total time of simulation (s)
n_steps = int(time_simulation / dt)  # Number of time-steps
