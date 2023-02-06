input_data = {

    # Material
    E : 275.0*10**9,     # Young's module (Pa)
    rho : 2750.0,        # Density (kg/m3)
    Gc : 100.0,          # Fracture energy (N/m)
    stress_critical : 300.0*10**6,   # Limit stress / critical stress (Pa)
    # Geometry
    A : 1*10**-3,        # Cross sectional area (m2)
    L : 50*10**-3,       # Lenght of the bar (m)
    x0 : -0.5*L,         # Left extremitiy x coordinate / 0-initial
    xf : 0.5*L,          # Rigth extremitiy x coordinate / f-final
    n_el : 4,         # Number of elements (n_el)
    hun : L/(n_el*0.5),  # Size of the elements (h) for a uniform mesh (un) 

    strain_rate : 10.0**4,
    alpha : (stress_critical**2 + 4.5 * strain_rate**(2/3) * E * Gc**(2/3) * rho**(1/3)) / (4.5 * Gc)



}