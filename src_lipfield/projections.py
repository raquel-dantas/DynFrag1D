from scipy.optimize import minimize
import numpy as np
import DFPlot
import DFMesh

# Lenght of the bar (m)
L = 1.05
x0 = 0.0
xf = L

# Number of linear elements
n_el = 21
# Number of points
n_points = n_el+1
# Element length(hel)
hel = L/n_el
# Coordenates of integration points
x = [hel*i + hel*0.5 for i in range(n_el)]
# Coordenates of the points in the mesh
coord = [x0 + hel*i for i in range(n_points)]


# Damage prediction field
ddash = [0.0, 0.03333333333333333, 0.06666666666666667,   0.08333333333333333, 0.11666666666666665, 0.13333333333333333,    0.16666666666666666, 0.19, 0.23999999999999997,   0.433333333333334, 0.9, 0.433333333333334, 0.23999999999999997,    0.19, 0.16666666666666666, 0.13333333333333333,    0.11666666666666665, 0.08333333333333333, 0.06666666666666667,    0.03333333333333333, 0]

# Lip-field regularization lenght
l = 0.2

lower = [minimize(
    lambda y: ddash[np.searchsorted(coord, y[0])-1] + abs(x[el]-y[0])/l,
    x0=0.5*L,
    method='SLSQP',
    bounds=[(x0, xf)],
    tol=1e-5
).fun
    for el in range(n_el)
]

upper = [-minimize(
    lambda y: -ddash[np.searchsorted(coord, y[0])-1] + abs(x[el]-y[0])/l,
    x0=0.5*L,
    method='SLSQP',
    bounds=[(x0, xf)],
    tol=1e-5
).fun
    for el in range(n_el)
]

# DFPlot.Plot(x,[ddash,lower,upper],"x","damage","Predicted damage field")


# Damage next time step

# Fabricated stress_c gives:

d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02204006, 0.06896462, 0.14971528, 0.10357966, 0.32905376, 0.70610548, 0.06212765, 0.2577371, 0.15601001, 0.14921157, 0.11367984, 0.04990276,
0.0, 0.0, 0.0, 0.0]

epsilon = [0.00336, 0.00336, 0.00336, 0.00335999, 0.003359, 0.00329387, 0.00305922, 
0.00255058, 0.00186622, 0.00111102, 0.00062619, 0.00102213, 0.00173492, 0.00212252, 
0.00236141, 0.00259228, 0.00291374, 0.00322845, 0.00336015, 0.00335999, 0.00336]


Yc = DFMesh.stress_c**2/(2*DFMesh.E)
lip_constraint = 2.21*10**-6
lamb = 2*Yc*lip_constraint/DFMesh.Gc
def h(d): return (2*d-d**2)/(1-d+lamb*d**2)**2

# Numerical integration on func
# Weight quadrature
w = 2.
def func(d): return w*sum([
    (0.5*(1-d)**2*DFMesh.E*epsilon**2 + Yc*h(d))*DFMesh.DetJac(el)
    for el in range(DFMesh.n_el)
])


dnext = [
    minimize(
        fun=func,
        x0=d[i],
        method='SLSQP',
        bounds=[(d[i], 1)],
        tol=1e-5,
        constraints=(
            {'type': 'ineq', 'fun': lambda d:
                -(d[i] - d[i+1] - DFMesh.ElemLength(i)/l)},
            {'type': 'ineq', 'fun': lambda d:
                -(d[i] - d[i-1] - DFMesh.ElemLength(i)/l)}
        )
    ).x
    for i in range(1, DFMesh.n_el-1)
]

DFPlot.Plot(x,[ddash,lower,upper,dnext],"x","damage","Predicted damage field")