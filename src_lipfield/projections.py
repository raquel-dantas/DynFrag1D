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

DFPlot.Plot(x,[ddash,lower,upper],"x","damage","Predicted damage field")

