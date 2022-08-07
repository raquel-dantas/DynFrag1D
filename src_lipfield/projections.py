import scipy
from scipy.optimize import minimize
import numpy as np
import DFPlot

# Lenght of the bar (m)
L = 1.05
x0 = 0
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
# DFPlot.Plot(x,ddash,"x","ddash","Predicted damage field")

# np.searchsorted(coord,y)
# Lip-field regularization lenght
l = 0.2

# Upper projection
func_upper = lambda y: ddash[np.searchsorted(coord,y)] + abs(x-y)/l
upper = np.zeros(n_el)
v = np.searchsorted(coord,0.724)-1
u = ddash[np.searchsorted(coord,0.724)]
print(v)
print(u)
# print(upper[0])
px = x[0]
func = lambda y:  ddash[np.searchsorted(coord,y)] + abs(px-y)/l
bnds = (x0,xf)
x0 = 0.0
upper = minimize(fun, x0,method='SLSQP', bounds=bnds,tol=1e-5)

# for el in range(n_el):
#     upper[el] = minimize(lambda y: ddash[np.searchsorted(coord,y)] + abs(x[el]-y)/l, x0=0,method='SLSQP',bounds=(x0,xf),tol=1e-5).fun



# upper = [minimize(
#     lambda y: ddash[np.searchsorted(coord,y)] + abs(x[el]-y)/l,
#     x0=0.5*L,
#     method='SLSQP', 
#     bounds=(x0,xf),
#     tol=1e-5).fun for el in range(n_el)]

print(upper)

