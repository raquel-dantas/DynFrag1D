from scipy.optimize import minimize
from matplotlib import pyplot as plt
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
# Predicted (before regularization) damage field on the next time-step
ddash = [0.0, 0.03333333333333333, 0.06666666666666667,   0.08333333333333333, 0.11666666666666665, 0.13333333333333333,    0.16666666666666666, 0.19, 0.23999999999999997,   0.433333333333334, 0.9, 0.433333333333334, 0.23999999999999997,    0.19, 0.16666666666666666, 0.13333333333333333,    0.11666666666666665, 0.08333333333333333, 0.06666666666666667,    0.03333333333333333, 0]

# Lip-field regularization lenght
l = 0.000001

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


# Damage next time step

# Fabricated stress_c gives:

d = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02204006, 0.06896462, 0.14971528, 0.10357966, 0.32905376, 0.70610548, 0.06212765, 0.2577371, 0.15601001, 0.14921157, 0.11367984, 0.04990276,
0.0, 0.0, 0.0])

epsilon = np.array([0.00336, 0.00336, 0.00336, 0.00335999, 0.003359, 0.00329387, 0.00305922, 
0.00255058, 0.00186622, 0.00111102, 0.00062619, 0.00102213, 0.00173492, 0.00212252, 
0.00236141, 0.00259228, 0.00291374, 0.00322845, 0.00336015, 0.00335999, 0.00336])

# v =  [  0.           5.           9.99999992  14.99999331  19.99903094
#   24.69527139  26.66771372  27.68443174  27.73113945  31.5214166
#   30.95955144  42.51067674  73.04842123  63.43573856  71.63618304
#   74.7608629   78.01139162  81.37245421  87.7256013   95.00048933
#  100.00000722 105.          79.73795409  71.43249504  42.53141791
#   75.2993981   77.00772053  33.35232426  41.93743217  81.95979983
#   37.86716403  86.56523663  32.1361159   90.68760069  92.61150386
#   25.11394139]
# u =  [0.00000000e+00 1.68000000e-07 3.36000000e-07 5.03999997e-07
#  6.71999635e-07 8.39949421e-07 1.00464293e-06 1.16317665e-06
#  1.30952261e-06 1.44779588e-06 1.53885150e-06 1.71209733e-06
#  2.06778078e-06 2.16734061e-06 2.35086937e-06 2.51150679e-06
#  2.67884880e-06 2.85151220e-06 3.02414746e-06 3.19200036e-06
#  3.36000000e-06 3.52800000e-06 2.01667407e-06 2.08059480e-06
#  1.68078797e-06 2.24474349e-06 2.39343617e-06 1.48330038e-06
#  1.35448468e-06 2.54923458e-06 1.18199340e-06 2.70582519e-06
#  1.01021571e-06 2.86272508e-06 3.02399276e-06 8.39949421e-07]




Yc = DFMesh.stress_c**2/(2*DFMesh.E)
# lip_constraint = 2.21*10**-6
lamb = 2*Yc*l/DFMesh.Gc
def h(d): return (2*d-d**2)/(1-d+lamb*d**2)**2

# Numerical integration on func
# Weight quadrature
w = 2.
def func(d): return w*sum([
    (0.5*(1-d[el])**2*DFMesh.E*epsilon[el]**2 + Yc*h(d[el]))*hel/2.
    for el in range(n_el)
])


const22a = [{'type': 'ineq', 'fun': lambda d:
         -(d[i] - d[i+1] - hel/l)} for i in range(0,n_el-2)]
const22b = [{'type': 'ineq', 'fun': lambda d:
         -(d[i] - d[i-1] - hel/l)} for i in range(1,n_el-1)]
const = tuple([*const22a,*const22b])



dlip = minimize(
    fun=func,
    x0=d,
    method='SLSQP',
    bounds=zip(d, [1]*n_el),
    tol=1e-5,
    constraints=const
).x


def PlotDamage(x, d, dlip):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')

    plt.title(str("Predicted damage field"))
    plt.xlabel(str("x"))
    plt.ylabel(str("D"))
    plt.plot(x, d, label='d')
    plt.plot(x, dlip, label='dlip')
    plt.legend()
    plt.show()


# DFPlot.Plot(x,[dnext,d],"x","damage","Predicted damage field")
PlotDamage(x,d,dlip)