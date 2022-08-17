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


# Fabricated critical stress field gives:
dn = [0.        , 0.        , 0.        , 0.        , 0.        , 0.        ,                  
      0.02204006, 0.06896462, 0.14971528, 0.10357966, 0.32905376, 0.70610548,
      0.06212765, 0.2577371 , 0.15601001, 0.14921157, 0.11367984, 0.04990276,
      0.        , 0.        , 0.        ]

dnext = [0.       , 0.        , 0.        , 0.        , 0.        , 0.00141144,
        0.03161133, 0.08465093, 0.16895514, 0.1049445 , 0.33939497, 0.74111494,
        0.06212765, 0.27381818, 0.16386678, 0.16067532, 0.12819722, 0.06647385,
        0.01825111, 0.        , 0.        ]


epsilon = [0.00336   , 0.00336   , 0.00336   , 0.00335999, 0.003359  , 0.00329387, 0.00305922,
           0.00255058, 0.00186622, 0.00111102, 0.00062619, 0.00102213,0.00173492,  0.00212252, 0.00236141, 0.00259228, 0.00291374, 0.00322845, 0.00336015, 0.00335999, 0.00336  ]

epsilon_next = [0.0034    , 0.0034    , 0.0034    , 0.00339998, 0.00339469, 0.00330176,
                0.00301833, 0.00246591, 0.00178436, 0.00109707, 0.00062842, 0.00096637,
                0.0016679 , 0.00209195, 0.00234327, 0.00256056, 0.00287197, 0.0031831 ,
                0.00335775, 0.00339999, 0.0034     ]


# Damage prediction field

# l = 0.000001 # Lip-field regularization lenght
l = 2.21*10**-6 # Lip-field regularization lenght
# l = 0.002 # Lip-field regularization lenght
Yc = DFMesh.stress_c**2/(2*DFMesh.E)
# lip_constraint = 2.21*10**-6

lamb = 2*Yc*l/DFMesh.Gc
def h(d): return (2*d-d**2)/(1-d+lamb*d**2)**2

# Numerical integration on func
# Weight quadrature
w = 2.
def func(d): return w*sum([
    (0.5*(1-d[el])**2*DFMesh.E*epsilon_next[el]**2 + Yc*h(d[el]))*hel/2.
    for el in range(n_el)
])


 
# Unregularized damage (ddash) for time step n+1 (Optimization 23)
ddash = minimize(
    fun=func,
    x0=dn,
    method='SLSQP',
    bounds=zip(dn, [1]*n_el),
    tol=1e-5,
).x

funproj = [lambda y: ddash[np.searchsorted(coord, y[0])-1] + abs(x[el]-y[0])/l for el in range(n_el)]

# Projections (Optimization problem 24a and 24b)
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




# Constraints for regularized damage for time step n+1 (dlip)
const22a = [{'type': 'ineq', 'fun': lambda d:
         -(d[i] - d[i+1] - hel/l)} for i in range(0,n_el-2)]
const22b = [{'type': 'ineq', 'fun': lambda d:
         -(d[i] - d[i-1] - hel/l)} for i in range(1,n_el-1)]
const = tuple([*const22a,*const22b])


# Regularized damage for time step n+1 (dlip) (Optimization 21b)
dlip = minimize(
    fun=func,
    x0=dn,
    method='SLSQP',
    bounds=zip(dn, [1]*n_el),
    tol=1e-5,
    constraints=const
).x







def PlotDamage(x, dn, dnext, ddash, upper, lower, dlip):
    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')

    plt.title(str("Damage field"))
    plt.xlabel(str("x"))
    plt.ylabel(str("D"))
    plt.plot(x, dn, label='dn')
    # plt.plot(x, dnext, label='dnext')
    plt.plot(x, ddash, label='ddash')
    # plt.plot(x, upper, label='upper')
    # plt.plot(x, lower, label='lower')
    plt.plot(x, dlip, label='dlip')
    plt.legend()
    plt.show()


PlotDamage(x, dn, dnext, ddash, upper, lower, dlip)