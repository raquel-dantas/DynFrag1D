import DFMesh
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint


# l = 0.00075     # Regularization length
l = 0.5*DFMesh.hun     # Regularization length
# l = 2.2*10**-6     # Regularization length
w = 2.          # Weight quadrature

# Coordenates of integration points
x = [DFMesh.hun*i + DFMesh.hun*0.5 for i in range(DFMesh.n_el)]


Yc = [DFMesh.sigmac[el]**2 / (2.*DFMesh.E) for el in range(DFMesh.n_el)] 
lamb = [2. * Yc[el] * l / DFMesh.Gc for el in range(DFMesh.n_el)]  
def h(lamb,d): return (2.*d - d**2) / (1. - d + lamb * d**2) **2

def DissipatedEnergy(d):
    Edis = 0.0
    for el in range(len(DFMesh.materials)):
        if DFMesh.materials[el] == 0:
            Edis += 0.5 * Yc[el] * h(lamb[el],d[el]) * DFMesh.A

    return Edis

# def func(d): return w*sum([
#             (0.5*(1. - d[el])** 2 * 
#             DFMesh.E*strain**2 + 
#             Yc[el] * 
#             h(lamb[el], d[el])) * 
#             DFMesh.hun/2.
#             for el in range(DFMesh.n_el)
#     ])

# A = scipy.sparse.eye(DFMesh.n_el-1, DFMesh.n_el) - scipy.sparse.eye(DFMesh.n_el-1, DFMesh.n_el, 1)
# b = DFMesh.hun/l
# const = LinearConstraint(A, -b * np.ones(DFMesh.n_el-1), b * np.ones(DFMesh.n_el-1))

# dlip_opt = minimize(
#     fun=func,
#     x0=dn,
#     method='SLSQP',
#     bounds=zip(dn, [1.]*DFMesh.n_el),
#     tol=1e-9,
#     constraints=const,
# )
# dlip = dlip_opt.x