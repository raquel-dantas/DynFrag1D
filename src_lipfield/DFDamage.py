import DFMesh
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from itertools import groupby
from operator import itemgetter


# Regularization inputs
n_element_reg = 10              # Nb elem in the reg length
l = n_element_reg*DFMesh.hun    # Regularization length
w = 2.                          # Weight quadrature

# Energy release rate
Yc = [DFMesh.sigmac[el]**2 / (2.*DFMesh.E) for el in range(DFMesh.n_el)] 

# Constant lambda
lamb_lip = [2. * Yc[el] * l / DFMesh.Gc for el in range(DFMesh.n_el)]  
lamb_czm = [DFMesh.sigmac[el] * DFMesh.hun / (DFMesh.E * DFMesh.deltac[el]) for el in range(DFMesh.n_el)]  
# lamb = lamb_czm
lamb = lamb_lip

# Softening function
def h_lip(lamb,d): return (2.*d - d**2) / (1. - d + lamb * d**2)**2
def h_czm(lamb,d): return 1./(1.-lamb)*((1./((1.-lamb)*(1.-d)**2 + lamb))-1.)
# h = h_czm
h = h_lip



# Functional given by Moes, Le and Stershic (2022)
def getFunctional(strain): 
    """Returns a lambda function of damage to be use in the optimization problem for the whole domain.\n
    Arguments: \n
    strain -- strain computed using the displacement at the next time-step(n+1)."""
    
    functional_whole_domain = lambda damage: w*sum([
        (0.5*(1. - damage[el])**2 * 
        DFMesh.E*strain[el]**2 + 
        Yc[el] * 
        h(lamb[el], damage[el])) * 
        DFMesh.hun*0.5
        for el in range(DFMesh.n_el)
        ])
    
    return functional_whole_domain


def getFunctionalSubdomain(strain, region_optimization): 
    """Returns a lambda function of damage to be use in the optimization problem for a specific region in the domain (subdomain).\n
    Arguments: \n
    strain -- strain computed using the displacement at the next time-step(n+1).\n
    region -- an array with the element indexes of the subdomain to compute damage at next time-step(n+1)."""

    functional_subdomain = lambda damage: w*sum([
        (0.5*(1. - damage[region_optimization[i]])**2 * 
        DFMesh.E*strain[region_optimization[i]]**2 + 
        Yc[region_optimization[i]] * 
        h(lamb[region_optimization[i]], damage[i])) * 
        DFMesh.hun*0.5
        for i in range(len(region_optimization))
        ])
    
    return functional_subdomain
