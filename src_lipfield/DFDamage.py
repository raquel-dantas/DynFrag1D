import DFMesh



# l = 2.21*10**-6      # Regularization length
l = 10*DFMesh.hun      # Regularization length
w = 2.                  # Weight quadrature

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
