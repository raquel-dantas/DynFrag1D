import DFMesh



# l = 2.21*10**-6      # Regularization length
l = 10*DFMesh.hun      # Regularization length
w = 2.                 # Weight quadrature

# Energy release rate
Yc = [DFMesh.sigmac[el]**2 / (2.*DFMesh.E) for el in range(DFMesh.n_el)] 
# Constant lambda
lamb = [2. * Yc[el] * l / DFMesh.Gc for el in range(DFMesh.n_el)]  
# Softening function
def h(lamb,d): return (2.*d - d**2) / (1. - d + lamb * d**2) **2
