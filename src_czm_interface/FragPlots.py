import numpy as np
import DFPlot
import DFFragmentation


# # Analytical estimations of fragment size
# points = 10
# strainrate = np.array([[10**x, 5*10**x] for x in range(2,points)])
# strainrate = strainrate.flatten()
# print(strainrate)

# norm_strainrate = np.array([epsilon/(((DFMesh.E / DFMesh.rho)**0.5 * DFMesh.stress_critical**3)/(DFMesh.E**2 * DFMesh.Gc)) for epsilon in strainrate])
# print(norm_strainrate)

# s_grady, norm_grady = DFFragmentation.GradyFragSize(strainrate)

# s_gc = DFFragmentation.GlenChudnoviskFragSize(strainrate)

# s_zmr = DFFragmentation.ZhouMolinariRameshFragSize(strainrate)


# DFPlot.PlotLogAnalyticals(norm_grady, s_gc, s_zmr, norm_strainrate)

# Convergence study

# meshes_un = np.array([500,1000,1500,2000,2500,3000,3500,4000])
# meshes_non_un = np.array([500,1000,1500,2000,2500,3000,3500,4000])

# # energies_un = np.array([19900.472340806784, 49840.19889092332, 66509.65502058562, 69595.74986314122, 69112.69004289804,  68360.64077513237, 68528.23425253623, 68814.43084842109, 67793.0109585675 ])
# nfrags_un = np.array([]) #nfrag src_akantu 
# nfrags_non_un = np.array([]) #nfrag src_czm_interface

# # DFPlot.PlotConvergenceEnergy(energies_un, energies_nun, meshes)
# # DFPlot.PlotLogConvergenceEnergy(energies_un, energies_nun, meshes)
# DFPlot.PlotConvergenceNumfrag(nfrags_un, nfrags_non_un, nfrags_non_un, meshes_un, meshes_non_un)



# Candidacy

meshes_un = np.array([250,500,750,1000,1250,1500,2000,2500,3000,3500,4000])
meshes_non_un = np.array([250,500,750,1000,1250,1500,2000,2500,3000,3500,4000])

nfrags_un = np.array([]) #nfrag src_akantu 
nfrags_non_un = np.array([]) #nfrag src_czm_interface

# DFPlot.PlotConvergenceEnergy(energies_un, energies_nun, meshes)
# DFPlot.PlotLogConvergenceEnergy(energies_un, energies_nun, meshes)
DFPlot.PlotConvergenceNumfrag(nfrags_un, nfrags_non_un, nfrags_non_un, meshes_un, meshes_non_un)