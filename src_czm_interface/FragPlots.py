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

meshes = np.array([1000,5000,7000,10000,14000])
meshes2 = np.array([500,1000,2000,2500,3000])

# energies_un = np.array([19900.472340806784, 49840.19889092332, 66509.65502058562, 69595.74986314122, 69112.69004289804,  68360.64077513237, 68528.23425253623, 68814.43084842109, 67793.0109585675 ])
nfrags_un = np.array([80,103,96,96,96]) #nfrag src_akantu 
nfrags_2 = np.array([216,165,146,142,136]) #nfrag src_czm_interface

# energies_nun = np.array([19899.13354102981, 42803.23452480667, 54405.863345645725, 56881.30579003729, 59454.90498136255, 62704.33683292938,63658.51791004801, 64110.87886908749, 65689.65727670692])
nfrags_nun = np.array([200.0, 339.0, 338.0, 351.0, 359.0, 378.0, 361.0, 359.0, 371.0])

# DFPlot.PlotConvergenceEnergy(energies_un, energies_nun, meshes)
# DFPlot.PlotLogConvergenceEnergy(energies_un, energies_nun, meshes)
DFPlot.PlotConvergenceNumfrag(nfrags_un, nfrags_2, nfrags_nun, meshes, meshes2)
