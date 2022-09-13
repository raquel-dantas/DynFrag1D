import numpy as np
import DFPlot2d



# Convergence study

meshes = np.array([200,500,1000,1250,1500,2000,2500, 3000, 4000])

energies_un = np.array([19900.472340806784, 49840.19889092332, 66509.65502058562, 69595.74986314122, 69112.69004289804,  68360.64077513237, 68528.23425253623, 68814.43084842109, 67793.0109585675 ])
nfrags_un = np.array([200.0, 498.0, 431.0, 406.0, 384.0, 380.0, 375.0, 387.0, 375.0])

energies_nun = np.array([19899.13354102981, 42803.23452480667, 54405.863345645725, 56881.30579003729, 59454.90498136255, 62704.33683292938,63658.51791004801, 64110.87886908749, 65689.65727670692])
nfrags_nun = np.array([200.0, 339.0, 338.0, 351.0, 359.0, 378.0, 361.0, 359.0, 371.0])

DFPlot2d.PlotConvergenceEnergy(energies_un, energies_nun, meshes)
DFPlot2d.PlotLogConvergenceEnergy(energies_un, energies_nun, meshes)
DFPlot2d.PlotConvergenceNumfrag(nfrags_un, nfrags_nun, meshes)