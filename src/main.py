from matplotlib.pyplot import connect
import DFMesh
import DFFem
import DFPostprocess
import DFNewmark
import DFInterface
import DFPlot
import numpy as np


# Global algorithm

u = DFMesh.u0
v = DFMesh.v0
acel = DFMesh.acel0

Epot = np.zeros((DFMesh.n_steps))
Ekin = np.zeros((DFMesh.n_steps))
Edis = np.zeros((DFMesh.n_steps))
Erev = np.zeros((DFMesh.n_steps))
Econ = np.zeros((DFMesh.n_steps))
Wext = np.zeros((DFMesh.n_steps))
work = 0.0


els_step = DFMesh.n_el

stress_evl = np.zeros((2*len(DFMesh.materials),DFMesh.n_steps))
av_stress_bar = np.zeros((DFMesh.n_steps))

for n in range(DFMesh.n_steps):

    up = u

    # Post process (stress, strain, energies)
    strain, stress, average_stress = DFPostprocess.PostProcess(u)
    stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
    av_stress_bar[n] = DFPostprocess.StressBar(stress, els_step)

    Epot[n], Ekin[n], Edis[n], Erev[n], Econ[n], Wext[n] = DFPostprocess.Energy(up, u, v, n, work)
    work =  Wext[n]

    # Get K, M and F
    K, M, F = DFFem.GlobalSystem()
    DFMesh.C = np.resize(DFMesh.C,K.shape)

    # Plots
    # DFPlot.PlotByDOF(acel)
    # DFPlot.PlotByElement(stress)

    # u,v,acel returns a vector for u,v and acel at every dof at the n step
    u, v, acel = DFNewmark.Newmark_exp(K, M, DFMesh.C, u, v, acel, F, DFMesh.dt, DFMesh.gamma)

    # DFPostprocess.VerifyFractureState(stress, u, v, up, vp)
    
    # Check limit stress
    for el in range(DFMesh.n_el-1):
        if average_stress[el] > DFMesh.stress_c:
            # Fracture happens: creat new interface element
            u, v, acel = DFInterface.InsertInterface(el, el+1, u, v, acel)
            els_step = els_step + 1
    

    # D returns a vector contained damage parameter for cohesive elements
    D = [DFInterface.DamageParameter(el) for el in range(len(DFMesh.materials))]
    # DFPlot.PlotByInterface(D)

# Variation of energy
varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPostprocess.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)

# Plots
# DFPlot.PlotStressByTime(stress_evl)
# DFPlot.PlotAverageStressBar(av_stress_bar)
DFPlot.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)
DFPlot.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot)
# DFPlot.PlotVarEnergy(varEtot)