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
usi = DFMesh.u0
vsi = DFMesh.v0
acelsi = DFMesh.acel0

Epot = np.zeros((DFMesh.n_steps))
Ekin = np.zeros((DFMesh.n_steps))
Edis = np.zeros((DFMesh.n_steps))
Erev = np.zeros((DFMesh.n_steps))
Econ = np.zeros((DFMesh.n_steps))
Wext = np.zeros((DFMesh.n_steps))
Eimp = np.zeros((DFMesh.n_steps))
work = 0.0
eimp = 0.0

uedge = np.zeros((DFMesh.n_steps))

# For the calculation of external work, it is needed to save the value of the displacement u at the previous time step 'up' for the boundary elements, and this is an input on the energy function.
up_bc_left = np.array([0,0])
up_bc_right = np.array([0,0])

els_step = DFMesh.n_el

stress_evl = np.zeros((2*len(DFMesh.materials),DFMesh.n_steps))
av_stress_bar = np.zeros((DFMesh.n_steps))

for n in range(DFMesh.n_steps):

    # Post process (stress, strain, energies)
    strain, stress, average_stress = DFPostprocess.PostProcess(u)
    stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
    av_stress_bar[n] = DFPostprocess.StressBar(stress, els_step)
    Epot[n], Ekin[n], Edis[n], Erev[n], Econ[n], Wext[n], Eimp[n] = DFPostprocess.Energy(up_bc_left,up_bc_right, n, u, v, acel,usi, vsi, acelsi, work, eimp)
    eimp =  Eimp[n]

    # Get K, M and F
    K, M, F = DFFem.GlobalSystem()
    # print(F)
    DFMesh.C = np.resize(DFMesh.C,K.shape)

    # Plots
    # DFPlot.PlotByDOF(u)
    # print(u)
    # DFPlot.PlotByElement(stress)

    DFPlot.PlotVTK('animation/chalk',n,u,stress)
    # u,v,acel returns a vector for u,v and acel at every dof at the n step
    p_next = DFMesh.p[n+1]
    u, v, acel = DFNewmark.Newmark_exp(n, K, M, DFMesh.C, u, v, acel, p_next, DFMesh.dt, DFMesh.gamma)

    usi, vsi, acelsi = DFNewmark.Newmark_exp(n, K, M, DFMesh.C, usi, vsi, acelsi, p_next, DFMesh.dt, DFMesh.gamma)

# DFPlot.PlotDispEdge(uedge)

# Variation of energy [Energy, time] return the difference between the energy value between the time t and t0 
varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEimp, varEtot = DFPostprocess.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, Eimp)
PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEimp, PEtot = DFPostprocess.Power(Epot, Ekin, Edis, Erev, Econ, Wext, Eimp)

# Plots
# DFPlot.PlotStressByTime(stress_evl)
# DFPlot.PlotAverageStressBar(av_stress_bar)
DFPlot.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, Eimp)
DFPlot.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEimp, varEtot)
DFPlot.PlotVarEnergy(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEimp, PEtot)
