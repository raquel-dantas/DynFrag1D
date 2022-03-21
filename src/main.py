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
Etot = np.zeros((DFMesh.n_steps))

els_step = DFMesh.n_el

stress_evl = np.zeros((2*len(DFMesh.materials),DFMesh.n_steps))


for n in range(DFMesh.n_steps):

    # Post process (stress, strain, energies)
    strain, stress, average_stress = DFPostprocess.PostProcess(u)
    stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
    Ekin[n], Epot[n], Edis[n], Erev[n], Etot[n] = DFPostprocess.Energy(u, v)

    # Get K, M and F
    K, M, F = DFFem.GlobalSystem()
    DFMesh.C = np.resize(DFMesh.C,K.shape)

    # Plots
    # DFPlot.PlotByDOF(v)
    # DFPlot.PlotByElement(stress)

    # u,v,acel returns a vector for u,v and acel at every dof at the n step
    u, v, acel = DFNewmark.Newmark_exp(K, M, DFMesh.C, u, v, acel, F, DFMesh.dt, DFMesh.gamma)

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
varEkin, varEpot, varEdis, varErev, varEtot = DFPostprocess.VarEnergy(Ekin, Epot, Edis, Erev, Etot)

# Plots
DFPlot.PlotStressByTime(DFMesh.n_steps, stress_evl)
DFPlot.PlotEnergy(DFMesh.n_steps, Epot, Ekin, Edis, Erev, Etot)
DFPlot.PlotVarEnergy(DFMesh.n_steps, varEpot, varEkin, varEdis, varErev, varEtot)
