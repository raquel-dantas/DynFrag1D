from matplotlib.pyplot import connect
import DFMesh
import DFFem
import DFPostprocess
import DFNewmark
import DFInterface
import DFPlot
import DFFragmentation
import numpy as np



def Run_simulation(strain_rate):

    # Initiation of variables
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
    up_bc_left = np.array([0,0])
    up_bc_right = np.array([0,0])

    nfrag = np.zeros((DFMesh.n_steps))
    # fraglen = np.zeros((DFMesh.n_steps))
    avg_fraglen = np.zeros((DFMesh.n_steps))


    for n in range(DFMesh.n_steps):

        # Plots at each time step
        # DFPlot.PlotByDOF(v)
        # DFPlot.PlotByElement(stress)
        # DFPlot.PlotByInterface(D)


        # Post process (stress, strain, energies)
        strain, stress, average_stress = DFPostprocess.PostProcess(u)

        # D returns a vector contained damage parameter for cohesive elements
        D = [DFInterface.DamageParameter(el) for el in range(len(DFMesh.materials))]
        # nfrag retuns a vector contained the number of fragments 
        nfrag[n] = DFFragmentation.NumberFragments(D)
        fraglen, avg_fraglen[n] = DFFragmentation.SizeFragments(D)

        stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
        av_stress_bar[n] = DFPostprocess.StressBar(stress, els_step)
        Epot[n], Ekin[n], Edis[n], Erev[n], Econ[n], Wext[n] = DFPostprocess.Energy(up_bc_left,up_bc_right, u, v, work)
        work =  Wext[n]

        # DFPlot.PlotVTK('animation/test',n,u,stress)
        # Get K, M and F
        K, M, F = DFFem.GlobalSystem()
        DFMesh.C = np.resize(DFMesh.C,K.shape)


        # up_bc is the previous displacement vector for the local dofs in the boundary elements (left and right)
        up_bc_left = np.array([0,0])
        up_bc_right = np.array([0,0])
        for bc in range(len(DFMesh.materials)):
            if DFMesh.materials[bc] == 4 or DFMesh.materials[bc] == 5:
                if DFMesh.materials[bc] == 4:
                    elbc = 0
                    up_bc_left = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])
                else:
                    elbc = DFMesh.n_el - 1
                    up_bc_right = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])

        # u,v,acel returns a vector for u,v and acel at every dof at the n step
        u, v, acel = DFNewmark.Newmark_exp(K, M, DFMesh.C, u, v, acel, F, DFMesh.dt)

        # Check limit stress for possible insertion of interface elements
        for el in range(DFMesh.n_el-1):
            if average_stress[el] > DFMesh.stress_c:
                # Fracture happens: creat new interface element
                u, v, acel = DFInterface.InsertInterface(el, el+1, u, v, acel)
                els_step = els_step + 1

        



    # Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
    varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPostprocess.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)

    # Power [Energy, time] returns the energy difference between consecutive time steps
    PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot = DFPostprocess.Power(Epot, Ekin, Edis, Erev, Econ, Wext)

    # Plots for the whole simulation
    # DFPlot.PlotStressByTime(stress_evl)
    # DFPlot.PlotAverageStressBar(av_stress_bar)
    # DFPlot.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)
    # DFPlot.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot)
    # DFPlot.PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot)
    DFPlot.PlotNumberFragments(nfrag)
    DFPlot.PlotAvgFragmentSize(avg_fraglen)



if __name__ == '__main__':
    Run_simulation(DFMesh.strain_rate)