import numpy as np
import progressbar
import DFMesh
import DFFem
import DFPostprocess
import DFNewmark
import DFPlot
import DFFragmentation
import pickle
import time

def Run_simulation(strain_rate):

    bar = progressbar.ProgressBar(maxval=50, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    # Initiation of variables
    u = DFMesh.u0
    v = DFMesh.v0
    acel = DFMesh.acel0
    d = DFMesh.d0
    Epot = np.zeros((DFMesh.n_steps))
    Ekin = np.zeros((DFMesh.n_steps))
    Edis = np.zeros((DFMesh.n_steps))
    Wext = np.zeros((DFMesh.n_steps))
    work = 0.0
    Edis_prev = 0.0
    els_step = DFMesh.n_el
    stress_evl = np.zeros((2*len(DFMesh.materials),DFMesh.n_steps))
    av_stress_bar = np.zeros((DFMesh.n_steps))
    up_bc_left = np.array([0,0])
    up_bc_right = np.array([0,0])
    dp_bc_left = 0.0
    dp_bc_right = 0.0

    nfrag = np.zeros((DFMesh.n_steps))
    avg_fraglen = np.zeros((DFMesh.n_steps))


    for n in range(DFMesh.n_steps):        

        progress = int(bar.maxval*float(n/DFMesh.n_steps))
        bar.update(progress)
        
        # DFPlot.PlotByDOF(v)

        # Post process (stress, strain, energies)
        strain, stress, average_stress = DFPostprocess.PostProcess(u, d)

        # nfrag retuns a vector contained the number of fragments 
        nfrag[n] = DFFragmentation.NumberFragments(d)
        fraglen, avg_fraglen[n] = DFFragmentation.SizeFragments(d)

        stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
        av_stress_bar[n] = DFPostprocess.StressBar(stress, els_step)
        Epot[n], Ekin[n], Edis[n], Wext[n] = DFPostprocess.Energy(up_bc_left,
        up_bc_right, u, v, stress, work, d, dp_bc_left, dp_bc_right, Edis_prev)
        work =  Wext[n]
        Edis_prev = Edis[n]

        # DFPlot.PlotVTK('LOG/animation_lipfield/test',n,u,stress)

        # Get K, M and F
        M, F = DFFem.GlobalSystem()


        # up_bc is the previous displacement vector for the local dofs in the boundary elements (left and right)
        up_bc_left = np.array([0,0])
        up_bc_right = np.array([0,0])
        for bc in range(len(DFMesh.materials)):
            if DFMesh.materials[bc] == 4 or DFMesh.materials[bc] == 5:
                if DFMesh.materials[bc] == 4:
                    elbc = 0
                    up_bc_left = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])
                    dp_bc_left = d[elbc]
                else:
                    elbc = DFMesh.n_el - 1
                    up_bc_right = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])
                    dp_bc_right = d[elbc]
        # u,v,acel returns a vector for u,v and acel at every dof at the n step
        u, v, acel, d = DFNewmark.Newmark_exp(M, u, v, acel, d, F, DFMesh.dt)
    
    bar.finish()
    print('\n')


    # Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
    varEkin, varEpot, varEdis, varWext, varEtot = DFPostprocess.VarEnergy(Epot, Ekin, Edis, Wext)
    # Power [Energy, time] returns the energy difference between consecutive time steps
    PEkin, PEpot, PEdis, PWext, PEtot = DFPostprocess.Power(Epot, Ekin, Edis, Wext)


    # Plots for the whole simulation
    DFPlot.PlotAverageStressBar(av_stress_bar)
    # DFPlot.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)
    DFPlot.PlotVarEnergy(varEpot, varEkin, varEdis, varWext, varEtot)
    DFPlot.PlotPower(PEpot, PEkin, PEdis, PWext, PEtot)
    DFPlot.PlotNumberFragments(nfrag)
    # DFPlot.PlotAvgFragmentSize(avg_fraglen)
    # DFPlot.PlotFragmentSizeHistogram(fraglen)    
    
    # Number of fragments
    # with open('LOG/src_czm_interface_number_fragments.pickle', 'wb') as handle:
    #     pickle.dump(nfrag, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # # Average fragment size
    # with open('LOG/src_czm_interface_average_fragment_size.pickle', 'wb') as handle:
    #     pickle.dump(avg_fraglen, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # # Average stress for the bar 
    # with open('LOG/src_czm_interface_avg_stress.pickle', 'wb') as handle:
    #     pickle.dump(av_stress_bar, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # # Energy
    # with open('LOG/src_czm_interface_epot.pickle', 'wb') as handle:
    #     pickle.dump(Epot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('LOG/src_czm_interface_ekin.pickle', 'wb') as handle:
    #     pickle.dump(Ekin, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('LOG/src_czm_interface_edis.pickle', 'wb') as handle:
    #     pickle.dump(Edis, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('LOG/src_czm_interface_wext.pickle', 'wb') as handle:
    #     pickle.dump(Wext, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # # Variation of energy
    # with open('LOG/src_czm_interface_var_epot.pickle', 'wb') as handle:
    #     pickle.dump(varEpot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('LOG/src_czm_interface_var_ekin.pickle', 'wb') as handle:
    #     pickle.dump(varEkin, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('LOG/src_czm_interface_var_edis.pickle', 'wb') as handle:
    #     pickle.dump(varEdis, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('LOG/src_czm_interface_var_wext.pickle', 'wb') as handle:
    #     pickle.dump(varWext, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    start_time = time.time()
    Run_simulation(DFMesh.strain_rate)
    print("--- %s seconds ---" % (time.time() - start_time))