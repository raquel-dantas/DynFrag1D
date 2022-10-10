from matplotlib.pyplot import connect
from matplotlib import pyplot as plt
import pickle
import DFMesh
import DFFem
import DFPostprocess
import DFNewmark
import DFInterface
import DFPlot
import DFFragmentation
import numpy as np
import progressbar


def Run_simulation(strain_rate):

    bar = progressbar.ProgressBar(maxval=50, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

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
    avg_stress = np.zeros((DFMesh.n_steps))
    up_bc_left = np.array([0,0])
    up_bc_right = np.array([0,0])

    nfrag = np.zeros((DFMesh.n_steps))
    avg_sfrag = np.zeros((DFMesh.n_steps))
    # fraglen = np.zeros((DFMesh.n_steps, DFMesh.n_el),dtype=float)
    datahist = []



    for n in range(DFMesh.n_steps):        

        progress = int(bar.maxval*float(n/DFMesh.n_steps))
        bar.update(progress)

        # Plots at each time step
        # DFPlot.PlotByDOF(v)
        # DFPlot.PlotByElement(stress)

        # Post process (stress, strain, energies)
        strain, stress, average_stress = DFPostprocess.PostProcess(u)

        # D returns a vector contained damage parameter for cohesive elements
        D = [DFInterface.DamageParameter(el) for el in range(len(DFMesh.materials))]
        # DFPlot.PlotByInterface(D)
        # nfrag retuns a vector contained the number of fragments 
        nfrag[n] = DFFragmentation.NumberFragments(D)
        fraglen, avg_sfrag[n] = DFFragmentation.SizeFragments(D)
        datahist = plt.hist(fraglen,10)

        stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
        avg_stress[n] = DFPostprocess.StressBar(stress, els_step)
        Epot[n], Ekin[n], Edis[n], Erev[n], Econ[n], Wext[n] = DFPostprocess.Energy(up_bc_left,
        up_bc_right, u, v, stress, work)
        work =  Wext[n]

        # DFPlot.PlotVTK('animation/test',n,u,stress)
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
                else:
                    elbc = DFMesh.n_el - 1
                    up_bc_right = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])

        # u,v,acel returns a vector for u,v and acel at every dof at the n step
        u, v, acel = DFNewmark.Newmark_exp(M, u, v, acel, F, DFMesh.dt)

        # Check limit stress for possible insertion of interface elements
        for el in range(DFMesh.n_el-1):
            if average_stress[el] > DFMesh.sigmac[el]:
                # Fracture happens: creat new interface element
                u, v, acel = DFInterface.InsertInterface(el, el+1, u, v, acel)
                els_step = els_step + 1
    

    
    bar.finish()
    print('\n')


    # Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
    varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPostprocess.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)

    # Power [Energy, time] returns the energy difference between consecutive time steps
    PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot = DFPostprocess.Power(Epot, Ekin, Edis, Erev, Econ, Wext)



    # Plots
    DFPlot.PlotAverageStressBar(avg_stress)
    DFPlot.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)
    DFPlot.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot)
    DFPlot.PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot)
    DFPlot.PlotNumberFragments(nfrag)
    DFPlot.PlotFragmentSizeHistogram(fraglen)

    # Save results 
    # Number of fragments
    with open('LOG/src_czm_interface_number_fragments.pickle', 'wb') as handle:
        pickle.dump(nfrag, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Average fragment size
    with open('LOG/src_czm_interface_average_fragment_size.pickle', 'wb') as handle:
        pickle.dump(avg_sfrag, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Histogram fragment size
    with open('LOG/src_czm_interface_datahist_fragment_size.pickle', 'wb') as handle:
        pickle.dump(datahist, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Average stress for the bar 
    with open('LOG/src_czm_interface_avg_stress.pickle', 'wb') as handle:
        pickle.dump(avg_stress, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Energy
    with open('LOG/src_czm_interface_epot.pickle', 'wb') as handle:
        pickle.dump(Epot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_ekin.pickle', 'wb') as handle:
        pickle.dump(Ekin, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_edis.pickle', 'wb') as handle:
        pickle.dump(Edis, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_erev.pickle', 'wb') as handle:
        pickle.dump(Erev, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_econ.pickle', 'wb') as handle:
        pickle.dump(Econ, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_wext.pickle', 'wb') as handle:
        pickle.dump(Wext, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Variation of energy
    with open('LOG/src_czm_interface_var_epot.pickle', 'wb') as handle:
        pickle.dump(varEpot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_var_ekin.pickle', 'wb') as handle:
        pickle.dump(varEkin, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_var_edis.pickle', 'wb') as handle:
        pickle.dump(varEdis, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_var_erev.pickle', 'wb') as handle:
        pickle.dump(varErev, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_var_econ.pickle', 'wb') as handle:
        pickle.dump(varEcon, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('LOG/src_czm_interface_var_wext.pickle', 'wb') as handle:
        pickle.dump(varWext, handle, protocol=pickle.HIGHEST_PROTOCOL)

    

if __name__ == '__main__':
    Run_simulation(DFMesh.strain_rate)