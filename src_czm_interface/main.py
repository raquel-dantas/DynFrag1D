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
    avg_stress_bar = np.zeros((DFMesh.n_steps))
    up_bc_left = np.array([0,0])
    up_bc_right = np.array([0,0])

    nfrag = np.zeros((DFMesh.n_steps))
    mean_sfrag = np.zeros((DFMesh.n_steps))
    # fraglen = np.zeros((DFMesh.n_steps, DFMesh.n_el),dtype=float)
    datahist = []


    # Create files to save outputs
    # Save mean size of fragment
    f = str(mean_sfrag)
    average_fragment_size = f
    with open('LOG/avgsize_fragments_czm_interface.txt','w') as f: 
        f.write("--\n")
    # Save data of histogram of size of fragments
    f = str(datahist)
    average_fragment_size = f
    with open('LOG/avgsize_fragments_czm_interface.txt','w') as f: 
        f.write("--\n")
    # Save number of fragments
    f = str(nfrag)
    number_fragments = f
    with open('LOG/number_fragments_czm_interface.txt','w') as f: 
        f.write("--\n")




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
        fraglen, mean_sfrag[n] = DFFragmentation.SizeFragments(D)
        datahist = plt.hist(fraglen,10)

        stress_evl = DFPostprocess.LogStress(n,stress_evl,stress)
        avg_stress_bar[n] = DFPostprocess.StressBar(stress, els_step)
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
    

    # write outputs in .txt
    f = str(datahist)
    datahist = f
    with open('LOG/datahist_czm_interface.txt','a') as f: 
        f.write(datahist + "\n")
    
    f = str(mean_sfrag[n])
    avgsize_fragments = f
    with open('LOG/avgsize_fragments_czm_interface.txt','a') as f: 
        f.write(avgsize_fragments + "\n")

    f = str(nfrag[n])
    number_fragments = f
    with open('LOG/number_fragments_czm_interface.txt','a') as f: 
        f.write(number_fragments + "\n")

    
    bar.finish()
    print('\n')


    # Variation of energy [Energy, time] returns the difference of energy value between time t and t0 
    varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot = DFPostprocess.VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)

    # Power [Energy, time] returns the energy difference between consecutive time steps
    PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot = DFPostprocess.Power(Epot, Ekin, Edis, Erev, Econ, Wext)



    # Plots and results

    # Average stress for the bar 
    DFPlot.PlotAverageStressBar(avg_stress_bar)



    # f = str(avg_stress_bar)
    # average_stress = f
    # with open('LOG/average_stress_czm_interface.txt','w') as f: 
    #     f.write(average_stress)

    # Energy and variation of energy
    DFPlot.PlotEnergy(Epot, Ekin, Edis, Erev, Econ, Wext)
    DFPlot.PlotVarEnergy(varEpot, varEkin, varEdis, varErev, varEcon, varWext, varEtot)
    DFPlot.PlotPower(PEpot, PEkin, PEdis, PErev, PEcon, PWext, PEtot)

    # f = str(Epot)
    # potential_energy = f
    # with open('LOG/energy_pot_czm_interface.txt','w') as f: 
    #     f.write(potential_energy)

    # f = str(varEpot)
    # potential_varenergy = f
    # with open('LOG/varenergy_pot_czm_interface.txt','w') as f: 
    #     f.write(potential_varenergy)

    # f = str(Ekin)
    # kinetic_energy = f
    # with open('LOG/energy_kin_czm_interface.txt','w') as f: 
    #     f.write(kinetic_energy)

    # f = str(varEkin)
    # kinetic_varenergy = f
    # with open('LOG/varenergy_kin_czm_interface.txt','w') as f: 
    #     f.write(kinetic_varenergy)

    # f = str(Edis)
    # dissipated_energy = f
    # with open('LOG/energy_diss_czm_interface.txt','w') as f: 
    #     f.write(dissipated_energy)

    # f = str(varEdis)
    # dissipated_varenergy = f
    # with open('LOG/varenergy_diss_czm_interface.txt','w') as f: 
    #     f.write(dissipated_varenergy)

    # f = str(Econ)
    # contact_energy = f
    # with open('LOG/energy_con_czm_interface.txt','w') as f: 
    #     f.write(contact_energy)

    # f = str(varEcon)
    # contact_varenergy = f
    # with open('LOG/varenergy_con_czm_interface.txt','w') as f: 
    #     f.write(contact_varenergy)

    # f = str(Wext)
    # external_energy = f
    # with open('LOG/energy_external_czm_interface.txt','w') as f: 
    #     f.write(external_energy)

    # f = str(varWext)
    # external_varenergy = f
    # with open('LOG/varenergy_external_czm_interface.txt','w') as f: 
    #     f.write(external_varenergy)



    # Plots and results fragmentation data

    # Number of fragments
    DFPlot.PlotNumberFragments(nfrag)
    with open('LOG/number_fragments_3000_czmint.pickle', 'wb') as handle:
        pickle.dump(nfrag, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Histogram final size fragments 
    DFPlot.PlotFragmentSizeHistogram(fraglen)
    

if __name__ == '__main__':
    Run_simulation(DFMesh.strain_rate)