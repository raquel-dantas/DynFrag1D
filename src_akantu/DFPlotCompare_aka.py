import numpy as np
from matplotlib import pyplot as plt

# czm_interface
# akantu


def PlotCompareAverageStressBar(avg_stress_bar_sim1, avg_stress_bar_sim2, time_simulation, n_steps_sim1, n_steps_sim2):
    """Plot average stress between all elements at eacth time step in the analysis for different simulatons.\n
    Arguments: \n
    avg_stress_bar_sim1 -- average stress for the whole bar from simulation 1 for all time-steps; \n
    avg_stress_bar_sim2 -- average stress for the whole bar from simulation 2 for all time-steps."""

    

    plt.title(str("Average stress bar comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Average stress (Pa)"))

    time_sim1 = np.linspace(0, time_simulation, n_steps_sim1)
    time_sim2 = np.linspace(0, time_simulation, n_steps_sim2)

    plt.plot(time_sim1, avg_stress_bar_sim1, label='sim1')
    plt.plot(time_sim2, avg_stress_bar_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_avg_stress.svg")
    plt.show()



def PlotCompareEnergies(Epot_sim1, Ekin_sim1, Edis_sim1, Erev_sim1, Econ_sim1, Wext_sim1,Epot_sim2, Ekin_sim2, Edis_sim2, Erev_sim2, Econ_sim2, Wext_sim2, time_simulation, n_steps):
    """Plot energies values at eacth time step in the analysis for different simulatons."""

    time = np.linspace(0, time_simulation, n_steps)

    # Potential energy comparison
    plt.title(str("Potential energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Epot"))
    plt.plot(time, Epot_sim1, label='sim1')
    plt.plot(time, Epot_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_epot.svg")
    plt.show()

    # Kinetic energy comparison
    plt.title(str("Kinetic energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Ekin"))
    plt.plot(time, Ekin_sim1, label='sim1')
    plt.plot(time, Ekin_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_ekin.svg")
    plt.show()

    # Dissipated energy comparison
    plt.title(str("Dissipated energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Edis"))
    plt.plot(time, Edis_sim1, label='sim1')
    plt.plot(time, Edis_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_edis.svg")
    plt.show()

    # Reversible energy comparison
    plt.title(str("Reversible energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Erev"))
    plt.plot(time, Erev_sim1, label='sim1')
    plt.plot(time, Erev_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_erev.svg")
    plt.show()

    # Contact energy comparison
    plt.title(str("Contact energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Econ"))
    plt.plot(time, Econ_sim1, label='sim1')
    plt.plot(time, Econ_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_econ.svg")
    plt.show()

    # External work comparison
    plt.title(str("External work comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Wext"))
    plt.plot(time, Wext_sim1, label='sim1')
    plt.plot(time, Wext_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_wext.svg")
    plt.show()

    # Potential and Kinetic energies comparison
    plt.title(str("Potential and Kinetic energies comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energy"))
    plt.plot(time, Epot_sim1, label='Epot')
    plt.plot(time, Epot_sim2, label='Epot')
    plt.plot(time, Ekin_sim1, label='Ekin')
    plt.plot(time, Ekin_sim2, label='Ekin')
    plt.legend()
    plt.savefig("LOG/comp_epot_ekin.svg")
    plt.show()

    # Cohesive energies (Dissipated, reversible, contact)
    plt.title(str("Cohesive energies comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, Edis_sim1, label='Edis')
    plt.plot(time, Edis_sim2, label='Edis')
    plt.plot(time, Erev_sim1, label='Erev')
    plt.plot(time, Erev_sim2, label='Erev')
    plt.plot(time, Econ_sim1, label='Econ')
    plt.plot(time, Econ_sim2, label='Econ')
    plt.legend()
    plt.savefig("LOG/comp_edis_erev_econ.svg")
    plt.show()

    # All energies 
    plt.title(str("Energies comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, Epot_sim1, label='Epot')
    plt.plot(time, Epot_sim2, label='Epot')
    plt.plot(time, Ekin_sim1, label='Ekin')
    plt.plot(time, Ekin_sim2, label='Ekin')
    plt.plot(time, Edis_sim1, label='Edis')
    plt.plot(time, Edis_sim2, label='Edis')
    plt.plot(time, Erev_sim1, label='Erev')
    plt.plot(time, Erev_sim2, label='Erev')
    plt.plot(time, Econ_sim1, label='Econ')
    plt.plot(time, Econ_sim2, label='Econ')
    plt.plot(time, Wext_sim1, label='sim1')
    plt.plot(time, Wext_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_energies.svg")
    plt.show()



def PlotCompareVarEnergies(varEpot_sim1, varEkin_sim1, varEdis_sim1, varErev_sim1, varEcon_sim1, varWext_sim1, varEpot_sim2, varEkin_sim2, varEdis_sim2, varErev_sim2, varEcon_sim2, varWext_sim2, time_simulation, n_steps):
    """Plot variation of energies values at eacth time step in the analysis for different simulatons."""

    time = np.linspace(0, time_simulation, n_steps)

    # Potential energy comparison
    plt.title(str("Potential energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Epot"))
    plt.plot(time, varEpot_sim1, label='sim1')
    plt.plot(time, varEpot_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_varepot.svg")
    plt.show()

    # Kinetic energy comparison
    plt.title(str("Kinetic energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Ekin"))
    plt.plot(time, varEkin_sim1, label='sim1')
    plt.plot(time, varEkin_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_varekin.svg")
    plt.show()

    # Dissipated energy comparison
    plt.title(str("Dissipated energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Edis"))
    plt.plot(time, varEdis_sim1, label='sim1')
    plt.plot(time, varEdis_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_varedis.svg")
    plt.show()

    # Reversible energy comparison
    plt.title(str("Reversible energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Erev"))
    plt.plot(time, varErev_sim1, label='sim1')
    plt.plot(time, varErev_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_varerev.svg")
    plt.show()

    # Contact energy comparison
    plt.title(str("Contact energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Econ"))
    plt.plot(time, varEcon_sim1, label='sim1')
    plt.plot(time, varEcon_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_varecon.svg")
    plt.show()

    # External work comparison
    plt.title(str("External work variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Wext"))
    plt.plot(time, varWext_sim1, label='sim1')
    plt.plot(time, varWext_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_varwext.svg")
    plt.show()

    # Potential and Kinetic energies comparison
    plt.title(str("Potential and Kinetic energies variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energy"))
    plt.plot(time, varEpot_sim1, label='varEpot')
    plt.plot(time, varEpot_sim2, label='varEpot')
    plt.plot(time, varEkin_sim1, label='varEkin')
    plt.plot(time, varEkin_sim2, label='varEkin')
    plt.legend()
    plt.savefig("LOG/comp_varepot_varekin.svg")
    plt.show()

    # Cohesive energies (Dissipated, reversible, contact)
    plt.title(str("Cohesive energies variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, varEdis_sim1, label='varEdis')
    plt.plot(time, varEdis_sim2, label='varEdis')
    plt.plot(time, varErev_sim1, label='varErev')
    plt.plot(time, varErev_sim2, label='varErev')
    plt.plot(time, varEcon_sim1, label='varEcon')
    plt.plot(time, varEcon_sim2, label='varEcon')
    plt.legend()
    plt.savefig("LOG/comp_varedis_varerev_varecon.svg")
    plt.show()

    # All energies 
    plt.title(str("Energies variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, varEpot_sim1, label='varEpot')
    plt.plot(time, varEpot_sim2, label='varEpot')
    plt.plot(time, varEkin_sim1, label='varEkin')
    plt.plot(time, varEkin_sim2, label='varEkin')
    plt.plot(time, varEdis_sim1, label='varEdis')
    plt.plot(time, varEdis_sim2, label='varEdis')
    plt.plot(time, varErev_sim1, label='varErev')
    plt.plot(time, varErev_sim2, label='varErev')
    plt.plot(time, varEcon_sim1, label='varEcon')
    plt.plot(time, varEcon_sim2, label='varEcon')
    plt.plot(time, varWext_sim1, label='varsim1')
    plt.plot(time, varWext_sim2, label='varsim2')
    plt.legend()
    plt.savefig("LOG/comp_varenergies.svg")
    plt.show()


def PlotCompareNumberFragments(nfrag_sim1, nfrag_sim2, nfrag_sim3, nfrag_sim4, nfrag_sim5, time_simulation, n_steps_sim1, n_steps_sim2, n_steps_sim3, n_steps_sim4, n_steps_sim5):
    """Plot the number of fragments at eacth time step in the analysis for different simulatons.\n
    Arguments: \n
    nfrag_sim1 -- number of fragments from simulation 1 for all time-steps; \n
    nfrag_sim2 -- number of fragments from simulation 2 for all time-steps."""

    plt.title(str("Number of fragments comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Number of fragments"))

    time_sim1 = np.linspace(0, time_simulation, n_steps_sim1)
    time_sim2 = np.linspace(0, time_simulation, n_steps_sim2)
    time_sim3 = np.linspace(0, time_simulation, n_steps_sim3)
    time_sim4 = np.linspace(0, time_simulation, n_steps_sim4)
    time_sim5 = np.linspace(0, time_simulation, n_steps_sim5)

    plt.plot(time_sim1, nfrag_sim1, label='sim1')
    plt.plot(time_sim2, nfrag_sim2, label='sim2')
    plt.plot(time_sim3, nfrag_sim3, label='sim3')
    plt.plot(time_sim4, nfrag_sim4, label='sim4')
    plt.plot(time_sim5, nfrag_sim5, label='sim5')
    plt.legend()
    plt.savefig("LOG/comp_nfrags.svg")
    plt.show()



def PlotCompareAvgFragmentSize(sfrag_sim1, sfrag_sim2, time_simulation, n_steps):
    """Plot the average size of fragments at eacth time step in the analysis for different simulatons.\n
    Arguments: \n
    sfrag_sim1 -- mean fragment size from simulation 1 for all time-steps; \n
    sfrag_sim2 -- mean fragment size from simulation 2 for all time-steps."""

    plt.title(str("Mean fragment size comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Fragment size (m)"))

    x = np.linspace(0, time_simulation, n_steps)

    plt.plot(x, sfrag_sim1, label='sim1')
    plt.plot(x, sfrag_sim2, label='sim2')
    plt.legend()
    plt.savefig("LOG/comp_sfrag.svg")
    plt.show()



def PlotCompareFragmentSizeHistogram(frag_sizes_sim1,frag_sizes_sim2):
    """Plot fragments size distribution for different simulatons.\n
    Arguments: \n
    frag_size_sim1 -- fragment size from simulation 1 for an specific time-step; \n
    frag_size_sim2 -- fragment size from simulation 1 for an specific time-step."""

    fig, axes = plt.subplots()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    plt.title("Fragment size distribution")
    plt.xlabel("Fragment size (m)")
    plt.ylabel("Number of fragments")

    plt.hist(frag_sizes_sim1,10)
    plt.hist(frag_sizes_sim2,10)
    plt.savefig("LOG/comp_fragment_size_distribution.svg")
    plt.show()

