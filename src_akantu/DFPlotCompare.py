import numpy as np
from matplotlib import pyplot as plt




def PlotCompareAverageStressBar(avg_stress_bar_sim1, avg_stress_bar_sim2, time_simulation, n_steps):
    """Plot average stress between all elements at eacth time step in the analysis for different simulatons.\n
    Arguments: \n
    avg_stress_bar_sim1 -- average stress for the whole bar from simulation 1 for all time-steps; \n
    avg_stress_bar_sim2 -- average stress for the whole bar from simulation 2 for all time-steps."""

    plt.title(str("Average stress bar comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Average stress (Pa)"))

    time = np.linspace(0, time_simulation, n_steps)

    plt.plot(time, avg_stress_bar_sim1, label='sim1 - czm_interface')
    plt.plot(time, avg_stress_bar_sim2, label='sim2 - akantu')
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
    plt.plot(time, Epot_sim1, label='sim1 - czm_interface')
    plt.plot(time, Epot_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_epot.svg")
    plt.show()

    # Kinetic energy comparison
    plt.title(str("Kinetic energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Ekin"))
    plt.plot(time, Ekin_sim1, label='sim1 - czm_interface')
    plt.plot(time, Ekin_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_ekin.svg")
    plt.show()

    # Dissipated energy comparison
    plt.title(str("Dissipated energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Edis"))
    plt.plot(time, Edis_sim1, label='sim1 - czm_interface')
    plt.plot(time, Edis_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_edis.svg")
    plt.show()

    # Reversible energy comparison
    plt.title(str("Reversible energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Erev"))
    plt.plot(time, Erev_sim1, label='sim1 - czm_interface')
    plt.plot(time, Erev_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_erev.svg")
    plt.show()

    # Contact energy comparison
    plt.title(str("Contact energy comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Econ"))
    plt.plot(time, Econ_sim1, label='sim1 - czm_interface')
    plt.plot(time, Econ_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_econ.svg")
    plt.show()

    # External work comparison
    plt.title(str("External work comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Wext"))
    plt.plot(time, Wext_sim1, label='sim1 - czm_interface')
    plt.plot(time, Wext_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_wext.svg")
    plt.show()

    # Potential and Kinetic energies comparison
    plt.title(str("Potential and Kinetic energies comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energy"))
    plt.plot(time, Epot_sim1, label='Epot - czm_interface')
    plt.plot(time, Epot_sim2, label='Epot - akantu')
    plt.plot(time, Ekin_sim1, label='Ekin - czm_interface')
    plt.plot(time, Ekin_sim2, label='Ekin - akantu')
    plt.legend()
    plt.savefig("LOG/comp_epot_ekin.svg")
    plt.show()

    # Cohesive energies (Dissipated, reversible, contact)
    plt.title(str("Cohesive energies comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, Edis_sim1, label='Edis - czm_interface')
    plt.plot(time, Edis_sim2, label='Edis - akantu')
    plt.plot(time, Erev_sim1, label='Erev - czm_interface')
    plt.plot(time, Erev_sim2, label='Erev - akantu')
    plt.plot(time, Econ_sim1, label='Econ - czm_interface')
    plt.plot(time, Econ_sim2, label='Econ - akantu')
    plt.legend()
    plt.savefig("LOG/comp_edis_erev_econ.svg")
    plt.show()

    # All energies 
    plt.title(str("Energies comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, Epot_sim1, label='Epot - czm_interface')
    plt.plot(time, Epot_sim2, label='Epot - akantu')
    plt.plot(time, Ekin_sim1, label='Ekin - czm_interface')
    plt.plot(time, Ekin_sim2, label='Ekin - akantu')
    plt.plot(time, Edis_sim1, label='Edis - czm_interface')
    plt.plot(time, Edis_sim2, label='Edis - akantu')
    plt.plot(time, Erev_sim1, label='Erev - czm_interface')
    plt.plot(time, Erev_sim2, label='Erev - akantu')
    plt.plot(time, Econ_sim1, label='Econ - czm_interface')
    plt.plot(time, Econ_sim2, label='Econ - akantu')
    plt.plot(time, Wext_sim1, label='sim1 - czm_interface')
    plt.plot(time, Wext_sim2, label='sim2 - akantu')
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
    plt.plot(time, varEpot_sim1, label='sim1 - czm_interface')
    plt.plot(time, varEpot_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varepot.svg")
    plt.show()

    # Kinetic energy comparison
    plt.title(str("Kinetic energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Ekin"))
    plt.plot(time, varEkin_sim1, label='sim1 - czm_interface')
    plt.plot(time, varEkin_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varekin.svg")
    plt.show()

    # Dissipated energy comparison
    plt.title(str("Dissipated energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Edis"))
    plt.plot(time, varEdis_sim1, label='sim1 - czm_interface')
    plt.plot(time, varEdis_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varedis.svg")
    plt.show()

    # Reversible energy comparison
    plt.title(str("Reversible energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Erev"))
    plt.plot(time, varErev_sim1, label='sim1 - czm_interface')
    plt.plot(time, varErev_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varerev.svg")
    plt.show()

    # Contact energy comparison
    plt.title(str("Contact energy variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Econ"))
    plt.plot(time, varEcon_sim1, label='sim1 - czm_interface')
    plt.plot(time, varEcon_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varecon.svg")
    plt.show()

    # External work comparison
    plt.title(str("External work variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Wext"))
    plt.plot(time, varWext_sim1, label='sim1 - czm_interface')
    plt.plot(time, varWext_sim2, label='sim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varwext.svg")
    plt.show()

    # Potential and Kinetic energies comparison
    plt.title(str("Potential and Kinetic energies variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energy"))
    plt.plot(time, varEpot_sim1, label='varEpot - czm_interface')
    plt.plot(time, varEpot_sim2, label='varEpot - akantu')
    plt.plot(time, varEkin_sim1, label='varEkin - czm_interface')
    plt.plot(time, varEkin_sim2, label='varEkin - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varepot_varekin.svg")
    plt.show()

    # Cohesive energies (Dissipated, reversible, contact)
    plt.title(str("Cohesive energies variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, varEdis_sim1, label='varEdis - czm_interface')
    plt.plot(time, varEdis_sim2, label='varEdis - akantu')
    plt.plot(time, varErev_sim1, label='varErev - czm_interface')
    plt.plot(time, varErev_sim2, label='varErev - akantu')
    plt.plot(time, varEcon_sim1, label='varEcon - czm_interface')
    plt.plot(time, varEcon_sim2, label='varEcon - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varedis_varerev_varecon.svg")
    plt.show()

    # All energies 
    plt.title(str("Energies variation comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energies"))
    plt.plot(time, varEpot_sim1, label='varEpot - czm_interface')
    plt.plot(time, varEpot_sim2, label='varEpot - akantu')
    plt.plot(time, varEkin_sim1, label='varEkin - czm_interface')
    plt.plot(time, varEkin_sim2, label='varEkin - akantu')
    plt.plot(time, varEdis_sim1, label='varEdis - czm_interface')
    plt.plot(time, varEdis_sim2, label='varEdis - akantu')
    plt.plot(time, varErev_sim1, label='varErev - czm_interface')
    plt.plot(time, varErev_sim2, label='varErev - akantu')
    plt.plot(time, varEcon_sim1, label='varEcon - czm_interface')
    plt.plot(time, varEcon_sim2, label='varEcon - akantu')
    plt.plot(time, varWext_sim1, label='varsim1 - czm_interface')
    plt.plot(time, varWext_sim2, label='varsim2 - akantu')
    plt.legend()
    plt.savefig("LOG/comp_varenergies.svg")
    plt.show()


def PlotCompareNumberFragments(nfrag_sim1, nfrag_sim2, time_simulation, n_steps):
    """Plot the number of fragments at eacth time step in the analysis for different simulatons.\n
    Arguments: \n
    nfrag_sim1 -- number of fragments from simulation 1 for all time-steps; \n
    nfrag_sim2 -- number of fragments from simulation 2 for all time-steps."""

    plt.title(str("Number of fragments comparison"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Number of fragments"))

    x = np.linspace(0, time_simulation, n_steps)

    plt.plot(x, nfrag_sim1, label='sim1 - czm_interface')
    plt.plot(x, nfrag_sim2, label='sim2 - akantu')
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

    plt.plot(x, sfrag_sim1, label='sim1 - czm_interface')
    plt.plot(x, sfrag_sim2, label='sim2 - akantu')
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

