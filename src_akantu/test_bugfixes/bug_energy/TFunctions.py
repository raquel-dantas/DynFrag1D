import akantu as aka
import numpy as np
from matplotlib import pyplot as plt


class FixedVelocity(aka.DirichletFunctor):
    """Fixed velocity at the boundaries."""

    def __init__(self, axis, vel):
        super().__init__(axis)
        self.axis = axis
        self.time = 0
        self.vel = vel

    def set_time(self, t):
        self.time = t

    def __call__(self, node, flags, disp, coord):
        flags[int(self.axis)] = True
        disp[int(self.axis)] = self.vel * self.time



def getEnergy(energies, energy_name):
    for i in range(len(energies)):
        if energies[i][0] == energy_name:
            energy = energies[i][1]
            return energy


def computeVariationEnergy(energies, n_steps):
    """Returns the variation of energies between the current time step and the time step 0."""

    energy_potential = np.zeros(n_steps)
    energy_kinetic = np.zeros(n_steps)
    energy_dissipated = np.zeros(n_steps)
    energy_contact = np.zeros(n_steps)
    energy_reversible = np.zeros(n_steps)
    external_work = np.zeros(n_steps)
    
    var_energy_potential = np.zeros(n_steps)
    var_energy_kinetic = np.zeros(n_steps)
    var_energy_dissipated = np.zeros(n_steps)
    var_energy_contact = np.zeros(n_steps)
    var_energy_reversible = np.zeros(n_steps)
    var_energy_total = np.zeros(n_steps)
    var_external_work = np.zeros(n_steps)

    for i in range(n_steps):

        energy_potential[i] = getEnergy(energies[i], "energy potential") 
        energy_kinetic[i] = getEnergy(energies[i], "energy kinetic") 
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated") 
        energy_reversible[i] = getEnergy(energies[i], "energy reversible") 
        energy_contact[i] = getEnergy(energies[i], "energy contact") 
        external_work[i] = getEnergy(energies[i], "external work") 

    for n in range(1, n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_energy_reversible[n] = energy_reversible[n] - energy_reversible[0]
        var_energy_contact[n] = energy_contact[n] - energy_contact[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (
            var_energy_potential[n]
            + var_energy_kinetic[n]
            + var_energy_dissipated[n]
            + var_energy_reversible[n]
            + var_energy_contact[n]
        )

    var_energies = [
        ["var energy potential", var_energy_potential],
        ["var energy kinetic", var_energy_kinetic],
        ["var energy dissipated", var_energy_dissipated],
        ["var energy reversible", var_energy_reversible],
        ["var energy contact", var_energy_contact],
        ["var external work", var_external_work],
        ["var energy total", var_energy_total],
    ]

    return var_energies


def computePower(energies, n_steps):
    """Returns the variation of energies between the current time step and the time step 0."""

    energy_potential = np.zeros(n_steps)
    energy_kinetic = np.zeros(n_steps)
    energy_dissipated = np.zeros(n_steps)
    energy_contact = np.zeros(n_steps)
    energy_reversible = np.zeros(n_steps)
    external_work = np.zeros(n_steps)
    
    power_energy_potential = np.zeros(n_steps)
    power_energy_kinetic = np.zeros(n_steps)
    power_energy_dissipated = np.zeros(n_steps)
    power_energy_contact = np.zeros(n_steps)
    power_energy_reversible = np.zeros(n_steps)
    power_energy_total = np.zeros(n_steps)
    power_external_work = np.zeros(n_steps)

    for i in range(n_steps):

        energy_potential[i] = getEnergy(energies[i], "energy potential") 
        energy_kinetic[i] = getEnergy(energies[i], "energy kinetic") 
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated") 
        energy_reversible[i] = getEnergy(energies[i], "energy reversible") 
        energy_contact[i] = getEnergy(energies[i], "energy contact") 
        external_work[i] = getEnergy(energies[i], "external work") 

    for n in range(1, n_steps):
        power_energy_potential[n] = energy_potential[n] - energy_potential[n-1]
        power_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[n-1]
        power_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[n-1]
        power_energy_reversible[n] = energy_reversible[n] - energy_reversible[n-1]
        power_energy_contact[n] = energy_contact[n] - energy_contact[n-1]
        power_external_work[n] = external_work[n] - external_work[n-1]
        power_energy_total[n] = power_external_work[n] - (
            power_energy_potential[n]
            + power_energy_kinetic[n]
            + power_energy_dissipated[n]
            + power_energy_reversible[n]
            + power_energy_contact[n]
        )

    power_energies = [
        ["power energy potential", power_energy_potential],
        ["power energy kinetic", power_energy_kinetic],
        ["power energy dissipated", power_energy_dissipated],
        ["power energy reversible", power_energy_reversible],
        ["power energy contact", power_energy_contact],
        ["power external work", power_external_work],
        ["power energy total", power_energy_total],
    ]

    return power_energies




def plotEnergies(energies, time_simulation, n_steps):


    energy_potential = np.zeros(n_steps)
    energy_kinetic = np.zeros(n_steps)
    energy_dissipated = np.zeros(n_steps)
    energy_contact = np.zeros(n_steps)
    energy_reversible = np.zeros(n_steps)
    external_work = np.zeros(n_steps)
    
    for i in range(n_steps):

        energy_potential[i] = getEnergy(energies[i], "energy potential") 
        energy_kinetic[i] = getEnergy(energies[i], "energy kinetic") 
        energy_dissipated[i] = getEnergy(energies[i], "energy dissipated") 
        energy_reversible[i] = getEnergy(energies[i], "energy reversible") 
        energy_contact[i] = getEnergy(energies[i], "energy contact") 
        external_work[i] = getEnergy(energies[i], "external work") 
    

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Energy"))

    x = np.linspace(0, time_simulation, n_steps)
    plt.plot(x, energy_potential, label="Epot")
    plt.plot(x, energy_kinetic, label="Ekin")
    plt.plot(x, energy_dissipated, label="Edis")
    plt.plot(x, energy_reversible, label="Erev")
    plt.plot(x, energy_contact, label="Econ")
    plt.plot(x, -external_work, label="-Wext")
    plt.legend()
    plt.show()



def plotVarEnergies(var_energies, time_simulation, n_steps):
    """Plot variation of energy from time t to t0."""

    var_energy_potential = getEnergy(var_energies, "var energy potential")
    var_energy_kinetic = getEnergy(var_energies, "var energy kinetic")
    var_energy_dissipated = getEnergy(
        var_energies, "var energy dissipated"
    )
    var_energy_reversible = getEnergy(
        var_energies, "var energy reversible"
    )
    var_energy_contact = getEnergy(var_energies, "var energy contact")
    var_external_work = getEnergy(var_energies, "var external work")
    var_energy_total = getEnergy(var_energies, "var energy total")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Variation of energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Variation of energy"))

    x = np.linspace(0, time_simulation, n_steps)
    plt.plot(x, var_energy_potential, label="varEpot")
    plt.plot(x, var_energy_kinetic, label="varEkin")
    plt.plot(x, var_energy_dissipated, label="varEdis")
    plt.plot(x, var_energy_reversible, label="varErev")
    plt.plot(x, var_energy_contact, label="varEcon")
    plt.plot(x, -var_external_work, label="-varWext")
    plt.plot(x, var_energy_total, label="varEtot")
    plt.legend()
    plt.show()



def plotPower(power, time_simulation, n_steps):
    """Plot variation of energy between time steps"""

    power_energy_potential = getEnergy(power, "power energy potential")
    power_energy_kinetic = getEnergy(power, "power energy kinetic")
    power_energy_dissipated = getEnergy(
        power, "power energy dissipated"
    )
    power_energy_reversible = getEnergy(
        power, "power energy reversible"
    )
    power_energy_contact = getEnergy(power, "power energy contact")
    power_external_work = getEnergy(power, "power external work")
    power_energy_total = getEnergy(power, "power energy total")

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title(str("Variation of energy"))
    plt.xlabel(str("Time (s)"))
    plt.ylabel(str("Variation of energy"))

    x = np.linspace(0, time_simulation, n_steps)
    plt.plot(x, power_energy_potential, label="powerEpot")
    plt.plot(x, power_energy_kinetic, label="powerEkin")
    plt.plot(x, power_energy_dissipated, label="powerEdis")
    plt.plot(x, power_energy_reversible, label="powerErev")
    plt.plot(x, power_energy_contact, label="powerEcon")
    plt.plot(x, -power_external_work, label="-powerWext")
    plt.plot(x, power_energy_total, label="powerEtot")
    plt.legend()
    plt.show()



def plotAverageStressBar(average_stress_bar, time_simulation, n_steps):
    """Plot a vector of values that corresponds to the average stress between all elements at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Average stress bar")
    plt.xlabel("Time (s)")
    plt.ylabel("Average Stress (Pa)")

    x = np.linspace(0, time_simulation, n_steps)
    y = average_stress_bar
    plt.plot(x, y)
    plt.show()



def plotNumberFragments(nfrag, time_simulation, n_steps):
    """Plot a vector of values that corresponds to the number of fragments at eacth time step in the analysis"""

    fig, axes = plt.subplots()
    axes.grid(True, which="both")
    axes.axhline(y=0, color="k")
    plt.title("Number of fragments")
    plt.xlabel("Time (s)")
    plt.ylabel("N")

    x = np.linspace(0, time_simulation, n_steps)
    y = nfrag
    plt.plot(x, y)
    plt.show()

