import numpy as np

def ExternalWork(mesh, fint, fp_left, fp_right, work_previous_step, vel, dt): 

    Wext = work_previous_step
    # External energy
    # Reaction force at the boundaries 
    nodes_left = mesh.getElementGroup('left').getNodeGroup().getNodes()
    nodes_right = mesh.getElementGroup('right').getNodeGroup().getNodes()
    # Current time step (fn)
    fn_left = -np.sum(fint[nodes_left])
    fn_right = -np.sum(fint[nodes_right])
    # The reaction force (fr) is taken as an average between fn (current time step) and fp (previous time step)
    fr_left = (fn_left + fp_left)*0.5
    fr_right = (fn_right + fp_right)*0.5
    # Stress at the boundary
    stress_bound_left = fr_left 
    stress_bound_right = fr_right 
    # External work (Wext)
    Wext = Wext + (stress_bound_left * -vel + stress_bound_right * vel)*dt

    return Wext, fn_left, fn_right
    


def VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps):
    """Returns the variation of energies between the current time step and the time step 0."""

    varEkin = np.zeros((n_steps))
    varEpot = np.zeros((n_steps))
    varEdis = np.zeros((n_steps))
    varWext = np.zeros((n_steps))
    varErev = np.zeros((n_steps))
    varEcon = np.zeros((n_steps))
    varEtot = np.zeros((n_steps))
    for n in range(1,n_steps):
        varEpot[n] = Epot[n] - Epot[0]
        varEkin[n] = Ekin[n] - Ekin[0]
        varEdis[n] = Edis[n] - Edis[0]
        varErev[n] = Erev[n] - Erev[0]
        varEcon[n] = Econ[n] - Econ[0]
        varWext[n] = Wext[n] - Wext[0]
        varEtot[n] = varWext[n] - (varEpot[n] + varEkin[n] + varEdis[n]  + varErev[n] + varEcon[n])

    return varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot


def Power(Epot, Ekin, Edis, Erev, Econ, Wext, n_steps):
    """Returns the variation of energies between two consecutives time steps."""

    PEkin = np.zeros((n_steps))
    PEpot = np.zeros((n_steps))
    PEdis = np.zeros((n_steps))
    PWext = np.zeros((n_steps))
    PErev = np.zeros((n_steps))
    PEcon = np.zeros((n_steps))
    PEtot = np.zeros((n_steps))
    for n in range(1,n_steps):
        PEpot[n] = Epot[n] - Epot[n-1]
        PEkin[n] = Ekin[n] - Ekin[n-1]
        PEdis[n] = Edis[n] - Edis[n-1]
        PErev[n] = Erev[n] - Erev[n-1]
        PEcon[n] = Econ[n] - Econ[n-1]
        PWext[n] = Wext[n] - Wext[n-1]
        PEtot[n] = PWext[n] - (PEpot[n] + PEkin[n] + PEdis[n]  + PErev[n] + PEcon[n])

    return PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot