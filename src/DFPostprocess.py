from matplotlib.pyplot import connect
import numpy as np
import DFMesh
import DFInterface

def Energy(u, v, nstep):
    """Returns Kinetic, potential, and total energy.\n
    Arguments:\n
    u -- displacements vector;\n
    v -- velocities vector."""

    # Element siffness and mass matrix
    h = DFMesh.L/DFMesh.n_el
    k_elem = DFMesh.E*DFMesh.A/h * np.array([[1.0, -1.0], [-1.0, 1.0]])
    m_elem = DFMesh.rho*DFMesh.A*h/6 * np.array([[2.0, 1.0], [1.0, 2.0]])

    Epot = 0.0
    Ekin = 0.0
    Eext = 0.0
    Edis = 0.0
    Erev = 0.0
    time = nstep * DFMesh.dt

    for el in range(DFMesh.n_el):
        # uloc,vloc[u,v,el] returns vectors contained u and v for local dof
        uloc = np.array([u[DFMesh.connect[el][0]], u[DFMesh.connect[el][1]]])
        vloc = np.array([v[DFMesh.connect[el][0]], v[DFMesh.connect[el][1]]])
        # Epot[kelem,uloc] returns the sum of strain energy values calculate per linear element
        Epot += (1.0/2.0*np.dot(np.matmul(k_elem, uloc), uloc))/DFMesh.A
        # Ekin[melem,vloc] returns the sum of kinetic energy values calculate per linear element
        Ekin += (1.0/2.0*np.dot(np.matmul(m_elem, vloc), vloc))/DFMesh.A

    for el in range(len(DFMesh.materials)):
        if DFMesh.materials[el] == 1:

            # Edis[stress_c, delta_max] returns the sum of dissipated energy caulate  per cohesive element
            Edis += 1.0/2.0*DFMesh.stress_c*DFMesh.delta_max[el]

            # jump_u returns the jump in the displacement between two consecutive linear elements 
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            stress_coh = DFInterface.CohesiveLaw(jump_u,el)

            if jump_u < DFMesh.delta_max[el]:
                # Erev[stress_coh,jump_u] returns the sum of reversible energy caulate per cohesive element for closing cracks (jump_u < delta_max) 
                Erev += 1.0/2.0*stress_coh*jump_u
        
        if DFMesh.materials[el] == 4:
            # Velocity on the boundary
            vo = np.array([-DFMesh.vel, 0])
            # DOF of the applied material type 4
            dofbc = 0
            uloc = np.array([u[DFMesh.connect[dofbc][0]], u[DFMesh.connect[dofbc][1]]])
            # External work
            Eext += (1.0/2.0*np.dot(np.matmul(k_elem, uloc), vo))/DFMesh.A*time
        
        if DFMesh.materials[el] == 5:
            # Velocity on the boundary
            vo = np.array([0, DFMesh.vel])
            # DOF of the applied material type 4
            dofbc = DFMesh.n_el - 1
            uloc = np.array([u[DFMesh.connect[dofbc][0]], u[DFMesh.connect[dofbc][1]]])
            # External work
            Eext += (1.0/2.0*np.dot(np.matmul(k_elem, uloc), vo))/DFMesh.A*time

    return Epot, Ekin, Edis, Eext, Erev


def VarEnergy(Ekin, Epot, Edis, Eext, Erev):
    """Returns the variation of energies between two consecutives time steps."""

    varEkin = np.zeros((DFMesh.n_steps))
    varEpot = np.zeros((DFMesh.n_steps))
    varEdis = np.zeros((DFMesh.n_steps))
    varEext = np.zeros((DFMesh.n_steps))
    varErev = np.zeros((DFMesh.n_steps))
    varEtot = np.zeros((DFMesh.n_steps))
    for n in range(DFMesh.n_steps - 1):
        varEpot[n+1] = Epot[n+1] - Epot[n]
        varEkin[n+1] = Ekin[n+1] - Ekin[n]
        varEdis[n+1] = Edis[n+1] - Edis[n]
        varEext[n+1] = Eext[n+1] - Eext[n]
        varErev[n+1] = Erev[n+1] - Erev[n]
        varEtot[n+1] = varEpot[n+1] + varEkin[n+1] + varEdis[n+1] + varEext[n+1] + varErev[n+1]

    return varEkin, varEpot, varEdis, varEext, varErev, varEtot
    

def PostProcess(u):
    """ Returns the strain and stress vectors, and the average stress vector between two consecutives line elements.\n
    Arguments:\n
    u -- displacemnt vector."""

    numel = len(DFMesh.materials)
    strain = np.zeros(numel)
    stress = np.zeros(numel)

    for el in range(numel):
        if DFMesh.materials[el] == 0:
            # Strain[u,L] returns the strain value at each linear element 'el' 
            strain[el] = (u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]) / DFMesh.h
            # Stress[strain] retuns the stress value at each linear element 'el' 
            stress[el] = DFMesh.E * strain[el]
        elif DFMesh.materials[el] == 1:
            # jump_u[u] returns the jump in the displacement between two consecutive linear elements 
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            # Stress[el retuns the stress value at each cohesive element 'el' 
            stress[el] = DFInterface.CohesiveLaw(jump_u,el)

    # average[stress] returns the average stress between two consecutive linear elements
    average = lambda el: (stress[el] + stress[el+1])/2.0 if DFMesh.connect[el][1] == DFMesh.connect[el+1][0] else 0       
    average_stress = [average(el) for el in range(DFMesh.n_el-1)]

    return strain, stress, average_stress


def LogStress(time_step,stress_evl,current_stress):
    """Returns a matrix that contains the stress for all elements (lin) at all time steps (cols).\n
    Arguments:\n
    time_step: current time step of the analysis;\n
    stress_evl: current stress evolution matrix;\n
    current_stress: the stress vector of the current time step."""

    numel = len(DFMesh.materials)
    for el in range(numel):
        stress_evl[el,time_step] = current_stress[el]

    return stress_evl

def StressBar(current_stress, els_step):
    """Returns the average stress at the whole bar at each time step.\n
    Arguments:\n
    current_stress: the stress vector of the current time step.\n
    els_step: number of elements (linear + cohesive) at the current time step.\n"""
    
    sumstress = 0.0
    numel = len(DFMesh.materials)
    for el in range(numel):
        sumstress += current_stress[el]
    av_stress_bar = sumstress/els_step

    return av_stress_bar
