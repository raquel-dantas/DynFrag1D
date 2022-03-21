from matplotlib.pyplot import connect
import numpy as np
import DFMesh
import DFInterface

def Energy(u, v):
    """Returns Kinetic, potential, and total energy.\n
    Arguments:\n
    u -- displacements vector;\n
    v -- velocities vector."""

    h = DFMesh.L/DFMesh.n_el
    k_elem = DFMesh.E*DFMesh.A/h * np.array([[1.0, -1.0], [-1.0, 1.0]])
    m_elem = DFMesh.rho*DFMesh.A*h/6 * np.array([[2.0, 1.0], [1.0, 2.0]])

    for el in range(DFMesh.n_el):
        Epot = 0.0
        Ekin = 0.0
        uloc = np.array([u[DFMesh.connect[el][0]], u[DFMesh.connect[el][1]]])
        vloc = np.array([v[DFMesh.connect[el][0]], v[DFMesh.connect[el][1]]])
        Epot += 1.0/2.0*np.dot(np.matmul(k_elem, uloc), uloc)
        Ekin += 1.0/2.0*np.dot(np.matmul(m_elem, vloc), vloc)
        Etot = Ekin + Epot
    
    return Epot, Ekin, Etot


def VarEnergy(Ekin, Epot, Etot):
    """Returns the variation of energies between two consecutives time steps."""

    varEkin = np.zeros((DFMesh.n_steps))
    varEpot = np.zeros((DFMesh.n_steps))
    varEtot = np.zeros((DFMesh.n_steps))
    for n in range(DFMesh.n_steps - 1):
        varEpot[n+1] = Epot[n+1] - Epot[n]
        varEkin[n+1] = Ekin[n+1] - Ekin[n]
        varEtot[n+1] = Etot[n+1] - Etot[n]

    return varEkin, varEpot, varEtot
    

def PostProcess(u):
    """ Returns the strain and stress vectors, and the average stress vector between two consecutives line elements.\n
    Arguments:\n
    u -- displacemnt vector."""

    numel = len(DFMesh.materials)
    strain = np.zeros(numel)
    stress = np.zeros(numel)

    for el in range(numel):
        if DFMesh.materials[el] == 0:
            # Strain[u,L] returns the strain value at each element 'el' for the step 'n'
            strain[el] = (u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]) / DFMesh.h
            # Stress[strain] retuns the stress value at each element 'el' for the step 'n'
            stress[el] = DFMesh.E * strain[el]
        elif DFMesh.materials[el] == 1:
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            stress[el] = DFInterface.CohesiveLaw(jump_u,el)

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