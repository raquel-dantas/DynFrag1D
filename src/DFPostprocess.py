import numpy as np
import DFMesh
import DFInterface

def Energy(K, M, u, v):
    """Returns Kinetic, potential, and total energy.\n
    Arguments:\n
    K -- global stiffness matrix;\n
    M -- global mass matrix;\n
    u -- displacements vector;\n
    v -- velocities vector.\n
    """
    # Kinetic energy
    Ekin = 1.0/2.0*np.dot(np.matmul(M, v), v)
    # Potential energy
    Epot = 1.0/2.0*np.dot(np.matmul(K, u), u)
    # Dissipated energy
    # External energy
    # Total
    Etot = Ekin + Epot



    return Ekin, Epot, Etot

# def EnergyBalance(Ekin,Epot,Etot):


def PostProcess(u):
    """ Returns the strain, stress vectors, and the average stress vector between two consecutives line elements.\n
    Arguments:\n
    u -- displacemnt vector
    """    
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
            stress[el] = DFInterface.CohesiveLaw(jump_u)

    average = lambda el: (stress[el] + stress[el+1])/2.0 if DFMesh.connect[el][1] == DFMesh.connect[el+1][0] else 0       
       
    average_stress = [average(el) for el in range(DFMesh.n_el-1)]

    return strain, stress, average_stress

def LogStress(time_step,stress_evl,current_stress):
    n_steps = len(stress_evl[0])
    numel = len(DFMesh.materials)

    for el in range(numel):
        stress_evl[el,time_step] = current_stress[el]

    return stress_evl