from matplotlib.pyplot import connect
import numpy as np
import DFMesh
import DFInterface



def Energy(up_bc_left, up_bc_right, u, v, work_previous_step):
    """Returns potential, kinetic, dissipated, reversible and external energies.\n
    Arguments:\n
    u -- displacements vector;\n
    v -- velocities vector;\n
    work_previous_step -- external energy from the previous time step."""
    
    # Element siffness and mass matrix
    h = DFMesh.L/DFMesh.n_el
    k_elem = DFMesh.E*DFMesh.A/h * np.array([[1.0, -1.0], [-1.0, 1.0]])
    m_elem = DFMesh.rho*DFMesh.A*h/6 * np.array([[2.0, 1.0], [1.0, 2.0]])
    # Initial values for the time step n
    Epot = 0.0
    Ekin = 0.0
    Edis = 0.0
    Erev = 0.0
    Econ = 0.0
    Wext = work_previous_step
    Ekin_bc = 0.0
    
    for el in range(len(DFMesh.materials)):
        if DFMesh.materials[el] == 0:
            # uloc,vloc[u,v,el] returns vectors contained u and v for local dof
            uloc = np.array([u[DFMesh.connect[el][0]], u[DFMesh.connect[el][1]]])
            vloc = np.array([v[DFMesh.connect[el][0]], v[DFMesh.connect[el][1]]])
            # Epot[kelem,uloc] returns the sum of strain energy values calculate per linear element
            Epot += (0.5*np.dot(np.matmul(k_elem, uloc), uloc))/DFMesh.A
            # Ekin[melem,vloc] returns the sum of kinetic energy values calculate per linear element
            Ekin += (0.5*np.dot(np.matmul(m_elem, vloc), vloc))/DFMesh.A

        if DFMesh.materials[el] == 1:
            # Edis[stress_c, delta_max] returns the sum of dissipated energy caulate  per cohesive element
            Edis += 0.5*DFMesh.stress_c*DFMesh.delta_max[el]
            # jump_u returns the jump in the displacement between two consecutive linear elements 
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            # stress_coh returns the stress in the cohesive elements give by an cohesive law 
            stress_coh = DFInterface.CohesiveLaw(jump_u,el)
            if jump_u < DFMesh.delta_max[el]:
                D = DFInterface.DamageParameter(el)
                # Erev[stress_coh,jump_u] returns the sum of reversible energy caulate per cohesive element for closing cracks (jump_u < delta_max) 
                Erev += 0.5*stress_coh*jump_u
            # Contact
            # if jump_u < 0:
            #     alpha = 10.0**15
            #     Econ += 0.5*alpha*jump_u**2

        if DFMesh.materials[el] == 4 or DFMesh.materials[el] == 5:
            # vo is the velocity applied on the extremity
            # elbc is the element index of the applied velocity
            # uploc is the local displacement of elbc from the previous time step 
            # c is the correct contribution of the local velocity vector of boundary elements on the calculus of kinectic energy
            if DFMesh.materials[el] == 4:
                vo = np.array([-DFMesh.vel, 0])
                elbc = 0
                uploc = up_bc_left
                c = np.array([0.0, v[DFMesh.connect[elbc][1]]])
            else:
                vo = np.array([0, DFMesh.vel])
                elbc = DFMesh.n_el - 1
                uploc = up_bc_right
                c = np.array([v[DFMesh.connect[elbc][0]], 0.0])
            # uloc,vloc[u,v,elbc] returns vectors contained u and v for local dofs of elbc
            uloc = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])
            # fn is the internal force on the current time step
            fn = np.matmul(k_elem, uloc)
            # fp is the internal force on the previous time step 
            fp = np.matmul(k_elem, uploc)
            # The reaction force is taken s an average between fn and fp
            fr = (fn+fp)*0.5
            # Stress on the boundary
            stress_bound = fr/DFMesh.A
            # Work is power integrated in time
            work = np.dot(stress_bound,vo)*DFMesh.dt
            Wext += work

            # Correction of the Kinectic energy: subtract the kinectic energy from the boundary
            Ekin_bc += (0.5*np.dot(np.matmul(m_elem, vo), (2.*c + vo)))/DFMesh.A
    
    Ekin = Ekin - Ekin_bc

    return Epot, Ekin, Edis, Erev, Econ, Wext


def VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext):
    """Returns the variation of energies between the current time step and the time step 0."""

    varEkin = np.zeros((DFMesh.n_steps))
    varEpot = np.zeros((DFMesh.n_steps))
    varEdis = np.zeros((DFMesh.n_steps))
    varWext = np.zeros((DFMesh.n_steps))
    varErev = np.zeros((DFMesh.n_steps))
    varEcon = np.zeros((DFMesh.n_steps))
    varEtot = np.zeros((DFMesh.n_steps))
    for n in range(1,DFMesh.n_steps):
        varEpot[n] = Epot[n] - Epot[0]
        varEkin[n] = Ekin[n] - Ekin[0]
        varEdis[n] = Edis[n] - Edis[0]
        varErev[n] = Erev[n] - Erev[0]
        varEcon[n] = Econ[n] - Econ[0]
        varWext[n] = Wext[n] - Wext[0]
        varEtot[n] = varWext[n] - (varEpot[n] + varEkin[n] + varEdis[n]  + varErev[n] + varEcon[n])

    return varEkin, varEpot, varEdis, varErev, varEcon, varWext, varEtot

def Power(Epot, Ekin, Edis, Erev, Econ, Wext):
    """Returns the variation of energies between two consecutives time steps."""

    PEkin = np.zeros((DFMesh.n_steps))
    PEpot = np.zeros((DFMesh.n_steps))
    PEdis = np.zeros((DFMesh.n_steps))
    PWext = np.zeros((DFMesh.n_steps))
    PErev = np.zeros((DFMesh.n_steps))
    PEcon = np.zeros((DFMesh.n_steps))
    PEtot = np.zeros((DFMesh.n_steps))
    for n in range(1,DFMesh.n_steps):
        PEpot[n] = Epot[n] - Epot[n-1]
        PEkin[n] = Ekin[n] - Ekin[n-1]
        PEdis[n] = Edis[n] - Edis[n-1]
        PErev[n] = Erev[n] - Erev[n-1]
        PEcon[n] = Econ[n] - Econ[n-1]
        PWext[n] = Wext[n] - Wext[n-1]
        PEtot[n] = PWext[n] - (PEpot[n] + PEkin[n] + PEdis[n]  + PErev[n] + PEcon[n])

    return PEkin, PEpot, PEdis, PErev, PEcon, PWext, PEtot
    

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
    els_step: number of elements (linear + cohesive) at the current time step."""

    sumstress = 0.0
    numel = len(DFMesh.materials)
    for el in range(numel):
        sumstress += current_stress[el]
    av_stress_bar = sumstress/els_step

    return av_stress_bar


# def FractureStateBoundary(u):
#     fint = DFInterface.ForceInt(u)
#     if DFMesh.materials.__contains__(1):
#         for el in range(len(DFMesh.materials)):
#             if DFMesh.materials[el] == 4 or DFMesh.materials[el] == 5:
#                 if DFMesh.materials[el] == 4:
#                     elbc = 0
#                     if fint[DFMesh.connect[elbc][1]] == 0.0:
#                         DFMesh.materials[el] = 6
#                 else:
#                     elbc = DFMesh.n_el - 1
#                     if fint[DFMesh.connect[elbc][0]] == 0.0:
#                         DFMesh.materials[el] = 6