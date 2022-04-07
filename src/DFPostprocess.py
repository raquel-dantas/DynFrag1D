from matplotlib.pyplot import connect
import numpy as np
import DFMesh
import DFInterface



def Energy(up, u, v, nstep, work_previous_step):
    """Returns potential, kinetic, dissipated, reversible and external energies.\n
    Arguments:\n
    u -- displacements vector;\n
    v -- velocities vector;\n
    work_previous_step -- external energy from the previous time step."""

    # Element siffness and mass matrix
    h = DFMesh.L/DFMesh.n_el
    k_elem = DFMesh.E*DFMesh.A/h * np.array([[1.0, -1.0], [-1.0, 1.0]])
    m_elem = DFMesh.rho*DFMesh.A*h/6 * np.array([[2.0, 1.0], [1.0, 2.0]])

    # Epot = ep
    # Ekin = ek
    Epot = 0.0
    Ekin = 0.0
    Edis = 0.0
    Erev = 0.0
    Econ = 0.0
    # Wext = 0.0
    Wext = work_previous_step

    # time = nstep * DFMesh.dt
    time = DFMesh.dt

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
            # jump_u returns the jump in the displacement between two consecutive linear elements 
            jump_u = u[DFMesh.connect[el][1]] - u[DFMesh.connect[el][0]]
            stress_coh = DFInterface.CohesiveLaw(jump_u,el)

            # Edis[stress_c, delta_max] returns the sum of dissipated energy caulate  per cohesive element
            Edis += 0.5*DFMesh.stress_c*DFMesh.delta_max[el]
            if jump_u < DFMesh.delta_max[el]:
                # Erev[stress_coh,jump_u] returns the sum of reversible energy caulate per cohesive element for closing cracks (jump_u < delta_max) 
                Erev += 0.5*stress_coh*jump_u
            # Contact
            # if jump_u < 0:
            #     alpha = 10.0**15
            #     Econ += 0.5*alpha*jump_u**2

        if DFMesh.materials[el] == 4 or DFMesh.materials[el] == 5:
            # Velocity on the boundary
            if DFMesh.materials[el] == 4:
                vo = np.array([-DFMesh.vel, 0])
                elbc = 0
            else:
                vo = np.array([0, DFMesh.vel])
                elbc = DFMesh.n_el - 1
            uloc = np.array([u[DFMesh.connect[elbc][0]], u[DFMesh.connect[elbc][1]]])
            uploc = np.array([up[DFMesh.connect[elbc][0]], up[DFMesh.connect[elbc][1]]])
            # External work
            f = np.matmul(k_elem, uloc)
            fp = np.matmul(k_elem, uploc)
            fr = (f+fp)*0.5
            stress_bound = fr/DFMesh.A
            work = np.dot(stress_bound,vo)*time
            Wext += work

    return Epot, Ekin, Edis, Erev, Econ, Wext


def VarEnergy(Epot, Ekin, Edis, Erev, Econ, Wext):
    """Returns the variation of energies between two consecutives time steps."""

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


# def VerifyFractureState(stress, u, v, up, vp):
#     for el in range(len(DFMesh.materials)):
#         if DFMesh.materials[el] == 1:
#             if stress[el] == 0.0:
#                 u[DFMesh.connect[el][0]] = 0.0
#                 u[DFMesh.connect[el][1]] = 0.0
#                 v[DFMesh.connect[el][0]] = 0.0
#                 v[DFMesh.connect[el][1]] = 0.0

#                 if DFMesh.connect[el][0] == DFMesh.connect[0][1]:
#                     u[DFMesh.connect[0][0]] = 0.0
#                     v[DFMesh.connect[0][0]] = 0.0

#                 if DFMesh.connect[el][1] == DFMesh.connect[DFMesh.n_el-1][0]:
#                     u[DFMesh.connect[DFMesh.n_el-1][1]] = 0.0
#                     v[DFMesh.connect[DFMesh.n_el-1][1]] = 0.0

#                 # u[DFMesh.connect[el][0]] = up[DFMesh.connect[el][0]]
#                 # u[DFMesh.connect[el][1]] = up[DFMesh.connect[el][1]]
#                 # v[DFMesh.connect[el][0]] = vp[DFMesh.connect[el][0]]
#                 # v[DFMesh.connect[el][1]] = vp[DFMesh.connect[el][1]]

#                 # if DFMesh.connect[el][0] == DFMesh.connect[0][1]:
#                 #     u[DFMesh.connect[0][0]] = up[DFMesh.connect[0][0]]
#                 # if DFMesh.connect[el][1] == DFMesh.connect[DFMesh.n_el-1][0]:
#                 #     u[DFMesh.connect[DFMesh.n_el-1][1]] = up[DFMesh.connect[DFMesh.n_el-1][1]]