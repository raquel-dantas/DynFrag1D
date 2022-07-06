import numpy as np
import DFMesh
import DFFem
import DFPlot


# Number and size of fragments

def NumberFragments(damage):
    return 1 + sum(1 for D in damage if D > 0.999)


def SizeFragments(damage):
    """Returns an array of fragment sizes and its average.    
    """

    nfrag = NumberFragments(damage)
    fraglen = np.zeros(nfrag, dtype=float)

    coord = DFMesh.ListDofCoord()
    lastrup = DFMesh.x0
    j = 0
    for i in range(len(DFMesh.materials)):
        if damage[i] > 0.999:
            fraglen[j] = coord[DFFem.Gl_index(i,0) , 0] - lastrup
            lastrup = coord[DFFem.Gl_index(i,1) , 0]
            j = j + 1
    
    return fraglen, np.average(fraglen)


# Analytical formulation for average fragment size

def GradyFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Grady (1982) for a given interval of strain rate."""

    s = np.array([((24.* DFMesh.Gc)/(DFMesh.rho * epsilon**2))**(1/3) for epsilon in values_strainrate])

    return s



def GlenChudnoviskFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Glen and Chudnovisk (1986) for a given interval of strain rate."""

    s = np.zeros(len(values_strainrate))
    for i in range(len(values_strainrate)):
        alpha = 3.*DFMesh.stress_c / (DFMesh.rho * DFMesh.E * values_strainrate[i]**2)
        beta = 3.*DFMesh.Gc / (2.* DFMesh.rho * DFMesh.E * values_strainrate[i]**2)
        phi = np.arcsinh(beta*(3./alpha)**(3/2))
        s[i] = 4.*(alpha/3.)*np.sinh(phi/3.)

    return s



def ZhouMolinariRameshFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Zhou, Molinari and Ramesh (2006) for a given interval of strain rate."""

    s = np.array([4.5 / (1 + 4.5*epsilon**(-2/3)) for epsilon in values_strainrate])

    return s

