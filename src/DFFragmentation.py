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


# Analytical formulation for average fragment size (non dimensional form)

def GradyFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Grady (1982) for a given interval of strain rate."""

    norm_strainrate = np.array([epsilon / (( (DFMesh.E / DFMesh.rho)**0.5 * DFMesh.stress_c**3) / (DFMesh.E**2 * DFMesh.Gc) ) for epsilon in values_strainrate])

    s = np.array([(24./(normepsilon**2))**(1/3) for normepsilon in norm_strainrate])
    norm_s = s / (DFMesh.E * DFMesh.Gc / DFMesh.stress_c**2)

    return norm_s



def GlenChudnoviskFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Glen and Chudnovisk (1986) for a given interval of strain rate."""

    norm_strainrate = np.array([epsilon / (( (DFMesh.E / DFMesh.rho)**0.5 * DFMesh.stress_c**3) / (DFMesh.E**2 * DFMesh.Gc) ) for epsilon in values_strainrate])

    s = np.array([4./normepsilon * np.sinh(1./3.*(np.sinh(3./2.*normepsilon))**-1) for normepsilon in norm_strainrate])
    norm_s = s / (DFMesh.E * DFMesh.Gc / DFMesh.stress_c**2)
  
    return norm_s



def ZhouMolinariRameshFragSize(values_strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Zhou, Molinari and Ramesh (2006) for a given interval of strain rate."""

    norm_strainrate = np.array([epsilon / (( (DFMesh.E / DFMesh.rho)**0.5 * DFMesh.stress_c**3) / (DFMesh.E**2 * DFMesh.Gc) ) for epsilon in values_strainrate])

    norm_s = np.array([4.5 / (1. + 4.5*normepsilon**-2/3) for normepsilon in norm_strainrate])

    return norm_s

