import numpy as np
from matplotlib import pyplot as plt

import DFMesh
import DFFem
import DFPlot



def getNumberFragments(damage):
    return 1 + sum(1 for D in damage if D > 0.999)


def getFragmentSizes(damage):
    """Returns an array of fragment sizes and its average."""

    n_fragments = getNumberFragments(damage)
    fragments_lengths = np.zeros(n_fragments, dtype=float)

    coord = DFMesh.listDofCoord()
    previous_crack_coord = DFMesh.x0
    j = 0
    for i in range(len(DFMesh.materials)):
        if damage[i] > 0.999:
            fragments_lengths[j] = coord[DFFem.getGlobalIndex(i, 0), 0] - previous_crack_coord
            previous_crack_coord = coord[DFFem.getGlobalIndex(i, 1), 0]
            j = j + 1

    return fragments_lengths, np.average(fragments_lengths)


def getFragmentSizeHistogramData(fragments_legths, hist_number_columns):
    data_histogram = plt.hist(fragments_legths,hist_number_columns)


# Analytical formulation for average fragment size


def GradyFragSize(strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Grady (1982) for a given interval of strain rate."""

    s = np.array(
        [
            ((24.0 * DFMesh.Gc) / (DFMesh.rho * strain**2)) ** (1 / 3)
            for strain in strainrate
        ]
    )

    # Non dimensional form
    norm_strainrate = np.array(
        [
            epsilon
            / (
                ((DFMesh.young_modulus / DFMesh.rho) ** 0.5 * DFMesh.stress_critical**3)
                / (DFMesh.young_modulus**2 * DFMesh.Gc)
            )
            for epsilon in strainrate
        ]
    )
    norm_s = np.array(
        [(24.0 / (normepsilon**2)) ** (1 / 3) for normepsilon in norm_strainrate]
    )
    # norm_s = s / (DFMesh.young_modulus * DFMesh.Gc / DFMesh.stress_critical**2)

    return s, norm_s


def GlenChudnoviskFragSize(strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Glen and Chudnovisk (1986) for a given interval of strain rate."""

    # Non dimensional form
    norm_strainrate = np.array(
        [
            epsilon
            / (
                ((DFMesh.young_modulus / DFMesh.rho) ** 0.5 * DFMesh.stress_critical**3)
                / (DFMesh.young_modulus**2 * DFMesh.Gc)
            )
            for epsilon in strainrate
        ]
    )
    norm_s = np.array(
        [
            4.0
            / normepsilon
            * np.sinh(1.0 / 3.0 * (np.arcsinh(3.0 / 2.0 * normepsilon)))
            for normepsilon in norm_strainrate
        ]
    )
    # norm_s = s / (DFMesh.young_modulus * DFMesh.Gc / DFMesh.stress_critical**2)

    return norm_s


def ZhouMolinariRameshFragSize(strainrate):
    """Returns the fragment size estimation through analytical formulation proposed by Zhou, Molinari and Ramesh (2006) for a given interval of strain rate."""

    # Non dimensional form
    norm_strainrate = np.array(
        [
            epsilon
            / (
                ((DFMesh.young_modulus / DFMesh.rho) ** 0.5 * DFMesh.stress_critical**3)
                / (DFMesh.young_modulus**2 * DFMesh.Gc)
            )
            for epsilon in strainrate
        ]
    )
    # norm_s = np.array([4.5 / (1.0 + 4.5*normepsilon**(-2/3)) for normepsilon in norm_strainrate])
    norm_s = np.array(
        [
            4.5
            * (1.0 + 0.77 * normepsilon ** (1 / 4) + 5.4 * normepsilon ** (3 / 4)) ** -1
            for normepsilon in norm_strainrate
        ]
    )

    return norm_s
