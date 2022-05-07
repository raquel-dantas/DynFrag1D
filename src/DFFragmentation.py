import numpy as np
import DFMesh


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
    for i in range(len(coord)):
        if damage[i] > 0.999:
            fraglen[j] = coord[i,0] - lastrup
            lastrup = coord[i,0]
            j = j + 1
    
    return fraglen, np.average(fraglen)


