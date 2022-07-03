import numpy as np
import DFMesh
import DFFem


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


