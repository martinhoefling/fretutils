cimport numpy as np
from numpy import random

def tryGetCythonPhoton(double p0,np.ndarray pvar,long start,long seed):
    """Try to generate a photon and returns 0 and 1 for donor and acceptor and -1 for failure"""
    cdef double rnd
    cdef long end,curndx
    
    end = len(pvar)  

    for curndx in range(start,end):
        rnd=random.random()
        if rnd < p0:
            return 0,curndx
        if rnd < pvar[curndx]:
            return 1,curndx
        
    return -1,len(pvar)