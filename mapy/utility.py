import numpy as np
import datatime

eigen = np.linalg.eig
savetxt = np.savetxt

def abs2(x):
    return np.abs(x*np.conj(x))

def savestate(title, x, vec, header=None):
    data = [x, abs2(vec), vec.real, vec.imag]
    np.savetxt(title, np.array(data).transpose(), header=header)
    
    
    