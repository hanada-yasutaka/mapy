import sys
import types

import numpy
import sympy
import mpmath
#import mpArray
import matplotlib
twopi=2*numpy.pi

savetxt = numpy.savetxt


def _getattr(x, funcname):
    if isinstance(x, sympy.Symbol):
        return getattr(sympy, funcname)
        
    elif isinstance(x, mpmath.mpf) or isinstance(x, mpmath.mpc):
        return getattr(mpmath, funcname)
        
#    elif isinstance(x,mpArray):
#        return getattr(mpArray, funcname)
        
    else:
        return getattr(numpy, funcname)


class sin(object):
    def __new__(cls, x):
        name = cls.__name__
        f = _getattr(x, name)
        return f(x)

class cos(object):
    def __new__(cls, x):
        name = cls.__name__
        f = _getattr(x, name)
        return f(x)
        
        
def abs2(x):
    return numpy.abs(x*numpy.conj(x))

def eigen(mat):
    evals, vecs = numpy.linalg.eig(mat)
    evecs = [vec for vec in vecs.transpose()]
    return evals, evecs

hsm_cdict = {'blue': ((0.0, 0, 1),
                  (0.25, 1, 1),                  
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0)),
        'green': ((0.0, 0, 1),
                  (0.325, 0.1, 0),                  
                  (0.375, 1, 1),
                  (0.64, 1, 1),
                  (0.91, 0, 0),
                  (1, 0, 0)),
        'red': ((0.0, 0, 1), 
                (0.2, 0, 0),
                (0.35, 0, 0), 
                (0.66, 1, 1), 
                (0.89, 1, 1), 
                (1, 0.5, 0.5))}
                  
hsm_cmap = matplotlib.colors.LinearSegmentedColormap('hsm_colormap',hsm_cdict,256)   