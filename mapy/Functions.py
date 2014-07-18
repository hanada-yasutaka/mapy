import sys
import types

import numpy
import sympy
import mpmath
#import mpArray
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
