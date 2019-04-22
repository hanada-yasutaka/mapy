import sys
import types

import numpy
import sympy
import mpmath
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

class cosh(object):
    def __new__(cls, x):
        name = cls.__name__
        f = _getattr(x, name)
        return f(x)
    
class tanh(object):
    def __new__(cls, x):
        name = cls.__name__
        f = _getattr(x, name)
        return f(x)    

class arctan(object):
    def __new__(cls, x):
        name = cls.__name__
        f = _getattr(x, name)
        return f(x)    
        
        
def abs2(x):
    return numpy.abs(x*numpy.conj(x))

def eigen(mat,sort=False):
    if isinstance(mat, numpy.matrix):
        mat = numpy.array(mat)
    evals, vecs = numpy.linalg.eig(mat)
    evecs = [vec for vec in vecs.transpose()]
#    return evals, evecs
    if sort:
        index = [i[0] for i in sorted(enumerate(evals.real), key=lambda x:x[1])]
        evals = numpy.array([evals[i] for i in index])
        evecs = [evecs[i] for i in index]
    return evals, evecs





