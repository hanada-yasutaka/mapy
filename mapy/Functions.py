import sys
import types

import numpy
import sympy
import mpmath
#import mpArray
twopi=2*numpy.pi


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