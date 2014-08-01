import sys
import types

import numpy
import sympy
import mpmath
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

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


hsm_cdict ={'blue': ((0.0, 1, 1),
                  (0.1, 1, 1),
                  (0.25, 1, 1),
                  (0.4, 1, 1),
                  (0.6, 0.56, 0.56),                  
                  (1, 0, 0)),
        'green': ((0.0, 1, 1),
                  (0.1, 1, 1),
                  (0.325, 0, 0),
                  (0.4, 1, 1),
                  (0.6, 1, 1),
                  (0.8, 1, 1),
                  (1, 0, 0)),
        'red': ((0.0, 1, 1),
                (0.1, 1, 1),
                (0.3, 0, 0),
                (0.4, 0, 0),
                (0.6, 0.56, 0.56),
                (0.7, 0.7, 0.7),
                (0.8, 1, 1), 
                (1, 1, 1))}
                 
hsm_cmap = matplotlib.colors.LinearSegmentedColormap('hsm_colormap',hsm_cdict,256)

def plot_hsm(evecs, traj=None,grid=[50,50], vrange=None, contour=None):
    nullfmt   = NullFormatter()         # no labels

    left, width = 0.1, 0.5
    bottom, height = 0.1, 0.5
    bottom_h = left_h = left+width+0.02

    rect_hsm= [left, bottom, width, height]
    rect_qrep = [left, bottom_h, width, 0.3]
    rect_prep = [left_h, bottom, 0.3, height]    
    
    plt.ion() 
    fig = plt.figure(figsize=(7,7))
    ax1 = fig.add_axes(rect_hsm)
    ax2 = fig.add_axes(rect_qrep)
    ax3 = fig.add_axes(rect_prep)


    for n, vec in enumerate(evecs):
        domain = vec.domain
        ax2.set_title("%d/%d-th, grid=[%d,%d]" % (n, len(evecs), grid[0],grid[1]))

        if traj != None: 
            ax1.plot(traj[0],traj[1],',',color='k',alpha=1)
        if contour != None:
            if len(contour) == 3:
                ax1.contour(contour[0],contour[1],contour[2],20, colors='k')
            if len(contour) == 4:
                ax1.contour(contour[0],contour[1],contour[2],contour[3], colors='k')
                                
        x,y,z = vec.hsmrep(grid,vrange)
        levels = numpy.linspace(0,z.max(),50)                    
        ax1.contourf(x,y,z,levels,cmap=hsm_cmap)        
        ax1.set_xlim(domain[0])
        ax1.set_ylim(domain[1])            
        
        ax2.plot(vec.x[0],numpy.log10(vec.abs2()), '-k',lw=3)
        ax2.set_xlim(domain[0])
        ax2.set_ylim(-20,0)

        pvec = vec.q2p()
        ax3.plot(numpy.log10(pvec.abs2()),pvec.x[1], '-k',lw=3)

        ax3.set_ylim(domain[1])
        ax3.set_xlim(-20,0)        
        tic = numpy.array(ax3.get_xticks(), dtype=numpy.int)
        ax3.set_xticklabels(tic,rotation=-90)        
#        ax3.semilogx()        
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.yaxis.set_major_formatter(nullfmt)

        plt.show()
        _ = input("Press [enter] to continue.")
        ax1.cla()
        ax2.cla()        
        ax3.cla()
        plt.cla()

    plt.ioff()
    plt.close()
