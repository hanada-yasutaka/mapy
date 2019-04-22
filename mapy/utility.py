import re
import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

def sort_index(l, ascending=True):
    index = [i[0] for i in sorted(enumerate(l), key=lambda x:x[1])]
    return index if ascending else index[::-1]

def nsort( l ): 
    """ 
    Sort the given list in the way that humans expect. For example,
    
    >>> l = ["1", "4", "2" , "0"]
    >>> natural(l)
    >>> l
    ['0', '1', '2', '4']
    
    Note that negative values (its argument takes only string!) is shifted to the backward  

    >>> l = ["-1","3","-2", "4"]
    >>> natural(l)
    >>> l
    ['3', '4', '-1', '-2']
     
    """ 
    ll = l[:]
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    ll.sort( key=alphanum_key )
    return ll

def hsm_axes(fig,nfm=True):
    nullfmt   = NullFormatter()         # no labels

    left, width = 0.1, 0.5
    bottom, height = 0.1, 0.5
    bottom_h = left_h = left+width+0.01

    rect_hsm= [left, bottom, width, height]
    rect_qrep = [left, bottom_h, width, 0.3]
    rect_prep = [left_h, bottom, 0.3, height]    
    
    ax1 = fig.add_axes(rect_hsm)
    ax2 = fig.add_axes(rect_qrep)
    ax3 = fig.add_axes(rect_prep)
    if nfm:
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.yaxis.set_major_formatter(nullfmt)

    return [ax1, ax2,ax3]

def hsmax(fig,nfm=True):
    nullfmt   = NullFormatter()         # no labels

    left, width = 0.1, 0.5
    bottom, height = 0.1, 0.5
    bottom_h = left_h = left+width+0.02

    rect_hsm= [left, bottom, width, height]
    rect_qrep = [left, bottom_h, width, 0.3]
    rect_prep = [left_h, bottom, 0.3, height]    
    
    ax1 = fig.add_axes(rect_hsm)
    ax2 = fig.add_axes(rect_qrep)
    ax3 = fig.add_axes(rect_prep)
    if nfm:
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.yaxis.set_major_formatter(nullfmt)

    return [ax1, ax2,ax3]


def plot_hsm(evecs, ref=None,traj=None,grid=[50,50], vrange=None, contour=None, vmin=[-30,-30]):
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
        if ref!=None:
            qvec = ref[n]
            ax2.plot(ref[n].x[0], numpy.log10(qvec.abs2()), '--k',lw=2)        
        
        ax2.set_xlim(domain[0])
        ax2.set_ylim(vmin[0],0)

        pvec = vec.q2p()
#        index = pvec.abs2() > 1e-30
#        ax3.plot(numpy.log10(pvec.abs2()[index]),pvec.x[1][index], '-k',lw=3)
        ax3.plot(numpy.log10(pvec.abs2()),pvec.x[1], '-k',lw=3)
        
        if ref!=None:
            pvec = ref[n].q2p()
            #index = pvec.abs2() > vmin[1]
            ax3.plot(numpy.log10(pvec.abs2()),pvec.x[1], '--k',lw=2)        
        

        ax3.set_ylim(domain[1])
        ax3.set_xlim(vmin[1],0)        
        tic = numpy.array(ax3.get_xticks(), dtype=numpy.int)
        ax3.set_xticklabels(tic,rotation=-90)        
#        ax3.semilogx()        
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.yaxis.set_major_formatter(nullfmt)

        plt.show()
        str_in = input("Press [enter] to continue.")
        if str_in in ["save", "s"]:
            fig.savefig("hsm_rep_%d.png" % n)
        ax1.cla()
        ax2.cla()        
        ax3.cla()
        plt.cla()

    plt.ioff()
    plt.close()
    
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
hsm_cmap._init()
#hsm_cmap._lut[:,-1] = numpy.linspace(0, 1, hsm_cmap.N+3)
