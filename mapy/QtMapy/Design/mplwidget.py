from PyQt4 import QtGui,QtCore
#from PyQt4 import QtGui,QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar2


import numpy as np
twopi=2*np.pi

class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        self.updateGeometry()

        self.position=[None, None]
"""        
    def onclick(self,event):
        try:
            print('(x,y) = (%f, %f) :button=%d, '%(event.xdata, event.ydata, event.button))
            self.position = (event.xdata, event.ydata)
        except TypeError:
            self.position = (None, None)            
            pass
"""            

class MplWidget(QtGui.QWidget):
    """Widget defined in Qt Designer"""
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.navigation_toolbar = NavigationToolbar2(self.canvas, self)

        """ connect click event """
#        self.canvas.mpl_connect('button_press_event', self.canvas.onclick)


        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.navigation_toolbar,0)
        self.setLayout(self.vbl)

