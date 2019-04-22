from PyQt5 import QtGui,QtCore,QtWidgets
#from PyQt4 import QtWidgets,QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar2


import numpy as np
twopi=2*np.pi

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)


class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(1,1,1)
        super(MplCanvas, self).__init__(self.fig)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
        self.updateGeometry()

 

class MplWidget(QtWidgets.QWidget):
    """Widget defined in Qt Designer"""
    def __init__(self, parent = None, navi=True):
        QtWidgets.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        if navi:
            self.navigation_toolbar = NavigationToolbar2(self.canvas, self)

        """ connect click event """

        self.vbl = QtWidgets.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        if navi:
            self.vbl.addWidget(self.navigation_toolbar,0)
        self.setLayout(self.vbl)

        ####
        self.setGeometry(QtCore.QRect(10, 20, 571, 561))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setObjectName(_fromUtf8("mpl"))



class MplDialog(QtWidgets.QDialog):
    def __init__(self,size=(450,450),navi=True,parent=None):
        QtWidgets.QDialog.__init__(self,parent)

        self.setObjectName(_fromUtf8("MPLDialog"))
        self.resize(size[0],size[1])
        self.mpl = MplWidget(self,navi)
        self.mpl.setGeometry(QtCore.QRect(0, 0, size[0],size[1]))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mpl.sizePolicy().hasHeightForWidth())
        self.mpl.setSizePolicy(sizePolicy)
        self.mpl.setObjectName(_fromUtf8("mpl"))
        self.setWindowTitle(_translate("MPLDialog", "canvas", None))
        QtCore.QMetaObject.connectSlotsByName(self)
        


