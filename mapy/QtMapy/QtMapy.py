
from __future__ import with_statement
from PyQt4 import QtCore, QtGui
import numpy as np
import sympy
import sys
import mapy

import inspect

version="0.2"
twopi = 2*np.pi

# import the MainWindow widget from the converted .ui files
from .Design.maindesign import Ui_MainWindow
from .Design.popupAbout import Ui_AboutDialog
from .Design.mplwidget import MplWidget
#from .Design.energydomain_dialog import Ui_EnergyDomainDialog
from . import QuantumMapping
from . import Mapping

MapNames=['Standard', 'Harper', 'Harmonic','Henon','HaradaNormal','JeremyNormal', 'NormannNormal']
#ClassicalAnalysisNames=['Mapping','TimeEvolve']
def decode_str2number(text, isint=False):
    strs = text.split(",")
    xstrs = [x.strip() for x in strs]
    x = []
    for y in xstrs:
        if isint:
            x.append(int(sympy.N(y)))
        else:
            try:
                x.append(float(sympy.N(y)))                
            except TypeError:
                x.append(complex(sympy.N(y)))                
    return x[0] if len(x) == 1 else x
        

class SetupDialog(QtGui.QDialog):
    def __init__(self, ui, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = ui()
        self.ui.setupUi(self)

class AboutDialog(QtGui.QDialog, Ui_AboutDialog):
    def __init__(self,parent=None):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)

class ParameterDialog(QtGui.QDialog):
    def __init__(self, mapsys, parent=None):
        QtGui.QDialog.__init__(self,parent)

        grid = QtGui.QGridLayout()
        args = inspect.getargspec(mapsys.__init__)[0][1:]
        parameters = [getattr(mapsys, p) for p in args]

        self.LineEdits=[]
        for i in range(len(args)):
            label = QtGui.QLabel("%s:" % args[i])
            le = QtGui.QLineEdit(self)
            le.setText("%f" % parameters[i])
            self.LineEdits.append(le)          
            
            grid.addWidget(label, i, 0)
            grid.addWidget(le, i, 1)
            
        label = QtGui.QLabel("Priodicity:")
        grid.addWidget(label,i+1,0)
        i = i + 2
        qp = ["q:","p:"]
        self.Combs=[]
        for j in range(2):
            i = i + j
            label = QtGui.QLabel("%s" % qp[j])
            grid.addWidget(label, i, 0)
                        
            combo = QtGui.QComboBox()
            combo.addItems(['True','False'])
            if mapsys.periodicity[j]:
                combo.setCurrentIndex(0)
            else:
                combo.setCurrentIndex(1)                
            self.Combs.append(combo)
            grid.addWidget(combo, i, 1)

        
        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        #self.buttonBox.clicked.connect(self.actClose)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.close)                

        grid.addWidget(self.buttonBox, i+1,1)
        self.setLayout(grid)
        self.setWindowTitle('Parameter Setting')
        
    def accept(self):
        self.paras = [decode_str2number(x.text()) for x in self.LineEdits]
        self.periodicity = [True if x.currentText() == 'True' else False for x in self.Combs]
        self.close()
        


class DesignerMainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None):
        super(DesignerMainWindow, self).__init__(parent)

        # setup the GUI --> function generated by pyuic4
        self.setupUi(self)
        #self.mpl.deleteLater()
        self.mpl = MplWidget(self.canvasBox)
        #self.mpl.show()
        self.creatActions()        
        #self.mpl.

        self.mapsys = mapy.Systems.Standard()
        self.mapsys = getattr(mapy.Systems, "Standard")()        
        
        self.SelectMapping()
        self.analysis = Mapping.MappingUI(self)

    def creatActions(self):
        #signal slot
        #self.plotButton.clicked.connect(self.actPlot)
        self.canvasClearButton.clicked.connect(self.actCanvasClear)
        self.ReplotButton.clicked.connect(self.actReplot)

        self.mpl.canvas.mpl_connect('button_press_event', self.actCanvasPress)
        self.mpl.canvas.mpl_connect('button_release_event', self.actCanvasRelease)        
        self.mpl.canvas.mpl_connect('mouse', self.actCanvasPress)        
        self.mpl.canvas.mpl_connect('mouse', self.actCanvasRelease)   
                  
#        self.actionAboutQtMapy.triggered.connect(self.on_AboutInfomation)
        """ menubar action setting """
        self.dialogTextBrowser = AboutDialog()        
        self.actionAboutQtMapy.triggered.connect(lambda :self.dialogTextBrowser.exec_())
        
        self.ParameterButton.clicked.connect(self.openParameterSettingDialog)        
        
        
        # grouping of definition of mapping 
        self.MappingGroup = QtGui.QActionGroup(self)
        self.MappingGroup.triggered.connect(self.SelectMapping)
        for name in MapNames:
            act = "action" + name + "_map"
            setattr(self, act,  QtGui.QAction(self))
            action = getattr(self, act)
            action.setObjectName(act)
            action.setText(name)
            action.setCheckable(True)
            self.menuMapping.addAction(action)
            
            
            self.MappingGroup.addAction(getattr(self, act))
        self.actionStandard_map.setChecked(True)            

        self.actionMapping.triggered.connect(lambda : self.SelectClassicalAnalysis("Mapping"))
        
        self.actionEnergyDomain.triggered.connect(lambda : self.SelectQuantumAnalysis("EnergyDomain"))
        self.actionEnergy_contour.triggered.connect(self.show_energyContour)
    
    def actCanvasPress(self,event):
        print("H")
    
    def actCanvasRelease(self,event):
        pass    

    def creatCanvas(self):
        try:
            self.mpl.deleteLater()
        except RuntimeError:
            pass

        self.mpl = MplWidget(self.canvasBox)
        self.mpl.show()

    def SelectMapping(self):
        mapping = self.MappingGroup.checkedAction()
        sysname = mapping.text().split(" ")[0]
        self.mapsys = getattr(mapy.Systems, sysname)()
        self.show_system_def(sysname, self.label_SystemInfo)
        self.show_system_parameter(self.label_Parameters)
        self.parameterDialog = ParameterDialog(self.mapsys)                          


    def SelectClassicalAnalysis(self, subject):
        if subject == "Mapping":
            self.analysis = Mapping.MappingUI(self)
        
    def SelectQuantumAnalysis(self, subject):
        if subject == "EnergyDomain":
            self.analysis = QuantumMapping.EnergyDomainUI(self)        
        self.analysis.show()



    """ show system infomation """
    def openParameterSettingDialog(self):
        self.parameterDialog.exec_()
        paras = self.parameterDialog.paras
        periodicity = self.parameterDialog.periodicity
        mapping = self.MappingGroup.checkedAction()        
        sysname = mapping.text().split(" ")[0]
        f = getattr(mapy.Systems, sysname)
        self.mapsys = f(*paras) # apply(f, paras)
        self.mapsys.periodicity = periodicity
        self.show_system_parameter(self.label_Parameters)
            
    def show_system_def(self, sysname,label):
        q,p = sympy.Symbol("q"),sympy.Symbol("p")
        args = inspect.getargspec(self.mapsys.__init__)[0][1:]
        paras = [sympy.Symbol(p) for p in args]
        f = getattr(mapy.Systems, sysname)
        mapsys = f(*paras)
        T, V = mapsys.get_Hamiltonian()
        funcT = T if len(inspect.getargspec(T)[0]) == 2 else lambda q,p: T(p)
        funcV = V if len(inspect.getargspec(V)[0]) == 2 else lambda q,p: V(q)        
        kinetic_info = "T(q,p) = %s" % funcT(q,p)
        potential_info = "V(q,p) = %s" % funcV(q,p)
        sysinfo = kinetic_info + "\n" + potential_info
        label.setText("%s" % sysinfo)

    def show_system_parameter(self, label):
        args = inspect.getargspec(self.mapsys.__init__)[0][1:]
        paras = [getattr(self.mapsys, p) for p in args]
        parainfo = ""
        for i in range(len(args)):
            parainfo += "%s = %f\n"  % (args[i], paras[i])
        label.setText("%s" % parainfo)
        self.show_system_periodicity(self.label_Periodicity)                

    def show_system_periodicity(self,label):
        try:
            p = self.mapsys.periodicity
        except AttributeError:
            p = [False, False]
            
        text = "q: %s\t p:%s" % (p[0],p[1])
        label.setText("%s"% text)
        

    """ mpl setting """
    def set_mpl_canvas_plot_range(self):
        x = self.decode_str2number(self.lineEdit_qrange.text())
        y = self.decode_str2number(self.lineEdit_prange.text())        
        self.mpl.canvas.ax.set_xlim(x[0],x[1])
        self.mpl.canvas.ax.set_ylim(y[0],y[1])
        return [x,y]
        
    def actReplot(self):
        self.set_mpl_canvas_plot_range()    
        self.mpl.canvas.draw()

    def actCanvasClear(self):
        #self.mpl.canvas.fig.clear()
        self.mpl.canvas.ax.clear()
        self.mpl.canvas.draw()

    """
    def plot_energy_contour(self):
        xr = self.set_mpl_canvas_plot_range()
        q = np.linspace(xr[0],xr[1])
        p = np.linspace(xr[2],xr[3])
        q, p = np.meshgrid(q,p)
        T,V = self.mapsys.get_Hamiltonian()
        H = T(p) + V(q)
        levels = np.linspace(H.min(), H.max(), 30)
        self.mpl.canvas.ax.contour(q,p,H,levels)
        self.mpl.canvas.draw()
    """ 

    """ others """
    def decode_str2number(self, text, isint=False):
        strs = text.split(",")
        xstrs = [x.strip() for x in strs]
        x = []
        for y in xstrs:
            if isint:
                x.append(int(sympy.N(y)))
            else:
                try:
                    x.append(float(sympy.N(y)))                
                except TypeError:
                    x.append(complex(sympy.N(y)))                
        return x
        

    def show_energyContour(self):
        T, V = self.mapsys.get_Hamiltonian()
        funcT = T if len(inspect.getargspec(T)[0]) == 2 else lambda q,p: T(p)
        funcV = V if len(inspect.getargspec(V)[0]) == 2 else lambda q,p: V(q)
        xr = self.decode_str2number(self.lineEdit_qrange.text())
        yr = self.decode_str2number(self.lineEdit_prange.text())
        x = np.linspace(xr[0],xr[1],100)
        y = np.linspace(yr[0],yr[1],100)
        x,y = np.meshgrid(x,y)
        H = funcT(x,y) + funcV(x,y)
        self.mpl.canvas.ax.contour(x,y,H,100)
        self.mpl.canvas.draw()
        
        
        
            
def run():
    app = QtGui.QApplication(sys.argv)
    dmw = DesignerMainWindow()
    dmw.show()
    sys.exit(app.exec_())

