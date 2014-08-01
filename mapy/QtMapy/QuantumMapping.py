#from __future__ import with_statement
# -*- coding:utf-8 -*-
from .Design.energydomain_dialog import Ui_EnergyDomainDialog
from .Design.mplwidget import MplDialog
from .Mapping import Mapping
from PyQt4 import QtCore, QtGui
import numpy as np
import mapy
import sympy
import inspect
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from .Design.mplwidget import MplWidget
nullfmt   = NullFormatter()   
                
class QuantumUI(QtGui.QDialog):
    def __init__(self, mainwindow):
        QtGui.QDialog.__init__(self, parent=mainwindow)
        self.mainwindow=mainwindow        

        
        
    def makeAxes(self):
        #self.mainwindow.creatCanvas()

        self.mainwindow.mpl.canvas.fig.clf()
        ax1 = self.mainwindow.mpl.canvas.fig.add_axes([0.1, 0.1, 0.5, 0.5])
        ax2 = self.mainwindow.mpl.canvas.fig.add_axes([0.1, 0.61, 0.5, 0.3])
        ax3 = self.mainwindow.mpl.canvas.fig.add_axes([0.61, 0.1, 0.3, 0.5])
        self.ax = [ax1, ax2, ax3]
        self.ax[1].xaxis.set_major_formatter(nullfmt)
        self.ax[2].yaxis.set_major_formatter(nullfmt)        
        self.mainwindow.mpl.canvas.draw()

        
        
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


class EnergyDomainUI(QuantumUI, Ui_EnergyDomainDialog):
    def __init__(self, mainwindow):
        QuantumUI.__init__(self, mainwindow)
        self.setupUi(self)
        self.mainwindow.creatCanvas()
        self.makeAxes()
        self.trajectory = []

        qr = self.mainwindow.lineEdit_qrange.text()
        pr = self.mainwindow.lineEdit_qrange.text()

        self.submpl = MplDialog(parent=self.mainwindow)

        self.lineEdit_qdomain.setText(qr)
        self.lineEdit_pdomain.setText(pr)
        self.lineEdit_hsmqrange.setText(self.lineEdit_qdomain.text())
        self.lineEdit_hsmprange.setText(self.lineEdit_pdomain.text())                

        self.creatActions()
        #self.get_parameters()

    def creatActions(self):
        self.GetEigenButton.clicked.connect(self.push_eigen_button)        
        self.drawButton.clicked.connect(self.push_draw_button)
        self.orderingButton.clicked.connect(self.ordering)
        self.eigenvaluesButton.clicked.connect(self.submpl.show)
        self.mainwindow.canvasClearButton.clicked.connect(self.actCanvasClear)
        self.mainwindow.mpl.canvas.mpl_connect('button_press_event', self.actCanvasPress)
        #self.mainwindow.actionEnergy_contour.triggered.connect(self.contour_plot)
                
    def get_parameters(self):
        self.dim = int(self.lineEdit_dim.text())        
        self.qr = self.decode_str2number(self.lineEdit_qdomain.text())
        self.pr = self.decode_str2number(self.lineEdit_pdomain.text())
        self.lineEdit_hsmqrange.setText(self.lineEdit_qdomain.text())
        self.lineEdit_hsmprange.setText(self.lineEdit_pdomain.text())    
        funcT, funcV = self.mainwindow.mapsys.get_Hamiltonian()
        self.func_T = funcT if len(inspect.getargspec(funcT)[0]) == 2 else lambda q,p: funcT(p)
        self.func_V = funcV if len(inspect.getargspec(funcV)[0]) == 2 else lambda q,p: funcV(q)
    
    def push_eigen_button(self):
        self.lineEdit_QuantumNumber.setText("0")
        self.get_parameters()        
        self.SelectAlgorism()

        self.evals, self.evecs = self.solve()        

        self.groupBox_Hsm.setEnabled(True)
        self.checkBox_withcontour.setEnabled(True)        
        self.groupBox_ordering.setEnabled(True)
        self.lineEdit_ordering.setEnabled(True)
                        
        self.push_draw_button()
    
                
    def push_draw_button(self, n=None):
        #if n==None:
        self.domain = [self.qr, self.pr]
        self.vqr = self.decode_str2number(self.lineEdit_hsmqrange.text())
        self.vpr = self.decode_str2number(self.lineEdit_hsmprange.text())
        self.hsmrange= [self.vqr, self.vpr]
        self.grid = self.decode_str2number(self.lineEdit_hsmgrid.text(), isint=True)
        
        n = int(self.lineEdit_QuantumNumber.text())
        self.get_parameters()
        self.hsmplot(n)
        self.plot_energy(n)
    
        if n == (self.dim-1):
            self.lineEdit_QuantumNumber.setText("%d" % 0)
        else:
            self.lineEdit_QuantumNumber.setText("%d" % (n+1))
        

    def solve(self):
        pass
    
    def solveQmapEigen(self):
        func_T, func_V = self.mainwindow.mapsys.get_Hamiltonian()
        dt = self.decode_str2number(self.mainwindow.lineEdit_SI_dt.text())[0]
        order = int(self.mainwindow.lineEdit_SI_order.text())
        qmap = self.Quantum.Qmap(func_T, func_V, dt, order)
        mat = qmap.UnitaryMatrix()
        evals, evecs = self.Quantum.eigen(mat, sort=False)
        
        #if self.checkBox_ordering.isChecked():
        #    if "integ." in self.comboBox_ordering.currentText().split(" "):
        #        evals, evecs = self.sortEigenOverlap(evals, evecs)
        return evals, evecs                

    def solveHamiltonPosition(self):
        funcT, funcV = self.mainwindow.mapsys.get_Hamiltonian()

        func_T = funcT if len(inspect.getargspec(funcT)[0]) == 2 else lambda q,p: funcT(p)
        func_V = funcV if len(inspect.getargspec(funcV)[0]) == 2 else lambda q,p: funcV(q)
    
        HSolver = self.Quantum.HamiltonSolver(base='position')
        if isinstance(self.mainwindow.mapsys, mapy.Systems.NormannNormal):
            matCosQ = HSolver.matrix_function(lambda x: mapy.cos(x), 'q')
            matP = HSolver.matrix_function(lambda x:x,'p')
            matHam = func_T(matCosQ, matP) + func_V(matCosQ, matP)
        
        elif isinstance(self.mainwindow.mapsys, mapy.Systems.Polynomial):
            matQ = HSolver.matrix_function(lambda x:x,'q')
            matP = HSolver.matrix_function(lambda x:x,'p')
            matHam = func_T(matQ, matP) + func_V(matQ, matP)
        else:
            matV = HSolver.matrix_function(funcV, 'q')            
            matT = HSolver.matrix_function(funcT, 'p')
            matHam = matT + matV
    
        return self.Quantum.eigen(matHam, sort=True)

    def solveHamiltonHarmonic(self):
        funcT, funcV = self.mainwindow.mapsys.get_Hamiltonian()        
        func_T = funcT if len(inspect.getargspec(funcT)[0]) == 2 else lambda q,p: funcT(p)
        func_V = funcV if len(inspect.getargspec(funcV)[0]) == 2 else lambda q,p: funcV(q)
        HSolver = self.Quantum.HamiltonSolver(base='harmonic')
        matQ = HSolver.matrixQ()
        matP = HSolver.matrixP()
        matHam = func_T(matQ, matP) + func_V(matQ, matP)

        evals, coeff = self.Quantum.eigen(matHam, sort=True)
        evecs = [HSolver.coeff2qrep(HSolver.x[0],coeff[n]) for n in range(self.dim)]
        return evals, evecs
    
        
    
    def hsmplot(self,n=0):
        self.makeAxes()
        self.get_parameters()
        #[ ax.cla() for ax in self.ax]
        vec = self.evecs[n]
        eval = self.evals[n]
        try:
            ene, gamma = self.qmap.quasienergy()
        except AttributeError:
            ene, gamma = np.nan, np.nan
        q,p = vec.x
        
        title = "n=%d,\nRe($u_n$)=%f,\nIm($u_n$)=%f,\n $\\varepsilon_n$=%f,\n $\\gamma_n$=%f" % (n, eval.real, eval.imag, ene, gamma)
        self.mainwindow.mpl.canvas.fig.suptitle(title,position=(0.8,0.9))


        try:
            x,y,z = vec.hsmrep(grid=self.grid, vrange=self.hsmrange)
            levels = np.linspace(0,z.max(),100)
          
            
            self.ax[0].contourf(x,y,z,levels,cmap=mapy.hsm_cmap)        
            self.ax[0].set_xlim(self.vqr[0],self.vqr[1])
            self.ax[0].set_ylim(self.vpr[0],self.vpr[1])
        except OSError:
            pass
        
        if self.checkBox_withcontour.isChecked():
            #x = np.linspace(self.vqr[0], self.vqr[1], 1000)
            #y = np.linspace(self.vqr[0], self.vqr[1], 1000)
            #x,y = np.meshgrid(x,y)            
            H = self.func_T(x,y) + self.func_V(x,y)
            variables = self.decode_str2number(self.lineEdit_contourlinspace.text())
            levels= int(variables[0]) if len(variables) ==1 else np.linspace(variables[0], variables[1], int(variables[2]))
            #levels = np.linspace(H.min(), H.max(), 100)                
            self.ax[0].contour(x,y,H,levels,colors='k') #,colors='k',vmin=variables[0],vmax=variables[1])          


        qvec = vec.abs2()
        pvec = vec.q2p().abs2()
        self.ax[1].plot(q, np.log10(qvec), '-k', lw=3)
        self.ax[2].plot(np.log10(pvec),p, '-k', lw=3)
        self.ax[1].set_xlim(self.vqr[0],self.vqr[1])        
        self.ax[1].set_ylim(-32,0)
        self.ax[2].set_xlim(-32,0)                                
        self.ax[2].set_ylim(self.vpr[0],self.vpr[1])
        tic = np.array(self.ax[2].get_xticks(),dtype=np.int)
        self.ax[2].set_xticklabels(tic, rotation=-90)

        if len(self.trajectory) != 0:
            [self.ax[0].plot(x[0], x[1],',k') for x in self.trajectory]
                                
        self.mainwindow.mpl.canvas.draw()


    def plot_energy(self,n):
        self.submpl.mpl.canvas.ax.clear()
        if self.solve_mode == "Qmap":
            theta = np.linspace(0, 2*np.pi,100)
            z = np.exp(-1.j*theta)
            self.submpl.mpl.canvas.ax.plot(z.real, z.imag,'-k')
            self.submpl.mpl.canvas.ax.plot(self.evals.real, self.evals.imag, 'ob',markersize=5)
            self.submpl.mpl.canvas.ax.plot(self.evals.real[n], self.evals.imag[n], 'ro',markersize=8)
            self.submpl.mpl.canvas.ax.set_xlim(-1.1,1.1)
            self.submpl.mpl.canvas.ax.set_ylim(-1.1,1.1)                        
        else:
            x = np.arange(self.dim)
            self.submpl.mpl.canvas.ax.plot(x, self.evals.real, 'ob',markersize=5)
            self.submpl.mpl.canvas.ax.plot(x[n],self.evals.real[n], 'ro',markersize=8)

        self.submpl.mpl.canvas.draw()
        

    def SelectAlgorism(self):
        self.Quantum = mapy.Quantum(self.dim, [self.qr, self.pr])
#        self.qmap = mapy.Qmap(self.dim, [self.qr, self.pr], dt=dt, order=order)        
        algo = self.comboBox_SelectAlgorism.currentText()
        if "(SplitOperator)" in algo.split(" "):
            self.solve = self.solveQmapEigen
            self.solve_mode = "Qmap"
        elif "(PositionBase)" in algo.split(" "):
            self.solve = self.solveHamiltonPosition
            self.solve_mode = "PosHam"
        elif "(HarmonicBase)" in algo.split(" "):
            self.solve = self.solveHamiltonHarmonic
            self.solve_mode = "HarHam"
        else:
            print("なんかえらー")        
            
    def SetSystem(self):
        func_T, func_V = self.mainwindow.mapsys.get_Hamiltonian()
        qmap = self.QuantumSystem(dim, domain, func_T, func_V, dt, order)         

    def actCanvasPress(self, event):
        if not self.mainwindow.checkBoxClickable.isChecked():
            return
        try:
            print('(x,y) = (%f, %f) :button=%d, '% (event.xdata, event.ydata, event.button))
            self.mainwindow.lineEdit_q.setText("%f" % event.xdata)
            self.mainwindow.lineEdit_p.setText("%f" % event.ydata)
            init = [np.array([event.xdata]), np.array([event.ydata])]
            
        except TypeError:
            init = None
            return 

        iteration = int(self.mainwindow.itaration_lineEdit.text())
        dt = float(self.mainwindow.lineEdit_SI_dt.text())
        order = int(self.mainwindow.lineEdit_SI_order.text())                        
        mapping = Mapping(self.mainwindow.mapsys, dt, order)
        x = mapping.iterate(init, iteration)
        self.trajectory.append(x)
        
        self.ax[0].plot(x[0],x[1],',k')
        self.mainwindow.mpl.canvas.draw()                       
    
    def actCanvasClear(self):
        #self.mainwindow.creatCanvas()
        self.makeAxes()
        #self.creatCanvas()
        self.trajectory = []

    def ordering(self):
        if "integ." in self.comboBox_ordering.currentText().split(" "):
            evals0, evecs0 = self.solveHamiltonPosition()
            index = self.Quantum.overlap_sort_index(self.evecs, evecs0)
            self.evals, self.evecs =  self.Quantum.sorteigen(self.evals, self.evecs, index)
        elif "q-vari" in self.comboBox_ordering.currentText().split(" "):
            x0 = self.decode_str2number(self.lineEdit_ordering.text())[0]
            index = self.Quantum.variance_sort_index(self.evecs, x0=x0,qp='q')
            self.evals, self.evecs = self.Quantum.sorteigen(self.evals, self.evecs, index)
        elif "p-vari" in self.comboBox_ordering.currentText().split(" "):
            x0 = self.decode_str2number(self.lineEdit_ordering.text())[0]
            index = self.Quantum.variance_sort_index(self.evecs, x0=x0,qp='p')
            self.evals, self.evecs = self.Quantum.sorteigen(self.evals, self.evecs, index)
        else:
            if np.all(self.evals.imag == 0.0):
                index = [i[0] for i in sorted(enumerate(self.evals), key=lambda x:x[1])]
            else:
                index = [i[0] for i in sorted(enumerate(np.angle(self.evals)), key=lambda x:x[1])]
            self.evals, self.evecs = self.Quantum.sorteigen(self.evals, self.evecs, index)
        self.lineEdit_QuantumNumber.setText("0")            
        self.push_draw_button()

