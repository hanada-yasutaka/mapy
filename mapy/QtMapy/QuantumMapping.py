#from __future__ import with_statement
# -*- coding:utf-8 -*-
from .Design.energydomain_dialog import Ui_EnergyDomainDialog
from .Design.mplwidget import MplDialog
from .Mapping import Mapping
from PyQt4 import QtCore, QtGui
import numpy as np
import mapy
import sympy
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
        self.lineEdit_hsmqrange.setText(qr)        
        self.lineEdit_hsmprange.setText(pr)

        self.creatActions()
        self.get_parameters()

    def creatActions(self):
        self.GetEigenButton.clicked.connect(self.push_eigen_button)        
        self.drawButton.clicked.connect(self.push_draw_button)
        self.eigenvaluesButton.clicked.connect(self.submpl.show)
        self.mainwindow.canvasClearButton.clicked.connect(self.actCanvasClear)
        self.mainwindow.mpl.canvas.mpl_connect('button_press_event', self.actCanvasPress)
                
    def get_parameters(self):
        self.dim = int(self.lineEdit_dim.text())        
        self.qr = self.decode_str2number(self.lineEdit_qdomain.text())
        self.pr = self.decode_str2number(self.lineEdit_pdomain.text())
        self.domain = [self.qr, self.pr]
        self.vqr = self.decode_str2number(self.lineEdit_hsmqrange.text())
        self.vpr = self.decode_str2number(self.lineEdit_hsmprange.text())
        self.hsmrange= [self.vqr, self.vpr]
        self.grid = self.decode_str2number(self.lineEdit_hsmgrid.text(), isint=True)
                
    def push_draw_button(self):
        n = int(self.lineEdit_QuantumNumber.text())
        self.get_parameters()
        self.hsmplot(n)
        self.plot_energy(n)
    
        if n == (self.dim-1):
            self.lineEdit_QuantumNumber.setText("%d" % 0)
        else:
            self.lineEdit_QuantumNumber.setText("%d" % (n+1))
        
    def push_eigen_button(self):
        self.lineEdit_QuantumNumber.setText("0")
        self.SelectAlgorism()
        
        self.evals, self.evecs = self.solve()

        self.groupBox_2.setEnabled(True)
        self.groupBox_3.setEnabled(True)        
        self.push_draw_button()
    
    def solve(self):
        pass
    
    def solveQmapEigen(self):
        func_T, func_V = self.mainwindow.mapsys.get_Hamiltonian()
        self.qmap.setFunctions(func_T, func_V)
        mat = self.qmap.UnitaryMatrix()
        evals, evecs = self.qmap.eigen(mat, sort=False)
        
        if self.checkBox_ordering.isChecked():
            if "integ." in self.comboBox_ordering.currentText().split(" "):
                evals, evecs = self.sortEigenOverlap(evals, evecs)
        return evals, evecs                

    def solveHamiltonPosition(self):
        func_T, func_V = self.mainwindow.mapsys.get_Hamiltonian()
        matHam = self.qmap.HamiltonMatrix(func_T, func_V)
        evals, evecs = self.qmap.eigen(matHam, sort=True)
        return evals, evecs
    
    def sortEigenOverlap(self, evals, evecs):
        evals0, evecs0 = self.solveHamiltonPosition()
        index = self.qmap.overlap_index(evecs, evecs0)
        return self.qmap.sorteigen(evals, evecs, index)                
        
    
    def hsmplot(self,n=0):
        self.makeAxes()
        #[ ax.cla() for ax in self.ax]
        vec = self.evecs[n]
        eval = self.evals[n]
        ene, gamma = self.qmap.quasienergy()
        q,p = self.qmap.x
        
        title = "n=%d,\nRe($u_n$)=%f,\nIm($u_n$)=%f,\n $\\varepsilon_n$=%f,\n $\\gamma_n$=%f" % (n, eval.real, eval.imag, ene, gamma)
        self.mainwindow.mpl.canvas.fig.suptitle(title,position=(0.8,0.9))
#        try(aaa):
        x,y,z = self.qmap.hsmrep(vec,grid=self.grid, vrange=self.hsmrange)
        levels = np.linspace(0,z.max(),50)

        self.ax[0].contourf(x,y,z,levels,cmap=mapy.hsm_cmap)        
        self.ax[0].set_xlim(self.vqr)
        self.ax[0].set_ylim(self.vpr)
#        except OSError:
#            pass

        qvec = self.qmap.abs2(vec)
        pvec = self.qmap.abs2(self.qmap.q2p(vec))
        self.ax[1].plot(q, np.log10(qvec), '-k', lw=3)
        self.ax[2].plot(np.log10(pvec),p, '-k', lw=3)
        self.ax[1].set_xlim(self.vqr)        
        self.ax[1].set_ylim(-25,0)
        self.ax[2].set_xlim(-25,0)                                
        self.ax[2].set_ylim(self.vpr)
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
        elif self.solve_mode == "Ham":
            x = np.arange(self.dim)
            self.submpl.mpl.canvas.ax.plot(x, self.evals.real, 'ob',markersize=5)
            self.submpl.mpl.canvas.ax.plot(x[n],self.evals.real[n], 'ro',markersize=8)

        self.submpl.mpl.canvas.draw()
        

    def SelectAlgorism(self):
        algo = self.comboBox_SelectAlgorism.currentText()
        dt = self.decode_str2number(self.mainwindow.lineEdit_SI_dt.text())[0]
        order = int(self.mainwindow.lineEdit_SI_order.text())
        self.qmap = mapy.Qmap(self.dim, self.domain, dt=dt, order=order)
        if "(SplitOperator)" in algo.split(" "):
            self.solve = self.solveQmapEigen
            self.solve_mode = "Qmap"
        elif "(PositionBase)" in algo.split(" "):
            self.solve = self.solveHamiltonPosition
            self.solve_mode = "Ham"
        elif "(HarmonicBase)" in algo.split(" "):
            print("まだじゅんびちゅう")
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



