from matplotlib.figure import Figure
from .Design.mplwidget import MplWidget
import mapy
import numpy as np
twopi=2*np.pi

class MappingUI(object):
    def __init__(self, mainwindow):
        self.mw=mainwindow

        self.mw.creatCanvas()
        self.createAction()


        
    def createAction(self):
        self.mw.plotButton.clicked.connect(self.actPlot)
        self.mw.mpl.canvas.mpl_connect('button_press_event', self.actCanvasPress)
        self.mw.mpl.canvas.mpl_connect('button_release_event', self.actCanvasRelease)        
#        self.mw.mpl.canvas.mpl_connect('mouse', self.actCanvasPress)        
#        self.mw.mpl.canvas.mpl_connect('mouse', self.actCanvasRelease)   
          
    def actCanvasRelease(self, event):
        if not self.mw.checkBoxClickable.isChecked():
            return
        elif self.mw.comboBox_Select_Init.currentText() in ['linear', 'box']:
            q0 = self.mw.decode_str2number(self.mw.lineEdit_q.text())[0]
            p0 = self.mw.decode_str2number(self.mw.lineEdit_p.text())[0]
            q1,p1 = event.xdata, event.ydata
            self.mw.lineEdit_q.setText("%f,%f" % (q0,q1))
            self.mw.lineEdit_p.setText("%f,%f" % (p0,p1))

    def actCanvasPress(self, event):
        if not self.mw.checkBoxClickable.isChecked():
            return
        try:
            print('(x,y) = (%f, %f) :button=%d, '% (event.xdata, event.ydata, event.button))
            self.mw.lineEdit_q.setText("%f" % event.xdata)
            self.mw.lineEdit_p.setText("%f" % event.ydata)
            init = [np.array([event.xdata]), np.array([event.ydata])]
            
        except TypeError:
            init = None
            return 

        if self.mw.comboBox_Select_Init.currentText() == 'point':  
            iteration = int(self.mw.itaration_lineEdit.text())
            dt = float(self.mw.lineEdit_SI_dt.text())
            order = int(self.mw.lineEdit_SI_order.text())                        
            mapping = Mapping(self.mw.mapsys, dt, order)
            x = mapping.iterate(init, iteration)
            self.mw.mpl.canvas.ax.plot(x[0],x[1],',')
            self.mw.mpl.canvas.draw()                       

    def actPlot(self):
        dt = self.mw.decode_str2number(self.mw.lineEdit_SI_dt.text())
        order = int(self.mw.lineEdit_SI_order.text())
        mapping = Mapping(self.mw.mapsys, dt, order)
        iteration = int(self.mw.itaration_lineEdit.text())
        samplenum = int(self.mw.samplenum_lineEdit.text())         
                   
        if self.mw.comboBox_Select_Init.currentText() == 'point':
            q0 = self.mw.decode_str2number(self.mw.lineEdit_q.text())
            p0 = self.mw.decode_str2number(self.mw.lineEdit_p.text())
            init = [np.array(q0), np.array(p0)]
            x = mapping.iterate(init, iteration)

        elif self.mw.comboBox_Select_Init.currentText() == 'random':        
            x = mapping.random_iterate(iteration, samplenum)

        elif self.mw.comboBox_Select_Init.currentText() == 'linear':
            xr = self.mw.decode_str2number(self.mw.lineEdit_q.text())
            yr = self.mw.decode_str2number(self.mw.lineEdit_p.text())
            if np.abs(xr[1] - xr[0])<=1e-16:
                x = np.array([0.0]*samplenum)
                y = np.linspace(yr[0], yr[1], samplenum)
            else:
                a = (yr[1] - yr[0])/(xr[1] - xr[0])
                x = np.linspace(xr[0],xr[1], samplenum)
                y = a*(x - xr[0]) + yr[0]
            init = [x,y]
            x = mapping.iterate(init, iteration)
        else:
            x = [None, None]
                        
        self.mw.mpl.canvas.ax.plot(x[0],x[1],',')
        self.mw.mpl.canvas.draw()            


class Mapping(object):
    def __init__(self, mapsys, dt=1, order=1):
        self.mapsys = mapsys
        self.dt = dt
        self.order = order
        self.init = []

    def iterate(self, init,iteration):
        funcT, funcV = self.mapsys.get_DHamiltonian()                        
        SI = mapy.SymplecticIntegrator(funcT,funcV,dt=self.dt, order=self.order)

        self.init.append(init)
        x = init[:]
        traj = init[:]

        for i in range(iteration):
            for j in range(2):
                if self.mapsys.periodicity[j]:
                    x[j] = x[j] - np.floor((x[j]-np.pi)/twopi)*twopi - twopi             
                traj[j] = np.append(traj[j], x[j])                
            x = SI.evolve(x)            
        return traj
        
        
    def random_iterate(self, iteration, samplenum):
        funcT, funcV = self.mapsys.get_DHamiltonian()                        
        SI = mapy.SymplecticIntegrator(funcT,funcV,dt=self.dt, order=self.order)

        q = (np.random.random(samplenum)-0.5)*twopi
        p = (np.random.random(samplenum)-0.5)*twopi
        x = [q,p]
        traj = [np.array([]), np.array([])]        

        for i in range(iteration):
            for j in range(2):
                if self.mapsys.periodicity[j]:
                    x[j] = x[j] - np.floor((x[j]-np.pi)/twopi)*twopi - twopi             
                traj[j] = np.append(traj[j], x[j])                
            x = SI.evolve(x)            
        return traj
        
                