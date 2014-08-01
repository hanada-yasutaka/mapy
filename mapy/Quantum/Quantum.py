import mapy
import numpy as np
import math
import sympy
import functools
import datetime
import mpmath
twopi=2*np.pi

class HilbertSpace(object):
    def __init__(self, dim, domain,dtype='float'):
        self.dim = self._setDim(dim)
        #self.dtype = dtype
        if domain !=None:        
            self.domain = self._setDomain(domain)
            self.h = self._setPlanck()
            self.x = [mapy.Array.linspace(self.domain[0][0], self.domain[0][1], self.dim, endpoint=False, dtype=dtype),
                  mapy.Array.linspace(self.domain[1][0], self.domain[1][1], self.dim, endpoint=False, dtype=dtype)]
            self._qrep = True                
            self._Shift()

    
    def _Shift(self):
        if self._qrep:
            if self.domain[1][0]*self.domain[1][1]<0:
                self.Shift = [False, True]
            else:
                self.Shift = [False, False]
        else:
            if self.domain[0][0]*self.domain[0][1]<0:
                self.Shift = [True, False]
            else:
                self.Shift = [False, False]        
        
    def _setDim(self, dim):
        if (dim <= 0):
            raise AttributeError("negative dimensions are not allowed")
        return dim
        
    def _setDomain(self, domain):
        qr, pr  = domain[0], domain[1]

        if (qr[0] >= qr[1]):
            raise AttributeError("qmin > qmax (%f,%f)" % (qr[0],qr[1]))
        if (pr[0] >= pr[1]):
            raise AttributeError("pmin > pmax (%f,%f)" % (pr[0],pr[1]))
        return domain
        
    def _setPlanck(self):
        qr, pr = self.domain
        return (pr[1] - pr[0])*(qr[1] - qr[0])/self.dim

class State(mapy.Array, HilbertSpace):
    def __init__(self, input, domain, dtype='float'):
        HilbertSpace.__init__(self, len(input), domain, dtype=dtype)
        
    def __new__(cls, input, domain,dtype='float'):
        if dtype == 'float':
            obj = np.asarray(input, dtype=getattr(np, 'complex')).view(cls)
        else:
            obj = np.asarray(input, dtype=getattr(np, dtype)).view(cls)            
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        
    def q2p(self):
        vec = np.fft.fft(self)/np.sqrt(self.dim)
        if self.Shift[1]:
            vec = np.fft.fftshift(vec)
        return State(vec, self.domain)
             
    def p2q(self):
        if self.Shift[1]:
            vec = np.fft.fftshift(vec)        
        vec = np.fft.ifft(self)*np.sqrt(self.dim)
        return State(vec, self.domain)
    
    def coherent(self, q_c, p_c, x=None):
        if x == None:
            x = self.x[0]
        res = mapy.Array([np.exp(-(xx-q_c)**2*np.pi/self.h + 1.j*(xx - q_c)*p_c*twopi/self.h) for xx in x])
        return State(res, self.domain).normalize() #.normalize()

    def cs(self, q_c, p_c):
        qrange = self.domain[0]
        d = qrange[1] - qrange[0]
        lqmin, lqmax = qrange[0] - 2*d, qrange[1] + 2*d
        long_q = np.linspace(lqmin, lqmax, 5*self.dim, endpoint=False)
        
        coh_state = self.coherent(q_c, p_c, long_q)

        vec = np.zeros(self.dim,dtype=np.complex128) 
        m = len(coh_state)//self.dim
        coh_state = coh_state.reshape(m,self.dim)
        
        for i in range(m):
            vec += coh_state[i][::1]
        return State(vec,self.domain).normalize()
    
    def qconst(self):
        pass
    def pconst(self):
        pass

    def normalize(self):
        return self/np.sqrt(self.norm())
    
    def norm(self):
        norm = np.abs(np.sum(self*np.conj(self)))
        return norm.tolist()

    def hsmrep(self, grid=[50,50], vrange=None):
        if vrange==None:
            vrange = self.domain

        from mapy.ctypes_wrapper import wrapper
        file_path = mapy.__file__.replace("__init__.py", "shared/libhsm.so")

        try:
            cw = wrapper.call_hsm_rep(file_path)
        except NameError:
            raise(RuntimeError,"libhsm.so not found.")
            
        hsm_imag = cw.husimi_rep(self.toarray(), self.dim, self.domain, vrange, grid)

        x = np.linspace(vrange[0][0], vrange[0][1], grid[0])
        y = np.linspace(vrange[1][0], vrange[1][1], grid[1])

        X,Y = np.meshgrid(x,y)
        return X,Y,hsm_imag    
    
    def conj(self):
        return np.conj(self)
    
    def abs2(self):
        vec = np.abs(self*self.conj())
        return vec
    
    def inner(self, vec):
        return np.sum(np.conj(vec)*self)    

class PositionBase(HilbertSpace):
    def __init__(self, dim, domain, dtype='float'):
        HilbertSpace.__init__(self, dim, domain, dtype=dtype)

    def _setDim(self, dim):
        if (dim <= 0) or (dim %2 !=0):
            raise AttributeError("dimension must be positive and even integer")
        return dim
    
                
    def matrix_function(self, f, qp):
        if qp == 'q':
            mat = self.matrix_funcV(f)
        elif qp == 'p':
            mat = self.matrix_funcT(f)
        else:
            raise TypeError("qp expects 'q' or 'p'")
        return np.matrix(mat)
        
    
    def HamiltonMatrix(self, funcT, funcV):
        return [self.matrix_funcT(funcT), self.matrix_funcV(funcV)]
    
    def matrix_funcV(self, func):
        if self._qrep:
            mat = func(self.x[0])*np.eye(self.dim, dtype=np.complex128)
        return mat

    def matrix_funcT(self, func):
        iden = np.eye(self.dim,dtype=np.complex128)
        if self._qrep:
            if self.Shift[1]:
                mat = np.array([np.fft.ifft(np.fft.fftshift(func(self.x[1]))*pvec) for pvec in map(np.fft.fft, iden)])
                #qvecs = np.array([np.fft.ifft(pvec) for pvec in pvecs])
#                mat = qvecs.reshape(self.dim, self.dim)
                #mat = np.array([np.fft.ifft(pvec*np.fft.fftshift(func(self.x[1])))/np.sqrt(self.dim) for pvec in map(np.fft.fft, iden)])
            else:
                mat = np.array([np.fft.ifft(pvec*func(self.x[1])) for pvec in map(np.fft.fft, iden)])
        return np.matrix(mat.transpose(),dtype=np.complex128)

class HarmonicBase(HilbertSpace):
    def __init__(self, dim, domain=None,planck=None, q=None,dtype='float'):
        HilbertSpace.__init__(self, dim, domain, dtype=dtype)
        if planck !=None:
            self.h = planck
        
                    
    def matrixQ(self):
        hbar = self.h/twopi
        m = np.sqrt(np.arange(1,self.dim))
        q = np.sqrt(hbar/2.0)*(np.diag(m,-1) + np.diag(m,1))
        return np.matrix(q,dtype=np.complex128)

    def matrixP(self):
        hbar = self.h/twopi
        m = np.sqrt(np.arange(1,self.dim))
        p = np.sqrt(hbar/2.0)*(np.diag(m,-1) - np.diag(m,1))*1.j
        return np.matrix(p,dtype=np.complex128)

    def function2matrix(self, f, qp='q'):
        if qp == 'q':
            mat  = np.matrix(self.matrixQ())
        elif qp == 'p':
            mat  = np.matrix(self.matrixP())
        #res = f(mat).tolist()            
        return mat #np.matrix(res,dtype=np.complex128)
    
    def qrep_base(self,x, n):
        x = x/np.sqrt(self.h/twopi)
        coef = np.zeros(self.dim)
        coef[n] = 1.0
        Herm = np.polynomial.hermite.Hermite(coef,domain=self.domain[0])#,window=self.domain[0])
        return Herm*np.exp(-x*x/2.0) + 0.j
    
    def qrep_basis(self, x, truncate=100):
        x = x/np.sqrt(self.h/twopi)    
        eye = np.eye(truncate)

        hermit= lambda x: np.polynomial.hermite.Hermite(x,domain=self.domain[0],window=self.domain[0])
        Her_poly= map(hermit, eye)
        def pref(n):
            return 1/math.sqrt(2**n*math.factorial(n))/(np.pi*self.h/twopi)**(1/4)
        
        vecs = [ pref(n)*Hn(x)*np.exp(-x*x/2.0) + 0.j for n, Hn in enumerate(Her_poly)]
        return  [self.normalize(vec) for vec in vecs]
    
    def eigen(self, mat,sort=True):
        self.evals, self.coeff = mapy.eigen(mat, sort)
        return self.evals, self.coeff
    
    def coeff2qrep(self, x, coef, truncate=100):
        basis = self.qrep_basis(x, truncate)
        index = truncate if self.dim > truncate else self.dim
        vec = np.array([coef[n]*basis[n] for n in range(index)]).sum(axis=0)
        return State(self.normalize(vec), self.domain)
        
    def hsmrep(self, vec, grid=[50,50], vrange=None):
        if vrange==None:
            vrange = self.domain

        from mapy.ctypes_wrapper import wrapper
        file_path = mapy.__file__.replace("__init__.py", "shared/libhsm.so")

        try:
            cw = wrapper.call_hsm_rep(file_path)
        except NameError:
            raise(RuntimeError,"libhsm.so not found.")
            
        hsm_imag = cw.husimi_rep(vec, self.dim, self.domain, vrange, grid)

        x = np.linspace(vrange[0][0], vrange[0][1], grid[0])
        y = np.linspace(vrange[1][0], vrange[1][1], grid[1])

        X,Y = np.meshgrid(x,y)
        return X,Y,hsm_imag
    
    def normalize(self,x):
        norm = np.abs(np.sum(x*np.conj(x)))
        return x/np.sqrt(norm)

                                
class SplitOperator(PositionBase):
    def __init__(self, dim, domain, funcT, funcV,dtype='float'):
        PositionBase.__init__(self, dim, domain,dtype=dtype)
        self.funcT = funcT
        self.funcV = funcV
 
    def operator_expV(self, x, dt=1,isShift=False):
        if self._qrep and self.Shift[0]:
            x = np.fft.fftshift(x)
        return lambda invec: np.exp(-1.j*self.funcV(x)*dt*twopi/self.h)*invec
        
    def operator_expT(self,x, dt=1, isShift=True):
        if self._qrep and self.Shift[1]:
            x = np.fft.fftshift(x)        
        return lambda invec: np.exp(-1.j*self.funcT(x)*dt*twopi/self.h)*invec
    
    def evolve(self):
        pass

    def UnitaryMatrix(self):
        mat = np.eye(self.dim, dtype=np.complex128)
        res = np.array([x for x in map(self.evolve, mat)])
        return res.transpose()        
             
    def __annotation(self):
        ann ="DATE: %s\n" % datetime.datetime.now()
        ann += "DIM %d\n" % self.dim
        ann += "QMIN %s\n" % self.domain[0][0]
        ann += "QMAX %s\n" % self.domain[0][1]
        ann += "PMIN %s\n" % self.domain[1][0]
        ann += "PMAX %s\n" % self.domain[1][1]
        return ann

    def savetxt(self, title, x, vec):
        data = [x, abs2(vec), vec.real, vec.imag]
        np.savetxt(title, np.array(data).transpose(), header=self.__annotation())

    def loadtxt(self, title):
        pass
        

class QIntegrator(SplitOperator):
    def __init__(self, dim, domain, func_T, func_V, dt=1, order=1,dtype='float'):
        SplitOperator.__init__(self, dim, domain, func_T, func_V, dtype=dtype)
        self.dt = dt
        self.order=order
        if order == 1:        
            self.evolve = self.Symplectic1
        elif order==2:
            self.evolve = self.Symplectic2
        elif order==4:
            self.evolve = self.Symplectic4
        elif order==6:
            self.evolve = self.Symplectic6            
        elif order==8:
            self.evolve = self.Symplectic8
        elif order==10:
            self.evolve = self.Symplectic10                
        elif order==12:
            self.evolve = self.Symplectic12                            
        elif order==14:
            self.evolve = self.Symplectic14                                                                
        else:
            raise AttributeError("")        
    
    def Symplectic1(self, invec, z=1):
        c = z*np.array([1,1])
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])

        pvec = np.fft.fft(expV(invec))
        qvec = np.fft.ifft(expT(pvec))
        return qvec
        
    def Symplectic2(self, invec, z=1):
        c = z*np.array([1/2,1])
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
    
        pvec = np.fft.fft(expV(invec))
        qvec = np.fft.ifft(expT(pvec))
        pvec = np.fft.fft(expV(qvec))
        qvec = np.fft.ifft(pvec)
        return qvec
    
    def Symplectic4(self, invec, z=1):
        beta = 2**(1/3)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
        qvec =self.Symplectic2(invec, c[0])
        qvec =self.Symplectic2(qvec, c[1])
        qvec =self.Symplectic2(qvec, c[0])        
        return qvec
    
    def Symplectic6(self, invec, z=1):
        beta = 2**(1/5)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
        qvec =self.Symplectic4(invec, c[0])
        qvec =self.Symplectic4(qvec, c[1])
        qvec =self.Symplectic4(qvec, c[0])        
        return qvec
    
    def Symplectic8(self, invec, z=1):
        beta = 2**(1/7)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
        qvec =self.Symplectic6(invec, c[0])
        qvec =self.Symplectic6(qvec, c[1])
        qvec =self.Symplectic6(qvec, c[0])        
        return qvec

    def Symplectic10(self, invec, z=1):
        beta = 2**(1/9)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
        qvec =self.Symplectic8(invec, c[0])
        qvec =self.Symplectic8(qvec, c[1])
        qvec =self.Symplectic8(qvec, c[0])        
        return qvec
    
    def Symplectic12(self, invec, z=1):
        beta = 2**(1/11)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
        qvec =self.Symplectic10(invec, c[0])
        qvec =self.Symplectic10(qvec, c[1])
        qvec =self.Symplectic10(qvec, c[0])        
        return qvec    

    def Symplectic14(self, invec, z=1):
        beta = 2**(1/13)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.dt*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.dt*c[1], self.Shift[1])
        qvec =self.Symplectic12(invec, c[0])
        qvec =self.Symplectic12(qvec, c[1])
        qvec =self.Symplectic12(qvec, c[0])        
        return qvec    
        
    
class Quantum(HilbertSpace):
    def __init__(self, dim, domain, dt=1, order=1, dtype='float'):
        HilbertSpace.__init__(self, dim, domain, dtype=dtype)
        self.dt = dt
        self.order = order 

    def Qmap(self, func_T, func_V, dt, order):
        return QIntegrator(self.dim, self.domain, func_T, func_V, dt, order)
    
    def HamiltonSolver(self, base='position',planck=None):
        if base in ['position', 'fourier', 'p', 'f']:
            return PositionBase(self.dim, self.domain)
        elif base in ['Harmonic','harmonic','h']:
            return HarmonicBase(self.dim, self.domain, planck=planck)
        else:
            raise TypeError("base expect 'position' or 'harmonic'")
    """
    def HamiltonMatrix(self, func_T, func_V):
        PB = PositionBase(self.dim, self.domain)
        matT, matV = PB.get_matrix_funcs(func_T, func_V)
        return matT + matV
    """
    def State(self):
        return State([0]*self.dim, self.domain)
    
    def eigen(self, mat, sort=True):
        evals, evecs = mapy.eigen(mat, sort)
        evecs = [State(vec, self.domain) for vec in evecs]
        return evals, evecs

    def sorteigen(self, evals, evecs, index):
        evals = np.array([evals[i] for i in index])
        evecs = [evecs[i] for i in index]
        return evals, evecs        

    def inner(self, vec0, vec1):
        return np.sum(np.conj(vec0)*vec1)
    
    def inners(self, vec0 ,basis):
        return [self.inner(vec0,base) for base in basis]

    def overlap_sort_index(self, evecs, basis):
        index = []
        for vec in basis:
            ovlp = self.inners(vec, evecs)
            ovlp2 = np.abs(np.conj(ovlp)*ovlp)
            i = np.where(ovlp2 == ovlp2.max())[0][0]
            index.append(i)
        return index
    
    def variance_sort_index(self, evecs, x0=0, qp='q'):
        #pvecs = [self.q2p(vec) for vec in evecs]
        if qp == 'q':
            variance = [np.abs(np.sum(np.conj(vec)*(self.x[0] - x0)**2*vec)) for vec in evecs]
        elif qp == 'p':
            pvec = [ evec.q2p() for evec in evecs]
            variance = [np.abs(np.sum(np.conj(vec)*(self.x[0] - x0)**2*vec)) for vec in pvec]            
        index = [i[0] for i in sorted(enumerate(variance), key=lambda x:x[1])]
        return index

    def quasienergy(self):
        return [np.nan, np.nan]
    """
    def norm(self, vec):
        return np.abs(np.sum(vec*np.conj(vec)))
    
    def q2p(self, vec):
        vec = np.fft.fft(vec)/np.sqrt(self.dim)
        if self.Shift[1]:
            vec = np.fft.fftshift(vec)
        return vec
             
    def p2q(self, vec):
        if self.Shift[1]:
            vec = np.fft.fftshift(vec)        
        vec = np.fft.ifft(vec)*np.sqrt(self.dim)
        return vec
    def abs2(self,vec):
        return np.abs(np.conj(vec)*vec)
    def quasienergy(self):
        return [np.nan, np.nan]
    
    def hsmrep(self, vec, grid=[50,50], vrange=None):
        if vrange==None:
            vrange = self.domain

        from mapy.ctypes_wrapper import wrapper
        file_path = mapy.__file__.replace("__init__.py", "shared/libhsm.so")

        try:
            cw = wrapper.call_hsm_rep(file_path)
        except NameError:
            raise(RuntimeError,"libhsm.so not found.")
            
        hsm_imag = cw.husimi_rep(vec, self.dim, self.domain, vrange, grid)

        x = np.linspace(vrange[0][0], vrange[0][1], grid[0])
        y = np.linspace(vrange[1][0], vrange[1][1], grid[1])

        X,Y = np.meshgrid(x,y)
        return X,Y,hsm_imag        
    """

        
        
    
        


        
        

        