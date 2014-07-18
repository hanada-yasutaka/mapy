import mapy
import numpy as np
import functools
import datetime
twopi=2*np.pi

class HilbertSpace(object):
    def __init__(self, dim, domain):
        self.dim = self.__setDim(dim)
        self.domain = self.__setDomain(domain)
        self.h = self.__setPlanck()
        
    def __setDim(self, dim):
        if (dim <= 0) or (dim %2 !=0):
            raise AttributeError("expects dim >0 and even (%d)" % (dim))
        return dim
        

    def __setDomain(self, domain):
        qr, pr  = domain[0], domain[1]

        if (qr[0] >= qr[1]):
            raise AttributeError("qmin > qmax (%f,%f)" % (qr[0],qr[1]))
        if (pr[0] >= pr[1]):
            raise AttributeError("pmin > pmax (%f,%f)" % (pr[0],pr[1]))
        return domain
        
    def __setPlanck(self):
        qr, pr = self.domain
        return (pr[1] - pr[0])*(qr[1] - qr[0])/self.dim

class State(HilbertSpace):
    def __init__(self, dim, domain):
        HilbertSpace.__init__(self, dim, domain)
        
    def zeros(self):
        return np.zeros(self.dim,dtype=np.complex128)
    def cs(self):
        pass
    def qconst(self):
        pass
    def pconst(self):
        pass

class SplitOperator(HilbertSpace):
    def __init__(self, dim, domain, funcT, funcV):
        HilbertSpace.__init__(self, dim, domain)
        self.funcT = funcT
        self.funcV = funcV
        self.x = [np.linspace(self.domain[0][0], self.domain[0][1], self.dim, endpoint=False),
                  np.linspace(self.domain[1][0], self.domain[1][1], self.dim, endpoint=False)]
        self._qrep = True
        self.__Shift()
    
    def __Shift(self):
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
    
    
    def operator_expV(self, x, tau=1,isShift=False):
        if self._qrep and self.Shift[0]:
            x = np.fft.fftshift(x)
        return lambda invec: np.exp(-1.j*self.funcV(x)*tau*twopi/self.h)*invec
        
    def operator_expT(self,x, tau=1, isShift=True):
        if self._qrep and self.Shift[1]:
            x = np.fft.fftshift(x)        
        return lambda invec: np.exp(-1.j*self.funcT(x)*tau*twopi/self.h)*invec
    
    def matrix_function(self, f, qp='q'):
        if self._qrep:
            if qp=='q':
                mat = np.eye(self.dim)*f(self.x[0])
            else:
                if self.Shift[1]:
                    mat = np.array([pvec*np.fft.fftshift(f(self.x[1])) for pvec in map(np.fft.fft, np.eye(self.dim))])                    
                else:
                    mat = np.array([pvec*f(self.x[1]) for pvec in map(np.fft.fft, np.eye(self.dim))])
        else:
            pass
        return mat.transpose()
    
    def get_matrix_funcs(self, funcT, funcV):
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
        return mat.transpose()
             
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
        

    
class Qmap(SplitOperator):
    def __init__(self, dim, domain, funcT, funcV, tau=1, order=1):
        SplitOperator.__init__(self,dim, domain, funcT, funcV)
        self.tau = tau
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

    def UnitaryMatrix(self):
        mat = np.eye(self.dim, dtype=np.complex128)
        res = np.array([x for x in map(self.evolve, mat)])
        return res.transpose()
        
    def Symplectic1(self, invec, z=1):
        c = z*np.array([1,1])
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])

        pvec = np.fft.fft(expV(invec))
        qvec = np.fft.ifft(expT(pvec))
        return qvec
        
    def Symplectic2(self, invec, z=1):
        c = z*np.array([1/2,1])
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
    
        pvec = np.fft.fft(expV(invec))
        qvec = np.fft.ifft(expT(pvec))
        pvec = np.fft.fft(expV(qvec))
        qvec = np.fft.ifft(pvec)
        return qvec
    
    def Symplectic4(self, invec, z=1):
        beta = 2**(1/3)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
        qvec =self.Symplectic2(invec, c[0])
        qvec =self.Symplectic2(qvec, c[1])
        qvec =self.Symplectic2(qvec, c[0])        
        return qvec
    
    def Symplectic6(self, invec, z=1):
        beta = 2**(1/5)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
        qvec =self.Symplectic4(invec, c[0])
        qvec =self.Symplectic4(qvec, c[1])
        qvec =self.Symplectic4(qvec, c[0])        
        return qvec
    
    def Symplectic8(self, invec, z=1):
        beta = 2**(1/7)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
        qvec =self.Symplectic6(invec, c[0])
        qvec =self.Symplectic6(qvec, c[1])
        qvec =self.Symplectic6(qvec, c[0])        
        return qvec

    def Symplectic10(self, invec, z=1):
        beta = 2**(1/9)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
        qvec =self.Symplectic8(invec, c[0])
        qvec =self.Symplectic8(qvec, c[1])
        qvec =self.Symplectic8(qvec, c[0])        
        return qvec
    
    def Symplectic12(self, invec, z=1):
        beta = 2**(1/11)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
        qvec =self.Symplectic10(invec, c[0])
        qvec =self.Symplectic10(qvec, c[1])
        qvec =self.Symplectic10(qvec, c[0])        
        return qvec    

    def Symplectic14(self, invec, z=1):
        beta = 2**(1/13)
        c = z*np.array([1/(2-beta),-beta/(2-beta)]) 
        q, p = self.x
        expV = self.operator_expV(q,self.tau*c[0], self.Shift[0])        
        expT = self.operator_expT(p,self.tau*c[1], self.Shift[1])
        qvec =self.Symplectic12(invec, c[0])
        qvec =self.Symplectic12(qvec, c[1])
        qvec =self.Symplectic12(qvec, c[0])        
        return qvec    


    

    def eigen(self, mat, sort=True):
        evals, evecs = mapy.eigen(mat)
        if sort:
            print(evals.real)
            index = [i[0] for i in sorted(enumerate(evals.real), key=lambda x:x[1])]
            print(index)
            evals = np.array([evals[i] for i in index])
            evecs = [evecs[i] for i in index]
        return evals, evecs

    def sorteigen(self, evals, evecs, index):
        evals = np.array([evals[i] for i in index])
        evecs = [evecs[i] for i in index]
        return evals, evecs        
    def inner(self, vec0, vec1):
        return np.sum(np.conj(vec0)*vec1)
    
    def inners(self, vec0 ,basis):
        return [self.inner(vec0,base) for base in basis]

    def overlap_index(self, basis, evecs):
        index = []
        for vec in evecs:
            ovlp = self.inners(vec, basis)
            ovlp2 = np.abs(np.conj(ovlp)*ovlp)
            i = np.where(ovlp2 == ovlp2.max())[0][0]
            index.append(i)
        return index
    
    def variance_index(self, evecs):
        #pvecs = [self.q2p(vec) for vec in evecs]
        variance = [np.abs(np.sum(np.conj(vec)*(self.x[0])**2*vec)) for vec in evecs]
        index = [i[0] for i in sorted(enumerate(variance), key=lambda x:x[1])]
        return index        
    
    def norm(self, vec):
        return np.abs(np.sum(vec*np.conj(vec)))

    def q2p(self, vec):
        if self._qrep:
            vec = np.fft.fft(vec)/np.sqrt(self.dim)
            if self.Shift[1]:
                vec = np.fft.fftshift(vec)
            return vec
        else:
            pass
             
    def p2q(self, vec):
        if self._qrep:
            if self.Shift[1]:
                vec = np.fft.fftshift(vec)        
            vec = np.fft.ifft(vec)*np.sqrt(self.dim)
            return vec
        else:
            pass
    
    def hsmrep(self, vec, col=50,row=50, vrange=None):
        if vrange==None:
            vrange = self.domain

        from mapy.ctypes_wrapper import wrapper
        file_path = mapy.__file__.replace("__init__.py", "shared/libhsm.so")

        try:
            cw = wrapper.call_hsm_rep(file_path)
        except NameError:
            raise(RuntimeError,"libhsm.so not found.")
            
        hsm_imag = cw.husimi_rep(vec, self.dim, self.domain, vrange, [row,col])

        x = np.linspace(vrange[0][0], vrange[0][1], row)
        y = np.linspace(vrange[1][0], vrange[1][1], col)

        X,Y = np.meshgrid(x,y)
        return X,Y,hsm_imag        
        


        
        
        
    
        


        
        

        