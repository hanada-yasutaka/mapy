import numpy 
twopi = 2.0*numpy.pi


class SymplecticIntegratorError(Exception):
    def __init__(self, value):
        self.value = value
    def __str___(self):
        return repr(self.value)
        
class SymplecticIntegrator(object):
    def __init__(self, funcT, funcV, dt=1, order=1):
        import inspect
        self.funcT = funcT if len(inspect.getargspec(funcT)[0]) == 2 else lambda q,p: funcT(p)
        self.funcV = funcV if len(inspect.getargspec(funcV)[0]) == 2 else lambda q,p: funcV(q)
        self.dt=dt
        self.order = order
        
        if order == 1:
            self.evolve = self.Symplectic1
        elif order == 2:
            self.evolve = self.Symplectic2
        elif order %2 != 0:
            raise SymplecticIntegratorError("Symplectic integrator except 1 or even value")  
        elif order == 4:
            self.evolve = self.Symplectic4

        #else:
        #    self.evolve = lambda x: self.Symplectic(x, 1, self.order)

        elif order == 6:
            self.evolve = self.Symplectic6
            
        elif order == 8:
            self.evolve = self.Symplectic8
        elif order == 10:
            self.evolve = self.Symplectic10
        elif order ==12:
            self.evolve = self.Symplectic12
        elif order ==14:
            self.evolve = self.Symplectic14
        elif order ==16:
            self.evolve = self.Symplectic16
        else:
            raise SymplecticIntegratorError("Symplectic integrator except 1 or even value")        

        
    def _SP(self, x, c=1):
        return numpy.array([x[0], x[1] - c*self.dt*self.funcV(x[0], x[1])])
    def _SQ(self, x, c=1):
        return numpy.array([x[0] + c*self.dt*self.funcT(x[0],x[1]), x[1]])        
        
    def Symplectic1(self, x):
        c = [1,1]
        x =self._SQ(self._SP(x,c[0]),c[1])
#        x =self._SP(self._SQ(x,c[0]),c[1])
        return x

    def Symplectic2(self,x, z=1):
        c = [0.5, 1, 0.5]
        x = self._SP(self._SQ(self._SP(x,z*c[0]),z*c[1]), z*c[2])
#        x = self._SQ(self._SP(self._SQ(x,z*c[0]),z*c[1]), z*c[2])
        return x
    """
    def Symplectic(self,x,z=1,n=2):
        if n == 2:
            return self.Symplectic2(x)
        else:
            beta = 2**(1/(n-1))
            print(n)
            c = [-1/(2-beta), -beta/(2-beta)]
            x = self.Symplectic(self.Symplectic(self.Symplectic(x, z*c[0], n-2), z*c[1],n-2), z*c[0],n-2) 
            return x

    def Symplectic4(self, x, z=1):
        beta=2**(1/3)
        c1=1/(2*(2-beta))
        c2=(1-beta)/(2*(2-beta))
        d1=1/(2-beta)
        d2=-beta/(2-beta)
        y=self._SQ(self._SP(x, c1),d1)
        y=self._SQ(self._SP(y, c2),d2)
        y=self._SQ(self._SP(y, c2),d1)
        y=self._SP(y, c1)
        return y
        #S
    """        
    def Symplectic4(self,x,z=1):
        beta = 2**(1/3)
        c = [1/(2-beta), -beta/(2-beta) ]
        #x = self.Symplectic2(self.Symplectic2(self.Symplectic2(x,z*c[0]),z*c[1]),z*c[0])
        x = self.Symplectic2(self.Symplectic2(self.Symplectic2(x,z*c[0]),z*c[1]),z*c[0])
        return x
        
    def Symplectic6(self,x, z=1):
        beta = 2**(1/5)
        c = [1/(2-beta), -beta/(2-beta)]
        x = self.Symplectic4(self.Symplectic4(self.Symplectic4(x,z*c[0]), z*c[1]),z*c[0])
        return x
        
    def Symplectic8(self, x, z=1):
        beta = 2**(1/7)
        c = [1/(2-beta), -beta/(2-beta)]        
        x = self.Symplectic6(self.Symplectic6(self.Symplectic6(x,z*c[0]), z*c[1]),z*c[0])        
        return x
        
    def Symplectic10(self, x, z=1):
        beta = 2**(1/9)
        c = [1/(2-beta), -beta/(2-beta)]        
        x = self.Symplectic8(self.Symplectic8(self.Symplectic8(x,z*c[0]), z*c[1]),z*c[0])        
        return x

    def Symplectic12(self, x, z=1):
        beta = 2**(1/11)
        c = [1/(2-beta), -beta/(2-beta)]        
        x = self.Symplectic10(self.Symplectic10(self.Symplectic10(x,z*c[0]), z*c[1]),z*c[0])        
        return x
    def Symplectic14(self, x, z=1):
        beta = 2**(1/13)
        c = [1/(2-beta), -beta/(2-beta)]        
        x = self.Symplectic12(self.Symplectic12(self.Symplectic12(x,z*c[0]), z*c[1]),z*c[0])        
        return x
    def Symplectic16(self, x, z=1):
        beta = 2**(1/15)
        c = [1/(2-beta), -beta/(2-beta)]        
        x = self.Symplectic14(self.Symplectic14(self.Symplectic14(x,z*c[0]), z*c[1]),z*c[0])        
        return x


