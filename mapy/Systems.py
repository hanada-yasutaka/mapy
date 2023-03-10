import numpy as np
#import mapy
import mapy

Cos = lambda x: mapy.cos(x)
Sin = lambda x: mapy.sin(x)
Abs = np.abs
class Symplectic(object):
    def get_Hamiltonian(self):
        return [self.funcT(), self.funcV()]

    def get_DHamiltonian(self):
        return [self.func_dHdp(), self.func_dHdq()]        

    def get_Jacobi(self):
        pass
class Polynomial(Symplectic):
    def __init__(self):
        self.periodicity = [False, False]         


class Standard(Symplectic):
    def __init__(self, k=1):
        self.periodicity = [True, True]
        self.k = k
    def func_dHdq(self):
        return lambda x: -self.k*mapy.sin(x)
    def func_dHdp(self):
        return lambda x: x
    def funcV(self):
        return lambda x: self.k*mapy.cos(x)
    def funcT(self):
        return lambda x: x*x/2

class Suris(Standard):
    def __init__(self, k=1/4):
        self.periodicity = [False,False]
        self.k = k
    def func_dHdq(self):
        return lambda x: 1/np.pi * mapy.arctan(self.k*mapy.sin(2*np.pi*x)/(1+self.k*mapy.cos(2*np.pi*x)))

    def func_dHdp(self):
        return lambda x: x

class PerturbedSuris(Suris):
    def __init__(self, k=1/4,eps = 1/10):
        self.periodicity = [False,False]
        self.k = k
        self.eps = eps 
    def func_dHdq(self):
        return lambda x: -1/np.pi * mapy.arctan(self.k*mapy.sin(2*np.pi*x)/(1-self.k*mapy.cos(2*np.pi*x))) - self.eps*np.pi*mapy.sin(2*np.pi*x)



class Standard_shift(Symplectic):
    def __init__(self, k=1):
        self.periodicity = [True, False]
        self.k = k
    def func_dHdq(self):
        return lambda x: -self.k*mapy.sin(x)
    def func_dHdp(self):
        return lambda x: (x-1.8849)
    def funcV(self):
        return lambda x: self.k*mapy.cos(x)
    def funcT(self):
        return lambda x: (x-3)**2/2

    
class LinearStandard(Symplectic):
    def __init__(self, k=1,omega=np.sqrt(2)):
        self.periodicity = [True, True]
        self.k = k
        self.omega = omega
    def func_dHdq(self):
        return lambda x: -self.k*mapy.sin(x)
    def func_dHdp(self):
        return lambda x: self.omega
    def funcV(self):
        return lambda x: self.k*mapy.cos(x)
    def funcT(self):
        return lambda x: self.omega*x
        
class IntegrableStandard(Symplectic):
    def __init__(self, k=1):
        self.periodicity = [True, True]
        self.k = k
    def func_dHdq(self):
        return lambda x: -1/(4*np.pi**2)*(self.k*mapy.sin(x)/(2+self.k*mapy.cos(x)))
    def func_dHdp(self):
        return lambda x: x
    def funcV(self):
        return lambda x: x**2
    def funcT(self):
        return lambda x: x**2
        
    
class ShudoIkedaMap(Standard):
    def __init__(self, k=1.2, pd=5.0,omega=1.0,nu=3.0):
        Standard.__init__(self, k)
        self.periodicity = [True, False]
        self.k = k
        self.pd = pd
        self.omega = omega
        self.nu = nu
    def funcT(self):
        h = lambda x: (x/self.pd)**(2*self.nu)
        return lambda x: x*x/2*h(x)/(h(x) +1) + self.omega*x
    def func_dHdp(self):
        h = lambda x: (x/self.pd)**(2*self.nu)
        #g=numpy.power(x/self.Para.para[1],2*self.Para.para[3])
        #return ( x*g*(1+self.Para.para[3]+g)/(1+g)/(1+g)+self.Para.para[2] )/twopi        
        return lambda x: x*h(x)*(h(x) + self.nu+1)/(h(x)+1)**2 + self.omega
    
class HypTanStandard(Standard):
    def __init__(self,k=2.0,omega=0.6418, s=4.173, beta=100.0, d=0.4):
        Standard.__init__(self,k)
        self.periodicity = [True, False]            
        self.k = k
        self.omega = omega
        self.s = s
        self.beta = beta
        self.d = d
        
    def func_dHdp(self):
        g = lambda x: (x - self.d)
        s = self.s
        beta = self.beta
        return lambda x: s/2*g(x)*(1+np.tanh(beta*g(x))) + s*g(x)**2*beta/4.0/(np.cosh(beta*g(x))**2) + self.omega

    def funcT(self):
        g = lambda x: (x - self.d)
        s = self.s
        beta = self.beta
        return lambda x: s*g(x)**2/2.0 *( 1 + np.tanh( beta*g(x) ) )/2.0 + self.omega*g(x)

    
        
class Harper(Symplectic):
    def __init__(self, a=1,b=1):
        self.periodicity = [True, True]        
        self.a = a
        self.b = b        
    def func_dHdq(self):
        return lambda x: -self.a*mapy.sin(x)
    def func_dHdp(self):
        return lambda x: -self.b*mapy.sin(x)        
        
    def funcV(self):
        return lambda x: self.a*mapy.cos(x)
    def funcT(self):
        return lambda x: self.b*mapy.cos(x)        
class HarperNormal(Symplectic):
    def __init__(self,a=-0.55):
        self.periodicity = [True, True]
        self.a = a
    def funcT(self):
        return lambda x: mapy.cos(x) + self.a*mapy.cos(x)**2
    def funcV(self):
        return lambda x: mapy.cos(x) + self.a*mapy.cos(x)**2
    def func_dHdq(self):
        return lambda x: - mapy.sin(x) * (1 + 2*self.a*mapy.cos(x)) 
    def func_dHdp(self):
        return lambda x: - mapy.sin(x) * (1 + 2*self.a*mapy.cos(x)) 



        
class Harmonic(Polynomial):
    def __init__(self, omega=1):
        Polynomial.__init__(self)
        self.omega = omega
    def func_dHdq(self):
        return lambda x: self.omega**2*x
    def func_dHdp(self):
        return lambda x: x
    def funcV(self):
        return lambda x: self.omega**2*x**2/2
    def funcT(self):
        return lambda x: x*x/2
        
class Henon(Polynomial):
    def __init__(self, epsilon=4/10):
        Polynomial.__init__(self)           
        self.epsilon = epsilon
    def func_dHdq(self):
        return lambda x: self.epsilon*(4*x + x**2)
    def func_dHdp(self):
        return lambda x: x
    def funcV(self):
        return lambda x: self.epsilon*(2*x*x + x**3/3)
    def funcT(self):
        return lambda x: x*x/2
    
class DoubleWell(Polynomial):
    def __init__(self):
        Polynomial.__init__(self)           
        #self.epsilon = epsilon
    def func_dHdq(self):
        return lambda x: 4*x*(-1+x**2)
    def func_dHdp(self):
        return lambda x: x
    def funcV(self):
        return lambda x: (x**2 - 1)**2
    def funcT(self):
        return lambda x: x*x/2


class Test(Polynomial):
    def __init__(self, a=-1):
        Polynomial.__init__(self)
        self.a = a
    def funcT(self):
        return lambda x: x**2/2 + self.a*x**4/4
    def funcV(self):
        return lambda x: x**2/2 + self.a*x**4/4        
    def func_dHdq(self):
        return lambda x:(x+self.a*x**3)
    def func_dHdp(self):
        return lambda x:(x+self.a*x**3)    
        
class JeremyNormal(Symplectic):
    def __init__(self, a1=1, a2=-55/100, b=5/100, phi=0):
        self.periodicity = [True, True]
        self.a1 = a1
        self.a2 = a2
        self.b = b
        self.phi = phi
        
    def funcT(self):
        return lambda q,p: self.a1/2*(q**2 + p**2) + self.a2*(q**4 + q**4 + q*q*p*p + p*p*q*q)
    
    def funcV(self):
        return lambda q,p: np.abs(self.b)*(q**4 + p**4 - 3*(q*q*p*p + q*q*p*p))*mapy.cos(self.phi) \
                - ((p**3*q + p*p*q*p + p*q*p*p + q*p**3) - (q**3*p + q*q*p*q + q*p*q*q + p*q**3))*mapy.sin(self.phi)

    def _funcT(self):
        return lambda q,p: (self.a1/2*(Cos(p)**2 + Cos(q)**2)) + self.a2*(Cos(p)**2 + Cos(q)**2)**2
        
    def _funcV(self):
	    return lambda q,p:\
            Abs(self.b)*(\
                (Cos(p)**4 - 6*Cos(p)**2*Cos(q)**2 + Cos(q)**4)*Cos(self.phi) \
                - 4*(Cos(p)**3*Cos(q) - Cos(p)*Cos(q)**3)*Sin(self.phi) \
            )
    def func_dHdq(self):
        return lambda q,p:\
            -self.a1*Cos(q)*Sin(q) - 4*self.a2*Cos(q)*Sin(q)*( Cos(p)**2 + Cos(q)**2 ) + \
            Abs(self.b)*(\
                Cos(self.phi)*Sin(q)*(12*Cos(p)**2*Cos(q)*Sin(q) - 4*Cos(q)**3) \
                - 4*Sin(self.phi)*(-(Cos(p)**3*Sin(q)) + 3*Cos(p)*Cos(q)**2*Sin(q)) \
            )
    def func_dHdp(self):
        return lambda q,p:\
            -self.a1*Cos(p)*Sin(p) - 4*self.a2*Cos(p)*Sin(p)*(Cos(p)**2 + Cos(q)**2) + \
            Abs(self.b)*( \
                Cos(self.phi)*(-4*Cos(p)**3*Sin(p) + 12*Cos(p)*Cos(q)**2*Sin(p)) \
                - 4*(-3*Cos(p)**2*Cos(q)*Sin(p) + Cos(q)**3*Sin(p))*Sin(self.phi)
            )            

class HaradaNormal(Polynomial):
    def __init__(self, a=-1,b=-1):
        Polynomial.__init__(self)
        self.a = a
        self.b = b
    def funcT(self):
        return lambda x,y: (x**2 + y**2)/2  + self.a * (x**4 + y**4 + x**2*y**2 + y**2*x**2)/4

    def funcV(self):
        return lambda x,y: self.b*(x**2*y**2 + y**2*x**2)/2
    def func_dHdq(self):
        return lambda x,y: x + 2*self.b*x*y**2 + 2*self.a*x*(x**2/2 + y**2/2) 
    def func_dHdp(self):
        return lambda x,y: y + 2*self.b*x**2*y + 2*self.a*y*(x**2/2 + y**2/2)
    
class NormannNormal(HaradaNormal):
    def __init__(self, a=-1,b=0.05):
        HaradaNormal.__init__(self, a, b)
        self.periodicity = [True, False]


class JeremyNormal(HaradaNormal):
    def __init__(self, a1=1, a2=-55/100, b=5/100, phi=0):
        self.periodicity = [True, True]
        self.a1 = a1
        self.a2 = a2
        self.b = b
        self.phi = phi
        
    def funcT(self):
        return lambda q,p: self.a1/2*(q**2 + p**2) + self.a2*(q**4 + p**4 + q**2*p**2 + p**2*q**2)
    
    def funcV(self):
        return lambda q,p: np.abs(self.b)*((q**4 + p**4 - 3*(q*q*p*p + p*p*q*q))*mapy.cos(self.phi)
                - ((p**3*q + p*p*q*p + p*q*p*p + q*p**3) - (q**3*p + q*q*p*q + q*p*q*q + p*q**3))*mapy.sin(self.phi))

    def _funcT(self):
        return lambda q,p: (self.a1/2*(Cos(p)**2 + Cos(q)**2)) + self.a2*(Cos(p)**2 + Cos(q)**2)**2
        
    def _funcV(self):
        return lambda q,p:\
            Abs(self.b)*(\
                (Cos(p)**4 - 6*Cos(p)**2*Cos(q)**2 + Cos(q)**4)*Cos(self.phi) \
                - 4*(Cos(p)**3*Cos(q) - Cos(p)*Cos(q)**3)*Sin(self.phi) \
            )
    def func_dHdq(self):
        return lambda q,p:\
            -self.a1*Cos(q)*Sin(q) - 4*self.a2*Cos(q)*Sin(q)*( Cos(p)**2 + Cos(q)**2 ) + \
            Abs(self.b)*(\
                Cos(self.phi)*Sin(q)*(12*Cos(p)**2*Cos(q)*Sin(q) - 4*Cos(q)**3) \
                - 4*Sin(self.phi)*(-(Cos(p)**3*Sin(q)) + 3*Cos(p)*Cos(q)**2*Sin(q)) \
            )
    def func_dHdp(self):
        return lambda q,p:\
            -self.a1*Cos(p)*Sin(p) - 4*self.a2*Cos(p)*Sin(p)*(Cos(p)**2 + Cos(q)**2) + \
            Abs(self.b)*( \
                Cos(self.phi)*(-4*Cos(p)**3*Sin(p) + 12*Cos(p)*Cos(q)**2*Sin(p)) \
                - 4*(-3*Cos(p)**2*Cos(q)*Sin(p) + Cos(q)**3*Sin(p))*Sin(self.phi)
            )            

class Standard_SBCH(Standard):
    def __init__(self, k,tau,bch=1):
        self.k = k
        self.tau = tau
        if bch==1:
            self.funcV = self.bch1_func0
            self.funcT = self.bch1_func1
        if bch==3:
            self.funcV = self.bch3_func0
            self.funcT = self.bch3_func1
        if bch==5:
            self.funcV = self.bch5_func0
            self.funcT = self.bch5_func1
        if bch==7:
            self.funcV = self.bch7_func0
            self.funcT = self.bch7_func1            
        if bch==9:
            self.funcV = self.bch9_func0
            self.funcT = self.bch9_func1            
        if bch==11:
            self.funcV = self.bch11_func0
            self.funcT = self.bch11_func1            
                
    def bch1_func0(self,q,p):
        return -self.k*np.sin(twopi*q)/twopi
    def bch1_func1(self,q,p):
        return p
        
    def bch3_func0(self,q,p):
        return ( self.bch1_func0(q,p) + 
            self.tau**2*(
                (self.k*p**2*Pi*Sin(2*Pi*q))/6. - 
                (self.k**2*Cos(2*Pi*q)*Sin(2*Pi*q))/(24.*Pi)
            )
        )
                
    def bch3_func1(self,q,p):
        return ( self.bch1_func1(q,p) + 
            -self.tau**2*(self.k*p*Cos(2*Pi*q))/6.
        )
        
    def bch5_func0(self,q,p):
        return (self.bch3_func0(q,p) +
            self.tau**4*(
                (self.k*p**4*Pi**3*Sin(2*Pi*q))/90. - 
                (self.k**2*p**2*Pi*Cos(2*Pi*q)*Sin(2*Pi*q))/10. + 
                (self.k**3*Cos(2*Pi*q)**2*Sin(2*Pi*q))/(240.*Pi) - 
                (self.k**3*Sin(2*Pi*q)**3)/(480.*Pi)
            )
        )
    def bch5_func1(self,q,p):
        return (self.bch3_func1(q,p) +
            self.tau**4*(
                -(self.k*p**3*Pi**2*Cos(2*Pi*q))/45. + 
                (self.k**2*p*Cos(2*Pi*q)**2)/30. - 
                (self.k**2*p*Sin(2*Pi*q)**2)/60.
            )
        )
    def bch7_func0(self,q,p):
        return (self.bch5_func0(q,p) +
            self.tau**6*(
                (self.k*p**6*Pi**5*Sin(2*Pi*q))/945. - 
                (3*self.k**2*p**4*Pi**3*Cos(2*Pi*q)*Sin(2*Pi*q))/70. + 
                (43*self.k**3*p**2*Pi*Cos(2*Pi*q)**2*Sin(2*Pi*q))/840. - 
                (self.k**4*Cos(2*Pi*q)**3*Sin(2*Pi*q))/(1680.*Pi) - 
                (5*self.k**3*p**2*Pi*Sin(2*Pi*q)**3)/336. + 
                (self.k**4*Cos(2*Pi*q)*Sin(2*Pi*q)**3)/(630.*Pi)
            )        
        )
    def bch7_func1(self,q,p):
        return (self.bch5_func1(q,p) + 
            self.tau**6*(
                -(self.k*p**5*Pi**4*Cos(2*Pi*q))/315. + 
                (8*self.k**2*p**3*Pi**2*Cos(2*Pi*q)**2)/315. - 
                (self.k**3*p*Cos(2*Pi*q)**3)/140. - 
                (11*self.k**2*p**3*Pi**2*Sin(2*Pi*q)**2)/630. + 
                (5*self.k**3*p*Cos(2*Pi*q)*Sin(2*Pi*q)**2)/336.
            )
        )
    def bch9_func0(self, q,p):
        return (self.bch7_func0(q,p) + 
            self.tau**8*(
                (self.k*p**8*Pi**7*Sin(2*Pi*q))/9450. - 
                (3*self.k**2*p**6*Pi**5*Cos(2*Pi*q)*Sin(2*Pi*q))/175. + 
                (191*self.k**3*p**4*Pi**3*Cos(2*Pi*q)**2*Sin(2*Pi*q))/2520. - 
                (8*self.k**4*p**2*Pi*Cos(2*Pi*q)**3*Sin(2*Pi*q))/315. + 
                (self.k**5*Cos(2*Pi*q)**4*Sin(2*Pi*q))/(10080.*Pi) - 
                (17*self.k**3*p**4*Pi**3*Sin(2*Pi*q)**3)/720. + 
                (73*self.k**4*p**2*Pi*Cos(2*Pi*q)*Sin(2*Pi*q)**3)/2520. - 
                (37*self.k**5*Cos(2*Pi*q)**2*Sin(2*Pi*q)**3)/(40320.*Pi) + 
                (31*self.k**5*Sin(2*Pi*q)**5)/(161280.*Pi)
            )        
        )
    def bch9_func1(self,q,p):
        return (self.bch7_func1(q,p) + 
            self.tau**8*(
                (-2*self.k*p**7*Pi**6*Cos(2*Pi*q))/4725. + 
                (22*self.k**2*p**5*Pi**4*Cos(2*Pi*q)**2)/1575. - 
                (2*self.k**3*p**3*Pi**2*Cos(2*Pi*q)**3)/105. + 
                (self.k**4*p*Cos(2*Pi*q)**4)/630. - 
                (37*self.k**2*p**5*Pi**4*Sin(2*Pi*q)**2)/3150. + 
                (17*self.k**3*p**3*Pi**2*Cos(2*Pi*q)*Sin(2*Pi*q)**2)/360. - 
                (self.k**4*p*Cos(2*Pi*q)**2*Sin(2*Pi*q)**2)/105. + 
                (5*self.k**4*p*Sin(2*Pi*q)**4)/2016.
            )            
        )
    def bch11_func0(self,q,p):
        return ( self.bch9_func0(q,p) +
            self.tau**10*(
                (self.k*p**10*Pi**9*Sin(2*Pi*q))/93555. - 
                (79*self.k**2*p**8*Pi**7*Cos(2*Pi*q)*Sin(2*Pi*q))/11550. + 
                (17183*self.k**3*p**6*Pi**5*Cos(2*Pi*q)**2*Sin(2*Pi*q))/207900. - 
                (4211*self.k**4*p**4*Pi**3*Cos(2*Pi*q)**3*Sin(2*Pi*q))/41580. + 
                (29*self.k**5*p**2*Pi*Cos(2*Pi*q)**4*Sin(2*Pi*q))/2310. - 
                (self.k**6*Cos(2*Pi*q)**5*Sin(2*Pi*q))/(55440.*Pi) - 
                (2227*self.k**3*p**6*Pi**5*Sin(2*Pi*q)**3)/83160. + 
                (9139*self.k**4*p**4*Pi**3*Cos(2*Pi*q)*Sin(2*Pi*q)**3)/83160. - 
                (137*self.k**5*p**2*Pi*Cos(2*Pi*q)**2*Sin(2*Pi*q)**3)/3520. + 
                (3*self.k**6*Cos(2*Pi*q)**3*Sin(2*Pi*q)**3)/(6160.*Pi) + 
                (337*self.k**5*p**2*Pi*Sin(2*Pi*q)**5)/59136. - 
                (155*self.k**6*Cos(2*Pi*q)*Sin(2*Pi*q)**5)/(354816.*Pi)
            )        
        )
    def bch11_func1(self,q,p):
        return ( self.bch9_func1(q,p) + 
            self.tau**10*(
                -(self.k*p**9*Pi**8*Cos(2*Pi*q))/18711. + 
                (368*self.k**2*p**7*Pi**6*Cos(2*Pi*q)**2)/51975. - 
                (8*self.k**3*p**5*Pi**4*Cos(2*Pi*q)**3)/275. + 
                (124*self.k**4*p**3*Pi**2*Cos(2*Pi*q)**4)/10395. - 
                (self.k**5*p*Cos(2*Pi*q)**5)/2772. - 
                (49*self.k**2*p**7*Pi**6*Sin(2*Pi*q)**2)/7425. + 
                (2227*self.k**3*p**5*Pi**4*Cos(2*Pi*q)*Sin(2*Pi*q)**2)/27720. - 
                (1073*self.k**4*p**3*Pi**2*Cos(2*Pi*q)**2*Sin(2*Pi*q)**2)/13860. + 
                (149*self.k**5*p*Cos(2*Pi*q)**3*Sin(2*Pi*q)**2)/27720. + 
                (2701*self.k**4*p**3*Pi**2*Sin(2*Pi*q)**4)/166320. - 
                (337*self.k**5*p*Cos(2*Pi*q)*Sin(2*Pi*q)**4)/59136.
            )
        )
        
        



