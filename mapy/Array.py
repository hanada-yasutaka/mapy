# -*- coding:utf-8 -*- 
import numpy
import mpmath

class Array(numpy.ndarray):    
    def __new__(cls, input=[],dtype='float'):
        return numpy.asarray(input,dtype=getattr(numpy, dtype)).view(cls)

    def __array_finalize__(self, obj):
        if obj is None: return
    
    @property
    def real(self):
        return numpy.asarray([x.real for x in self]).view(type(self))
    @property
    def imag(self):
        return numpy.asarray([x.imag for x in self]).view(type(self))
    
    def toarray(self):
        # to numpy array, dtype is numpy.complex128
        return numpy.array(self.tolist())
    
    def insert(self, i, x):
        self[i] = x
    def append(self, x):
        l = self.tolist()
        l.append(x)
        return Array(l)
    
    def fftshift(self):
        return Array(numpy.fft.fftshift(self))

    def conj(self):
        return Array(numpy.conj(self))
    
    @classmethod
    def append(arr,x):
        if not isinstance(arr, list):
            arr = arr.tolist()
            arr.append(x)
        return Array(arr)    
    
    @classmethod
    def zeros(cls, len, dtype='int'):
        return Array([0]*len, dtype=dtype)
    
    @classmethod
    def linspace(cls, min, max, num, endpoint=True,dtype=None):
        if dtype!="object":
            return Array(numpy.linspace(min, max,num, endpoint=endpoint))
        else:
            return Array(mpmath.linspace(min, max, num, endpoint=endpoint),dtype=dtype)
    
    @classmethod
    def savetxt(cls, fname, data):
        numpy.savetxt(fname, data)
        #with open(fname,"wb") as f:
        #    for x in data:
        #        [f.write("%s " % xx) for xx in x]
        #        f.write("\n")
    @classmethod
    def loadtxt(cls,fname,dtype='float',row=False):
        if dtype == 'object':
            data = numpy.loadtxt(fname,dtype=numpy.bytes_)
            if len(data.shape) == 1:
                data = [mpmath.mpf(x.decode()) for x in data]
            else:
                data = [[mpmath.mpf(data[i][j].decode()) for i in range(data.shape[0])] for j in range(data.shape[1])]
            return Array(data, dtype=dtype).transpose()
        else:            
            return numpy.loadtxt(fname,dtype=getattr(numpy, dtype))

class mpMatrix(numpy.ndarray):
    def __new__(cls, input, dtype='object'):
        if isinstance(input, int):
            data  = [[mpmath.mpc("0","0")]*input]*input
            obj = numpy.asarray(data, dtype=object).view(cls)
        else:
            obj = numpy.asarray(input, dtype=dtype).view(cls)            
            if len(obj.shape) != 2 and obj.shape[0] != obj.shape[1]:
                raise TypeError("excepted: n by n matrix, but input data shape",obj.shape )
        return obj
    
    def eigen(self,left=False, solver='qeispack', verbose=False):
        solver_list = {'qeispack':self.qeispack, 'lapack':self.lapack}
        if solver not in solver_list:
            raise TypeError("excepted solver:", solver_list.keys())
        self.solver = solver_list[solver]
        return self.solver(verbose)
        

    def qeispack(self,verbose):
        from utility import ctypes_wrapper as cw

        start=time.time()
        rfname = 'matrix_real.dat'
        ifname = 'matrix_imag.dat'
        self._save_matrix(rfname, ifname)
        end = time.time()
        t = time.time() - start
        if verbose:
            print("save matrix size", self.shape, "in", t,'sec.')            

        start= time.time()
        qeispack = cw.qEispack()
        qeispack.file2call_eig(len(self), rfname, ifname)
        end = time.time()
        t = time.time() -start
        if verbose:
            print("(QEISPACK) get eigen-values and -vectors of size", len(self), "in", t,'sec.')
        return self._load_file(verbose)
