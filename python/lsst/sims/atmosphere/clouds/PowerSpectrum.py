# -*- coding: utf-8 -*-
from numpy import *
from scipy import interpolate

class PowerSpectrum:
    
    def __init__(self, ws, s):
        self.windowsize = ws
        self.sampling = s

    def ComputeStructureFunction(self):
        """create structure function based on an analytical model (fitted on data)"""
        # define the x range
        self.xstep = self.windowsize/float(self.sampling)
        x = []
        for i in range(self.sampling):
            x.append(i*self.xstep)
        self.x = array(x)
        # useful variables
        x0 = 0.
        x1 = 1.75
        y1 = 0.04
        ymm = 0.08
        al = -(1/(x1-x0))*log(1-(y1/ymm))
        # structure function
        self.SF = ymm*(1.-exp(-al*self.x))   
        
    def ComputeCorrelationFunction(self):
        """compute the correlation function from the structure function"""
        # following the definition astro-ph/0703157v1
        tmp = -.5*self.SF**2
        self.correl = -min(tmp)+tmp

    def getCorrel2D(self):
        """return a 2D correlation function"""
        # 1D windowssize correspond to 2D diagonal
        localxstep = self.xstep/sqrt(2.)
        tmp = -.5*self.SF**2
        correltmp = -min(tmp)+tmp
        # interpolation of the correlation function on 2D real space
        interpfunc = interpolate.interp1d(self.x, correltmp, bounds_error=False, fill_value = 0.)
        # 2D sampling = 1D sampling (arbitrary)
        correl2D = zeros(self.sampling*self.sampling).reshape(self.sampling, self.sampling)
        # symetry for the 2D correlation function
        # to define a periodic function
        # correl[size-i][size-j] = correl[i][j]
        # correl[i][size-j] = correl[i][j]
        # correl[size-i][j] = correl[i][j]
        for i in range(self.sampling/2+self.sampling%2):
            for j in range(self.sampling/2+self.sampling%2):
                co = interpfunc(sqrt((i*i+j*j))*localxstep)
                correl2D[i,j] = co
                correl2D[self.sampling-1-i,self.sampling-1-j] = co
                correl2D[i,self.sampling-1-j] = co
                correl2D[self.sampling-1-i,j] = co
        return correl2D
       
    def writeSF(self, filename):
        """write structure function in a text file"""
        f = open(filename, 'w')
        for i in range(len(self.x)):
            f.write(str(self.x[i])+'\t'+str(self.SF[i])+'\n')
        f.close()
        
    def writeCorrel(self, filename):
        """write correlation function in a text file"""
        f = open(filename, 'w')
        for i in range(len(self.x)):
            f.write(str(self.x[i])+'\t'+str(self.correl[i])+'\n')
        f.close()
        
