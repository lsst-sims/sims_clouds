# -*- coding: utf-8 -*-
from pylab import *
from numpy import *

from scipy import interpolate, fftpack
#from pylab import imshow

## get clean ...

class Clouds:
    
    def __init__(self, ws, s):
        # set the half-diagonal size of the 2D Fourier space equal to the maximum frequency of the 1D power spectrum
        self.windowsize = ws/sqrt(2.)
        # set the sampling directly in the real space : same sampling than the correlation function (arbitrary)
        self.sampling = s
        # and the spacing
        self.xstep = self.windowsize/float(self.sampling)
        
    def DirectClouds(self, correl2D):
        """Directly compute clouds using random noise in real space"""
        # nothing to do with symetry, quarters shift, etc.
        # compute 2D power spectrum from 2D symetric correlation function
        PowerSpec2D = abs(fftpack.fft2(correl2D))
        # get a 2D random gaussien noise with rms = 1
        noise2D = random.normal(zeros(self.sampling*self.sampling), 1.).reshape(self.sampling, self.sampling)
        # a realization is given in Fourier space calculating tf(noise)*sqrt(powerspectum)
        fourierclouds = fftpack.fft2(noise2D)*sqrt(PowerSpec2D)
        # then, inverse Fourier transform to get clouds
        self.clouds = real(fftpack.ifft2(fourierclouds))


    def WriteClouds(self, filename):
        """write clouds in a text file x y z"""
        f = open(filename, 'w')
        # only take 1/4 of the 2D real space (other 3/4 is only use for symetry of the correlation function)
        # but 4 fields of view can be taken separatly from one realization
        # !! DO NOT take the whole field of view because, clouds are correlated over larger distances due to the symetry !!
        for i in range(self.sampling/2):
            for j in range(self.sampling/2):
                f.write(str(i*self.xstep)+'\t'+str(j*self.xstep)+'\t'+str(real(self.clouds[i,j]))+'\n')
        f.close()

## essai
    def FitClouds(self):
        """fit clouds with a 2D spline"""
        # build vectors
        x=[]
        y=[]
        z=[]
        for i in range(self.sampling/2):
            for j in range(self.sampling/2):
                x=append(x,i*self.xstep)
                y=append(y,j*self.xstep)
                z=append(z,self.clouds[i,j])
                
#        print len(x), len(y), len(self.clouds), len(z)
       # print self.sampling/2
        figure(1,(6,4))

        im = imshow(self.clouds, cmap=cm.jet)

        im.set_interpolation('bicubic')
        colorbar()
        show()
        ## tck = interpolate.bisplrep(x,y,z)

##         print  'nombre noeuds en x : ', len(tck[0])
        
##         print 'deg : ', tck[3], tck[4] 
##         print  'coefs : ', tck[2] 
        
##   #      print len(x), len(y), len(self.clouds), len(z), len(tck[3])

##         ## plot
##         ## normalement il faut mettre le nombre de noeuds en x/2 ? 
##         figure(2,(6,4))
##         xnew,ynew = mgrid[-1:1:34j,-1:1:34j] 
##         znew = interpolate.bisplev(xnew[:,0],ynew[0,:],tck)
        
##         print znew.shape
##         print self.clouds.shape
##         imbis = imshow(znew, cmap=cm.jet)
##         colorbar()
##         show()
        stop=raw_input()
