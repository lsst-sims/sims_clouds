#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from PowerSpectrum import *
from Clouds import *

# call one realization of the cloud map 

# set dimensions
windowsize = 20 # degrees, diagonal
sampling = 128 # bins

# for fft symetry, double the size
windowsize *= 2
# get power spectrum from structure function
ps = PowerSpectrum(windowsize, sampling)
ps.ComputeStructureFunction()
ps.ComputeCorrelationFunction()
ps.writeCorrel('correlfunc.dat')

# loop on experiences
for i in range(nbexp):
    print 'Exp', i
    # clouds construction
    c = Clouds(windowsize, sampling)
    
    # new method: faster and more simple
    c.DirectClouds(ps.getCorrel2D())
    
    # write clouds to a text file
#    c.WriteClouds('clouds'+str(i)+'.dat')

    c.WriteClouds(cloudfile+str(i)+'.dat')

