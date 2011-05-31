# -*- coding: cp1252 -*-
#==============================================================================
#
#               module makesequence.py  (alias : Pointing.py):
#
#-------------------------------------------------------------------------
# Read description of a sequence of LSST pointings  from Opsim361
# Call MakeCouds.py to generate gray absorption factor
#_______________________________________________________________________________
# common items

## nettoyer les appels : place mem deg

from numpy import *
from Astrotools import * # utilitaire perso 
from copy import * 
from PowerSpectrum import *
from Clouds import *

d2r = pi/180.
r2d = 1./d2r
mlsstpath = '/Users/cecile/LSST/Nuages/'
expname = 'opsim3.61.list'
LSSTdatafiles = 'LSSToutputs/opsim361/'

#_______________________________________________________________________________
#
#               Observing History acquisition
#               including MJD start and MJD end
#

mydata = mlsstpath + LSSTdatafiles+expname
#
# open obs history file and convert to list
#
mydatf = open(mydata,'r')
fullhlines= list(mydatf)
mydatf.close()


# get list of retrievable parameters in file
headers = fullhlines[0]
head = headers.split('\t')
for i in range(0,len(head)):
	head[i]=head[i].strip()
#  list of obs.history parameters to be extracted
headextra= ['min(obshistid)', 'expmjd','fieldra','fielddec']

# identify column index of parameters to be retrieved
hecolist = []
for hh in headextra:
##	hecolist = hecolist + [head.index(hh)]
	if hh == 'min(obshistid)' :
		colID = head.index(hh)			
	if hh == 'expmjd' :
		colmjd = head.index(hh)
	if hh == 'fieldra' :
		colRA = head.index(hh)
	if hh == 'fielddec' :
		colDec = head.index(hh)

## get the input file line by line and extract arguments
		
## print 'Succesfully read ', len(fullhlines),' pointings in file : ' , mydata
## print  'from  istart ',  1,' ex : ', fullhlines[1]
## print  'to iend   ', len(fullhlines), 'ex : ',fullhlines[len(fullhlines)-1]

inlen = len(fullhlines)
spl_line = []
for xx in range (1,inlen):
## test sur un seul pointage 
##for xx in range (1,2):
	line = fullhlines[xx]
	spl_line = line.split()
	obs = spl_line[colID] 
	mjd = spl_line[colmjd] 	
	ra = spl_line[colRA]
	dec = spl_line[colDec]
#	print obs, mjd, ra, dec  

# Call the generation of clouds simulation for each

# call one realization of the cloud map 

# set dimensions
windowsize = 20 # degrees, diagonal
sampling = 128 # bins

# for fft symetry, double the size
windowsize *= 2
# get power spectrum from structure function
ps = PowerSpectrum(windowsize, sampling)
ps.ComputeStructureFunction()
#ps.writeSF('SF.dat')
ps.ComputeCorrelationFunction()
ps.writeCorrel('correlfunc.dat')

## generate clouds map for each pointing
## test on 3 iterations 
for i in range(0,3):
  # clouds construction
    c = Clouds(windowsize, sampling)
    
    # new method: faster and simpler
    c.DirectClouds(ps.getCorrel2D())

    # fit the clouds with a 2d spline
    c.FitClouds()
    
    # write clouds to a text file
    c.WriteClouds('clouds'+str(i)+'.dat')


