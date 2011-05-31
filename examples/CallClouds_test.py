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
from numpy import *
import Astrotools
import copy
d2r = pi/180.
r2d = 1./d2r

#________________________________________________________________________________
#
# pointing file info (opsim)

## on apccal01 mlsstpath = '/groups/LSST/cecile/'
## on Uchuu
mlsstpath = '/Users/cecile/LSST/Nuages/'
LSSTdatafile = 'LSSToutputs/opsim361/'
expname = 'Opsim3.61.list'
#____________
#
# atmos parameters history file
## independant des scenar d'obs 
Atmos_parmfile= mlsstpath+'Atmosphereoutputs/Opsim129.02Atm_history.dat'


aparmf = open(Atmos_parmfile,'r')
aparmlis = list(aparmf)
aparmf.close()
#
#               
#               Adding atmosphere parameters at appropriate MJD time 
#
print 'succesfully read ', len(aparmlis),'atmos parms in file : ', Atmos_parmfile

#
# extract MJD start 
#
MJDS = float( aparmlis[0].split('$')[0].split('=')[1])
mjdstep = 0.01
#
# extract Aerosol drivers
#
natmos = len(aparmlis)
avpl = [[]]* natmos # list of aerosol parameters
atpl = [[]]* natmos # list of strings containing other atmosphere parameters
for iat in range(0,natmos):
    for parv in aparmlis[iat].split('$') :
        if parv.split('=')[0].strip()=='VIS0':
            vis0 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='VISAMP':
            visamp = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='VISAZ':
            visaz = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='ZAER11':
            zaer11 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='ZAER12':
            zaer12 = float(parv.split('=')[1])
    avpl[iat] = [vis0,visamp,visaz,zaer11,zaer12]
    atpl[iat] = [aparmlis[iat].split('VIS0')[0]]
print 'Aerosol drivers set for ' , iat, 'atmospheres'
print 'last are :' , avpl[iat],atpl[iat]
#

#_______________________________________________________________________________
#



#_______________________________________________________________________________
#
#               Observing History acquisition
#               including MJD start and MJD end
#

mydata = mlsstpath + LSSTdatafile+expname
#
# open obs history file and convert to list
#
mydatf = open(mydata,'r')
fullhlines= list(mydatf)
mydatf.close()

#
# get list of retrievable parameters in file
headers = fullhlines[0]
head = headers.split('\t')
for i in range(0,len(head)):
	head[i]=head[i].strip()
#  list of obs.history parameters to be extracted
## 
headextra= ['obshistid', 'expmjd', \
	    'fieldra','fielddec']

# identify column index of parameters to be retrieved

print headextra

hecolist = []
for hh in headextra:
    hecolist = hecolist + [head.index(hh)]
    if hh == 'expMJD' :
        colmjd = head.index(hh)
    if hh == 'fieldRA' :
        colRA = head.index(hh)
    if hh == 'fieldDec' :
        colDec = head.index(hh)
# 
#
#
# extract usefull data and format in single string
#      field data are stored in table :
#                 pplis
#     as single strings in the format :
#           parm_0 = val_0 $ parm_1 = val_1 $ ....
#      numerical values of mjd and and zenith angle
#      are stored in eponym tables
#
batsize = 30000
batextens = []
#
#building batch name extension list
#
for j in range(1,40):
    icar = '.0'
    if j > 9 :
        icar='.'
    batextens = batextens +[icar+str(j)]

for ibat in range(0,39):
    isb = 1 + ibat*batsize
    ieb = isb + batsize
    obshlines = fullhlines[isb: ieb]
#
    print 'Succesfully read ', len(obshlines),' pointings in file : '\
          ,LSSTdatafile 
    print  'from  istart ',isb , obshlines[0]
    print  '   to iend   ',ieb , obshlines[batsize-1]


      
    nobs = len(obshlines)
    pplis = [[]]*nobs
    mjdlis = [[]]*nobs
    mjdlis[0] = MJDS
    RAlis = [[]]*nobs
    Declis = [[]]*nobs
    zanglis = [[]]*nobs
    for ipl in range(0,nobs):
        pl = obshlines[ipl][:-1]
        pli=''
        pline=pl.split('\t')
        for ic in hecolist:
            pli = pli + head[ic] + ' = ' + pline[ic] + ' $ '
            pplis[ipl] = pli
            mjdlis[ipl] = float(pline[colmjd])
            RAlis[ipl] = float(pline[colRA])
            Declis[ipl] = float(pline[colDec])
        nend = ipl+1
        #print nend
    print 'pointing parms extracted for  ', nend,'  pointings'
# pplis contains one parameter string per LSST pointing
# mjdlis contains MJD times of LSST pointings 
#_______________________________________________________________________________
#=====================================================
# Search for atmospheric parameters at pointing times
#   derive zangle and azimuth from RA, Dec, MJD
#   derive aerosol parameters from list
#   format pointing sequence
#
    for ipoint in range(0,nend) :
        Ra = RAlis[ipoint]*r2d
        Dec = Declis[ipoint]*r2d
        mjd = mjdlis[ipoint]
        azza = Astrotools.eq2loc([Ra,Dec,mjd])
        azim = azza[0] # field azimuth at mjd(degrees)
        zang = azza[1] # field zenith angle at mjd (degrees)
        idmjd = int((mjd-MJDS)/mjdstep)
        atmoparm = atpl[idmjd][0]
        vis0 = avpl[idmjd][0]
        vamp = avpl[idmjd][1] # origine des azimuth pour les aerosols
        vphi = avpl[idmjd][2]
        zaer11 = avpl[idmjd][3]
        zaer12 = avpl[idmjd][4]
        scale1 = avpl[idmjd][5]
        zaer21 = avpl[idmjd][6]
        zaer22 = avpl[idmjd][7]
        scale2 = avpl[idmjd][8]
        si2z = 1.
        if zang < 45. :
            si2z = sin(2 *zang * d2r)
        vis = vis0 + vamp * sin(azim * d2r - vphi)* si2z
        if vis < 4.0:
            vis = 4.0
        pplis[ipoint] = 'ID' + pplis[ipoint][12:]\
                    + ' ZANGLE = % 8.3f' %zang\
                    + ' $ '+ atmoparm \
                    + ' VIS = % 8.1f'%vis\
                    + ' $ ZAER11 = % 6.0f'%zaer11\
                    + ' $ ZAER12 = % 6.0f'%zaer12\
                    + ' $ SCALE1 = % 6.0f'%scale1\
                    + ' $ ZAER21 = % 6.0f'%zaer21\
                    + ' $ ZAER22 = % 6.0f'%zaer22\
                    + ' $ SCALE2 = % 6.0f'%scale2\

##########################################################
## Call cloud function 

## definition of characteristics for clouds simulation
## Pour le moment on definit ici ce qu'on va utiliser

	windowsize = 5  ## window size for simulation (diagonal length in deg) - roughly LSST FoV [side 3.5 deg]
	sampling = 20   ## nb of bins to sample the field - so that we get bins of roughly 15' on diag 

	## characteristics of the SF

	ext_m = 0.08 # long range gray extinction - to be played with for photometric nights or not 
	# position du point d'inflexion pour la fonction de structure (euh ?)
	xi = 1.75
	yi =0.04

	cloudfilebs='clouds'+ pplis[ipoint][12:]
	
	# for fft symetry, double the size
	windowsize *= 2
	# get power spectrum from structure function
	ps = PowerSpectrum(windowsize, sampling)
	ps.ComputeStructureFunction()
	# ps.writeSF('SF.dat')
	ps.ComputeCorrelationFunction()
	ps.writeCorrel('correlfunc.dat')

## get information : photometric nights 	
## modifies long range gray extinction based on opsim info 

## where do I get this info ? at this stage...


## generate clouds map for each pointing
## test on 3 iterations
	for i in range(0,3):
		# clouds construction
	   c = Clouds(windowsize, sampling)    
	   # new method: faster and more simple
	   c.DirectClouds(ps.getCorrel2D())
	   # write clouds to a text file
	   c.WriteClouds(cloudfilebs+str(i)+'.dat')

## keep file info in parameters  file 
    pplis[ipoint] + ' GRAY_EXT_NORM = %6.0 f' %ext_m  \
    + ' GRAY_EXT_FILE = %s' %cloudfile  + ' $ \n'

## remplacer par fonction de structure
    ## ajout d'un module qui reproduit write clouds 

    

######################################################
#Write pointing sequence with atmosphere parameters
#  
    sequence_Parmfile = mlsstpath +'Runparm/'\
                        + expname + batextens[ibat] + '_parmlist.dat'
    seqparmfil = open(sequence_Parmfile,'w')
    for ic in range(0,nend):
        seqparmfil.write(pplis[ic])
    seqparmfil.close()
    print nend,'parameter lines  written in batch ',expname+batextens[ibat]
    print '       from ',pplis[0][:20], '  to ',pplis[ipoint][:20]
