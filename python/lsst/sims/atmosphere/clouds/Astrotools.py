#               Module : Astrotools
#
#    This tool box provide basic tools needed for processing 
#       astronomical coordinates;
#       currently (02-2009)it contains only
#
#            tsid  : to get local sidereal time at  
#               given MJD and longitude
#            eq2loc : convert equatorial coordinates (RA,DEC of an astronomical
#               target to local coordinates (azimuth, zenith angle)
#               given Universal time, observers longitude and latitude
#
#       Times are not guaranteed to better than 0.1 second between 1980 and 2030
#       Zenith angles ignore atmospheric refraction.
#       Refer to document "Anisotropic Transmission" for details  
#

from numpy import *

# conversion from degree to radian
d2r = pi/180
r2d = 1./d2r
#_____________________________________________________________
#
#                       tsid
#
#    get greenwhich apparent sidereal time from given time (MJD)
#
def tsid(mjd):
        # time (UT  )ellapsed since January 1 2000

        delta = mjd - 51544.5
        
        # greenwhich mean sidereal time (hours)
        #        (approximation according to USNO)

        Tgms =  18.697374558 + 24.06570982441908 * delta

        # nutation in longitude (dpsi)
        omega = 125.04 - 0.052954*delta    #longitude of ascending node(degrees)
        L = 280.47 + 0.98565*delta         # mean longitude of the Sun (  "   )
        epsilon = 23.4393 - 0.0000004*delta# obliquity of equator      (  "   )
        dpsi = -0.000319 * sin(d2r*omega) - 0.000024 * sin(2.*d2r*L )
        # correct Tgms to apparent by adding
        #      equation of equinox (dpsi*cos(epsilon))
        Tgs = Tgms + dpsi * cos( d2r * epsilon )
        Tgs = mod(Tgs,24.)

        # return greenwhich apparent sidereal time

        return Tgs
#________________________________________________________________
#
#                      eq2loc
#
#     convert equatorial coordinates (RA,DEC in degrees) of an astronomical
#     target to local coordinates (azimuth, zenith angle)
#     given observation time (mjd), observer's longitude and latitude
#     
#      Input : aldeti=[alfa_target,delta_target,mjd_observation]
#      Output : locco = [azimuth, z_angle] in decimal degrees
#
def eq2loc(aldeti):
        # fix observation coordinates (decimal degrees)
        obslong = -70.749389 # LSST longitude 
        obslat  = -30.244333 # LSST latitude 
        # get observation data
        alfa = aldeti[0] # Right Ascencion (degrees)
        delta = aldeti[1]# Declination (degrees)
        mjd = aldeti[2]   # modified Julian date (days)
        
        # get local sidereal time
        Tls = tsid(mjd) + obslong/15.
        
        # get hour angle (degrees)
        H = (Tls*15. - alfa)

        # convert angles to radians
        Hr = H * d2r
        der = delta * d2r
        fir = obslat * d2r
        
        # get local coordinate projections ( azimuth and zenith angle)
        cz = sin(fir)*sin(der)+cos(fir)*cos(der)*cos(Hr)
        szsa = cos(der)*sin(Hr)
        szca = -cos(fir)*sin(der)+sin(fir)*cos(der)*cos(Hr)
        # separate trigon fun
        szi = 1./sqrt(1.-cz*cz)
        sa = szsa * szi
        ca = szca * szi
        # extract angles
        z_angle = arccos(cz)*r2d
        azimuth = arctan2(sa,ca)* r2d
        
        locco = [azimuth,z_angle]

        return locco
  #________________________________________________________________
#
#                      airmass
#
#  compute airmass at given zenith angle (degrees) using the approximation
#   by Hardie (1962) 
#     sz1 = secz - 1
#    airmass = secz - 0.0018167*sz1 - 0.002875*sz1^2 - 0.0008083*sz1^3
#    Warning : This model is for sea level
#
#
def airmass(z_angle) :
        secz = 1./cos(z_angle*d2r)
        sz1 = secz-1.
        airmass = secz - 0.0018167*sz1 - 0.002875*sz1**2 - 0.0008083*sz1**3

        return airmass

# Modtran airmass
def mdtram(z_angle):
        secz = 1./cos(z_angle*d2r)
        chi = log(secz)
        X =1.00003873*secz-0.00548117*chi**2\
            + 0.00832316*chi**3 - 0.00711221*chi**3
        return X
        
