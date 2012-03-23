import numpy
import pylab
import lsst.sims.atmosphere.clouds.Clouds as Clouds
import lsst.sims.atmosphere.clouds.PowerSpectrum as PowerSpectrum
#import Clouds
#import PowerSpectrum 

def make_clouds(sampling=256, seed=42):
    # Set up the power spectrum for the clouds.
    rad_fov = 1.8  #degrees? or radians?
    windowsize = numpy.sqrt(2.*(2*rad_fov)**2)
    # and 'for fft symmetry double the windowsize'
    #  (not sure why this is needed)
    windowsize *= 2
    # Really only have to set up PS once, since we're not changing
    # any of the parameters.  (so save correlfun & reuse)
    ps = PowerSpectrum.PowerSpectrum(windowsize, sampling)
    ps.ComputeStructureFunction()
    correlfunc = ps.getCorrel2D()
    clouds = Clouds.Clouds(windowsize, sampling)
    # All randomness in the clouds is set up in clouds.DirectClouds
    numpy.random.seed(seed)
    clouds.DirectClouds(correlfunc)
    cloudimage = clouds.clouds
    # and plot
    pylab.figure()
    pylab.imshow(cloudimage)
    return ps, clouds, cloudimage
    

if __name__=="__main__":
    make_clouds(256, 42)
    make_clouds(256, 43)
    make_clouds(256, 44)
    pylab.show()
