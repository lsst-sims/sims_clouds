import numpy
import pylab
from lsst.sims.atmosphere.clouds.Clouds import Clouds
from scipy import interpolate
#import Clouds
#import PowerSpectrum 

def make_clouds(sampling=256, seed=42):
    # Set up the power spectrum for the clouds.
    rad_fov = 1.8  #degrees? or radians?
    windowsize = numpy.sqrt(2.*(2*rad_fov)**2)

    clouds = Clouds(windowsize, sampling)
    # Set up the power spectrum for the clouds.                
    clouds.setupPowerSpectrum()
    # Generate the clouds in 'real space'. 
    clouds.DirectClouds(randomSeed=seed)
    cloudimage = clouds.clouds
    # The cloud extinction image is just an 'image' numbered from 0/0 to samplesize/samplesize
    #  So we need to add x/y locations appropriate to stellar x/y locations.
    # PROBABLY NEED TO CHANGE THIS SCALING, APPROPRIATELY FOR TANGENT PLANE
    x, y = numpy.mgrid[0:sampling, 0:sampling]
    x = x / (numpy.max(x)/2.0) - 1.0  # if focal plane goes from -1 to 1
    y = y / (numpy.max(y)/2.0) - 1.0   # probably need to change this.
    cloudxy = numpy.column_stack([numpy.ravel(x), numpy.ravel(y)])
    # Create an object we could use for interpolating cloud density at a point.
    print numpy.shape(cloudxy), numpy.shape(numpy.ravel(cloudimage))
    cloud_interp = interpolate.NearestNDInterpolator(cloudxy, cloudimage.ravel())
    # and plot
    pylab.figure()
    pylab.imshow(cloudimage)
    return cloud_interp
    

if __name__=="__main__":
    xy = numpy.array([[0.5, 0.5], [0.2, 0.3], [0.1, -.3]])
    cloud_interp = make_clouds(256, 42)
    c = cloud_interp(xy)
    print c
    cloud_interp = make_clouds(256, 43)
    c=  cloud_interp(xy)
    print c
    cloud_interp =make_clouds(256, 44)
    c=  cloud_interp(xy)
    print c

    pylab.show()
