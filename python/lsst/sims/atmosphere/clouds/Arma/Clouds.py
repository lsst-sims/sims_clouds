"""
Clouds: a class to compute a cloud image, given a structure function (theta / SF(theta)) -- where the units
of the SF are in mags (and angles are in degrees).
Requires numpy and PImage (and optionally, pylab and PImagePlots). [note that PImage itself requires scipy.fftpack]. 

Typical usage - 
 instantiate the class (cloud = Clouds())
 make an actual cloud image using SF for a given image (with particular sigma/kappa)
      (cloud.makeCloudImage(SFtheta, SFsf, kappa, [seed])
 use cloud image (a 2-d numpy array, cloud.cloudimage) in your own program.

The structure function is not generated here - any function which is described in theta (degrees) / SF (mag) will
 work. See the associated class, ArmaSf, for Tim's structure function generated from an ARMA process (and described in his cloud 
 extinction minipapers, uploaded to confluence http://www.lsstcorp.org:8090/display/CAL/Cloud+Simulation+Documentation ).

"""

import numpy
from pImage import PImage

class Clouds():
    def __init__(self):
        return

    def makeCloudImage(self, SFtheta, SFsf, kappa, fov=3.0, pixscale=9.0, oversample=1.7,
                       seed=None, verbose=False):
        """Create cloud image from SFtheta/SFsf, adding random phase spectrum and scaling cloud extinction
        to kappa and desired sigma across fov.

        SFtheta = numpy array of angular 'theta' values (in degrees) for the structure function. 
        SFsf = numpy array of actual structure function values (in magnitudes, not mags**2). 
               Obviously, the SF is desired to cover (at least) the entire range between 0 and fov/2.0*sqrt(2).
        fov = the desired (approximate) final field of view (diameter) for the cloud image
                - typical LSST value would be 3.0. 
              This may be modified slightly within the class (but a change be reported), as imsize
              must be divisible by 2.  Pixscale will remain the unchanged.
        pixscale = the desired final resolution for the cloud image, per pixel (in arcseconds/pixel)
                   - while LSST pixscale is 0.02"/pix, setting the pixscale this high here (with a large fov)
                   would mean working with a *very* large image. Setting a pixscale of 9.0 means
                   image sizes on the order of 1200x1200 pixels (for a fov=3.0 deg), which is reasonable.
                   NOTE - check timing of creating images using examples/arma_clouds.py, as it seems that some combinations of
                   fov/pixscale can be particularly bad (i.e. take > 2 times longer to create a cloud image)
        oversample = (>=1) .. if this is >1 then an image larger by this factor will first be created
                     from the SF, and then trimmed down to the final fov.
                     This is an option that I don't think really matters and could be safely set to 1
                     to reduce resource use ... I was originally playing with it to see if it altered
                     the final structure function calculated from the cloudimage, but it's a very minor effect.
        seed = random number seed for phase spectrum, if desired.
        """
        # Translate pixscale to degrees/pixel for internal use (from arcseconds/pixel).
        self.pixscale = pixscale
        self.dpixscale = pixscale / 3600.0
        # Generate desired final image size, and adjust (if necessary) the final_fov value
        #  so that imsize can be divisible by 2.        
        self.final_imsize = (fov / self.dpixscale)
        if (self.final_imsize%2 != 0):
            self.final_imsize = self.final_imsize + 1
        self.final_fov = self.final_imsize * self.dpixscale        
        # Calculate, from the overall structure function, the desired standard deviation in opacity
        #  for the clouds over the entire final fov.
        goal_rad_fov = self.final_fov / 2.0 * numpy.sqrt(2) 
        rms_goal = numpy.interp(goal_rad_fov, SFtheta, SFsf)
        # Determine size of image want to generate from SF (may want to generate an image
        #  larger than final desired fov if oversample > 1)
        fov = self.final_fov * oversample
        self.imsize = int(fov/self.dpixscale)
        if (self.imsize%2 != 0):
            self.imsize = self.imsize + 1
        fov = self.imsize * self.dpixscale
        # Translate SFtheta into pixels
        SFx = SFtheta / self.dpixscale         
        # Actually make cloud image.
        im = PImage(shift=True, nx=self.imsize, ny=self.imsize)
        im.makeImageFromSf(sfx=SFx, sf=SFsf)
        if verbose:
            print 'Image: fov', self.final_fov, 'Imsize (Final)', self.final_imsize, 'Imsize (oversampled)', self.imsize, 'Pixscale', self.dpixscale, \
                'deg/pix', '(', self.pixscale, 'arcsec/pix)'
        # Trim image to desired final size (in case of oversampling earlier in setupCloudImage()). 
        trim = round(self.imsize - self.final_imsize)/2.0
        self.cloudimage = im.imageI.real[trim:trim+self.final_imsize, trim:trim+self.final_imsize]
        # Rescale image to have proper mean extinction and variation in cloud extinction. 
        self.rescaleImage(rms_goal, kappa, verbose=verbose)
        if verbose:
            print 'Image: fov', self.final_fov, 'Imsize', self.final_imsize, 'Pixscale', self.dpixscale, \
                'deg/pix', '(', self.pixscale, 'arcsec/pix)'
        return
            
    def rescaleImage(self, sigma_goal, kappa, verbose=False):
        """Rescale a reconstructed image to have the desired variation and overall opacity."""
        if verbose:
            print '# Before rescaling:'
            self._cloudStats()
        # Subtract mean, multiply by new stdev, then add new mean opacity.
        self.cloudimage = self.cloudimage - self.cloudimage.mean()
        self.cloudimage *= sigma_goal / self.cloudimage.std()
        self.cloudimage += kappa
        # Make sure no 'negative' clouds 
        self.cloudimage = numpy.where(self.cloudimage<0, 0, self.cloudimage)
        if verbose:
            print '# After rescaling:'
            self._cloudStats()
        return 

    def _cloudStats(self):
        print '# CloudImage: mean/sigma/min/max:', self.cloudimage.std(), self.cloudimage.mean(), \
            self.cloudimage.min(), self.cloudimage.max()
        return

    def plotCloudImage(self):
        """Generate some additional information and plots about the cloud image, if desired."""
        from pImagePlots import PImagePlots
        import pylab
        im = PImagePlots()
        im.setImage(self.cloudimage)
        im.showImage(copy=True)
        im.hanningFilter()
        im.calcAll()
        im.showPsd2d()
        im.showAcovf2d()
        im.showAcovf1d()
        im.showSf(linear=True)
        #pylab.show()
        return
