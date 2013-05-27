import numpy
from scipy import fftpack

class PImage():
    def __init__(self, shift=True, nx=1000, ny=1000):
        """Init. Does nothing."""
        if (nx%2 != 0) | (ny%2 !=0):
            raise Exception('Make nx and ny divisible by 2')
        # Note that numpy array translate to images in [y][x] order!
        self.shift = shift
        self.nx = nx
        self.ny = ny
        self.xcen = round(self.nx/2.0)
        self.ycen = round(self.ny/2.0)
        self.image = numpy.zeros([self.ny, self.nx], 'float')
        self.padx = 0.0
        self.pady = 0.0
        self.yy, self.xx = numpy.indices(self.image.shape)
        return

    def setImage(self, imdat, copy=False):
        """Set the image using an external numpy array. If 'copy' is true, makes a copy. """
        if copy:
            self.image = numpy.copy(imdat)
        else:
            self.image = imdat    
        self.ny, self.nx = self.image.shape
        if (self.nx%2 != 0) | (self.ny%2 !=0):
            raise Exception('Make image size in x and y divisible by 2')
        self.yy, self.xx = numpy.indices(self.image.shape)
        self.padx = 0.0
        self.pady = 0.0
        self.xcen = round(self.nx/2.0)
        self.ycen = round(self.ny/2.0)
        return

    def zeroPad(self, width=None):
        """Add padding to the outside of the image data. Default width = 1/4 of image."""
        if width != None:
            self.padx += numpy.floor(width)
            self.pady += numpy.floor(width)
        else:
            self.padx += numpy.floor(self.nx / 4.0)
            self.pady += numpy.floor(self.ny / 4.0)
        newy = self.ny + int(2*self.padx)
        newx = self.nx + int(2*self.pady)
        newimage = numpy.zeros((newy, newx), 'float')
        newyy, newxx = numpy.indices(newimage.shape)
        condition = ((newxx >= self.padx) & (newxx < self.padx + self.nx) &
                     (newyy >= self.pady) & (newyy < self.pady + self.ny))
        newimage[condition] = self.image.flatten()
        self.nx = newx
        self.ny = newy
        self.xx = newxx
        self.yy = newyy
        self.image = newimage
        self.xcen = round(self.nx/2.0)
        self.ycen = round(self.ny/2.0)
        return

    def hanningFilter(self, rmax=None):
        """Apply a radial hanning filter to the image data.
        This removes noise in the FFT resulting from discontinuities in the edges of the images."""
        if rmax==None:
            rmax = min(self.nx, self.ny) / 2.0
        rvals = numpy.hypot((self.xx-self.xcen), (self.yy-self.ycen)) 
        self.hanning = numpy.where(rvals<=rmax, (0.5-0.5*numpy.cos(numpy.pi*(1-rvals/rmax))), 0.0)
        self.image *= self.hanning
        self.hanningFilter = True
        return

    def reflectEdges(self, width=None):
        """Extend the edges of the image by reflection.
        The corners aren't dealt with properly, but this might give some help when applying a hanningFilter after."""
        # Extend the size of the image and do some bookkeeping.
        if width == None:
            width = min(self.nx, self.ny) / 4.0            
        self.zeroPad(width)
        # And then put reflected copy of data into the boundaries.        
        #  Reflect/flip left edge.
        xmin = self.padx
        xmax = self.padx * 2 
        ymin = self.pady
        ymax = self.ny - self.pady
        self.image[ymin:ymax, 0:xmin] = numpy.fliplr(self.image[ymin:ymax, xmin:xmax])
        # Reflect/flip right edge
        xmin = self.nx - self.padx*2
        xmax = self.nx - self.padx
        self.image[ymin:ymax, (self.nx-self.padx):self.nx] = numpy.fliplr(self.image[ymin:ymax, xmin:xmax])
        # Reflect/flip bottom edge
        xmin = self.padx
        xmax = self.nx - self.padx
        ymin = self.padx
        ymax = self.padx * 2
        self.image[0:self.pady, xmin:xmax] = numpy.flipud(self.image[ymin:ymax, xmin:xmax])
        # Reflect/flip top edge
        ymin = self.ny - self.pady*2
        ymax = self.ny - self.pady
        self.image[(self.ny - self.pady):self.ny, xmin:xmax] = numpy.flipud(self.image[ymin:ymax, xmin:xmax])
        # I should interpolate over the corners, but .. todo.         
        return
        
        
    def calcFft(self):
        """Calculate the 2d FFT of the image (self.fimage).
        If 'shift', adds a shift to move the small spatial scales to the
        center of the FFT image to self.fimage. Also calculates the frequencies. """
        # Generate the FFT (note, scipy places 0 - largest spatial scale frequencies - at corners)
        self.fimage = fftpack.fft2(self.image)
        if self.shift:
            # Shift the FFT to put the largest spatial scale frequencies at center
            self.fimage = fftpack.fftshift(self.fimage)
        # Note, these frequencies follow unshifted order (0= first, largest spatial scale (with positive freq)). 
        self.xfreq = fftpack.fftfreq(self.nx, 1.0)
        self.yfreq = fftpack.fftfreq(self.ny, 1.0)
        self.xfreqscale = self.xfreq[1] - self.xfreq[0]
        self.yfreqscale = self.yfreq[1] - self.yfreq[0]
        return

    def calcPsd2d(self):
        """Calculate the 2d power and phase spectrum of the image.
        If 'shift', shifts small frequencies to the edges
        (just carries through for PSD, but changes calculation of phase)."""
        try:
            self.fimage
        except AttributeError:
            self.calcFft()
        # Calculate 2d power spectrum.
        # psd = <| R(u,v)^2 + I(u,v)^2| >
        self.psd2d = numpy.absolute(self.fimage)**2.0
        # phase spectrum
        # phase = arctan(I(u,v) / R(u,v))
        if self.shift:
            self.phasespec = numpy.arctan2(fftpack.ifftshift(self.fimage).imag,
                                           fftpack.ifftshift(self.fimage).real)
        else:
            self.phasespec = numpy.arctan2(self.fimage.imag, self.fimage.real)
        return

    def calcPsd1d(self, min_npix=3., min_dr=1.0):
        """Calculate the 1-D power spectrum. The 'tricky' part here is determining spatial scaling.
        At least npix pixels will be included in each radial bin, in frequency space, and the minimum
        spacing between bins will be minthresh pixels. This means the 'central' bins (large pixel scales)
        could have larger steps in the 1dPSD. """
        # Calculate 1d power spectrum                
        #  - uses shifted PSD so that can create radial bins from center, with largest scales at center.
        # Calculate all the radius values for all pixels. These are still in frequency space. 
        rvals = numpy.hypot((self.yy-self.ycen), (self.xx-self.xcen)) + 0.5
        # Sort the PSD2d by the radius values and make flattened representations of these.
        idx = numpy.argsort(rvals.flatten())
        if self.shift:
            dvals = self.psd2d.flatten()[idx]
        else:
            dvals = (fftpack.fftshift(self.psd2d)).flatten()[idx]
        rvals = rvals.flatten()[idx]
        # Set up bins uniform in min_dr pix per bin, but don't want to subdivide center with too
        #  few pixels, so rebin if needed if fall below min_npix. 
        self.min_dr = min_dr
        self.min_npix = min_npix
        rbins = numpy.arange(0, rvals.max()+min_dr*2., min_dr)
        pixperbin = 2.0*numpy.pi*(rbins[1:]**2 - rbins[:-1]**2)
        update = [rbins[0],]
        r0 = rbins[0]
        while r0 < rbins[numpy.where(pixperbin>min_npix)].min():
            r1 = min_npix / 2.0 / numpy.pi + r0**2
            update.append(r1)        
            r0 = r1
        update = numpy.array(update, 'float')
        rbins = numpy.concatenate((update, rbins[numpy.where(rbins>update.max())]))
        # Calculate how many data points are actually present in each radius bin (for weighting)
        nvals = numpy.histogram(rvals, bins=rbins)[0]
        # Calculate the value of the image in each radial bin (weighted by the # of pixels in each bin)
        rprof = numpy.histogram(rvals, bins=rbins, weights=dvals)[0]
        rprof = numpy.where(nvals>0, rprof/nvals, numpy.nan)
        # Calculate the central radius values used in the histograms (note this is still in frequency). 
        rcenters =  (rbins[1:] + rbins[:-1])/2.0
        # And calculate the relevant frequencies (using self.xfreq/yfreq & the corresponding pixel values).
        self.rfreq = rcenters * self.xfreqscale
        # Set the value of the 1d psd, and interpolate over any nans (where there were no pixels).
        self.psd1d = numpy.interp(rcenters, rcenters[rprof==rprof], rprof[rprof==rprof])
        # Scale rcenters to 'original pixels' scale (ie. in the FFT space, pixels are scaled x_fft = 1/(x_pix*2pi)
        #   but must also account for overall size of image 
        self.psdx = 1/(rcenters*2.0*numpy.pi) * numpy.sqrt(self.nx*self.ny)
        return

    def calcAcovf2d(self):
        """Calculate the 2d auto covariance function. """
        # See Wiener-Kinchine theorem
        if self.shift:
            # Note, the ACovF needs the unshifted 2d PSD for inverse FFT, so unshift.
            #  Then shift back again. 
            self.acovf = fftpack.fftshift(fftpack.ifft2(fftpack.ifftshift(self.psd2d)))
        else:
            self.acovf = fftpack.ifft2(self.psd2d)
        return

    def calcAcovf1d(self, min_npix=3, min_dr=1.0):
        """Calculate the 1d average of the ACovF. """
        # Calculate all the radius values for all pixels. These are actually in 'pixel' space
        #    (as ACovF is FFT of PSD). 
        rvals = numpy.hypot((self.xx-self.xcen), (self.yy-self.ycen)) + 0.5
        # Sort the ACovF2d by the radius values and make flattened representations of these.
        idx = numpy.argsort(rvals.flatten())
        if self.shift:
            dvals = self.acovf.flatten()[idx].real
        else:
            dvals = (fftpack.fftshift(self.acovf)).flatten()[idx].real
        rvals = rvals.flatten()[idx]
        # Set up bins uniform in npix per bin, for the 1d ACovF calculation (like with PSD case). 
        #  but want to subdivide outer parts of image too much, so use a minimum of min_dr pix
        self.min_dr = min_dr
        self.min_npix = min_npix
        rbins = numpy.arange(0, rvals.max()+min_dr, min_dr)
        pixperbin = 2.0*numpy.pi*(rbins[1:]**2 - rbins[:-1]**2)
        update = [rbins[0],]
        r0 = rbins[0]
        while r0 < rbins[numpy.where(pixperbin>min_npix)].min():
            r1 = npix / 2.0 / numpy.pi + r0**2
            update.append(r1)        
            r0 = r1
        update = numpy.array(update, 'float')
        rbins = numpy.concatenate((update, rbins[numpy.where(rbins>update.max())]))
        # Calculate how many data points are actually present in each radius bin (for weighting)
        nvals = numpy.histogram(rvals, bins=rbins)[0]
        # Calculate the value of the image in each radial bin (weighted by the # of pixels in each bin)
        rprof = numpy.histogram(rvals, bins=rbins, weights=dvals)[0] / nvals
        # Calculate the central radius values used in the histograms
        rcenters =  (rbins[1:] + rbins[:-1])/2.0
        # Set the value of the 1d ACovF, and interpolate over any nans (where there were no pixels)
        self.acovf1d = numpy.interp(rcenters, rcenters[rprof==rprof], rprof[rprof==rprof])
        self.acovfx = rcenters
        return

    def calcSf(self):
        """Calculate the structure function from the 1d ACovF, discounting zeroPadding region.
        Structure function is calculated as SF = sqrt(ACovF(0) - ACovF), as described in powers_qs doc."""
        self.sfx = numpy.arange(0, (numpy.sqrt((self.nx/2.0 - self.padx)**2 + (self.ny/2.0 - self.pady)**2)), 1.0)
        self.sf = numpy.interp(self.sfx, self.acovfx, numpy.sqrt(self.acovf1d[0] - self.acovf1d))        
        return

    def calcAll(self, min_npix=2, min_dr=1.):
        self.calcFft()
        self.calcPsd2d()
        self.calcPsd1d(min_npix=min_npix, min_dr=min_dr)
        self.calcAcovf2d()
        self.calcAcovf1d(min_npix=min_npix, min_dr=min_dr)
        self.calcSf()
        return

    def _makeRandomPhases(self, seed=None):
        # Generate random phases for construction of real image. 
        # There should be some symmetry here - real images should have negative and positive frequency conjugate symmetry. 
        # It's pretty close to upper right quarter / lower left quarter are negative, mirror copies (same for upper left/lower right quarters).
        # This is true for unshifted images, not for shifted images, so by copying/mirroring  these quarters, must reconstruct image accounting for this.
        # The first column/row are not generally part of this symmetry, as these are the 0 frequency (in some component) values.
        # To preserve the mean of the image, the phase for the 0,0 frequency should be 0 or pi. 
        if seed != None:
            numpy.random.seed(seed)
        # Generate random phases (uniform -360 to 360) for whole image -- will write over some sections. 
        self.phasespecI = numpy.random.uniform(low=-numpy.pi, high=numpy.pi, size=[self.ny, self.nx])
        # Set 0,0 to 0 to preserve overall mean value of image. (= 0 means no value into imaginary component of FFT of image). 
        self.phasespecI[0][0] = 0
        # Mirror upper left corner into lower right corner. 
        self.phasespecI[self.xcen+1:self.nx, self.ycen+1:self.ny] = -1.*numpy.fliplr(numpy.flipud(self.phasespecI[1:self.xcen, 1:self.ycen]))
        # Mirror upper right corner into lower left corner. 
        self.phasespecI[1:self.xcen, self.ycen+1:self.ny] = -1.*numpy.fliplr(numpy.flipud(self.phasespecI[self.xcen+1:self.nx, 1:self.ycen]))
        return

    def invertFft(self, useI=False, verbose=False):
        """Convert the 2d FFT into an image (imageI)."""
        # Checking this process with a simple (non-noisy) image shows that it will result in errors on the
        #  level of 1e-15 counts (in an original image with min/max scale of 1.0).
        if useI:
            fimage = self.fimageI
        else:
            fimage = self.fimage
        if self.shift:
            self.imageI = fftpack.ifft2(fftpack.ifftshift(fimage))
        else:
            self.imageI = fftpack.ifft2(fimage)
        if self.imageI.imag.max() < 1e-14:
            if verbose:  print "Inverse FFT created only small imaginary portion - discarding."
            self.imageI = self.imageI.real
        return

    def invertPsd2d(self, useI=False, usePhasespec=True, seed=None):
        """Convert the 2d PSD and phase spec into an FFT image (FftI). """
        # The PHASEs of the FFT are encoded in the phasespec ('where things are')
        # The AMPLITUDE of the FFT is encoded in the PSD ('how bright things are' .. also radial scale)
        # amp = sqrt(| R(uv)^2 + I(u,v)^2|) == length of 'z' (z = x + iy, in fourier image)
        # phase = arctan(y/x)
        if useI:
            psd2d = self.psd2dI
            phasespec = self.phasespecI
        else:
            psd2d = self.psd2d
            phasespec = self.phasespec
        if self.shift:
            amp = numpy.sqrt(fftpack.ifftshift(psd2d))
        else:
            amp = numpy.sqrt(psd2d)
        # Can override using own phasespec for random (if coming in here to use 2d PSD for image)
        if not(usePhasespec):
            self._makeRandomPhases(seed=seed)
            phasespec = self.phasespecI
        # Shift doesn't matter for phases, because 'unshifted' it above, before calculating phase.
        x = numpy.cos(phasespec) * amp
        y = numpy.sin(phasespec) * amp
        self.fimageI = x + 1j*y
        if self.shift:
            self.fimageI = fftpack.fftshift(self.fimageI)
        return
    
    def invertPsd1d(self, psdx=None, psd1d=None, phasespec=None, seed=None):
        """Convert a 1d PSD, generate a phase spectrum (or use user-supplied values) into a 2d PSD (psd2dI)."""
        # Converting the 2d PSD into a 1d PSD is a lossy process, and then there is additional randomness
        #  added when the phase spectrum is not the same as the original phase spectrum, so this may or may not
        #  look that much like the original image (unless you keep the phases, doesn't look like image). 
        # 'Swing' the 1d PSD across the whole fov.         
        if psd1d == None:
            psd1d = self.psd1d
            xr = self.rfreq / self.xfreqscale
        else:
            if psdx == None:
                xr = numpy.arange(0, len(psd1d), 1.0)
        # Resample into even bins (definitely necessarily if using psd1d).
        xrange = numpy.arange(0, numpy.sqrt(self.xcen**2 + self.ycen**2)+1.0, 1.0)
        psd1d = numpy.interp(xrange, xr, psd1d, right=psd1d[len(psd1d)-1])
        # Calculate radii - distance from center.
        rad = numpy.hypot((self.yy-self.ycen), (self.xx-self.xcen))
        # Calculate the PSD2D from the 1d value.
        self.psd2dI = numpy.interp(rad.flatten(), xrange, psd1d)
        self.psd2dI = self.psd2dI.reshape(self.ny, self.nx)
        if phasespec == None:
            self._makeRandomPhases(seed=seed)
        else:
            self.phasespecI = phasespec
        if not(self.shift):
            # The 1d PSD is centered, so the 2d PSD will be 'shifted' here. 
            self.psd2dI = fftpack.ifftshift(self.psd2dI)
        return

    def invertAcovf2d(self, useI=False, usePhasespec=True, seed=None):
        """Convert the 2d ACovF into a 2d PSD (psd2dI). """
        if useI:
            acovf = self.acovfI
        else:
            acovf = self.acovf
            self.phasespecI = self.phasespec
        # Let user override phasespec if want to use real 2d ACovF but random phases.
        if not(usePhasespec):
            self._makeRandomPhases(seed=seed) 
        # Calculate the 2dPSD from the ACovF. 
        # Note that if the ACovF has values reaching all the way to the edges, this will
        # induce noise into the PSD2d (just as the discontinuity would induce noise into the FFT of an image).
        if self.shift:
            self.psd2dI = fftpack.ifftshift(fftpack.fft2(fftpack.fftshift(acovf)))
        else:
            self.psd2dI = fftpack.fft2(acovf)
        # PSD2d should be entirely real and positive (PSD2d = |R(uv,)**2 + I(u,v)**2|
        #print 'PSD real limits', self.psd2dI.real.min(), self.psd2dI.real.max()
        #print 'PSD imaginary limits', self.psd2dI.imag.min(), self.psd2dI.imag.max()
        # Okay, I admit - this next line is a bit of a hack, but it does seem to work. 
        self.psd2dI = numpy.sqrt(numpy.abs(self.psd2dI)**2)
        return

    def invertAcovf1d(self, acovfx=None, acovf1d=None, phasespec=None, seed=None, reduceNoise=True):
        """Convert a 1d ACovF into a 2d ACovF (acovfI). """
        # 'Swing' the 1d ACovF across the whole fov.         
        if acovf1d == None:
            # Use self-set values.
            acovf1d = self.acovf1d
            xr = self.acovfx            
        else:
            if acovfx == None:
                # If not given r values for acovf1d, assume they are even, 1 pixel spacing.
                xr = numpy.arange(0, len(acovf1d), 1.0)
        # It turns out that if acovf1d extends past min(xcen/ycen) -- i.e., extends into the 'corners' of the image -- then
        #  we get noise in the calculation of the PSD2d, which is an FFT of the ACovF2d. This noise results in linear features in the
        #  reconstructed image, which is bad. 
        # So, to avoid this, we must make the acovf1d a single value after rmax = min(xcen/ycen). Could be zero, with a rolloff .. but it 
        #  seems like this is unnecessary. So, unless reduceNoise is set to False, just truncate acovf1d at rmax and extend with constant value.
        if (xr.max() > (min(self.xcen, self.ycen))) & reduceNoise:
            rmax = min(self.xcen, self.ycen)
            condition = (xr < rmax)
            acovf1d = numpy.interp(xr[condition], xr, acovf1d)
            xr = xr[condition]
            #rin = rmax - rmax/10.0
            #filter = numpy.where(xr > rmax, 0, numpy.where(xr < rin, 1, 0.5 - 0.5*numpy.cos(numpy.pi*(1-(xr-rin)/(rmax-rin)))))
            #acovf1d *= filter
        # Calculate radii - distance from center for all pixels in image.
        rad = numpy.hypot((self.yy-self.ycen), (self.xx-self.xcen))
        # Calculate the ACovF2D from the 1d value.
        self.acovfI = numpy.interp(rad.flatten(), xr, acovf1d)
        self.acovfI = self.acovfI.reshape(self.ny, self.nx)
        if phasespec == None:
            self._makeRandomPhases(seed=seed)
        else:
            self.phasespecI = phasespec
        if not(self.shift):
            self.acovfI = fftpack.fftshift(self.acovfI)
        return

    def invertSf(self, sfx, sf):
        """Convert a structure function (sf) with r coordinates sfx into the 1-d ACovF.
        Note definition of structure function here: SF = \sqrt(ACovF(0) - ACovF), as in powers_qs.pdf doc."""
        self.sfx = sfx
        self.sf = sf
        # Calcluate 1-d ACovF
        self.acovf1d = sf[len(sfx)-1]**2 - sf**2
        self.acovfx = numpy.copy(sfx)
        # Now you can call invertAcovf1d directly.  (invertAcovf1d(phasespec/seed))
        return

    def makeImageFromSf(self, sfx, sf, reduceNoise=True):
        """Invert SF all the way back to ImageI."""
        self.invertSf(sfx, sf)
        self.invertAcovf1d(reduceNoise=reduceNoise)
        self.invertAcovf2d(useI=True)
        self.invertPsd2d(useI=True)
        self.invertFft(useI=True)
        self.xfreq = fftpack.fftfreq(self.nx, 1.0)
        self.yfreq = fftpack.fftfreq(self.ny, 1.0)
        return
