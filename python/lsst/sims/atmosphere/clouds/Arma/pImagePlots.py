import numpy
import pylab
from scipy import fftpack
# For plotting an image in log scale
from matplotlib.colors import LogNorm
from pImage import PImage

class PImagePlots(PImage):
    """Inherits from PImage, which does all the 'math' work, but adds plotting capabilities.
    These are separate to make it easier to use PImage from a non-interactive terminal. """
    
    def showXSlice(self, y=None, source='image'):
        """Plot a 1d slice through the image (or fft or psd), at y (if None, defaults to center)."""
        if y == None:
            y = round(self.ny / 2.0)
        if source == 'image':
            x = numpy.arange(0, self.nx)
            pylab.plot(x, self.image[y][:])
        elif source == 'fft':
            # have to adjust y for the fact that fimage has 'shifted' axis versus fimage2
            #  (but xfreq goes with fimage)
            y = self.ny/2.0 - y
            pylab.plot(self.xfreq, fftpack.ifftshift(self.fimage)[y][:].real, 'k.')
        elif source == 'psd':
            y = self.ny/2.0 - y
            pylab.plot(self.xfreq, fftpack.ifftshift(self.psd2d)[y][:], 'k.')
        elif source == 'acovf':
            x = numpy.arange(0, self.nx)
            pylab.plot(x, self.acovf[y][:])
        else:
            raise Exception('Source must be one of image/fft/psd/acovf')
        return

    def showYSlice(self, x=None):
        """Plot a 1d slice through the image, at x (if None, defaults to center)."""
        if x == None:
            x = round(self.nx/2.0)
        if source == 'image':
            y = numpy.arange(0, self.ny)
            pylab.plot(y, self.image[:][x])
        elif source == 'fft':
            pylab.plot(self.yfreq, fftpack.ifftshift(self.fimage)[:][x].real, 'k.')
        elif source == 'psd':
            pylab.plot(self.yfreq, fftpack.ifftshift(self.psd2d)[:][x], 'k.')
        elif source == 'acovf':
            y = numpy.arange(0, self.ny)
            pylab.plot(y, self.acovf[:][x])
        else:
            raise Exception('Source must be one of image/fft/psd/acovf')
        return

    def showImage(self, xlim=None, ylim=None, clims=None, cmap=None, copy=False, stats=True):
        if copy:
            # useful if you're going to apply a hanning filter or something later, but
            # want to keep an imge of the original
            image = numpy.copy(self.image)
        else:
            image = self.image
        pylab.figure()
        pylab.title('Image')
        if xlim == None:
            x0 = 0
            x1 = self.nx
        else:
            x0 = xlim[0]
            x1 = xlim[1]
        if ylim == None:
            y0 = 0
            y1 = self.ny
        else:
            y0 = ylim[0]
            y1 = ylim[1]
        if clims == None:
            pylab.imshow(image, origin='lower', cmap=cmap)
        else:
            pylab.imshow(image, origin='lower', vmin=clims[0], vmax=clims[1], cmap=cmap)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        cb = pylab.colorbar()
        clims = cb.get_clim()
        pylab.xlim(x0, x1)
        pylab.ylim(y0, y1)
        if stats:
            statstxt = 'Mean/Stdev/Min/Max:\n %.2f/%.2f/%.2f/%.2f' %(numpy.mean(image), numpy.std(image), image.min(), image.max())
            pylab.figtext(0.75, 0.03, statstxt)
        return clims

    def showImageI(self, xlim=None, ylim=None, clims=None, cmap=None):
        pylab.figure()
        pylab.title('Reconstructed Image')
        if xlim == None:
            x0 = 0
            x1 = self.nx
        else:
            x0 = xlim[0]
            x1 = xlim[1]
        if ylim == None:
            y0 = 0
            y1 = self.ny
        else:
            y0 = ylim[0]
            y1 = ylim[1]
        if clims == None:
            pylab.imshow(self.imageI.real, origin='lower', cmap=cmap)
        else:
            pylab.imshow(self.imageI.real, origin='lower', vmin=clims[0], vmax=clims[1], cmap=cmap)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        cb = pylab.colorbar()
        pylab.xlim(x0, x1)
        pylab.ylim(y0, y1)
        return

    def showImageAndImageI(self, clims=None, climsI=None, title=None):
        shrinkratio=0.6
        pylab.figure()
        if title!=None:
            pylab.suptitle(title)
        pylab.subplot(121)
        pylab.title('Image')
        if clims == None:
            pylab.imshow(self.image, origin='lower')
        else:
            pylab.imshow(self.image, origin='lower', vmin=clims[0], vmax=clims[1])
        pylab.colorbar(shrink=shrinkratio)
        pylab.xticks(rotation=45)
        pylab.subplot(122)
        pylab.title('Reconstructed image')
        if climsI == None:
            pylab.imshow(self.imageI.real, origin='lower')
        else:
            pylab.imshow(self.imageI.real, origin='lower', vmin=climsI[0], vmax=climsI[1])
        pylab.colorbar(shrink=shrinkratio)
        pylab.xticks(rotation=45)
        return
        
    def showFft(self, real=True, imag=True, clims=None, log=False):
        if ((real == True) & (imag==True)):
            p = 2
            shrinkratio=0.7
        else:
            p = 1
            shrinkratio=0.9
        pylab.figure()        
        if log:
            from matplotlib.colors import LogNorm
            norml = LogNorm()
            if clims!=None:
                if (clims[0] <0) | (clims[1]<0):
                    print 'Using a log scale and clims < 0 does not work ...'
                norml = LogNorm(vmin=clims[0], vmax=clims[1])
        else:
            from matplotlib.colors import Normalize
            norml = Normalize()
            if clims!=None:
                norml = Normalize(vmin=clims[0], vmax=clims[1])
        if real:
            pylab.subplot(1,p,1)
            pylab.title('Real FFT')
            pylab.imshow(self.fimage.real, origin='lower', norm=norml,
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
            pylab.xlabel('u')
            pylab.ylabel('v')
            cb = pylab.colorbar(shrink=shrinkratio)
        if imag:
            pylab.subplot(1,p,p)
            pylab.title('Imaginary FFT')
            pylab.imshow(self.fimage.imag, origin='lower', norm=norml,
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
            pylab.xlabel('u')
            pylab.ylabel('v')
            cb = pylab.colorbar(shrink=shrinkratio)
        return

    def showFftI(self, real=True, imag=True, clims=None, log=False):
        if ((real == True) & (imag==True)):
            p = 2
            shrinkratio=0.7
        else:
            p = 1
            shrinkratio=0.9
        pylab.figure()        
        if log:
            from matplotlib.colors import LogNorm
            norml = LogNorm()
            if clims!=None:
                if (clims[0] <0) | (clims[1]<0):
                    print 'Using a log scale and clims < 0 does not work ...'
                norml = LogNorm(vmin=clims[0], vmax=clims[1])
        else:
            from matplotlib.colors import Normalize
            norml = Normalize()
            if clims!=None:
                norml = Normalize(vmin=clims[0], vmax=clims[1])
        if real:
            pylab.subplot(1,p,1)
            pylab.title('Reconstructed Real FFT')
            pylab.imshow(self.fimageI.real, origin='lower',
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
            pylab.xlabel('u')
            pylab.ylabel('v')
            cb = pylab.colorbar(shrink=shrinkratio)
        if imag:
            pylab.subplot(1,p,p)
            pylab.title('Reconstructed Imaginary FFT')
            if clims == None:
                pylab.imshow(self.fimageI.imag, origin='lower',
                             extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
            else:
                pylab.imshow(self.fimageI.imag, origin='lower', vmin=clims[0], vmax=clims[1],
                             extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
            pylab.xlabel('u')
            pylab.ylabel('v')
            cb = pylab.colorbar(shrink=shrinkratio)
        return

    def showPsd2d(self, log=True):
        pylab.figure()
        pylab.title('2D Power Spectrum')
        vmax = self.psd2d.max()
        vmin = max(vmax-10e20, self.psd2d.min())
        vmin = max(vmin, 0.0001)
        if log==True:
            from matplotlib.colors import LogNorm
            norml = LogNorm(vmin=vmin, vmax=vmax)
            pylab.imshow(self.psd2d, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(), \
                                                              self.yfreq.min(), self.yfreq.max()], \
                                                              norm=norml)
        else:
            pylab.imshow(self.psd2d, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(),
                                                              self.yfreq.min(), self.yfreq.max()])
        cb = pylab.colorbar()
        pylab.xlabel('u')
        pylab.ylabel('v')
        return

    def showPhasesI(self):
        pylab.figure()
        pylab.title('PSD Phases')
        pylab.imshow(self.phasespecI, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(),
                                                             self.yfreq.min(), self.yfreq.max()])
        cb = pylab.colorbar()
        pylab.xlabel('u')
        pylab.ylabel('v')
        return

    def showPhases(self):
        pylab.figure()
        pylab.title('PSD Phases')
        pylab.imshow(self.phasespec, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(),
                                                             self.yfreq.min(), self.yfreq.max()])
        cb = pylab.colorbar()
        pylab.xlabel('u')
        pylab.ylabel('v')
        return

    def showPsd2dI(self, log=True):
        pylab.figure()
        pylab.title('Reconstructed 2D Power Spectrum')
        vmax = self.psd2dI.max()
        vmin = max(vmax-10e20, self.psd2dI.min())
        vmin = max(vmin, 0.0001)
        if log==True:
            from matplotlib.colors import LogNorm
            norml = LogNorm(vmin=vmin, vmax=vmax)
            pylab.imshow(self.psd2dI, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(), \
                                                              self.yfreq.min(), self.yfreq.max()], norm=norml)
        else:
            pylab.imshow(self.psd2dI, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(),
                                                              self.yfreq.min(), self.yfreq.max()])
        cb = pylab.colorbar()
        pylab.xlabel('u')
        pylab.ylabel('v')
        return


    def showPsd1d(self, comparison=None, linear=False, legendlabels=['Image', 'Comparison']):
        pylab.figure()
        # limit y scale
        ymax = self.psd1d.max()
        if comparison != None:
            ymax = max(ymax, comparison.psd1d.max())
        ymin = max(ymax - 10e15, self.psd1d.min())
        if comparison != None:
            ymin = min(self.psd1d.min(), comparison.psd1d.min())
            ymin = max(ymax - 10e15, ymin)            
        ymin = max(ymin, 1e-12)
        xmax = max(self.xcen, self.ycen) - min(self.padx, self.pady)
        xmin = 0.5
        pylab.subplots_adjust(left=0.15, right=0.9, wspace=0.35, hspace=0.2)
        pylab.subplot(121)
        pylab.plot(self.rfreq, self.psd1d, 'b-')
        if comparison != None:
            pylab.plot(comparison.rfreq, comparison.psd1d, 'r-')            
        pylab.ylim(ymin, ymax)
        pylab.yscale('log', subsy=[2,3,4,5,6,7,8,9])
        if linear:
            pylab.yscale('linear')
        pylab.yticks(fontsize='smaller')
        pylab.xlabel('Frequency')
        pylab.ylabel('1-D Power Spectrum')
        pylab.subplot(122)
        pylab.plot(self.psdx, self.psd1d, 'b-', label=legendlabels[0])
        if comparison != None:
            pylab.plot(comparison.psdx, comparison.psd1d, 'r-', label=legendlabels[1])
        pylab.xlim(xmin, xmax)
        pylab.ylim(ymin, ymax)
        pylab.yscale('log', subsy=[2,3,4,5,6,7,8,9])
        pylab.xscale('log', subsx=[2,3,4,5,6,7,8,9])
        if linear:
            pylab.xscale('linear')
            pylab.yscale('linear')
        pylab.yticks(fontsize='smaller')
        pylab.xlabel('Spatial scale (pix)')        
        pylab.grid()
        if comparison!=None:
            pylab.legend(numpoints=1, fancybox=True, loc='lower right', fontsize='smaller')
        pylab.suptitle('1D Power Spectrum: min_npix %.0f, min_dr %.2f' %(self.min_npix, self.min_dr))
        return

    def showAcovf2d(self, real=True, imag=False, clims=None, log=False):
        if ((real == True) & (imag==True)):
            p = 2
            shrinkratio=0.7
        else:
            p = 1
            shrinkratio=0.9
        pylab.figure()
        if clims != None:
            vmin = clims[0]
            vmax = clims[1]
        if real:
            if clims == None:
                vmax = self.acovf.max().real
                vmin = max(vmax-10e20, self.acovf.min().real)
                if log:
                    vmax = max(vmax, 0.00001)
                    vmin = max(vmin, 0.000001)            
            pylab.subplot(1,p,1)
            pylab.title('Real ACovF')
            if log==True:
                from matplotlib.colors import LogNorm
                norml = LogNorm(vmin=vmin, vmax=vmax)
                pylab.imshow(self.acovf.real, origin='lower', norm=norml, extent=[0-self.xcen, self.nx-self.xcen,
                                                                                0-self.ycen, self.ny-self.ycen])
            else:
                pylab.imshow(self.acovf.real, origin='lower', vmin=vmin, vmax=vmax, 
                             extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
            pylab.xlabel('X')
            pylab.ylabel('Y')
            cb = pylab.colorbar(shrink=shrinkratio)
        if imag:
            if clims == None:
                vmax = self.acovf.max().imag
                vmin = max(vmax-10e20, self.acovf.min().imag)
                if log:
                    vmax = max(vmax, 0.00001)
                    vmin = max(vmin, 0.000001)            
            pylab.subplot(1,p,p)
            pylab.title('Imaginary ACovF')
            if log==True:
                from matplotlib.colors import LogNorm
                norml = LogNorm(vmin=vmin, vmax=vmax)
                pylab.imshow(self.acovf.imag, origin='lower', norm=norml, extent=[-self.xcen, self.nx-self.xcen, 
                                                                                 -self.ycen, self.ny-self.ycen])
            else:
                pylab.imshow(self.acovf.imag, origin='lower', vmin=vmin, vmax=vmax, 
                             extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
            pylab.xlabel('X')
            pylab.ylabel('Y')
            cb = pylab.colorbar(shrink=shrinkratio)
        return

    def showAcovf2dI(self, real=True, imag=True, clims=None, log=False):
        if ((real == True) & (imag==True)):
            p = 2
            shrinkratio=0.7
        else:
            p = 1
            shrinkratio=0.9
        pylab.figure()
        if clims != None:
            vmin = clims[0]
            vmax = clims[1]
        if real:
            if clims == None:
                vmax = self.acovfI.max().real
                vmin = max(vmax-10e20, self.acovfI.min().real)
                if log:
                    vmax = max(vmax, 0.00001)
                    vmin = max(vmin, 0.000001)            
            pylab.subplot(1,p,1)
            pylab.title('Reconstructed Real ACovF')
            if log==True:
                from matplotlib.colors import LogNorm
                norml = LogNorm(vmin=vmin, vmax=vmax)
                pylab.imshow(self.acovfI.real, origin='lower', norm=norml,
                             extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
            else:
                pylab.imshow(self.acovfI.real, origin='lower', vmin=vmin, vmax=vmax,
                             extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
            pylab.xlabel('X')
            pylab.ylabel('Y')
            cb = pylab.colorbar(shrink=shrinkratio)
        if imag:
            if clims == None:
                vmax = self.acovfI.max().imag
                vmin = max(vmax-10e20, self.acovfI.min().imag)
                if log:
                    vmax = max(vmax, 0.00001)
                    vmin = max(vmin, 0.000001)            
            pylab.subplot(1,p,p)
            pylab.title('Reconstructed Imag. ACovF')
            if log==True:
                from matplotlib.colors import LogNorm
                norml = LogNorm(vmin=vmin, vmax=vmax)
                pylab.imshow(self.acovfI.imag, origin='lower', norm=norml,
                             extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
            else:
                pylab.imshow(self.acovfI.imag, origin='lower', vmin=vmin, vmax=vmax,
                             extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
            pylab.xlabel('X')
            pylab.ylabel('Y')
            cb = pylab.colorbar(shrink=shrinkratio)
        return

    def showAcovf1d(self, linear=False, comparison=None, legendlabels=['Image', 'Comparison']):        
        pylab.figure()
        try:
            pylab.title('1D ACovF: min_npix %.0f, min_dr %.2f' %(self.min_npix, self.min_dr))
        except:
            pylab.title('1D ACovF')
        maxscale_image = (numpy.sqrt((self.nx/2.0 - self.padx)**2 + (self.ny/2.0 - self.pady)**2))
        condition = (self.acovfx <= maxscale_image)
        pylab.plot(self.acovfx[condition], self.acovf1d[condition], 'b-', label=legendlabels[0])
        if comparison != None:
            pylab.plot(comparison.acovfx[condition], comparison.acovf1d[condition], 'r-', label=legendlabels[1])
        if linear:
            pylab.yscale('linear')
        else:
            pylab.yscale('log', subsy=[2,3,4,5,6,7,8,9])
        #pylab.xscale('log', subsx=[2,3,4,5,6,7,8,9])
        pylab.xscale('linear')
        pylab.xlabel('Spatial Scale (pix)')
        pylab.ylabel('1d ACovF')
        pylab.grid()
        if comparison!=None:
            pylab.legend(numpoints=1, fancybox=True, loc='upper right', fontsize='smaller')
        return

    def showSf(self, linear=False, comparison=None, legendlabels=['Image', 'Comparison']):        
        pylab.figure()
        pylab.title('1-d SF: min_npix %.0f, min_dr %.2f' %(self.min_npix, self.min_dr))
        pylab.plot(self.sfx, self.sf, 'b-', label=legendlabels[0])
        if comparison != None:
            pylab.plot(comparison.sfx, comparison.sf, 'r-', label=legendlabels[1])
        if linear:
            pylab.yscale('linear')
        else:
            pylab.yscale('log', subsy=[2,3,4,5,6,7,8,9])
        #pylab.xscale('log', subsx=[2,3,4,5,6,7,8,9])
        pylab.xscale('linear')
        pylab.xlabel('Spatial Scale (pix)')
        pylab.ylabel('1-d SF')
        pylab.grid()
        if comparison!=None:
            pylab.legend(numpoints=1, fancybox=True, loc='lower right', fontsize='smaller')
        return

    def plotAll(self, title=None):
        pylab.figure()        
        pylab.subplots_adjust(left=0.1, right=0.97, wspace=0.45, hspace=0.2)
        shrinkratio = 0.7
        #
        ax1 = pylab.subplot2grid((2,3),(0,0))
        pylab.imshow(self.image, origin='lower')
        pylab.xticks(rotation=45)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        cb = pylab.colorbar(shrink=shrinkratio)
        clims = cb.get_clim()
        pylab.title('Image', fontsize=12)
        #
        ax2 = pylab.subplot2grid((2,3), (1,0))
        pylab.imshow(self.fimage.real, origin='lower', vmin=clims[0], vmax=clims[1],
                     extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
        pylab.xticks(rotation=45)
        pylab.xlabel('u')
        pylab.ylabel('v')
        cb = pylab.colorbar(shrink=shrinkratio)
        pylab.title('Real FFT', fontsize=12)
        #
        ax3 = pylab.subplot2grid((2,3), (0, 1))
        pylab.imshow(self.acovf.real, origin='lower',
                     extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
        pylab.xticks(rotation=45)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        cb = pylab.colorbar(shrink=shrinkratio)
        pylab.title('2D ACovF', fontsize=12)
        #
        ax4 = pylab.subplot2grid((2,3), (1,1))
        vmax = self.psd2d.max()
        vmin = max(vmax - 10.e20, self.psd2d.min())
        vmin = max(vmin, 0.01)
        norml = LogNorm(vmin=vmin, vmax=vmax)
        pylab.imshow(self.psd2d, origin='lower', extent=[self.xfreq.min(), self.xfreq.max(), \
                                                          self.yfreq.min(), self.yfreq.max()], \
                                                          norm=norml)
        cb = pylab.colorbar(shrink=shrinkratio)
        pylab.xticks(rotation=45)
        pylab.xlabel('u')
        pylab.ylabel('v')
        pylab.title('2D PSD', fontsize=12)
        #
        ax5 = pylab.subplot2grid((2,3),(0,2))
        pylab.plot(self.acovfx, self.acovf1d, 'k-')
        #pylab.xscale('linear')
        pylab.xscale('log', subsx=[2,3,4,5,6,7,8,9])
        pylab.yscale('log', subsy=[2,3,4,5,6,7,8,9])
        pylab.xticks(rotation=45)
        pylab.title('1D ACovF', fontsize=12)
        #
        ax6 = pylab.subplot2grid((2,3), (1,2))
        pylab.plot(self.psdx, self.psd1d, 'k-')
        pylab.yscale('log', subsy=[2,3,4,5,6,7,8,9])
        pylab.xscale('log', subsx=[2,3,4,5,6,7,8,9])
        #pylab.xticks(rotation=45)
        #pylab.xlabel('Pix')
        pylab.xlabel('1D PSD', fontsize=12)
        if title!=None:
            pylab.suptitle(title, fontsize=14)
        return


    def _colorbarTicks(self, vmin, vmax, steps=4., log=False):
        if not(log):
            stepsize = (vmax-vmin)/float(steps)
            ticks = numpy.arange(vmin, vmax+stepsize, stepsize)            
            cbformat = '%.2f'
            if (vmax-vmin) < 0.01:
                cbformat = '%.1e'
        if log:
            vmintmp = numpy.log10(vmin)
            vmaxtmp = numpy.log10(vmax)
            stepsize = (vmaxtmp - vmintmp)/float(steps)
            ticks = numpy.arange(vmintmp, vmaxtmp+stepsize/2.0, stepsize)
            ticks = numpy.floor(ticks)
            ticks = 10**ticks
            cbformat='%.0e'
        return ticks, cbformat

    def plotMore(self, title=None, useClims=True):
        pylab.figure()        
        pylab.subplots_adjust(left=0.1, right=0.95, bottom=0.1, wspace=0.45, hspace=0.45)
        shrinkratio = 0.7
        #
        ax1 = pylab.subplot2grid((3,3),(0,0))
        pylab.imshow(self.image, origin='lower')
        pylab.xticks(rotation=45, fontsize='smaller')
        pylab.yticks(fontsize='smaller')
        pylab.ylabel('X / Y', fontsize='smaller')
        if useClims:
            clims = [self.image.min(), self.image.max()]
            ticks, cbformat = self._colorbarTicks(clims[0], clims[1])
            cb = pylab.colorbar(shrink=shrinkratio, ticks=ticks, format=cbformat)
        else:
            pylab.colorbar(shrink=shrinkratio)
        pylab.title('Image', fontsize=12)
        #
        ax2 = pylab.subplot2grid((3,3), (0,1))
        if useClims:
            pylab.imshow(self.fimage.real, origin='lower', vmin=clims[0], vmax=clims[1],
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
        else:
            pylab.imshow(self.fimage.real, origin='lower',
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])            
        pylab.xticks(rotation=45, fontsize='smaller')
        pylab.yticks(fontsize='smaller')
        pylab.ylabel('u / v', fontsize='smaller')
        if useClims:
            cb = pylab.colorbar(shrink=shrinkratio, ticks=ticks, format=cbformat)
        else:
            pylab.colorbar(shrink=shrinkratio)
        pylab.title('Real FFT', fontsize=12)
        # 
        ax3 = pylab.subplot2grid((3,3), (0,2))
        if useClims:
            pylab.imshow(self.fimage.imag, origin='lower', vmin=clims[0], vmax=clims[1],
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
        else:
            pylab.imshow(self.fimage.imag, origin='lower', 
                         extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
        pylab.xticks(rotation=45, fontsize='smaller')
        pylab.yticks(fontsize='smaller')
        pylab.ylabel('u / v', fontsize='smaller')
        if useClims:
            cb = pylab.colorbar(shrink=shrinkratio, ticks=ticks, format=cbformat)
        else:
            pylab.colorbar(shrink=shrinkratio)
        pylab.title('Imag FFT', fontsize=12)
        #
        ax4 = pylab.subplot2grid((3,3), (1,0))
        pylab.imshow(self.phasespec, origin='lower', 
                     extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
        pylab.xticks(rotation=45, fontsize='smaller')
        pylab.yticks(fontsize='smaller')
        pylab.ylabel('u / v', fontsize='smaller')
        if useClims:
            ticks, cbformat = self._colorbarTicks(-numpy.pi, numpy.pi)
            cb = pylab.colorbar(shrink=shrinkratio, ticks=ticks, format=cbformat)
        else:
            pylab.colorbar(shrink=shrinkratio)
        pylab.title('Phase Spec.', fontsize=12)
        # 
        ax5 = pylab.subplot2grid((3,3), (1,1))
        vmax = self.psd2d.max()
        vmin = max(vmax - 10.e20, self.psd2d.min())
        vmin = max(vmin, 0.0001)
        norml = LogNorm(vmin=vmin, vmax=vmax)
        pylab.imshow(self.psd2d, origin='lower', norm=norml, 
                     extent=[self.xfreq.min(), self.xfreq.max(), self.yfreq.min(), self.yfreq.max()])
        ticks, cbformat = self._colorbarTicks(vmin, vmax, log=True)
        cb = pylab.colorbar(shrink=shrinkratio, ticks=ticks)
        pylab.xticks(rotation=45, fontsize='smaller')
        pylab.yticks(fontsize='smaller')
        pylab.ylabel('u / v', fontsize='smaller')
        pylab.title('2D PSD', fontsize=12)
        #
        ax6 = pylab.subplot2grid((3,3), (1, 2))
        pylab.imshow(self.acovf.real, origin='lower',
                     extent=[-self.xcen, self.nx-self.xcen, -self.ycen, self.ny-self.ycen])
        pylab.xticks(rotation=45, fontsize='smaller')
        pylab.yticks(fontsize='smaller')
        pylab.ylabel('X / Y', fontsize='smaller')
        if useClims:
            ticks, cbformat = self._colorbarTicks(0, self.acovf.real.max())
            cb = pylab.colorbar(shrink=shrinkratio, ticks=ticks, format=cbformat)
        else:
            pylab.colorbar(shrink=shrinkratio)
        pylab.title('2D ACovF', fontsize=12)
        #
        ax7 = pylab.subplot2grid((3,2), (2,0))
        pylab.plot(self.rfreq, self.psd1d, 'k-')
        pylab.yscale('log')#, subsy=[2,3,4,5,6,7,8,9])
        pylab.xscale('linear')
        pylab.xticks(rotation=45)
        ylim = pylab.ylim()
        ticks, cbformat = self._colorbarTicks(ylim[0], ylim[1], log=True)
        pylab.yticks(ticks)
        pylab.grid()
        pylab.xlabel('u')
        pylab.ylabel('1D PSD', fontsize=12)
        #
        ax8 = pylab.subplot2grid((3,2),(2,1))
        pylab.plot(self.acovfx, self.acovf1d, 'k-')
        pylab.xscale('linear')
        #pylab.xscale('log', subsx=[2,3,4,5,6,7,8,9])
        pylab.yscale('log')#, subsy=[2,3,4,5,6,7,8,9])
        #pylab.xticks(rotation=45)
        ylim = pylab.ylim()
        try:
            self.hanningFilter = True
            vmin  = self.acovf1d[self.acovfx > (self.acovfx.max()*4./5.)].mean() 
            vmin = max(vmin, ylim[0])
        except AttributeError:
            vmin = ylim[0]
        ticks, cbformat = self._colorbarTicks(vmin, ylim[1], log=True)
        xmax = max(self.xcen, self.ycen) - min(self.padx, self.pady)
        xmin = 0.5
        pylab.xlim(xmin, xmax)
        pylab.ylim(vmin, ylim[1])
        pylab.yticks(ticks)
        pylab.grid()
        pylab.xlabel('X')
        pylab.ylabel('1D ACovF', fontsize=12)
        if title!=None:
            pylab.suptitle(title, fontsize=14)
        return
