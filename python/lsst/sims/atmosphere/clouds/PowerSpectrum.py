import numpy
from scipy import interpolate

class PowerSpectrum:
    
    def __init__(self, ws, s):
        """Instantiate the PowerSpectrum object to hold the structure function for cloud variance.
        
        The windowsize (ws) and sampling rate (s) for the structure function are set here.
        The window size details the total size of the structure function (so set to be slightly
        larger than or equivalent to the fov). The sampling size details how often the power
        spectrum is sampled to create the structure function .. higher # means sharper variations."""
        self.windowsize = ws
        self.sampling = s

    def ComputeStructureFunction(self, x0=0.0, x1=1.75, y1=0.04, ymm=0.08):
        """Create the structure function based on an analytical model.

        This analytical model is determined by fitting data (done elsewhere).
        The variables x0, x1, y1, and ymm describe the structure function.
        x0=0, x1=1.75,y1=0.04, ymm=0.08 are defaults fit to the SDSS (?) structure function.
        ymm is related to the long-range grayscale extinction, while x1 and y2 are related to
        the inflection points in the structure function. """
        # define the x range
        self.xstep = self.windowsize/numpy.float(self.sampling)
        self.x = numpy.arange(0, self.sampling-self.xstep, self.xstep)
        # Calculate variable for the structure function.
        al = -(1/(x1-x0))*numpy.log(1-(y1/ymm))
        # Calculate the actual structure function.
        self.SF = ymm*(1.-numpy.exp(-al*self.x))   
        
    def ComputeCorrelationFunction(self):
        """Compute the correlation function from the structure function."""
        # following the definition astro-ph/0703157v1
        tmp = -0.5*self.SF**2
        self.correl = -numpy.min(tmp)+tmp
        # Is this actually used anywhere other than in writing out the correlation function?

    def getCorrel2D(self):
        """return a 2D correlation function"""
        # .. not sure why this method & ComputeStructureFunction should be separate if no separate arguments
        try:
            self.xstep
        except AttributeError:
            self.ComputeStructureFunction()            
        # 1D windowssize correspond to 2D diagonal
        localxstep = self.xstep/numpy.sqrt(2.)
        tmp = -.5*self.SF**2
        correltmp = -numpy.min(tmp)+tmp
        # interpolation of the correlation function on 2D real space
        interpfunc = interpolate.interp1d(self.x, correltmp, bounds_error=False, fill_value = 0.)
        # 2D sampling = 1D sampling (arbitrary)
        correl2D = numpy.zeros(self.sampling*self.sampling).reshape(self.sampling, self.sampling)
        # symmetry for the 2D correlation function
        # to define a periodic function
        # correl[size-i][size-j] = correl[i][j]
        # correl[i][size-j] = correl[i][j]
        # correl[size-i][j] = correl[i][j]
        for i in range(self.sampling/2+self.sampling%2):
            for j in range(self.sampling/2+self.sampling%2):
                co = interpfunc(numpy.sqrt((i*i+j*j))*localxstep)
                correl2D[i,j] = co
                correl2D[self.sampling-1-i,self.sampling-1-j] = co
                correl2D[i,self.sampling-1-j] = co
                correl2D[self.sampling-1-i,j] = co
        return correl2D
       
    def writeSF(self, filename):
        """write structure function in a text file"""
        f = open(filename, 'w')
        for i in range(len(self.x)):
            f.write(str(self.x[i])+'\t'+str(self.SF[i])+'\n')
        f.close()
        
    def writeCorrel(self, filename):
        """write correlation function in a text file"""
        f = open(filename, 'w')
        for i in range(len(self.x)):
            f.write(str(self.x[i])+'\t'+str(self.correl[i])+'\n')
        f.close()
        
