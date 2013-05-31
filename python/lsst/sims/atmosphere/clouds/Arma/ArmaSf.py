"""
ArmaSf - a class to calculate an ARMA-based structure function for the clouds.
For more information on the basis of this structure function for clouds, please see 
http://www.lsstcorp.org:8090/display/CAL/Cloud+Simulation+Documentation, in particular Tim's two mini-papers uploaded there.

Requires statsmodels, which in turns requires pandas and cython. 

Usage:
  Instantiate class (a = ArmaSf())
  Generate structure function using CloudSf (a.CloudSf(lambda_p, lambda_avg, lambda_s, kappa, c))
     (for each visit). 

"""

from math import *
import numpy
from statsmodels.tsa.arima_process import arma_acf

class ArmaSf():
    def __init__(self):
        return

    def CloudSf(self, lambda_p, lambda_avg, lambda_s, kappa, c, acfFileName=None, x_max=1000.0):
        """
        Set cloud ACF and the variance in the ACF that depend on lambda_p/lambda_avg (using ARMA process). 

        lambda_p - physical length scale of AR process
           the characteristic physical scale of the clouds (m)
           a typical value might be around 500
        lambda_avg - physical length over which AR process is to be integrated
           the length scale over which the structure is averaged by a combination of finite aperture size plus wind (m)
           a typical value might be around 300
        lambda_s - physical sampling interval of AR process
           a nuisance value, should be small compard to lambda_p 

        kappa is the average cloud layer extinction in mags
        c is the multiplier that gives the sigma for the cloud extinction (0.3 < c < 0.8)
          (sigma = c * kappa - this scales the overall structure function)

        If using 'acFileName' to read in an ACF calculated from SF_ARMA.R (using R), then 
        the acf from acfFileName must have been written by the ACF_Write function of SF_ARMA.R
        with the same values of lambda_p, lambda_s, lambda_avg
        x_max is the maximum size to calculate the SF over, and only relevant for calculating the 
        acf when acfFileName (where the ARMA ACF is calculated by R) is used. (1000 ~ 6 degrees)
        """    
        self.lambda_p = float(lambda_p)
        self.lambda_avg = float(lambda_avg)
        self.lambda_s = float(lambda_s)
        # Calculate (or read) ACF
        if acfFileName:
            self.acf = numpy.loadtxt(acfFileName)
        else:
            self.acf = self._ACF_ARMA(x_max)
        # Calculate cloud variance in ARMA model
        self.var = self._CloudVar()
        # Scale ACF for desired cloud variance. 
        SFtheta = self.acf[:,0]
        var = (c*kappa)**2 * self.var
        SFsf = var * (1 - self.acf[:,1])
        SFsf = numpy.sqrt(SFsf)
        return SFtheta, SFsf


    def _ACF_ARMA(self, x_max):
        """Calculate the ACF using an ARMA function (from statsmodels)."""
        nma = round(self.lambda_avg/self.lambda_s)
        ma_coeff = numpy.ones((nma))
        lag_max = ceil(x_max / self.lambda_s)
        a1 = exp(-self.lambda_s / self.lambda_p)
        # Use the statsmodels ARMA function
        armaAcf = arma_acf([1,-a1], ma_coeff, nobs=lag_max)
        x = numpy.arange(0, x_max, self.lambda_s)
        # Set self.acf - a 2-d array of theta/ACF
        acf = numpy.empty((len(x), 2))
        # convert distance x (meters) in the acf to angular scale in degrees - clouds assumed
        # at 10 km (thus dividing by 10,000)
        acf[:,0] = x / (10000.0*numpy.pi/180.0)
        acf[:,1] = armaAcf
        return acf

    def _CloudVar(self):
        """
        Map a CAR(1) process to an AR(1), then average over lambda_avg
        and calculate variance
        
        lambda_p - physical length scale of AR process
        lambda_s - physical sampling interval of AR process
        lambda_avg - physical length over which AR process is to be integrated
        """    
        # q is MA order of ARMA(1,q)
        q = int(round(self.lambda_avg/self.lambda_s))
        a = exp(-self.lambda_s / self.lambda_p)    
        (var, var_ratio) = self._ARMAvar(q, a)
        # This variance is a multiple of the variance of the noise driving the
        # AR(1) model.   This variance, in turn, is a multiple of the underlying
        # measurement variance, with the relationship given in Gillespie 96
        var = var * (1. - exp(-2*self.lambda_s / self.lambda_p))/2
        #    print q, a
        return var

    def _ARMAvar(self, q, a):
        """
        Calculate the variance of an ARMA(1,q) model, where the theta terms are
        all = 1/(q+1)
        Based on Brockwell & Davis 96 eqn 3.2.3
        Modified to have theta_0 = 1/(q+1) instead of = 1
        """    
        if q>0:
            theta = 1./(q+1.)
            psi_jm1 = 1./(q+1.)
            var = psi_jm1**2
            for j in range(q):
                psi_j = a*psi_jm1 + theta
                var += psi_j**2
                psi_jm1 = psi_j
        else:
            var = 1
        # Now, sum up terms from q to infinity
        var += 1./(1.-a**2) - (1. - a**(2*(q+1)))/(1.-a**2)
        # var_ratio is ratio of variance to that of AR(1) model
        var_ratio = var*(1.-a**2)
        return (var, var_ratio)

    
    
