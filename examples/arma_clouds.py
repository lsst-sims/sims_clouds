import numpy as np
import pylab
from lsst.sims.atmosphere.clouds.Arma.ArmaSf import ArmaSf
from lsst.sims.atmosphere.clouds.Arma.Clouds import Clouds

# PS - you need statsmodels (it's a dependency for atmosphere_clouds now..)
# PPS - apparenty some values of kappa/c can give nans for the cloud. This is bad, but not the time to fix.

# Snippet to support test timing
import time
def dtime(time_prev):
   return (time.time() - time_prev, time.time())


# Instantiate cloud and arma sf objects
cloud = Clouds()
armasf = ArmaSf()

# TIMING EVALUATED ON AN OLD MAC NOTEBOOK!! You should test this on the machines you're using, speed may likely be better. 
#  python version 2.7.3 and numpy v 1.7.1, scipy v 0.12.0, statsmodels v 0.4.3

# play with pixscale (arcseconds/pix) to test relative timings for creating cloud images 
# With fov = 3.0 degrees ... pixscale X -> time to generate successive cloud images 

# pixscale = 60 "/pix (180x180 images) -> 0.1 s/visit       [for a 2 year run, with 800 visits/night avg, this is about 16.2 hours of cloud generation]
# pixscsale = 40 "/pix (270x270 images) -> 0.13 s/visit
# pixscale = 25 "/pix (432x432 images) -> 0.22 s/visit
# pixscale = 17 "/pix (636x636 images) -> 0.44 s/visit      [for a 2 year run, with 800 visits/night avg, this is about 35.6 hrs of cloud generation]
# pixscale = 16 "/pix (676x676 images) -> 0.44 s/visit
# pixscale = 15 "/pix (720x720 images) -> 0.45 s/visit
# pixscale = 14 "/pix (772x772 images) -> 1.17 s/visit (?? another bad spot)
# pixscale = 13 "/pix (831x831 images) -> 0.59 s/visit     [for a 2 year run, with 800 visits/night avg, this is about 97.3 hours of cloud generation]
# pixscale = 12 "/pix (900x900 images) -> 0.64 s/visit
# pixscale = 11 "/pix (982x982 images) -> 4.7 s/visit (??? some fov/pixscale values are just BAD)
# pixscale = 10 "/pix (1080x1080 images) -> 0.9 s/visit
# pixscale = 9 "/pix (1200x1200 images) -> 1.1 s/visit
# pixscale = 8 "/pix (1350x1350 images) -> 1.4 s/visit     [for a 2 year run, with avg 800 visits/night, this is about 227 hours of cloud generation]

pixscale = 9.

# Set up lambda_p / lambda_avg / lambda_s (lambda_s = anything, but just much smaller than lambda_*)
lambda_p = 500.
lambda_avg = 300.
lambda_s = 2.

# Set kappa and c values
num_clouds = 10
cs = np.random.uniform(low=3, high=8, size=num_clouds)
kappas = np.random.normal(loc=0.5, scale=0.3, size=num_clouds)
kappas = np.where(kappas<0, 0, kappas)

# Start timer 
t = time.time()

#fig1 = pylab.figure(1)

outRoot = 'cloudTest_'
for i, (kappa, c) in enumerate(zip(kappas, cs)):
    sftheta, sfsf = armasf.CloudSf(lambda_p, lambda_avg, lambda_s, kappa, c)
    #pylab.figure(1)
    #pylab.plot(sftheta, sfsf, label='kappa %.2f c %.2f' %(kappa, c))
    cloud.makeCloudImage(sftheta, sfsf, kappa, fov=3.0, pixscale=pixscale, oversample=1.0)    
    print i, cloud.cloudimage.mean()
    np.save(outRoot+'%02d'%i, cloud.cloudimage)
    #cloud.plotCloudImage()
    dt, t = dtime(t)
    print '# To generate next cloud image: %f seconds' %(dt)

#pylab.show()

# At each point, the cloud image itself is cloud.cloudimage and has size as below: 
print '# Size of cloud images being generated:', np.shape(cloud.cloudimage)

