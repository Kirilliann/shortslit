import numpy as np
from astropy.io import fits
import math
from scipy import signal
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import pylab

datafile = fits.open('order_places.fits')
flat = datafile[0].data
datafile = fits.open('slit2_gyjx_g_4s_flat-1.fits')
data = datafile[0].data

hdul = []
hdu_ = fits.PrimaryHDU(data=data+flat)
hdul.append(hdu_)
hdulist = fits.HDUList(hdul)
hdulist.writeto('residuals.fits', clobber=True)