import numpy as np
from astropy.io import fits
import math
from scipy import signal
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import pylab


flatfile = fits.open('slit2_gyjx_g_4s_flat-1.fits')
flat = flatfile[0].data
datafile = fits.open('hip90967-GYJX-G-20170717203713-2.fits')
data = datafile[0].data
calib = fits.open('calib_disp.fits')
N_poly = calib[1].header['N_poly']
coeff = np.zeros(N_poly, dtype=np.float)
for i in range(N_poly):
    coeff[i] = calib[1].header['c'+str(i)]
p = np.poly1d(coeff)
window = 18
N_pix = data.shape[0]
extracted_spectra = np.zeros((data.shape[0],36),dtype=np.float)
flat[flat==0]=-1
data = data /flat
for i in range(N_pix):
    x_ = p(i)
    extracted_spectra[i,:] = data[i,int(x_-window):int(x_+window)]


pylab.plot(range(N_pix), extracted_spectra[:,14], 'r')
pylab.plot(range(N_pix), extracted_spectra[:,19], 'g')
pylab.plot(range(N_pix), extracted_spectra[:,24], 'b')

pylab.show()
hdu = fits.PrimaryHDU(data=extracted_spectra)
hdul = [hdu]
hdulist = fits.HDUList(hdul)
hdulist.writeto('spectrum_w_wl.fits', clobber=True)