import numpy as np
from astropy.io import fits
import math
from scipy import signal
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
from scipy import ndimage
import pylab
datafile = fits.open('slit2_gyjx_g_4s_flat-1.fits')
data = datafile[0].data
data_in = data
#corr_func = np.ones((1,34), dtype=np.float)
data = -data
#corr = signal.correlate2d(data,corr_func , boundary='symm', mode='same')



mask = np.ones((data.shape[0],data.shape[1]),dtype=np.float)
#mask[1391:1874,694:908] = 1.
k = [0.3587142,0.30862998,0.42492089]
b = [762.8660,531.0,1110.35108]
x = [[],[],[]]
y = [[],[],[]]
y2 = []
for i in range(data.shape[0]):
    for j in [0,1,2]:
        x_ = k[j]*i + b[j]
        inpt = np.ones(data.shape[1],dtype=np.float)
        window = 60
        inpt[int(x_-window):int(x_+window)] = data[i,int(x_-window):int(x_+window)]**4
        cm = ndimage.measurements.center_of_mass(inpt)[0]
        x[j].append(i)
        y[j].append(cm)
        y2.append(sum(inpt)/80.)
    
res = []    
x_poly = x[0][400:1600]
y_poly = y[0][400:1600]
z = np.polyfit(np.array(x_poly), np.array(y_poly), 3)
res.append(z)
print z 
p = np.poly1d(z)
pylab.plot(x[0], y[0]-p(x[0]), 'k')
x_poly = x[1][700:1900]
y_poly = y[1][700:1900]
z = np.polyfit(np.array(x_poly), np.array(y_poly), 3)
res.append(z)
print z 
p = np.poly1d(z)
pylab.plot(x[1], y[1]-p(x[1]), 'r')
x_poly = x[2][140:660]
y_poly = y[2][140:660]
z = np.polyfit(np.array(x_poly), np.array(y_poly), 3)
res.append(z)
print z 
p = np.poly1d(z)
pylab.plot(x[2], y[2]-p(x[2]), 'g')
pylab.show()



flat = data_in
flat[flat==0]=-1

hdu = fits.PrimaryHDU(data=[])
hdul = [hdu]

datafile = fits.open('slit2_gyjx_g_4s-1.fits')
data = datafile[0].data
data = data / flat
model = np.zeros((data.shape[0],data.shape[1]),dtype=np.float)
for i in range(data.shape[0]):
    x_ = p(i)
    model[i,int(x_-19):int(x_+19)] = 300.
hdu = fits.PrimaryHDU(data=data-model)
hdulist = fits.HDUList([hdu])
hdulist.writeto('order_places_mc.fits', clobber=True)
hdu = fits.PrimaryHDU(data=[])
hdul = [hdu]

for j in [0,1,2]:
    reduced_spectra = np.zeros((data.shape[0],36),dtype=np.float)
    for i in range(data.shape[0]):
        p = np.poly1d(res[j])
        x = p(i)
        reduced_spectra[i,:] = data[i,int(x-18):int(x+18)]
    hdu_ = fits.ImageHDU(data=reduced_spectra)
    hdu_.header.append(card=('c',res[j][0]))
    hdu_.header.append(card=('k',res[j][1]))
    hdu_.header.append(card=('b',res[j][2]))
    hdul.append(hdu_)
hdulist = fits.HDUList(hdul)
hdulist.writeto('extracted_spectra.fits', clobber=True)