import numpy as np
from astropy.io import fits
import math
from scipy import signal
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit


datafile = fits.open('slit2_gyjx_g_4s_flat-1.fits')
data = datafile[0].data
data_in = data
#corr_func = np.ones((1,34), dtype=np.float)
data = -data
#corr = signal.correlate2d(data,corr_func , boundary='symm', mode='same')



mask = np.ones((data.shape[0],data.shape[1]),dtype=np.float)
#mask[1391:1874,694:908] = 1.

def fcn2min(params, x, data):
    """ model decaying sine wave, subtract data"""
    k1 = params['k1']
    b1 = params['b1']
    c1 = params['c1']
    Q1 = params['Q1']
    k2 = params['k2']
    b2 = params['b2']
    Q2 = params['Q2']
    k3 = params['k3']
    b3 = params['b3']
    Q3 = params['Q3']
    k4 = params['k4']
    b4 = params['b4']
    Q4 = params['Q4']
    model = np.zeros((data.shape[0],data.shape[1]),dtype=np.float)
    print params.pretty_print()
    for i in range(data.shape[0]):
        x_ = c1*i**2+k1 *i + b1
        model[i,int(x_-19):int(x_+19)] = Q1*1
        x2_ = k2 *i + b2
        model[i,int(x2_-19):int(x2_+19)] = Q2*1
        x3_ = k3 *i + b3
        model[i,int(x3_-19):int(x3_+19)] = Q3*1
        x4_ = k4 *i + b4
        model[i,int(x4_-19):int(x4_+19)] = Q4*1
    res = (model-data)**2*mask
    return res


params = Parameters()
params.add('k1',   value= 0.360, vary=True)
params.add('b1', value= 764.01, vary=True)
params.add('c1', value= 4.3e-6, vary=True)
params.add('Q1', value= 345.01, vary=True)
params.add('k2',   value= 0.30862998, vary=False)
params.add('b2', value= 531., vary=False)
params.add('c2', value= 0, vary=False)
params.add('Q2', value= 132., vary=False)
params.add('k3',   value= 0.42492089, vary=False)
params.add('b3', value= 1110.35108, vary=False)
params.add('c3', value= 0, vary=False)
params.add('Q3', value= 205.646982, vary=False)
params.add('k4',   value= 0.27, vary=False)
params.add('b4', value= 350., vary=False)
params.add('c4', value= 0, vary=False)
params.add('Q4', value= 66., vary=False)
#print get_lambd(x,7000,-0.3,0)



# do fit, here with leastsq model
minner = Minimizer(fcn2min, params, fcn_args=( [], data))
result = minner.minimize(method='nelder')

report_fit(result)
params = result.params

'''
params['k1'].vary = False
params['b1'].vary = False
params['c1'].vary = False
params['Q1'].vary = False
params['k2'].vary = True
params['b2'].vary = True
params['Q2'].vary = True

minner = Minimizer(fcn2min, params, fcn_args=( [], data))
result = minner.minimize(method='nelder')

report_fit(result)
params = result.params

params['k2'].vary = False
params['b2'].vary = False
params['Q2'].vary = False
params['k3'].vary = True
params['b3'].vary = True
params['Q3'].vary = True

minner = Minimizer(fcn2min, params, fcn_args=( [], data))
result = minner.minimize(method='nelder')

report_fit(result)
params = result.params


mask = np.zeros((data.shape[0],data.shape[1]),dtype=np.float)
mask[1391:1874,694:908] = 1.

params['k3'].vary = False
params['b3'].vary = False
params['Q3'].vary = False
params['k4'].vary = True
params['b4'].vary = True
params['Q4'].vary = True

minner = Minimizer(fcn2min, params, fcn_args=( [], data))
result = minner.minimize(method='nelder')

report_fit(result)
params = result.params
'''
model = np.zeros((data.shape[0],data.shape[1]),dtype=np.float)
k1 = result.params['k1'].value
b1 = result.params['b1'].value
c1 = result.params['c1'].value
Q1 = result.params['Q1'].value
k2 = result.params['k2'].value
b2 = result.params['b2'].value
Q2 = result.params['Q2'].value
k3 = result.params['k3'].value
b3 = result.params['b3'].value
Q3 = result.params['Q3'].value
k4 = result.params['k4'].value
b4 = result.params['b4'].value
Q4 = result.params['Q4'].value
for i in range(data.shape[0]):
    x_ = c1*i**2+k1 *i + b1
    model[i,int(x_-19):int(x_+19)] = Q1*1
    x2_ = k2 *i + b2
    model[i,int(x2_-19):int(x2_+19)] = Q2*1
    x3_ = k3 *i + b3
    model[i,int(x3_-19):int(x3_+19)] = Q3*1
    x4_ = k4 *i + b4
    model[i,int(x4_-19):int(x4_+19)] = Q4*1
hdu = fits.PrimaryHDU(data=model)
hdulist = fits.HDUList([hdu])
hdulist.writeto('order_places.fits', clobber=True)


flat = data_in
flat[flat==0]=-1

hdu = fits.PrimaryHDU(data=[])
hdul = [hdu]

datafile = fits.open('slit2_gyjx_g_4s-1.fits')
data = datafile[0].data
data = data / flat
for num_order in [1,2,3,4]:
    reduced_spectra = np.zeros((data.shape[0],36),dtype=np.float)
    k = params['k'+str(num_order)].value
    b = params['b'+str(num_order)].value
    c = params['c'+str(num_order)].value
    Q = params['Q'+str(num_order)].value
    for i in range(data.shape[0]):
        x =c*i**2+ k *i + b
        reduced_spectra[i,:] = data[i,int(x-18):int(x+18)]
    hdu_ = fits.ImageHDU(data=reduced_spectra)
    hdu_.header.append(card=('k',k))
    hdu_.header.append(card=('b',b))
    hdu_.header.append(card=('Q',Q))
    hdul.append(hdu_)
hdulist = fits.HDUList(hdul)
hdulist.writeto('extracted_spectra.fits', clobber=True)
    
