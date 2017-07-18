import numpy as np
from astropy.io import fits
import math
from scipy import signal
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import pylab


datafile = fits.open('extracted_spectra.fits')

hdu = fits.PrimaryHDU(data=[])
hdul = [hdu]
hdul2 = [hdu]


for n_hdr in [1,2,3]:
    data = datafile[n_hdr].data
    data_extracted = data[:,5:25]
    summed_spectra = np.zeros(data_extracted.shape[0],dtype=np.float)
    data_poly = []
    x_poly = []
    x = np.linspace(0,data_extracted.shape[0], data_extracted.shape[0])
    for i in range(data_extracted.shape[0]):
        summed_spectra[i] = sum(data_extracted[i])
        if i>2:
            diff = summed_spectra[i] - summed_spectra[i-1]
            diff2 = summed_spectra[i] - summed_spectra[i-2]
            diff3 = summed_spectra[i] - summed_spectra[i-3]
            if (abs(diff)+abs(diff2)+abs(diff3))<30:
                data_poly.append(summed_spectra[i])
                x_poly.append(i)
    
    hdu_ = fits.ImageHDU(data=summed_spectra,header=datafile[n_hdr].header)
    hdul.append(hdu_)
    z = np.polyfit(x_poly, data_poly, 10)
    p = np.poly1d(z)
    hdu_2 = fits.ImageHDU(data=summed_spectra-p(x),header=datafile[n_hdr].header)
    hdul2.append(hdu_2)

hdulist = fits.HDUList(hdul)
hdulist.writeto('summed_spectra.fits', clobber=True)

hdulist2 = fits.HDUList(hdul2)
hdulist2.writeto('summed_spectra_no_cont.fits', clobber=True)


def get_lambd(x,a,b,c,d,e):
    return a + b*x + c*x**2 + d*x**3 + e*x**4
    
'''    
linelist=[
    [5051,1],
    [5370,0.8],
    [5590,0.6],
    [6240,0.7]
    ]
'''
linelist=[
    [14254.351,0.017],
    [14097.5,0.13],
    [13829.51,0.04],
    [13914.4,0.05],
    [13722.34,1],
    [13682.29,0.20],
    [13626.38,0.3],
    [13603.05,0.05],
    [13507.88,0.87],
    [13410.26,0.05],
    [13370.77,0.6],
    [13316.85,0.3],
    [13276.27,0.4],
    [13231.72,0.23],
    [13217.16,0.27],
    [13012.0,0.15],
    [12960.1,0.65],
    [12936.9,0.08],
    [12806.2,0.2],
    [12705.9,0.2],
    [12491.08,0.4],
    [12459.53,0.17],
    [12442.73,0.58],
    [12406.22,0.36],
    [12359.67,0.02],
    [12346.77,0.09],
    [12143.06,0.10],
    [12115.64,0.15]
    ]


linelist2=[
    
    
    
    ]

data = hdul[1].data

N_pix = data.shape[0]

x = np.linspace(0,N_pix, N_pix)
mask = np.zeros(N_pix,dtype=np.float)
mask[500:1400] = 1.

def H (x,n):
	if n==3:
		return(-12*x+8*x**3)
	else:
		return(12-48*x**2+16*x**4)

def profile_function(l_c, sigma,flux,x,h3,h4):
    return flux * np.exp(-np.power(x - l_c, 2.) / (2 * np.power(sigma, 2.)))*(1+H((x-l_c)/(sigma*(2**0.5)),3)*h3+h4*H((x-l_c)/(sigma*(2**0.5)),4))



model = np.zeros(N_pix,dtype=np.float)
def fcn2min(params, x,data):
    """ model decaying sine wave, subtract data"""
    a = params['a']
    b = params['b']
    c = params['c']
    d = params['d']
    e = params['e']
    Q = params['Q']
    sigma = params['sigma']
    h3 = params['h3']
    h4 = params['h4']
    model = np.zeros(N_pix,dtype=np.float)
    print params.pretty_print()
    print sum(mask)
    for line in linelist:
        if get_lambd(110,a,b,c,d,e)<line[0]:
            model = model + Q*profile_function(line[0],sigma,line[1],get_lambd(x,a,b,c,d,e),h3,h4)
        #print get_lambd(x,a,b,c)
    res = np.nan_to_num(mask*(model-data)**2)
    return mask*(model-data)





params = Parameters()
params.add('a',   value= 10305.3960, vary=True)
params.add('b', value= 2.47028877, vary=True)
params.add('c', value= 0.0002842, vary=True)
params.add('d', value= -7.8345e-08, vary=True)
params.add('e', value= 0., vary=False)
params.add('Q', value= 350.0, vary=True)
params.add('sigma', value= 4.5, vary=True, min=3.)
params.add('h3',value= 0,vary=False, min=-0.0005, max=0.0005)
params.add('h4',value= 0,vary=False, min=-0.0005, max=0.0005)


minner = Minimizer(fcn2min, params, fcn_args=(x, data))
result = minner.minimize()

datafile = fits.open('extracted_spectra.fits')
datafile[1].header.append(card=('a',result.params['a'].value))
datafile[1].header.append(card=('b',result.params['b'].value))
datafile[1].header.append(card=('c',result.params['c'].value))
datafile[1].header.append(card=('d',result.params['d'].value))
datafile[1].header.append(card=('e',result.params['e'].value))



report_fit(result)

model = np.zeros(N_pix,dtype=np.float)
print sum(model)
for line in linelist:
    print result.params['Q'].value
    model = model + result.params['Q'].value*profile_function(line[0],result.params['sigma'].value,line[1],get_lambd(x,result.params['a'].value,result.params['b'].value,result.params['c'].value,result.params['d'].value,result.params['e'].value),0,0)



pylab.plot(x, data, 'k')
pylab.plot(x, model, 'r')
pylab.show()






params['a'].value=12800
params['b'].value=3.0
params['b'].vary=True
params['sigma'].vary=True
params['c'].vary=True
params['d'].vary=True
params['Q'].vary=True
data = hdul[3].data

N_pix = data.shape[0]
x = np.linspace(0,N_pix, N_pix)
mask = np.zeros(N_pix,dtype=np.float)
mask[110:400] = 1.

minner = Minimizer(fcn2min, params, fcn_args=(x, data))
result = minner.minimize()


report_fit(result)


datafile[3].header.append(card=('a',result.params['a'].value))
datafile[3].header.append(card=('b',result.params['b'].value))
datafile[3].header.append(card=('c',result.params['c'].value))
datafile[3].header.append(card=('d',result.params['d'].value))
datafile[3].header.append(card=('e',result.params['e'].value))




hdulist = fits.HDUList(datafile)
hdulist.writeto('calib_disp.fits', clobber=True)

model = np.zeros(N_pix,dtype=np.float)
print sum(model)
for line in linelist:
    print result.params['Q'].value
    model = model + result.params['Q'].value*profile_function(line[0],result.params['sigma'].value,line[1],get_lambd(x,result.params['a'].value,result.params['b'].value,result.params['c'].value,result.params['d'].value,result.params['e'].value),0,0)



pylab.plot(x, data, 'k')
pylab.plot(x, model, 'r')
pylab.show()



