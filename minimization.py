from scipy.optimize import minimize, rosen, rosen_der
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
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

datafile = fits.open('summed_spectra_no_cont.fits')
data_in = datafile[1].data
data = np.zeros(data_in.shape[0],dtype=np.float)
mask = np.zeros(data_in.shape[0],dtype=np.float)
for j in range(data_in.shape[0]):
    data[j] = data_in[j]
    if (j>400) & (j<1500):
        mask[j]=1.
    else:
        data[j] = 0.
a = 10462.5174
b = 2.25271097
c = 0.00038127
d = -8.8837e-08
data[data<0] = 0.
def x_to_l(x):
    return a+b*x+c*x**2+d*x**3
    
Q = 337.192328
def profile_function(lsf, l_c,flux):
    model = np.zeros(data.shape[0],dtype=np.float)
    for i in range(model.shape[0]):
        x_ = x_to_l(i)
        step = b + 2*c*i+3*d*i**2
        arg = int((3./step)*(x_ - l_c) + len(lsf)/2.)
        if (arg>0) & (arg<len(lsf)):
            model[i] = flux * lsf[int(arg)]
    return model





model = np.zeros(data.shape[0],dtype=np.float)
def func2min(a):
    model = np.zeros(data.shape[0],dtype=np.float)
    for line in linelist:
        model = model + Q*profile_function(a,line[0],line[1])
    print max(a)
    return sum((model-data)**2/model) + 0*(sum(a**2)+sum(np.diff(a)**2))


x0 = np.ones(60,dtype=np.float)

#x = np.linspace(0,60, 60)
#x0 = np.exp(-np.power(x - len(x)/2, 2.) / (2 * np.power(4.5, 2.)))

'''
print func2min(x0)
import pylab
pylab.plot(range(1400), func2min(x0), 'r')
#pylab.plot(x, fcn2min(params, x,data), 'r')
pylab.show()
'''


model = np.zeros(data.shape[0],dtype=np.float)
for line in linelist:
    model = model + Q*profile_function(x0,line[0],line[1])


pylab.step(range(len(model)), data, 'r')
pylab.step(range(len(model)), model, 'k')

pylab.show()
''''''
res = minimize(func2min, x0, method='Nelder-Mead', tol=1e-4)
print res.x

pylab.step(range(len(res.x)), res.x, 'k')
pylab.show()
#pylab.plot(range(len(res.x)), res.x, 'r')




'''

'''