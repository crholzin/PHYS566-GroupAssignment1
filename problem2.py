import matplotlib.pyplot as pyplot
import numpy as np
from math import *
from pylab import *
from scipy.optimize import curve_fit

# Initial variables
dx = 0.1
dt = 0.0001       # dt <= (dx)^2/(2D)
n = 101
rho = [0.0 for col in range(n)]
D = 2.0
countstep = 0
s = [100, 500, 1000, 3000, 5000]
x = np.linspace(-5,5,101)
t = []
sigma2 = []
for i in range(48,53):
    rho[i] = 2


# 5 different time
for i in range(5):
    step = s[i]
    # Finite difference formular
    while countstep < step:
        for i in range(1,n-1):
                rho[i] = rho[i]+D*dt/((dx)**2)*(rho[i+1]+rho[i-1]-2*rho[i])
        countstep += 1 
    
                
    # Curve fit
    def func(x, sigma):
        return 1.0/((2*pi*(sigma**2))**0.5)*np.exp(-(x**2)/(2*(sigma**2)))
    popt, pocv = curve_fit(func, x, rho)
    t += [dt*step]
    sigma2 += [popt[0]**2]
    
    
    # Plot the density profile after certain time step
    pyplot.figure()
    pyplot.plot(x,rho,'bo', label='diffusion after %d time steps' %step)
    pyplot.legend(prop={'size':9})
    pyplot.ylabel('density')
    pyplot.xlabel('position')
    pyplot.hold(True)
    y = 1.0/((2*pi*(popt[0]**2))**0.5)*np.exp(-(x**2)/(2*(popt[0]**2)))
    pyplot.plot(x, y, label='Gaussian fit with sigma = ' + str(popt[0]))
    pyplot.legend(prop={'size':9})
    pyplot.savefig('2b%d.png' %(i+1))
    pyplot.show()

# relation between sigma square and t
fit = polyfit(t,sigma2,1)
fit_fn = poly1d(fit)
pyplot.figure()
pyplot.plot(t,sigma2,'kx',label='experimental sigma square')
pyplot.legend(prop={'size':9})
pyplot.ylabel('sigma square')
pyplot.xlabel('t')
pyplot.xlim(0,0.6)
pyplot.ylim(0,2.5)
hold(True)
pyplot.plot(t,fit_fn(t), label = 'linear fitting: y = %.2f x' %fit[0] )
pyplot.legend(prop={'size':9})
pyplot.savefig('2b6.png')
pyplot.show()
