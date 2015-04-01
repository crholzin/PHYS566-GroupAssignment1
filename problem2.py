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
x = np.linspace(-5,5,n)
t = []
sigma2 = []
for i in range(48,53):
    rho[i] = 2

# Define a function to fit the diffusion curve
def func(x, sigma):
    return 1.0/((2*pi*(sigma**2))**0.5)*np.exp(-(x**2)/(2*(sigma**2)))
    
# 5 different time
for j in range(5):
    step = s[j]
    # Finite difference formular
    while countstep < step:
        for k in range(1,n-1):
                rho[k] = rho[k]+D*dt/((dx)**2)*(rho[k+1]+rho[k-1]-2*rho[k])
        countstep += 1 
        
                                        
    # Curve fit
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
    pyplot.savefig('2b'+str(j+1)+'.eps')
    pyplot.show()
    pyplot.close()
    

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
pyplot.savefig('2b6.eps')
pyplot.show()
