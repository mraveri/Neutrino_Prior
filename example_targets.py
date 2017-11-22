"""

This file contains an example on how to use the neutrino fitting formula
to produce compute the confidence bound appearing in the conclusions of the
paper.

Developed by: Marco Raveri

"""

import neutrino_fitting_form as nup
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy import optimize
import math

# define cdf:
def cdf( x, func ):
    return 1.0-integrate.quad( func, 0, x )[0]/integrate.quad( func, 0, np.inf )[0]

cdf = np.vectorize(cdf);

# define the target:
target_p = 1.0-0.997

# plot the cdf in NO:
outdir = './results/'
x = np.linspace( 0.0001, 0.3, 1000 )
plt.plot( x, cdf( x, nup.fitting_function_M_NO ), label='Majorana NO' );
plt.plot( x, cdf( x, nup.fitting_function_D_NO ), label='Dirac NO' );
plt.plot( x, cdf( x, nup.fitting_function_S_NO ), label='Seesaw NO' );
plt.xlabel('$\\Sigma m_\\nu \,\,[\\, {\\rm eV}\\,]$');
plt.ylabel('$cdf$');
plt.legend();
plt.yscale('log')
plt.savefig(outdir+'2_cdf_NO.pdf')
plt.clf()

# plot the cdf in IO:
outdir = './results/'
x = np.linspace( 0.0001, 0.3, 1000 )
plt.plot( x, cdf( x, nup.fitting_function_M_IO ), label='Majorana IO' );
plt.plot( x, cdf( x, nup.fitting_function_D_IO ), label='Dirac IO' );
plt.plot( x, cdf( x, nup.fitting_function_S_IO ), label='Seesaw IO' );
plt.xlabel('$\\Sigma m_\\nu \,\,[\\, {\\rm eV}\\,]$');
plt.ylabel('$cdf$');
plt.legend();
plt.yscale('log')
plt.savefig(outdir+'2_cdf_IO.pdf')
plt.clf()

# get the 3 sigma values:
print 'Three sigma values:'
MNO = optimize.brentq( lambda x: cdf( x, nup.fitting_function_M_NO )-target_p, 0.0001, 0.5)
MIO = optimize.brentq( lambda x: cdf( x, nup.fitting_function_M_IO )-target_p, 0.0001, 0.5)

print 'Majorana NO : ', round(MNO,2)
print 'Majorana IO : ', round(MIO,2)

DNO = optimize.brentq( lambda x: cdf( x, nup.fitting_function_D_NO )-target_p, 0.0001, 0.5)
DIO = optimize.brentq( lambda x: cdf( x, nup.fitting_function_D_IO )-target_p, 0.0001, 0.5)

print 'Dirac NO    : ', round(DNO,2)
print 'Dirac IO    : ', round(DIO,2)

SNO = optimize.brentq( lambda x: cdf( x, nup.fitting_function_S_NO )-target_p, 0.0001, 0.5)
SIO = optimize.brentq( lambda x: cdf( x, nup.fitting_function_S_IO )-target_p, 0.0001, 0.5)

print 'Seesaw NO   : ', round(SNO,2)
print 'Seesaw IO   : ', round(SIO,2)
