"""

This file contains the definition of the fitting functions for neutrino mass priors.

Developed by: Marco Raveri

"""

# import modules:
import scipy.special as spec
import numpy         as np
import math

# measurements from PDG 2016:
Delta_m21         = 7.37*1.e-5 # eV**2
Delta_m21_sigma   = 0.17*1.e-5 # eV**2
Delta_m           = 2.50*1.e-3 # eV**2
Delta_m_sigma     = 0.04*1.e-3 # eV**2
Delta_m_inv       = 2.46*1.e-3 # eV**2
Delta_m_sigma_inv = 0.04*1.e-3 # eV**2

# minimal mass:
minimal_mass_NO = np.sqrt(Delta_m21)+np.sqrt(Delta_m+0.5*Delta_m21)
minimal_mass_IO = np.sqrt(Delta_m_inv-0.5*Delta_m21)+np.sqrt(Delta_m_inv+0.5*Delta_m21)

# minimal mass error (with error propagation):
sigma_minimal_mass_NO = np.sqrt( 1.0/(2.0*np.sqrt(Delta_m +0.5*Delta_m21))**2*Delta_m_sigma**2 +
                                (1.0/(4.0*np.sqrt(Delta_m +0.5*Delta_m21))+1.0/(2.0*np.sqrt(Delta_m21)))**2*Delta_m21_sigma**2 )
sigma_minimal_mass_IO = np.sqrt( ( 1.0/(2.0*np.sqrt(Delta_m_inv -Delta_m21/2.0)) +1.0/(2.0*np.sqrt(Delta_m_inv +Delta_m21/2.0)) )**2*Delta_m_sigma_inv**2 +
                                 ( -1.0/(4.0*np.sqrt(Delta_m_inv -Delta_m21/2.0))+1.0/(4.0*np.sqrt(Delta_m_inv +Delta_m21/2.0)) )**2*Delta_m21_sigma**2 )

# fit parameters: Normalization Mean Sigma Skewness PowerLawAmplitude PowerLawCenter PowerLawSigma, PowerLawExponent
fitting_params_M_NO = [ 0.9766, 0.06164, 0.02416, 0.003215, 1.565, 0.1336, 0.01251 , -10, minimal_mass_NO, sigma_minimal_mass_NO]
fitting_params_M_IO = [ 0.9389, 0.1039, 0.0249, 0.003634, 2.781, 0.1705, 0.01359 , -10, minimal_mass_IO, sigma_minimal_mass_IO]
fitting_params_D_NO = [ 0.9862, 0.06109, 0.017, 0.002484, 0.2275, 0.1315, 0.01603 , -16, minimal_mass_NO, sigma_minimal_mass_NO]
fitting_params_D_IO = [ 0.9695, 0.1033, 0.01792, 0.002936, 0.5281, 0.1713, 0.01876 , -16, minimal_mass_IO, sigma_minimal_mass_IO]
fitting_params_S_NO = [ 8.235e+26, -0.4099, 0.04301, 0.01624, 5.5e-27, 0.0779, 0.007016, -12, minimal_mass_NO, sigma_minimal_mass_NO]
fitting_params_S_IO = [ 1.186e+65, 5.791e-06, 0.01851, -0.006123, 1.124e-63, 0.09837, 0.003094 , -12, minimal_mass_IO, sigma_minimal_mass_IO]

# evidence results:
evidence_M_NO = 10.0**3.5756
evidence_M_IO = 10.0**2.3055
evidence_D_NO = 10.0**3.0760
evidence_D_IO = 10.0**0.9593
evidence_S_NO = 10.0**5.9995
evidence_S_IO = 10.0**3.3313

# definition of the general fitting function:

sqrt_two = np.sqrt(2.0)

def smooth_step_function( x, mu, sigma ):
    return 0.5*( 1.0 +spec.erf( (x-mu)/sqrt_two/sigma ) )

def neutrino_fitting_function( x, Nn, mean, sigma, skew, pamp, pmean, psigma, pexp, minimal, minimal_sigma ):
    # protect against zero:
    if x <= 0.0: return 0.0
    # get the minimal mass cut:
    min_cut    = smooth_step_function( x, minimal, minimal_sigma )
    # get the skewed Gaussian component:
    skew_gauss = sqrt_two/np.sqrt(math.pi)/sigma*smooth_step_function( x, mean, skew )*np.exp( -0.5*(x-mean)**2/sigma**2 )
    # get the power law:
    power_law  = pamp*smooth_step_function( x, pmean, psigma )*(x/pmean)**pexp
    # combine:
    y = Nn*min_cut*( skew_gauss+power_law )
    return y

# auxiliary definitions:
def fitting_function_M_NO( x ):
    return neutrino_fitting_function( x, *fitting_params_M_NO )
fitting_function_M_NO = np.vectorize( fitting_function_M_NO )
def fitting_function_M_IO( x ):
    return neutrino_fitting_function( x, *fitting_params_M_IO )
fitting_function_M_IO = np.vectorize( fitting_function_M_IO )
def fitting_function_D_NO( x ):
    return neutrino_fitting_function( x, *fitting_params_D_NO )
fitting_function_D_NO = np.vectorize( fitting_function_D_NO )
def fitting_function_D_IO( x ):
    return neutrino_fitting_function( x, *fitting_params_D_IO )
fitting_function_D_IO = np.vectorize( fitting_function_D_IO )
def fitting_function_S_NO( x ):
    return neutrino_fitting_function( x, *fitting_params_S_NO )
fitting_function_S_NO = np.vectorize( fitting_function_S_NO )
def fitting_function_S_IO( x ):
    return neutrino_fitting_function( x, *fitting_params_S_IO )
fitting_function_S_IO = np.vectorize( fitting_function_S_IO )

# exit if the file is directly called:
if __name__ == "__main__":

    # check the normalization:
    import scipy.integrate as integrate

    print 'Test of the normalization of the fitting functions.'
    print 'Integrating from zero to one:'
    print 'Normalization of M NO: ', integrate.quad( fitting_function_M_NO, 0, 1.0 )
    print 'Normalization of M IO: ', integrate.quad( fitting_function_M_IO, 0, 1.0 )
    print 'Normalization of D NO: ', integrate.quad( fitting_function_D_NO, 0, 1.0 )
    print 'Normalization of D IO: ', integrate.quad( fitting_function_D_IO, 0, 1.0 )
    print 'Normalization of S NO: ', integrate.quad( fitting_function_S_NO, 0, 1.0 )
    print 'Normalization of S IO: ', integrate.quad( fitting_function_S_IO, 0, 1.0 )
    print 'Integrating from zero to infinity:'
    print 'Normalization of M NO: ', integrate.quad( fitting_function_M_NO, 0, np.inf )
    print 'Normalization of M IO: ', integrate.quad( fitting_function_M_IO, 0, np.inf )
    print 'Normalization of D NO: ', integrate.quad( fitting_function_D_NO, 0, np.inf )
    print 'Normalization of D IO: ', integrate.quad( fitting_function_D_IO, 0, np.inf )
    print 'Normalization of S NO: ', integrate.quad( fitting_function_S_NO, 0, np.inf )
    print 'Normalization of S IO: ', integrate.quad( fitting_function_S_IO, 0, np.inf )
    exit(0)
