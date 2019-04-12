###################################################################################
# GM_spectrum
###################################################################################
"""
Compute GM spectrum

based on GarrettMunkMatlab by J. Klymak

http://jklymak.github.io/GarrettMunkMatlab/


"""

###################################################################################
#Load modules
###################################################################################

#for numeric functions
import numpy as np

# gamma function
from scipy.special import gamma

###################################################################################

def GMfreq(omega,f):
    '''
    Frequency part of the GM spectrum

    '''
    B = 2./np.pi * f/omega * (omega**2 - f**2)**(-0.5)

    return B

###################################################################################

def GMvert(K,f,jstar,jp,N,b,N0,I,s,t,Ef=0):
    '''
    vertical wavenumber part of the GM spectrum

    '''

    delta = jp * N / (2*b*N0)
    kstar = jstar * N / (2*b*N0)

    H = I / kstar * ( 1. + ((K - delta) / kstar)**s )**(-t/s)

    if Ef > 0:
        A10 = I / kstar * ( 1. + ((0.1 - delta) / kstar)**s )**(-t/s)
        Aa = A10 * ((K/0.1)**(-3))
        H = np.minimum(Aa,Ef*H)

    return H

#################################################

def GM_kf(omega,k,f,N,var='vel'):
    '''
    return 2D horizontal-frequency wavenumber spectra
    '''

    #param from GM76
    t = 2.
    s = 2.
    jp = 0.
    jstar = 3.

    # dimensional parameters (based on GM72)
    b = 1300.
    N0 = 5.2e-3
    E0 = 6.3e-5

    ###########

    Omega, K = np.meshgrid(omega,k)

    ###########
    
    delta = jp *2 * N * b * np.sqrt(Omega**2-f**2)
    kstar = jstar / (2*N*b) * np.sqrt(Omega**2-f**2)

    I = s * gamma(t/s) / gamma(1./s) / gamma((t-1)/s)
    H = I / kstar / ( 1. + ((K - delta) / kstar)**s )**(t/s)

    #print np.trapz(H,axis=0,dx=k[1]-k[0])

    ###########

    B = GMfreq(Omega,f)

    ###########

    if var=='vel':
        X = b**2 * N0 *N * (Omega**2 + f**2)/Omega**2

    S = B * H * X * E0 #* 2./np.pi

    # the 2./np.pi is in Klymak's code but seems to be a mistake

    return S


#################################################



def GM_mf(omega,m,f,N,var='vel',Ef=0):
    '''
    return 2D vertical-frequency wavenumber spectra
    '''

    #param from GM76
    t = 2.
    s = 2.
    jp = 0.
    jstar = 3.

    # dimensional parameters (based on GM72)
    b = 1300.
    N0 = 5.2e-3
    E0 = 6.3e-5

    ###########

    Omega, K = np.meshgrid(omega,m)

    ###########

    I = s * gamma(t/s) / gamma(1./s) / gamma((t-1)/s)
    H = GMvert(K,f,jstar,jp,N,b,N0,I,s,t,Ef=Ef)

    ###########

    B = GMfreq(Omega,f)

    ###########

    if var=='vel':
        X = b**2 * N0 *N * (Omega**2 + f**2)/Omega**2

    S = B * H * X * E0 

    return S


#################################################



   
