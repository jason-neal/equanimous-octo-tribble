#SpectralTools.py


# Collection of useful tools for dealing with spectra:



def BERVcorr(W1, Berv):
 """Barycentric Earth Radial Velocity correction from tapas 

A wavelength W0 of the original spectrum is seen at wavelength W1 in the spectrometer stretched by Doppler effect:

W1 / W0 = (1+ Vr/c)	=(1- Berv/c)

Therefore, if there is a BERV correction in the pipeline of the spectrometer to recover W0, the measure wavelength W1 has be transformed in W0 with the formula: 

W0 = W1 (1+ Berv/c) 

Inputs:
W1 - the spctrum in spectrograph frame.
Berv - The Barycentric Earth Radial Velocity value for this observation in km/s

Output:
 W0. - BERV corrected wavelength values 
 """
 c = 299792.458 # km/s
return W1 * (1 + Berv/c)  

#    return None
def inverse_BERV(W0, Berv):
   """Obtain un-BERV corrected wavelengths """
    c = 299792.458 # km/s
    return W1 * (1 - Berv/c)
    
#def dopplershift():  # standard doppler correction
#  return None


def Convolution():  # IP convolution ?
    return None