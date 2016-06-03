#SpectralTools.py


# Collection of useful tools for dealing with spectra:
from __future__ import division
import numpy as np

def BERVcorr(W1, Berv):
    """Barycentric Earth Radial Velocity correction from tapas 

    A wavelength W0 of the original spectrum is seen at wavelength W1 in the spectrometer stretched by Doppler effect:

    W1 / W0 = (1+ Vr/c)    =(1- Berv/c)

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



def air2vac(air):
    """ Conversion of air wavelenghts to vacuum wavelenghts.
    Adapted from wcal from Pat Hall's IRAF tasks for displaying SDSS spectra

    Input: Air wavelengths in nm

    Output: Vacuum wavelengths in nm

    """
    print("Probably best to use pyastronomy versions !!!!!!!!!!!!!!!!!!!!!!!!")
    air_in_angstroms = air * 10
   
    sigma2 = (10**8)/air**2
    n = 1 + 0.000064328 + 0.0294981/(146-sigma2) + 0.0002554/(41 - sigma2)

    vacuum_in_angstroms = air * n
    if (min(air_in_angstroms)<1600):
        print ("# WARNING! formula intended for use only at >1600 Ang!")

    return vacuum_in_angstroms/10

def vac2air(vac): 
    """ Conversion of vacuum wavelenghts to air wavelenghts.

    From http://classic.sdss.org/dr7/products/spectra/vacwavelength.html
    AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4) given in Morton (1991, ApJS, 77, 119)
    Better correction for infrared may be found in 
    http://adsabs.harvard.edu/abs/1966Metro...2...71E
    and 
    http://adsabs.harvard.edu/abs/1972JOSA...62..958P """
    print("Probably best to use pyastronomy versions !!!!!!!!!!!!!!!!!!!!!!!!")
    vac_in_angstroms = vac * 10
    air_in_angstroms = vac_in_angstroms / (1.0 + 2.735182E-4 + 131.4182 / vac_in_angstroms**2 + 2.76249E8 / vac_in_angstroms**4)
    ##Need to look at these for nir compatbile formular
    return air_in_ang/10








def fast_wav_selector(wav, flux, wav_min, wav_max, verbose=False):
    """ Faster Wavelenght selector
    
    If passed lists it will return lists.
    If passed np arrays it will return arrays
    
    Fastest is using np.ndarrays
    fast_wav_selector ~1000-2000 * quicker than wav_selector
    """

    if isinstance(wav, list): # if passed lists
        wav_sel = [wav_val for wav_val in wav if (wav_min < wav_val < wav_max)]
        flux_sel = [flux_val for wav_val, flux_val in zip(wav,flux) if (wav_min < wav_val < wav_max)]
    elif isinstance(wav, np.ndarray):
        # Super Fast masking with numpy
        mask = (wav > wav_min) & (wav < wav_max)
        if verbose:
            print("mask=", mask)
            print("len(mask)", len(mask))
            print("wav", wav)
            print("flux", flux)
        wav_sel = wav[mask]
        flux_sel = flux[mask]
    else:
          raise TypeError("Unsupported input wav type")
    return [wav_sel, flux_sel]