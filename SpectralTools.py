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

    Note:
    pyasl.dopplerShift is much smoother than this function so use that instead. 
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



def wav_selector(wav, flux, wav_min, wav_max, verbose=False):
    """ Fast Wavelength selector between wav_min and wav_max values 
    If passed lists it will return lists.
    If passed np arrays it will return arrays
    
    """
    if isinstance(wav, list): # If passed lists
        wav_sel = [wav_val for wav_val in wav if (wav_min < wav_val < wav_max)]
        flux_sel = [flux_val for wav_val, flux_val in zip(wav,flux) if (wav_min < wav_val < wav_max)]
    elif isinstance(wav, np.ndarray):
        # Super Fast masking with numpy
        mask = (wav > wav_min) & (wav < wav_max)
        wav_sel = wav[mask]
        flux_sel = flux[mask]
        if verbose:
            print("mask=", mask)
            print("len(mask)", len(mask))
            print("wav", wav)
            print("flux", flux)
    else:
          raise TypeError("Unsupported input wav type")
    return [wav_sel, flux_sel]




# Wavelength Interplation from telluric correct
def wl_interpolation(wl, spec, ref_wl, method="scipy", kind="linear", verbose=False):
    """Interpolate Wavelengths of spectra to common WL
    Most likely convert telluric to observed spectra wl after wl mapping performed"""
    starttime = time.time()
    if method == "scipy":
        print(kind + " scipy interpolation")
        linear_interp = interp1d(wl, spec, kind=kind)
        new_spec = linear_interp(ref_wl)
    elif method == "numpy":
        if kind.lower() is not "linear":
            print("Warning: Cannot do " + kind + " interpolation with numpy, switching to linear" )
        print("Linear numpy interpolation")
        new_spec = np.interp(ref_wl, wl, spec)  # 1-d peicewise linear interpolat
    else:
        print("Method was given as " + method)
        raise("Not correct interpolation method specified")
    print("Interpolation Time = " + str(time.time() - starttime) + " seconds")

    return new_spec  # test inperpolations 







###################################################################
#                            Convolution
###################################################################


def unitary_Gauss(x, center, FWHM):
    """
    Gaussian_function of area=1
    
    p[0] = A;
    p[1] = mean;
    p[2] = FWHM;
    """
    
    sigma = np.abs(FWHM) /( 2 * np.sqrt(2 * np.log(2)) );
    Amp = 1.0 / (sigma*np.sqrt(2*np.pi))
    tau = -((x - center)**2) / (2*(sigma**2))
    result = Amp * np.exp( tau );
    
    return result

def fast_convolve(wav_val, R, wav_extended, flux_extended, fwhm_lim):
    """IP convolution multiplication step for a single wavelength value"""
    FWHM = wav_val/R
    
    index_mask = (wav_extended > (wav_val - fwhm_lim*FWHM)) &  (wav_extended < (wav_val + fwhm_lim*FWHM))
    
    flux_2convolve = flux_extended[index_mask]
    IP = unitary_Gauss(wav_extended[index_mask], wav_val, FWHM)
    
    sum_val = np.sum(IP*flux_2convolve) 
    unitary_val = np.sum(IP*np.ones_like(flux_2convolve))  # Effect of convolution onUnitary. For changing number of points
        
    return sum_val/unitary_val

def instrument_convolution(wav, flux, chip_limits, R, fwhm_lim=5.0, plot=True, verbose=True):
    """Convolution code adapted from pedros code and speed up with np mask logic"""
    
    #print("types", type(wav), type(flux), type(chip))
    #print("lengths", len(wav), len(flux), len(chip))
    
    # CRIRES HDR vals for chip limits don't match well with calibrated values (get interpolation out of range error)
    # So will use limits from the obs data instead 
    #wav_chip, flux_chip = chip_selector(wav, flux, chip)
    wav_chip, flux_chip = fast_wav_selector(wav, flux, chip_limits[0], chip_limits[1])
    #we need to calculate the FWHM at this value in order to set the starting point for the convolution
    
    FWHM_min = wav_chip[0]/R    #FWHM at the extremes of vector
    FWHM_max = wav_chip[-1]/R       
    
    #wide wavelength bin for the resolution_convolution
    wav_extended, flux_extended = fast_wav_selector(wav, flux, wav_chip[0]-fwhm_lim*FWHM_min, wav_chip[-1]+fwhm_lim*FWHM_max, verbose=False) 
    # isinstance check is ~100*faster then arraying the array again.
    if not isinstance(wav_extended, np.ndarray):
        wav_extended = np.array(wav_extended, dtype="float64") 
    if not isinstance(flux_extended, np.ndarray):
        flux_extended = np.array(flux_extended, dtype="float64")
    
    print("Starting the Resolution convolution...")
    # Predefine np array space
    flux_conv_res = np.empty_like(wav_chip, dtype="float64")
    counter = 0 
    base_val = len(wav_chip)//20   # Adjust here to change % between reports
    
    for n, wav in enumerate(wav_chip):
        # put value directly into the array
        flux_conv_res[n] = fast_convolve(wav, R, wav_extended, flux_extended, fwhm_lim)
        if(n%base_val== 0) and verbose:
            counter = counter+5
            print("Resolution Convolution at {}%%...".format(counter))
    
    #if not isinstance(flux_conv_res, np.ndarray):
    #    flux_conv_res = np.array(flux_conv_res, dtype="float64")
        
    print("Done.\n")
    
    if(plot):
        fig=plt.figure(1)
        plt.xlabel(r"wavelength [ $\mu$m ])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wav_chip, flux_chip/np.max(flux_chip), color ='k', linestyle="-", label="Original spectra")
        plt.plot(wav_chip, flux_conv_res/np.max(flux_conv_res), color ='b', linestyle="-", label="Spectrum observed at and R=%d ." % (R))
        plt.legend(loc='best')
        plt.show() 
    return [wav_chip, flux_conv_res]


