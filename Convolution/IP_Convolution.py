
# Convolution of spectra to a Instrument profile of resolution R.
#
# The spectra does not have to be equidistant in wavelength.

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt


def wav_selector(wav, flux, wav_min, wav_max):
    """ Faster Wavelenght selector

    If passed lists it will return lists.
    If passed np arrays it will return arrays

    """
    wav = np.asarray(wav)
    flux = np.asarray(flux)
    # Super Fast masking with numpy
    mask = (wav > wav_min) & (wav < wav_max)
    wav_sel = wav[mask]
    flux_sel = flux[mask]
    return [wav_sel, flux_sel]


def unitary_Gauss(x, center, FWHM):
    """
    Gaussian_function of area=1

    p[0] = A;
    p[1] = mean;
    p[2] = FWHM;
    """

    sigma = np.abs(FWHM) / (2 * np.sqrt(2 * np.log(2)))
    Amp = 1.0 / (sigma*np.sqrt(2*np.pi))
    tau = -((x - center)**2) / (2*(sigma**2))
    result = Amp * np.exp(tau)

    return result


def fast_convolve(wav_val, R, wav_extended, flux_extended, FWHM_lim):
    """IP convolution multiplication step for a single wavelength value"""
    FWHM = wav_val/R
    # Mask of wavelength range within 5 FWHM of wav
    index_mask = ((wav_extended > (wav_val - FWHM_lim*FWHM)) &
                  (wav_extended < (wav_val + FWHM_lim*FWHM)))

    flux_2convolve = flux_extended[index_mask]
    # Gausian Instrument Profile for given resolution and wavelength
    IP = unitary_Gauss(wav_extended[index_mask], wav_val, FWHM)

    sum_val = np.sum(IP*flux_2convolve)
    # Correct for the effect of convolution with non-equidistant postions
    unitary_val = np.sum(IP*np.ones_like(flux_2convolve))

    return sum_val / unitary_val


def IPconvolution(wav, flux, chip_limits, R, FWHM_lim=5.0, plot=True,
                  verbose=True):
    """Spectral convolution which allows non-equidistance step values"""

    # Make sure they are numpy arrays
    wav = np.asarray(wav, dtype='float64')
    flux = np.asarray(flux, dtype='float64')

    timeInit = dt.now()
    wav_chip, flux_chip = wav_selector(wav, flux, chip_limits[0],
                                       chip_limits[1])
    # We need to calculate the FWHM at this value in order to set the starting
    # point for the convolution
    FWHM_min = wav_chip[0] / R    # FWHM at the extremes of vector
    FWHM_max = wav_chip[-1] / R

    # Wide wavelength bin for the resolution_convolution
    wav_min = wav_chip[0] - FWHM_lim * FWHM_min
    wav_max = wav_chip[-1] + FWHM_lim * FWHM_max
    wav_ext, flux_ext = wav_selector(wav, flux, wav_min, wav_max)

    print("Starting the Resolution convolution...")
    # Predefine array space
    flux_conv_res = np.empty_like(wav_chip, dtype="float64")
    counter = 0
    base_val = len(wav_chip)//20   # Adjust here to change % between reports

    for n, wav in enumerate(wav_chip):
        # Put convolution value directly into the array
        flux_conv_res[n] = fast_convolve(wav, R, wav_ext, flux_ext,
                                         FWHM_lim)
        if (n % base_val == 0) and verbose:
            counter = counter + 5  # And ajust here to change % between reports
            print("Resolution Convolution at {0}%%...".format(counter))

    timeEnd = dt.now()
    print("Single-Process convolution has been completed in"
          " {}.\n".format(timeEnd-timeInit))

    if plot:
        plt.figure(1)
        plt.xlabel(r"wavelength [ nm ])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wav_chip, flux_chip / np.max(flux_chip), color='k',
                 linestyle="-", label="Original spectra")
        plt.plot(wav_chip, flux_conv_res / np.max(flux_conv_res), color='r',
                 linestyle="-", label="Spectrum observed at R={0}.".format(R))
        plt.legend(loc='best')
        plt.title(r"Convolution by an Instrument Profile ")
        plt.show()
    return [wav_chip, flux_conv_res]


if __name__ == "__main__":
    # Example useage of this convolution
    wav = np.linspace(2040, 2050, 30000)
    flux = (np.ones_like(wav) - unitary_Gauss(wav, 2045, .6) -
            unitary_Gauss(wav, 2047, .9))
    # Range in which to have the convoled values. Be careful of the edges!
    chip_limits = [2042, 2049]

    R = 1000
    convolved_wav, convolved_flux = IPconvolution(wav, flux, chip_limits, R,
                                                  FWHM_lim=5.0, plot=True,
                                                  verbose=True)
