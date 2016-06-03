#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
# Testing for a snr of my spectra

import SpectralTools as s_tools
def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Telluric Removal')
    parser.add_argument('fname', help='Input fits file')
    #parser.add_argument('-x', '--export', default=False,
    #                    help='Export result to fits file True/False')
    parser.add_argument('-i', '--interactive_mode',  action='store_true',
                        help='Interactive mode lets you click the range '
                        + 'to calculate the signal to noise of.')
    parser.add_argument('-b', '--binsize', default=5,
                        help='SNR binsize')
    parser.add_argument('-l', '--limits', type=float, nargs=2,
                        help='Bounds for the wavelength range to '
                            +'calculate snr for in format [min, max]')
    parser.add_argument('-p', '--plot', action='store_true',
                        help="Display plots")
    args = parser.parse_args()
    return args

def snr(spectra, ddof=0):
    snr = np.mean(spectra) / np.std(spectra, ddof=ddof)
    return snr

def snrscan(Spec, width):
    snr_list = []
    for i in range(len(Spec)-width):
        #snr = np.mean(Spec[i:i+width]) / np.std(Spec[i:i+width], ddof=0)
        SNR = snr(Spec[i:i+width])
        snr_list.append(SNR)

    return snr_list

def click_snr(wl, Spec):
    """ Calculate snr in a specific range given by clicks on a plot """
    plt.figure()
    plt.plot(wl, Spec)
    plt.show(block=False)
    # points from click

    # temp values untill implement above
    point2 = np.max(wl) 
    point1 = np.min(wl) 
    map2 = wl < point2
    map1 = wl > point1
    wl_slice = wl[map1*map2]
    Spec_slice = Spec[map1*map2]
 
    # Calculate SNR on the slice
    SNR = snr(Spec_slice)

    return SNR

def simple_plot(wav, flux, label="Spectra", xscale="nm"):
    plt.plot(wav, flux, label=label)
    plt.ylabel("FLux")
    if xscale== "nm":
        plt.xlabel("Wavelenght (nm)")
    else:
        plt.xlabel("Pixel number")
    plt.show()

def main(fname, binsize=5, interactive_mode=False, limits=None, plot=False):
    # load the fits file 
    data = fits.getdata(fname)
    hdr = fits.getheader(fname)

    try:
        # Gasgano reduction
        I = data["Extracted_OPT"]
    except:
        # Dracs reduction
        I = data["Extracted_DRACS"]
   
    try:
        wl = data["Wavelength"]
        xscale = "nm"
    except:
        # no wavelength values so use pixel positions
        wl = range(len(I))
        xscale = "pixel"
        
    if interactive_mode:
        # do interactve mode stuff
        pass

    elif limits: 
        limited_wl, limited_flux = s_tools.fast_wav_selector(wl, I, limits[0], limits[1])
        limit_snr = snr(limited_flux)
        print("SNR of the spectra between the given limits {0}-{1} is {2}".format(limits[0], limits[1], limit_snr))
        if plot:
            simple_plot(wl, I, label="Spectra", xscale=xscale)
    else:
       # 
       binsize = int(binsize)

       # fancy colour plot if snr for different bin sizes verse central wavelenght value
       #
       #  b   |
       #  i   |
       #  n   |
       #      |    SNR Ratio
       #  s   |
       #  i   |
       #  z   |_________________
       #  e         Position
       #
       #

#test = [1,2,3,4,2,1,2,1,1,2,3,2,1,4,2,1,2,3,4,1,2,3,4,1,2,4,2,1,2,3,2,2.5,2,2,2,1,1,2,3,4]
#tests = [t / 100 for t in test]

# old code
    #x = range(len(test))

    #xbin = [xi + binsize/2 for xi in x]  # center snr values
    #xbin = xbin[:-binsize]
    
    #snrlist = snrscan(test, binsize)

    #plt.figure
    #plt.plot(x, test, label="data")
    #plt.plot(xbin, snrlist, label="snr sze"+ str(binsize))
    #plt.legend()
    #plt.show()

def fast_wav_selector(wav, flux, wav_min, wav_max, verbose=False):
    """ Faster Wavelenght selector
    
    If passed lists it will return lists.
    If passed np arrays it will return arrays
    
    Fastest is using np.ndarrays
    fast_wav_selector ~1000-2000 * quicker than wav_selector
    """

    if isinstance(wav, list): # if passed lists
        wav_sel = [wav_val for wav_val in wav if (wav_min < wav_val < wav_max)]
        flux_sel = [flux_val for wav_val, flux_val in zip(wav, flux) if (wav_min < wav_val < wav_max)]
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


if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    main(fname, **opts)
